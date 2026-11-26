#include "pubchem.hpp"

#include <cstdio>
#include <string>

#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/process.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/core/text.hpp"
#include "mdbridge/core/value.hpp"
#include "mdbridge/io/cif.hpp"

namespace mdbridge::app {

namespace {

using io::DictionaryError;

// f"{x:g}"
std::string formatG(double x) {
    char text[64];
    std::snprintf(text, sizeof text, "%g", x);
    return text;
}

// str.isdigit() of the whole (non-empty) text
bool allDigits(std::string_view text) {
    const std::u32string codes = core::decodeUtf8Surrogateescape(text);
    if (codes.empty()) {
        return false;
    }
    for (const char32_t c : codes) {
        if (!core::isPythonDigit(c)) {
            return false;
        }
    }
    return true;
}

// Python's truth value of a JSON value (entry.get(key))
bool truthy(const core::Value& v) {
    if (v.isNull()) return false;
    if (const bool* b = v.get<bool>()) return *b;
    if (const std::int64_t* i = v.get<std::int64_t>()) return *i != 0;
    if (const double* d = v.get<double>()) return *d != 0.0;
    if (const std::string* s = v.get<std::string>()) return !s->empty();
    if (const core::Array* a = v.get<core::Array>()) return !a->empty();
    if (const core::Table* t = v.get<core::Table>()) return !t->empty();
    return true;
}

}  // namespace

std::string pubchemRoot() {
    std::string url = core::process::environmentVariable(kPubchemUrlVariable);
    while (!url.empty() && url.back() == '/') {
        url.pop_back();
    }
    return url.empty() ? std::string(kPubchemRest) : url;
}

std::string percentEncoded(std::string_view text) {
    static const char* const hex = "0123456789ABCDEF";
    std::string out;
    for (const char ch : text) {
        const auto c = static_cast<unsigned char>(ch);
        if ((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9') || c == '_' || c == '.' ||
            c == '-' || c == '~') {
            out += static_cast<char>(c);
        } else {
            out += '%';
            out += hex[c >> 4];
            out += hex[c & 15];
        }
    }
    return out;
}

std::string fetchPubchemSmiles(std::string_view cidOrName, double timeoutSeconds, const HttpGet& get,
                               const std::string& root) {
    const std::string ident = core::pythonStrip(cidOrName);
    if (ident.empty()) {
        throw DictionaryError("no PubChem CID or name given");
    }
    if (!get) {
        throw DictionaryError("PubChem lookups need the mdbridge program (this one was built without network access)");
    }
    const std::string path = allDigits(ident) ? "cid/" + ident : "name/" + percentEncoded(ident);
    const std::string quoted = core::pythonStringRepr(ident);
    std::string problem;
    for (const std::string_view property : kPubchemSmilesProperties) {
        const std::string url = root + "/compound/" + path + "/property/" + std::string(property) + "/JSON";
        const HttpReply reply =
            get(url, {{"User-Agent", "MDBridge dict pubchem"}, {"Accept", "application/json"}}, timeoutSeconds);
        if (reply.outcome == HttpReply::Outcome::TimedOut) {
            throw DictionaryError("PubChem did not answer within " + formatG(timeoutSeconds) + " s");
        }
        if (reply.outcome == HttpReply::Outcome::Failed) {
            throw DictionaryError("could not reach PubChem (" + reply.error +
                                  "); check the network connection / proxy");
        }
        if (reply.status >= 400) {  // urllib's HTTPError (2xx and redirects are followed answers)
            if (reply.status == 404) {
                throw DictionaryError("PubChem has no compound " + quoted);
            }
            if (reply.status == 400) {  // unknown property name on this server: try the older one
                problem = "PubChem rejected the request for " + std::string(property) + " (HTTP 400)";
                continue;
            }
            throw DictionaryError("PubChem request failed: HTTP " + std::to_string(reply.status) + " " + reply.reason);
        }
        const core::Array* props = nullptr;
        try {
            if (core::encodeUtf8(core::decodeUtf8(reply.body)) != reply.body) {
                throw core::ValueError("not UTF-8");
            }
            const core::Value data = core::parseJson(reply.body, "PubChem");
            const core::Value* table = data.find("PropertyTable");
            const core::Value* list = table ? table->find("Properties") : nullptr;
            props = list ? list->get<core::Array>() : nullptr;
            if (props == nullptr) {
                throw core::ValueError("no properties");
            }
            for (const core::Value& entry : *props) {
                if (entry.get<core::Table>() == nullptr) {
                    throw core::ValueError("no property table");
                }
            }
            for (const core::Value& entry : *props) {
                for (const std::string_view key : kPubchemSmilesProperties) {
                    const core::Value* value = entry.find(key);
                    if (value != nullptr && truthy(*value)) {
                        return core::pythonStr(*value);
                    }
                }
            }
        } catch (const std::exception&) {
            throw DictionaryError("PubChem returned an unreadable answer");
        }
        problem = "PubChem's answer for " + quoted + " holds no SMILES";
    }
    throw DictionaryError(problem.empty() ? "no SMILES for " + quoted : problem);
}

}  // namespace mdbridge::app
