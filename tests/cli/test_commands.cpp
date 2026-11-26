// The MDATK commands of the command line (init, info, list, analyze, export-pdb, density, and the
// chem commands pdb, build, dict with a PubChem stand-in), run
// in-process on the inputs of tests/data/cli/cli_cases.json and compared with what MDATK v0.0.4.1
// printed and wrote there (tests/parity/cli_parity.py --write; normalised the same way: {work},
// {time}, {stamp}, the analysis cache's {manifest} and {key}, its files compared by meaning).
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <map>
#include <optional>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <CLI/CLI.hpp>

#include "commands.hpp"
#include "mdbridge/analysis/cache.hpp"
#include "mdbridge/analysis/pipeline.hpp"
#include "mdbridge/core/clock.hpp"
#include "mdbridge/core/hash.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/panels.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/project.hpp"
#include "mdbridge/core/project_file.hpp"
#include "mdbridge/core/text.hpp"
#include "pubchem.hpp"
#include "test_support.hpp"

namespace mdbridge {
namespace {

namespace fs = std::filesystem;
using core::Value;
using test::kData;

const Value& field(const Value& v, std::string_view key) {
    const Value* f = v.find(key);
    if (f == nullptr) {
        throw std::runtime_error("missing key " + std::string(key));
    }
    return *f;
}

std::string text(const Value& v, std::string_view key) { return field(v, key).stringOr(); }

std::string replaceAll(std::string s, const std::string& from, const std::string& to) {
    for (std::size_t pos = 0; (pos = s.find(from, pos)) != std::string::npos; pos += to.size()) {
        s.replace(pos, from.size(), to);
    }
    return s;
}

struct CommandRun {
    int exit = 0;
    std::string out;
    std::string err;
};

// ---- the PubChem stand-in (cli_parity.py PubchemStandIn, PUBCHEM_ANSWERS; owner decision D7) ----

// The root URL the cases' {pubchem} stands for; nothing listens there: the stand-in transport answers.
constexpr const char* kPubchemStandIn = "http://pubchem.stand-in/rest/pug";

std::string jsonText(const Value& v) { return core::toJson(v, 0); }

// GET {root}/compound/<cid/N | name/NAME>/property/<PROP>/JSON as cli_parity.py's server answers it
app::HttpReply pubchemStandIn(const std::string& url, const app::HttpHeaders& headers, double timeoutSeconds) {
    bool acceptsJson = false;
    for (const auto& [name, value] : headers) acceptsJson = acceptsJson || (name == "Accept" && value == "application/json");
    EXPECT_TRUE(acceptsJson) << url;
    app::HttpReply reply;
    if (!url.starts_with(kPubchemStandIn)) {
        reply.outcome = app::HttpReply::Outcome::Failed;
        reply.error = "Connection refused";
        return reply;
    }
    const std::string adenosineL = "C1=NC(=C2C(=N1)N(C=N2)[C@@H]3[C@H]([C@H]([C@@H](O3)CO)O)O)N";
    const std::string adenosineD = "C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)CO)O)O)N";
    const auto properties = [](const Value& entry) {
        return jsonText(Value(core::Table{{"PropertyTable", core::Table{{"Properties", core::Array{entry}}}}}));
    };
    struct Answer {
        int status;
        std::string body;
        double wait;
    };
    const std::map<std::pair<std::string, std::string>, Answer> answers{
        {{"cid/1001", "SMILES"}, {200, properties(core::Table{{"CID", 1001}, {"SMILES", adenosineL}}), 0}},
        {{"name/L-adenosine", "SMILES"}, {400, jsonText(core::Table{{"Fault", core::Table{{"Code", "PUGREST.BadRequest"}}}}), 0}},
        {{"name/L-adenosine", "IsomericSMILES"},
         {200, properties(core::Table{{"CID", 1001}, {"IsomericSMILES", adenosineL}}), 0}},
        {{"name/D-adenosine", "SMILES"}, {200, properties(core::Table{{"CID", 60961}, {"SMILES", adenosineD}}), 0}},
        {{"name/ethanol", "SMILES"}, {200, properties(core::Table{{"CID", 702}, {"SMILES", "CCO"}}), 0}},
        {{"cid/500", "SMILES"}, {500, "server failure", 0}},
        {{"name/garbage", "SMILES"}, {200, "not json", 0}},
        {{"name/empty", "SMILES"}, {200, properties(core::Table{{"CID", 1}}), 0}},
        {{"name/empty", "IsomericSMILES"}, {200, properties(core::Table{{"CID", 1}}), 0}},
        {{"name/slow", "SMILES"}, {200, properties(core::Table{{"CID", 2}, {"SMILES", "C"}}), 3.0}},
    };
    static const std::regex path{R"(^/compound/((?:cid|name)/[^/]+)/property/([^/]+)/JSON$)"};
    std::smatch m;
    const std::string rest = url.substr(std::string(kPubchemStandIn).size());
    Answer answer{404, jsonText(core::Table{{"Fault", core::Table{{"Code", "PUGREST.NotFound"}}}}), 0};
    if (std::regex_match(rest, m, path)) {
        std::string compound = m[1].str();  // urllib.parse.unquote
        std::string unquoted;
        for (std::size_t i = 0; i < compound.size(); ++i) {
            if (compound[i] == '%' && i + 2 < compound.size()) {
                unquoted += static_cast<char>(std::stoi(compound.substr(i + 1, 2), nullptr, 16));
                i += 2;
            } else {
                unquoted += compound[i];
            }
        }
        if (const auto it = answers.find({unquoted, m[2].str()}); it != answers.end()) answer = it->second;
    }
    if (answer.wait > timeoutSeconds) {
        reply.outcome = app::HttpReply::Outcome::TimedOut;
        return reply;
    }
    static const std::map<int, std::string> reasons{{200, "OK"}, {400, "Bad Request"}, {404, "Not Found"},
                                                    {500, "Internal Server Error"}};
    reply.outcome = app::HttpReply::Outcome::Answered;
    reply.status = answer.status;
    reply.reason = reasons.at(answer.status);
    reply.body = answer.body;
    return reply;
}

const app::NetworkBackend& pubchemNetwork() {
    static const app::NetworkBackend network{&pubchemStandIn};
    return network;
}

// The command line in `folder` (the working directory meanwhile).
CommandRun runIn(const fs::path& folder, const std::vector<std::string>& args, bool terminal = false) {
    std::vector<const char*> argv{"mdbridge"};
    for (const std::string& a : args) {
        argv.push_back(a.c_str());
    }
    std::ostringstream out;
    std::ostringstream err;
    app::Console console{out, err, terminal};
    const fs::path previous = fs::current_path();
    fs::current_path(folder);
    CommandRun run;
    try {
        run.exit = app::runCoreCommandLine(static_cast<int>(argv.size()), argv.data(), console, nullptr,
                                           &pubchemNetwork());
    } catch (...) {
        fs::current_path(previous);
        throw;
    }
    fs::current_path(previous);
    run.out = out.str();
    run.err = err.str();
    return run;
}

// Environment variables set for a while ("" is set, but empty, as Python's os.environ sets it).
class ScopedEnvironment {
public:
    ScopedEnvironment() = default;
    ScopedEnvironment(const ScopedEnvironment&) = delete;
    ScopedEnvironment& operator=(const ScopedEnvironment&) = delete;
    ~ScopedEnvironment() {
        for (auto it = m_saved.rbegin(); it != m_saved.rend(); ++it) {
            put(it->first, it->second);
        }
    }
    void set(const std::string& name, const std::string& value) {
        m_saved.emplace_back(name, get(name));
        put(name, value);
    }

private:
    static std::optional<std::string> get(const std::string& name) {
#ifdef _WIN32
        char* value = nullptr;
        std::size_t size = 0;
        if (_dupenv_s(&value, &size, name.c_str()) != 0 || value == nullptr) {
            return std::nullopt;
        }
        std::string out(value);
        std::free(value);
        return out;
#else
        const char* value = std::getenv(name.c_str());
        return value ? std::optional<std::string>(value) : std::nullopt;
#endif
    }
    static void put(const std::string& name, const std::optional<std::string>& value) {
#ifdef _WIN32
        _putenv_s(name.c_str(), value ? value->c_str() : "");
#else
        if (value) {
            setenv(name.c_str(), value->c_str(), 1);
        } else {
            unsetenv(name.c_str());
        }
#endif
    }
    std::vector<std::pair<std::string, std::optional<std::string>>> m_saved;
};

class Normaliser {
public:
    // the work folder's real path and the path as given (on Windows a temporary folder can be given
    // with 8.3 names, which the current directory keeps), longest first
    explicit Normaliser(const fs::path& work) {
        for (const fs::path& p : {fs::canonical(work), fs::absolute(work).lexically_normal()}) {
            const std::string native = core::pathToUtf8(p);
            m_work.push_back(native);
#ifdef _WIN32
            m_work.push_back(replaceAll(native, "\\", "/"));  // the cache's paths are in generic form
#endif
        }
        std::sort(m_work.begin(), m_work.end(),
                  [](const std::string& a, const std::string& b) { return a.size() > b.size(); });
    }

    // the work folder, time stamps, cache keys and the cache manifest's name
    std::string operator()(std::string s) const {
        for (const std::string& w : m_work) {
            s = replaceAll(std::move(s), w, "{work}");
        }
        s = std::regex_replace(s, m_stamp, "{stamp}");
        s = std::regex_replace(s, m_cleanedFolder, "$1{key}");
        s = std::regex_replace(s, m_generated, "# Generated by MDBridge {version} on {time}");
        s = std::regex_replace(s, m_key, "-{key}");
        return std::regex_replace(s, m_manifest, "{manifest}");
    }
    std::string file(std::string s) const {
#ifdef _WIN32
        s = replaceAll(std::move(s), "\r\n", "\n");
#endif
        return std::regex_replace((*this)(std::move(s)), m_saved, "$1saved = {time}");
    }
    // a file's relative path as compared: also <name>-{key}.mdbridge.npz -> <name>-{key}.npz
    std::string name(std::string rel) const {
        return std::regex_replace((*this)(std::move(rel)), m_entry, "-{key}.$1");
    }
    // a path stored in the cache: '\' -> '/' after {work} on Windows (the project's paths are native there)
    std::string cachePath(std::string s) const {
        s = (*this)(std::move(s));
#ifdef _WIN32
        if (s.starts_with("{work}")) {
            s = replaceAll(std::move(s), "\\", "/");
        }
#endif
        return s;
    }

private:
    std::vector<std::string> m_work;
    std::regex m_stamp{R"(\d{4}-\d{2}-\d{2}_\d{6})"};
    std::regex m_key{R"(-[0-9a-f]{16}(?=\)|(\.mdbridge)?\.(npz|json)\b))"};
    // the cleaning (cli_parity.py CLEANED_KEY_RE, GENERATED_RE; the marker is compared as JSON);
    // the script names MDBridge (CMD-08): MDATK's "MDATK 0.0.4.1" is not accepted
    std::regex m_cleanedFolder{R"(([/\\]cleaned[/\\])[0-9a-f]{16}\b)"};
    std::regex m_generated{R"(# Generated by MDBridge \S+ on \d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2})"};
    std::regex m_manifest{R"(\b(MDATK|MDBridge)_cache\.json\b)"};
    std::regex m_entry{R"(-\{key\}\.mdbridge\.(npz|json)$)"};
    std::regex m_saved{R"((^|\n)saved = \d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(\.\d+)?(?=\n|$))"};
};

bool validUtf8(const std::string& s) {
    return core::encodeUtf8(core::decodeUtf8(s)) == s;
}


// ---- the analysis cache, compared by meaning (cli_parity.py masked_cache_json, npz_differences) --

// A cache manifest or entry description with what legitimately differs masked.
Value maskedCacheJson(const Value& v, const Normaliser& normalise, bool top = true) {
    if (const core::Table* table = v.get<core::Table>()) {
        core::Table out;
        for (const auto& [k, x] : *table) {
            if (top && (k == "mdatk_version" || k == "mdbridge_version")) {
                out.emplace_back("version", "{version}");
            } else if (top && k == "project") {
                out.emplace_back(k, "{project}");
            } else if (top && (k == "updated" || k == "created")) {
                out.emplace_back(k, "{time}");
            } else if ((k == "key" || k == "cleaning_key") && x.isString()) {
                out.emplace_back(k, "{key}");
            } else if (k == "mtime_ns") {
                out.emplace_back(k, 0);
            } else if (k == "file" && x.isString()) {
                out.emplace_back(k, normalise.name(x.stringOr()));
            } else {
                out.emplace_back(k, maskedCacheJson(x, normalise, false));
            }
        }
        return out;
    }
    if (const core::Array* array = v.get<core::Array>()) {
        core::Array out;
        for (const Value& x : *array) {
            out.push_back(maskedCacheJson(x, normalise, false));
        }
        return out;
    }
    if (v.isString()) {
        return normalise.cachePath(v.stringOr());
    }
    return v;
}

// all_entries in one order (the generator sorts it its way)
Value sortedEntries(Value v) {
    if (core::Table* table = v.get<core::Table>()) {
        for (auto& [k, x] : *table) {
            if (core::Array* entries = x.get<core::Array>(); entries && k == "all_entries") {
                std::sort(entries->begin(), entries->end(),
                          [](const Value& a, const Value& b) { return core::toJson(a, 0) < core::toJson(b, 0); });
            }
        }
    }
    return v;
}

// numpy's dtype.str
std::string dtypeText(const core::NumArray& a) {
    switch (a.dtype()) {
        case core::NumArray::DType::Bool: return "|b1";
        case core::NumArray::DType::Int8: return "|i1";
        case core::NumArray::DType::Int32: return "<i4";
        case core::NumArray::DType::Int64: return "<i8";
        case core::NumArray::DType::Float32: return "<f4";
        case core::NumArray::DType::Float64: return "<f8";
        case core::NumArray::DType::String: return "<U" + std::to_string(a.stringWidth());
    }
    return "?";
}

// The flat values as the generator writes them (floats as numbers, non-finite ones as text).
std::vector<Value> flatValues(const core::NumArray& a) {
    std::vector<Value> out;
    const auto real = [&out](double x) {
        out.emplace_back(std::isfinite(x) ? Value(x) : Value(std::isnan(x) ? "nan" : (x > 0 ? "inf" : "-inf")));
    };
    switch (a.dtype()) {
        case core::NumArray::DType::Bool:
            for (const std::uint8_t x : a.values<std::uint8_t>()) {
                out.emplace_back(x != 0);
            }
            break;
        case core::NumArray::DType::Int8:
            for (const std::int8_t x : a.values<std::int8_t>()) {
                out.emplace_back(std::int64_t{x});
            }
            break;
        case core::NumArray::DType::Int32:
            for (const std::int32_t x : a.values<std::int32_t>()) {
                out.emplace_back(std::int64_t{x});
            }
            break;
        case core::NumArray::DType::Int64:
            for (const std::int64_t x : a.values<std::int64_t>()) {
                out.emplace_back(x);
            }
            break;
        case core::NumArray::DType::Float32:
            for (const float x : a.values<float>()) {
                real(x);
            }
            break;
        case core::NumArray::DType::Float64:
            for (const double x : a.values<double>()) {
                real(x);
            }
            break;
        case core::NumArray::DType::String:
            for (const std::string& x : a.values<std::string>()) {
                out.emplace_back(x);
            }
            break;
    }
    return out;
}

// Floats within the tolerance of cli_parity.py (FLOAT_TOLERANCE), everything else exactly.
bool close(const Value& e, const Value& a) {
    const double* x = e.get<double>();
    const double* y = a.get<double>();
    if (x == nullptr || y == nullptr) {
        return e == a;
    }
    return std::abs(*x - *y) <= 1e-9 * std::max({1.0, std::abs(*x), std::abs(*y)});
}

std::int64_t integer(const Value& v) {
    const std::int64_t* i = v.get<std::int64_t>();
    return i ? *i : -1;
}

// The differences between MDATK's arrays (expected) and an .npz MDBridge wrote; [] when they agree.
std::vector<std::string> npzDifferences(const Value& expected, const fs::path& path) {
    const core::NamedArrays got = analysis::readNpz(path);
    const core::Array& e = *expected.get<core::Array>();
    std::vector<std::string> problems;
    std::vector<std::string> eNames;
    std::vector<std::string> gNames;
    for (const Value& entry : e) {
        eNames.push_back(entry.get<core::Array>()->at(0).stringOr());
    }
    for (const auto& [name, array] : got) {
        gNames.push_back(name);
    }
    if (eNames != gNames) {
        return {"other arrays"};
    }
    std::map<std::string, std::vector<Value>> eValues;
    std::map<std::string, std::vector<Value>> gValues;
    for (std::size_t n = 0; n < e.size(); ++n) {
        const core::Array& entry = *e[n].get<core::Array>();
        eValues[eNames[n]] = *core::parseJson(entry.at(3).stringOr(), "values").get<core::Array>();
        gValues[eNames[n]] = flatValues(got[n].second);
    }
    const auto has = [&eValues](const std::string& name) { return eValues.count(name) > 0; };
    std::vector<std::string> pairmaps;
    for (const std::string& name : eNames) {
        if (name.ends_with("_offsets")) {
            const std::string p = name.substr(0, name.size() - 8);
            if (has(p + "_i") && has(p + "_j") && has(p + "_v")) {
                pairmaps.push_back(p);
            }
        }
    }
    const auto inPairmap = [&pairmaps](const std::string& name) {
        return std::any_of(pairmaps.begin(), pairmaps.end(), [&name](const std::string& p) {
            return name == p + "_i" || name == p + "_j" || name == p + "_v";
        });
    };
    for (std::size_t n = 0; n < e.size(); ++n) {
        const core::Array& entry = *e[n].get<core::Array>();
        const std::string& name = eNames[n];
        core::Array shape;
        for (const std::size_t s : got[n].second.shape()) {
            shape.emplace_back(static_cast<std::int64_t>(s));
        }
        if (entry.at(1).stringOr() != dtypeText(got[n].second) || entry.at(2) != Value(shape)) {
            problems.push_back(name + ": dtype or shape " + dtypeText(got[n].second));
            continue;
        }
        if (inPairmap(name)) {
            continue;  // compared with the offsets
        }
        const std::vector<Value>& ev = eValues[name];
        const std::vector<Value>& gv = gValues[name];
        if (name.ends_with("_offsets") &&
            std::find(pairmaps.begin(), pairmaps.end(), name.substr(0, name.size() - 8)) != pairmaps.end()) {
            if (ev != gv) {
                problems.push_back(name + " differs");
                continue;
            }
            const std::string p = name.substr(0, name.size() - 8);
            using Entry = std::tuple<std::int64_t, std::int64_t, double>;
            using Values = std::map<std::string, std::vector<Value>>;
            const auto rows = [&p](Values& values, std::int64_t lo, std::int64_t hi) {
                std::vector<Entry> out;
                for (std::int64_t k = lo; k < hi; ++k) {
                    const auto at = [&](const char* s) -> const Value& {
                        return values[p + s].at(static_cast<std::size_t>(k));
                    };
                    const double* v = at("_v").get<double>();
                    out.emplace_back(integer(at("_i")), integer(at("_j")), v ? *v : std::nan(""));
                }
                std::sort(out.begin(), out.end());
                return out;
            };
            for (std::size_t f = 0; f + 1 < ev.size(); ++f) {
                const std::vector<Entry> a = rows(eValues, integer(ev[f]), integer(ev[f + 1]));
                const std::vector<Entry> b = rows(gValues, integer(ev[f]), integer(ev[f + 1]));
                for (std::size_t k = 0; k < a.size(); ++k) {
                    if (std::get<0>(a[k]) != std::get<0>(b[k]) || std::get<1>(a[k]) != std::get<1>(b[k]) ||
                        !close(Value(std::get<2>(a[k])), Value(std::get<2>(b[k])))) {
                        problems.push_back("pair map " + p + ", frame " + std::to_string(f) + " differs");
                        break;
                    }
                }
            }
            continue;
        }
        for (std::size_t k = 0; k < ev.size(); ++k) {
            if (!close(ev[k], gv[k])) {
                problems.push_back(name + "[" + std::to_string(k) + "]: " + core::toJson(ev[k], 0) + " / " +
                                   core::toJson(gv[k], 0));
                break;
            }
        }
    }
    return problems;
}

// cli_parity.py texts_close (FITTED): equal but for their decimal numbers, which differ by at most
// `tolerance`, and for white space between them.
bool textsClose(const std::string& expected, const std::string& actual, double tolerance) {
    static const std::regex number{R"(-?\d+(?:\.\d+(?:e[-+]?\d+)?|e[-+]?\d+))"};  // 0.25, 1.9e-05, 7e-06
    struct Split {
        std::vector<std::string> parts;
        std::vector<double> numbers;
    };
    const auto split = [](const std::string& text) {
        Split out;
        std::size_t last = 0;
        for (auto it = std::sregex_iterator(text.begin(), text.end(), number); it != std::sregex_iterator(); ++it) {
            const std::size_t at = static_cast<std::size_t>(it->position());
            out.parts.push_back(text.substr(last, at - last));
            out.numbers.push_back(std::stod(it->str()));
            last = at + static_cast<std::size_t>(it->length());
        }
        out.parts.push_back(text.substr(last));
        return out;
    };
    Split a = split(expected);
    Split b = split(actual);
    if (a.parts.size() != b.parts.size()) {
        return false;
    }
    // white space between the numbers is not compared, nor a trailing comma before "]" (cli_parity.py
    // texts_close: PDB columns, a TOML list on one line or one number per line)
    const auto canon = [](const std::string& t) {
        std::string out;
        for (const char ch : t) {
            if (ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r') continue;
            if (ch == ']' && !out.empty() && out.back() == ',') out.pop_back();
            out += ch;
        }
        return out;
    };
    for (std::size_t i = 0; i < a.parts.size(); ++i) {
        a.parts[i] = canon(a.parts[i]);
        b.parts[i] = canon(b.parts[i]);
    }
    if (a.parts != b.parts) {
        return false;
    }
    for (std::size_t i = 0; i < a.numbers.size(); ++i) {
        if (!(std::abs(a.numbers[i] - b.numbers[i]) <= tolerance)) {
            return false;
        }
    }
    return true;
}

// cli_parity.py step_log_mask: the step log's time and duration (the chem cases, "needs")
std::string stepLogMasked(const std::string& text) {
    static const std::regex time{R"((^|\n)time = \d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(\.\d+)?(?=\n|$))"};
    static const std::regex duration{R"((^|\n)duration_s = [0-9.e+-]+(?=\n|$))"};
    return std::regex_replace(std::regex_replace(text, time, "$1time = {time}"), duration, "$1duration_s = 0.0");
}

// How a case's texts are compared: the step log masked (chem cases), numbers within a tolerance
// (fitted cases).
struct TextRules {
    bool logMask = false;
    std::optional<double> tolerance;
};

// cli_parity.py map_summary: an MRC file's header as text
std::string mapSummary(const std::string& data) {
    const auto i32 = [&data](std::size_t word) {
        std::int32_t v = 0;
        std::memcpy(&v, data.data() + 4 * word, 4);
        return static_cast<int>(v);
    };
    const auto f32 = [&data](std::size_t word) {
        float v = 0;
        std::memcpy(&v, data.data() + 4 * word, 4);
        return static_cast<double>(v);
    };
    char text[512];
    std::snprintf(text, sizeof text,
                  "grid %d %d %d mode %d start %d %d %d sampling %d %d %d cell %.4f %.4f %.4f axes %d %d %d "
                  "min %.6f max %.6f mean %.6f origin %.4f %.4f %.4f",
                  i32(0), i32(1), i32(2), i32(3), i32(4), i32(5), i32(6), i32(7), i32(8), i32(9), f32(10), f32(11),
                  f32(12), i32(16), i32(17), i32(18), f32(19), f32(20), f32(21), f32(49), f32(50), f32(51));
    return text;
}

// The differences between one expected file entry (text, sha256, json, npz, map) and a file; [] when equal.
std::vector<std::string> fileDifferences(const Value& expected, const fs::path& path, const Normaliser& normalise,
                                         const TextRules& rules = {}) {
    if (const Value* arrays = expected.find("npz")) {
        return npzDifferences(*arrays, path);
    }
    const std::string content = test::readFile(path);
    if (const Value* map = expected.find("map")) {  // a fitted case's map: its header within the tolerance
        if (content.size() < 1024) {
            return {"is no MRC map"};
        }
        const std::string got = mapSummary(content);
        const bool same = got == map->stringOr() || (rules.tolerance && textsClose(map->stringOr(), got, *rules.tolerance));
        return same ? std::vector<std::string>{} : std::vector<std::string>{"differs: " + got};
    }
    if (const Value* json = expected.find("json")) {
        const Value got = sortedEntries(maskedCacheJson(core::parseJson(content, core::pathToUtf8(path)), normalise));
        if (got != sortedEntries(*json)) {
            return {"differs:\n" + core::toJson(got, 1)};
        }
        return {};
    }
    if (validUtf8(content)) {
        const std::string got = rules.logMask ? stepLogMasked(normalise.file(content)) : normalise.file(content);
        if (const Value* t = expected.find("text")) {
            const bool same = got == t->stringOr() || (rules.tolerance && textsClose(t->stringOr(), got, *rules.tolerance));
            return same ? std::vector<std::string>{} : std::vector<std::string>{"differs:\n" + got};
        }
        if (core::sha256Hex(got) != text(expected, "sha256") ||
            static_cast<std::int64_t>(got.size()) != *field(expected, "size").get<std::int64_t>()) {
            return {"differs (SHA-256 / size)"};
        }
        return {};
    }
    if (core::sha256Hex(content) != text(expected, "sha256")) {
        return {"differs (SHA-256)"};
    }
    return {};
}

// The render, charts and preview cases ("render": true) that draw a figure or a video (their
// files compared by metadata, tests/parity/cli_parity.py render_file_entry) need the program's Qt
// backend: cli_parity.py runs them against the `mdbridge` program (the CI parity job; locally
// against a Qt build). This test runs the others in-process: those stopped by check() or by the
// format checks (RN-11 ... RN-15) and those that draw nothing.
bool drawsFigures(const Value& c) {
    const Value* render = c.find("render");
    if (render == nullptr || render->get<bool>() == nullptr || !*render->get<bool>()) {
        return false;
    }
    for (const auto& [rel, entry] : *field(c, "written").get<core::Table>()) {
        for (const char* kind : {"png", "jpeg", "tiff", "pdf", "ps", "svg", "mp4", "params_log"}) {
            if (entry.find(kind) != nullptr) {
                return true;
            }
        }
    }
    return false;
}

TEST(Commands, MatchMdatk) {
    const Value data = test::readJson(kData / "cli" / "cli_cases.json");
    std::size_t run = 0;
    for (const Value& c : *field(data, "cases").get<core::Array>()) {
        const std::string name = text(c, "name");
        SCOPED_TRACE(name);
#ifdef _WIN32
        // absolute paths, OS error texts and environment variables as on Linux, where parity is
        // defined (posix_only)
        if (*field(c, "posix_only").get<bool>()) {
            continue;
        }
#endif
        if (drawsFigures(c)) {
            continue;  // render / charts / preview cases that draw: see drawsFigures
        }
        if (const Value* pending = c.find("pending"); pending != nullptr && pending->get<bool>() != nullptr &&
                                                      *pending->get<bool>()) {
            continue;  // a command MDBridge does not have yet (cli_parity.py PENDING)
        }
        test::TempDir work;
        const std::string workText = core::pathToUtf8(fs::canonical(work.path()));
        for (const auto& [dest, source] : *field(c, "inputs").get<core::Table>()) {
            const fs::path target = work.path() / core::pathFromUtf8(dest);
            fs::create_directories(target.parent_path());
            fs::copy_file(kData / core::pathFromUtf8(source.stringOr()), target);
        }
        for (const auto& [dest, content] : *field(c, "setup").get<core::Table>()) {
            const fs::path target = work.path() / core::pathFromUtf8(dest);
            fs::create_directories(target.parent_path());
            test::writeFile(target, replaceAll(content.stringOr(), "{work}", workText));
        }
        for (const Value& dest : *field(c, "executable").get<core::Array>()) {
            fs::permissions(work.path() / core::pathFromUtf8(dest.stringOr()),
                            fs::perms::owner_exec | fs::perms::group_exec | fs::perms::others_exec,
                            fs::perm_options::add);
        }
        ScopedEnvironment environment;
        for (const auto& [variable, value] : *field(c, "env").get<core::Table>()) {
            environment.set(variable, replaceAll(replaceAll(value.stringOr(), "{work}", workText), "{pubchem}",
                                                 kPubchemStandIn));
        }
        TextRules rules;
        rules.logMask = c.find("needs") != nullptr;
        if (const Value* fitted = c.find("fitted")) {
            rules.tolerance = fitted->get<double>() ? *fitted->get<double>()
                                                    : static_cast<double>(*fitted->get<std::int64_t>());
        }
        const auto sameText = [&rules](const std::string& got, const std::string& expected) {
            return got == expected || (rules.tolerance && textsClose(expected, got, *rules.tolerance));
        };
        const auto arguments = [](const Value& list) {
            std::vector<std::string> args;
            for (const Value& a : *list.get<core::Array>()) {
                args.push_back(a.stringOr());
            }
            return args;
        };
        std::size_t index = 0;
        for (const Value& before : *field(c, "before").get<core::Array>()) {
            const CommandRun first = runIn(work.path(), arguments(before));
            EXPECT_EQ(first.exit, integer(field(c, "before_exits").get<core::Array>()->at(index++))) << first.err;
        }
        const CommandRun result = runIn(work.path(), arguments(field(c, "argv")));
        ++run;
        const Normaliser normalise(work.path());
        if (field(c, "exit_only").get<bool>() != nullptr && *field(c, "exit_only").get<bool>()) {
            EXPECT_EQ(result.exit, *field(c, "exit").get<std::int64_t>());
            continue;
        }
        EXPECT_EQ(result.exit, *field(c, "exit").get<std::int64_t>()) << result.err;
        if (!sameText(normalise(result.out), text(c, "stdout"))) {
            EXPECT_EQ(normalise(result.out), text(c, "stdout"));
        }
        if (field(c, "stderr").isNull()) {
            // MDATK crashed with a traceback (mdatk_crash): MDBridge reports the error instead
            EXPECT_TRUE(result.err.starts_with("error: ")) << text(c, "mdatk_crash") << " / " << result.err;
        } else if (!sameText(normalise(result.err), text(c, "stderr"))) {
            EXPECT_EQ(normalise(result.err), text(c, "stderr"));
        }

        // the files left in the folder: everything but the unchanged inputs; several cache entries
        // of one analysis (different keys) are matched in any order
        core::Table expectedFiles = *field(c, "written").get<core::Table>();
#ifdef _WIN32
        // MDATK's project lock file (.NAME.lock) needs fcntl, which Windows lacks: there neither MDATK
        // nor MDBridge creates one (projectfile._locked, detail::ProjectLock)
        static const std::regex projectLock{R"((^|/)\.[^/]+\.lock$)"};
        std::erase_if(expectedFiles,
                      [](const auto& entry) { return std::regex_search(entry.first, projectLock); });
#endif
        std::map<std::string, std::vector<fs::path>> seen;
        for (const auto& entry : fs::recursive_directory_iterator(work.path())) {
            if (!entry.is_regular_file()) {
                continue;
            }
            std::string rel = core::pathToUtf8(entry.path().lexically_relative(work.path()));
            rel = replaceAll(rel, "\\", "/");
            if (field(c, "inputs").find(rel) != nullptr) {
                continue;
            }
            if (const Value* given = field(c, "setup").find(rel);
                given != nullptr && replaceAll(given->stringOr(), "{work}", workText) == test::readFile(entry.path())) {
                continue;
            }
            seen[normalise.name(rel)].push_back(entry.path());
        }
        for (const auto& [rel, paths] : seen) {
            const Value* expected = nullptr;
            for (const auto& [key, value] : expectedFiles) {
                if (key == rel) {
                    expected = &value;
                }
            }
            if (expected == nullptr) {
                ADD_FAILURE() << "unexpected file " << rel;
                continue;
            }
            std::vector<Value> entries{*expected};
            if (const Value* several = expected->find("several")) {
                entries = *several->get<core::Array>();
            }
            if (entries.size() != paths.size()) {
                ADD_FAILURE() << rel << ": " << paths.size() << " files, MDATK wrote " << entries.size();
                continue;
            }
            std::vector<fs::path> unused = paths;
            for (const Value& e : entries) {
                std::vector<std::string> problems;
                const auto match = std::find_if(unused.begin(), unused.end(), [&](const fs::path& p) {
                    problems = fileDifferences(e, p, normalise, rules);
                    return problems.empty();
                });
                if (match == unused.end()) {
                    ADD_FAILURE() << rel << ": " << (problems.empty() ? std::string("no match") : problems.front());
                } else {
                    unused.erase(match);
                }
            }
        }
        for (const auto& [key, value] : expectedFiles) {
            if (seen.count(key) == 0) {
                ADD_FAILURE() << "missing file " << key;
            }
        }
    }
    EXPECT_GT(run, 0u);
}

// ---- pieces ------------------------------------------------------------------------------------

TEST(Commands, OverridesFollowMdatk) {
    core::Project p;
    app::applyOverride(p, "analyses.contacts.cutoff=5");
    EXPECT_EQ(p.analyses.contacts.cutoff, 5.0);
    app::applyOverride(p, " groups.-1.name =Last");  // the key is stripped, indices may be negative
    EXPECT_EQ(p.groups.back().name, "Last");
    app::applyOverride(p, "notes=[1, 2]");  // a text setting takes the text
    EXPECT_EQ(p.notes, "[1, 2]");
    app::applyOverride(p, "notes=\"quoted\"");
    EXPECT_EQ(p.notes, "quoted");
    app::applyOverride(p, "inputs.stop=7");
    EXPECT_EQ(p.inputs.stop, 7);
    app::applyOverride(p, "inputs.stop=null");
    EXPECT_EQ(p.inputs.stop, std::nullopt);
    app::applyOverride(p, "groups=[{\"name\": \"X\"}]");
    ASSERT_EQ(p.groups.size(), 1u);
    EXPECT_EQ(p.groups[0].residues, "");

    const auto message = [&p](const std::string& assignment) {
        try {
            app::applyOverride(p, assignment);
        } catch (const core::ValueError& e) {
            return std::string(e.what());
        }
        return std::string("no error");
    };
    EXPECT_EQ(message("x"), "--set expects key=value, got 'x'");
    EXPECT_EQ(message(" . =1"), "empty key in ' . =1'");
    EXPECT_EQ(message("inputs.bogus=1"), "unknown parameter 'inputs.bogus'");
    EXPECT_EQ(message("layout.panels.3.options.mode=x"), "unknown parameter 'layout.panels.3.options.mode'");
    EXPECT_EQ(message("validate.x=1"), "unknown parameter 'validate.x'");
    EXPECT_EQ(message("cleaning.enabled=1"), "cleaning.enabled expects true/false");
    EXPECT_EQ(message("groups.1.name=x"), "groups.1.name: no list element '1'");
    EXPECT_EQ(message("groups.99999999999999999999999.name=x"),
              "groups.99999999999999999999999.name: no list element '99999999999999999999999'");
    // MDATK: "invalid literal for int() with base 10: '--1'" (PRJ-04)
    EXPECT_EQ(message("groups.--1.name=x"), "groups.--1.name: no list element '--1'");
    EXPECT_EQ(message("groups.-.name=x"), "groups.-.name: no list element '-'");
    EXPECT_EQ(message("analyses.contacts.cutoff=far"),
              "project.analyses.contacts.cutoff: expected a number, got 'far'");
}

TEST(Commands, JsonIsAscii) {
    EXPECT_EQ(app::jsonForConsole(Value(core::Table{{"k\xc3\xa9", "\xf0\x9f\x98\x80\x7f"}})),
              "{\n  \"k\\u00e9\": \"\\ud83d\\ude00\\u007f\"\n}");
    EXPECT_EQ(app::jsonForConsole(Value(std::string("a\xff"))), "\"a\\udcff\"");
}

TEST(Commands, OnlyTheseCommandsAreOffered) {
    CLI::App root("test", "mdbridge");
    app::CoreCommands commands(root);
    std::vector<std::string> names;
    for (const CLI::App* sub : root.get_subcommands({})) {
        names.push_back(sub->get_name());
    }
    EXPECT_EQ(names, (std::vector<std::string>{"demo", "init", "info", "list", "clean", "analyze", "export-pdb",
                                               "render", "charts", "preview", "citations", "density", "pdb",
                                               "build", "dict"}));
    const std::string help = root.help();
    EXPECT_EQ(help.find("workflow"), std::string::npos);
}

TEST(Commands, UsageErrorsExitWith2) {
    test::TempDir work;
    EXPECT_EQ(runIn(work.path(), {"init"}).exit, 2);
    EXPECT_EQ(runIn(work.path(), {"density"}).exit, 2);
    EXPECT_EQ(runIn(work.path(), {"pdb"}).exit, 2);
    EXPECT_EQ(runIn(work.path(), {"build"}).exit, 2);
    EXPECT_EQ(runIn(work.path(), {"dict"}).exit, 2);
    EXPECT_EQ(runIn(work.path(), {"density", "fit", "--resolution", "x"}).exit, 2);
    EXPECT_EQ(runIn(work.path(), {"demo"}).exit, 2);
    EXPECT_EQ(runIn(work.path(), {"demo", "d", "--n-frames", "x"}).exit, 2);
    EXPECT_EQ(runIn(work.path(), {"clean", "--force=x"}).exit, 2);
    EXPECT_EQ(runIn(work.path(), {"citations", "--bogus"}).exit, 2);
    EXPECT_EQ(runIn(work.path(), {"preview", "--frame", "x"}).exit, 2);
    EXPECT_EQ(runIn(work.path(), {"list", "x"}).exit, 2);
    EXPECT_EQ(runIn(work.path(), {"analyze", "--force=x"}).exit, 2);
    EXPECT_EQ(runIn(work.path(), {}).exit, 2);
    EXPECT_EQ(runIn(work.path(), {"--help"}).exit, 0);
}

TEST(Commands, InitWritesTheProjectFileWithTheClock) {
    test::TempDir work;
    const core::Clock previous = core::setClock({[] {
                                                     core::DateTime t;
                                                     t.year = 2026, t.month = 10, t.day = 15, t.hour = 14, t.minute = 32;
                                                     return t;
                                                 },
                                                 {}});
    const CommandRun run = runIn(work.path(), {"init-params", "a.mdatk.toml", "-s", "notes=hello"});
    core::setClock(previous);
    EXPECT_EQ(run.exit, 0);
    EXPECT_EQ(run.out, "Wrote a project file with default settings to a.mdatk.toml\n");
    const std::string written = test::readFile(work.path() / "a.mdatk.toml");
    EXPECT_NE(written.find("saved = 2026-10-15T14:32:00"), std::string::npos);
    EXPECT_NE(written.find("notes = \"hello\""), std::string::npos);
    EXPECT_EQ(runIn(work.path(), {"init", "a.mdatk.toml"}).err, "a.mdatk.toml exists (use --overwrite)\n");
}

TEST(Commands, ListPrintsTheTables) {
    std::ostringstream out;
    app::printList(out);
    const std::string printed = out.str();
    EXPECT_TRUE(printed.starts_with("Analyses:\n  contacts   Pairwise contacts\n"));
    EXPECT_NE(printed.find("\n  stacking   \xcf\x80\xe2\x80\x93\xcf\x80 stacking\n"), std::string::npos);  // π–π
    EXPECT_NE(printed.find("\n  trace      Backbone trace (internal, used by the structure panel)\n"),
              std::string::npos);
    EXPECT_NE(printed.find("\nConda environments with cpptraj:\n"), std::string::npos);
    for (const core::PanelType& t : core::panelTypes()) {
        EXPECT_NE(printed.find(std::string(t.label) + " - " + std::string(t.description) + "\n"), std::string::npos)
            << t.name;
    }
    // names longer than the column push the label one space further (format(name, "14s"))
    EXPECT_NE(printed.find("\n  large_top_3_strip 1 large (top) + 3 small + strip\n"), std::string::npos);
    EXPECT_NE(printed.find("\n  single         1 large panel\n"), std::string::npos);
    EXPECT_TRUE(printed.ends_with("Run 'mdbridge analyze' to list the keys of a system.\n"));
}

TEST(Commands, ErrorsAreReportedAsMdatkReportsThem) {
    // MDATK's main() prints str(exc) of OSError, ValueError, RuntimeError and KeyError (the last with
    // quotes; MDBridge prints the message: CLI-01)
    EXPECT_EQ(app::pythonErrorText(analysis::PipelineError("Invalid parameters:\n  - x")),
              "Invalid parameters:\n  - x");
    EXPECT_EQ(app::pythonErrorText(core::NamedError("FileNotFoundError", "no MDBridge_cache.json in /x")),
              "no MDBridge_cache.json in /x");
    EXPECT_EQ(app::pythonErrorText(core::NamedError("CpptrajError", "topology not found: a.parm7")),
              "topology not found: a.parm7");
    EXPECT_EQ(app::pythonErrorText(std::out_of_range("unknown analysis 'x'")), "unknown analysis 'x'");

    test::TempDir work;  // cleaning with a cpptraj that does not exist: CpptrajError, exit 1
    fs::copy_file(kData / "topology" / "demo_dry.parm7", work.path() / "demo_dry.parm7");
    fs::copy_file(kData / "io" / "trajectories" / "demo_dry.nc", work.path() / "demo_dry.nc");
    const CommandRun run = runIn(work.path(), {"analyze", "--top", "demo_dry.parm7", "--traj", "demo_dry.nc",
                                               "--clean", "--cpptraj", "/nonexistent/cpptraj", "-q"});
    EXPECT_EQ(run.exit, 1);
    EXPECT_TRUE(run.err.starts_with("error: cpptraj")) << run.err;
    EXPECT_EQ(run.out, "");
}

TEST(Commands, AnalyzeRedrawsProgressOnATerminal) {
    test::TempDir work;
    fs::copy_file(kData / "topology" / "demo_dry.parm7", work.path() / "demo_dry.parm7");
    fs::copy_file(kData / "io" / "trajectories" / "demo_dry.nc", work.path() / "demo_dry.nc");
    const std::vector<std::string> args{"analyze", "--top", "demo_dry.parm7", "--traj",
                                        "demo_dry.nc", "--only", "rmsd"};
    const CommandRun first = runIn(work.path(), args, true);
    EXPECT_EQ(first.exit, 0) << first.err;
    // log lines and progress lines start by clearing the line; only progress at 100 % ends it
    EXPECT_TRUE(first.err.starts_with(
        "\r\033[K[  2.0%] inputs ready\r\033[KPreparing RMSD ...\n\r\033[KAnalysing 24 frame(s): RMSD\n"))
        << first.err;
    EXPECT_NE(first.err.find("\r\033[K[  "), std::string::npos);  // an intermediate fraction
    EXPECT_TRUE(first.err.ends_with("\r\033[K[100.0%] analysis done\n")) << first.err;
    const CommandRun second = runIn(work.path(), args, true);
    const std::regex loaded("\r\033\\[KRMSD: loaded from cache \\(rmsd-[0-9a-f]{16}\\)\n");
    EXPECT_TRUE(std::regex_search(second.err, loaded)) << second.err;
    const CommandRun quiet = runIn(work.path(), {"analyze", "--top", "demo_dry.parm7", "--traj", "demo_dry.nc",
                                                 "--only", "rmsd", "--force", "-q"},
                                   true);
    EXPECT_EQ(quiet.err, "");
    EXPECT_NE(quiet.out.find("\"n_frames\": 24"), std::string::npos);
}

}  // namespace
}  // namespace mdbridge
