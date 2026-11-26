#pragma once
// `dict pubchem`'s PubChem lookup (MDATK's dictionaries.fetch_pubchem_smiles, owner decision D7):
// PUG REST over HTTP(S), only when the command asks for it. The request logic, the answers and the
// messages are here, without Qt; the transport is an HttpGet the program passes (the `mdbridge`
// program's Qt Network one, app/pubchem_qt.cpp; the tests' stand-ins). `chem` never touches the
// network.
//
// The root of the service is MDATK's PUBCHEM_REST unless the environment variable
// MDBRIDGE_PUBCHEM_URL names another (a local stand-in server in the tests; D7).

#include <functional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mdbridge::app {

// What a GET gave.
struct HttpReply {
    enum class Outcome {
        Answered,  // the server answered (any HTTP status)
        TimedOut,  // no answer within the timeout
        Failed,    // no connection (refused, unknown host, proxy, TLS ...)
    };
    Outcome outcome = Outcome::Failed;
    int status = 0;      // Answered: the HTTP status
    std::string reason;  // Answered: the reason phrase ("Internal Server Error")
    std::string body;    // Answered: the body's bytes
    std::string error;   // Failed: what went wrong, for the message
};

using HttpHeaders = std::vector<std::pair<std::string, std::string>>;
// GET `url` with `headers`, giving up after `timeoutSeconds` (> 0).
using HttpGet = std::function<HttpReply(const std::string& url, const HttpHeaders& headers, double timeoutSeconds)>;

// What the commands may use beyond the file system.
struct NetworkBackend {
    HttpGet get;  // empty: `dict pubchem` fails with "PubChem lookups need the mdbridge program ..."
};

inline constexpr std::string_view kPubchemRest = "https://pubchem.ncbi.nlm.nih.gov/rest/pug";
inline constexpr std::string_view kPubchemUrlVariable = "MDBRIDGE_PUBCHEM_URL";
// PUBCHEM_SMILES_PROPERTIES: the current property name first, then the one older servers know
inline constexpr std::string_view kPubchemSmilesProperties[] = {"SMILES", "IsomericSMILES"};

// The PUG REST root: MDBRIDGE_PUBCHEM_URL without trailing '/' when set and not empty, else
// kPubchemRest.
std::string pubchemRoot();

// urllib.parse.quote(text, safe=""): every byte but A-Z a-z 0-9 _ . - ~ as %XX.
std::string percentEncoded(std::string_view text);

// fetch_pubchem_smiles(cid_or_name, timeout): the SMILES of a PubChem compound (all digits: a CID,
// else a name), the property SMILES first and IsomericSMILES when the server rejects that name
// (HTTP 400). io::DictionaryError with MDATK's messages: "no PubChem CID or name given", "PubChem
// has no compound 'x'" (404), "PubChem request failed: HTTP 500 Internal Server Error", "PubChem did
// not answer within 0.5 s", "could not reach PubChem (<reason>); check the network connection /
// proxy", "PubChem returned an unreadable answer" (also for an answer whose properties are no list
// of objects, where MDATK crashed), "PubChem's answer for 'x' holds no SMILES", "PubChem rejected
// the request for IsomericSMILES (HTTP 400)"; "PubChem lookups need the mdbridge program (this one
// was built without network access)" without a transport.
std::string fetchPubchemSmiles(std::string_view cidOrName, double timeoutSeconds, const HttpGet& get,
                               const std::string& root = pubchemRoot());

}  // namespace mdbridge::app
