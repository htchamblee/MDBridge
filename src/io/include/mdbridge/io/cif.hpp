#pragma once
// mmCIF documents as MDATK reads and writes them (the mmCIF part of mdatk/dictionaries.py): the
// restraint dictionaries of acedrg / eLBOW and the CCD, parsed into blocks of loops and tag/value
// pairs and written back in MDATK's layout. Text is UTF-8 and handled as Python handles str
// (characters are code points; bytes that are not UTF-8 read as U+FFFD).

#include <cstddef>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "mdbridge/core/errors.hpp"

namespace mdbridge::io {

// mdatk.dictionaries.DictionaryError: an unreadable dictionary (and, in later milestones, a
// SMILES or structure that does not match its code, or a failed lookup).
class DictionaryError : public core::ValueError {
public:
    using core::ValueError::ValueError;
    const char* pythonName() const noexcept override { return "DictionaryError"; }
};

// A category: a loop_ or a run of "_tag value" pairs (pairs = true, one row).
struct CifLoop {
    std::vector<std::string> tags;
    std::vector<std::vector<std::string>> rows;
    bool pairs = false;

    // The lower-case category of the first tag ("_chem_comp_atom"), "" without tags.
    std::string category() const;
    // Index of the item ("atom_id" of "_chem_comp_atom.atom_id"; case-insensitive), -1 if absent.
    std::ptrdiff_t col(std::string_view item) const;
    // The item's values (empty if absent). Throws std::out_of_range for a row that is too short.
    std::vector<std::string> column(std::string_view item) const;
};

struct CifBlock {
    std::string name;
    std::vector<CifLoop> items;

    // The first item of the category (case-insensitive), or nullptr.
    const CifLoop* get(std::string_view category) const;
    CifLoop* get(std::string_view category);
};

struct CifDocument {
    std::vector<std::string> header;  // comment lines before the first block
    std::vector<CifBlock> blocks;

    // The first block of that name (case-insensitive), or nullptr.
    const CifBlock* block(std::string_view name) const;
    CifBlock* block(std::string_view name);
};

// The tokens of _cif_tokens: text fields (";" lines) and quoted or bare words as values,
// "#..." comments (rest of the line), data_<name> (text: the name), loop_, save_* / global_ /
// stop_ (Other, text: the word) and _tags.
struct CifToken {
    enum class Kind { Value, Comment, Data, Loop, Other, Tag };
    Kind kind = Kind::Value;
    std::string text;

    bool operator==(const CifToken&) const = default;
};

// All tokens of the text. Throws DictionaryError("mmCIF: unterminated text field") or
// ("mmCIF line N: unterminated quote").
std::vector<CifToken> cifTokens(std::string_view text);

// parse_cif: blocks, loops and tag/value pairs (comments before the first block become the
// header, later ones are dropped; anything before the first data_ block is dropped). Throws
// DictionaryError with MDATK's messages, at the first problem in the order MDATK meets them:
// "mmCIF: loop <category or ?> has N values for M columns", "mmCIF: <tag> has no value",
// "mmCIF: value <repr> without a tag" and the token errors.
CifDocument parseCif(std::string_view text);

// cif_value: a value as an mmCIF token ('' for empty text, double quotes for whitespace, a
// leading _#$'"[]; or a reserved word prefix, single quotes when it holds ", a text field
// "\n;...\n;\n" when it holds a line break or both quotes).
std::string cifValue(std::string_view value);

// write_cif: the document in MDATK's layout (byte for byte).
std::string writeCif(const CifDocument& doc);

// is_dictionary_file: the first 2,000,000 characters of the file (read as UTF-8 with
// replacement, universal newlines) contain "_chem_comp_atom." and "_chem_comp_bond."; false
// when it cannot be read.
bool isDictionaryFile(const std::string& path);

// find_dictionary: the restraint dictionary a program wrote, from files, folders (their *.cif)
// or glob patterns (Python's glob: hidden names only matched by patterns starting with "."),
// without the `exclude` paths (compared as absolute paths); files named after `code` first
// (case-insensitive), then paths without "tmp", the newest, the path. Throws
// DictionaryError("no restraint dictionary (.cif with _chem_comp_atom and _chem_comp_bond loops)
// in <sources joined with ', '>"). Paths are UTF-8, split and joined with POSIX rules (os.path
// on Linux, where MDATK is the reference); names match case-sensitively on every platform.
std::string findDictionary(std::span<const std::string> sources, std::string_view code = {},
                           std::span<const std::string> exclude = {});

}  // namespace mdbridge::io
