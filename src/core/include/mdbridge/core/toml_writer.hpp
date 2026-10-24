#pragma once
// Commented TOML writing, as MDATK's `tomlio` writes its project files: every entry carries the
// comment its documentation callback gives, and the text matches MDATK's byte for byte.

#include <functional>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "mdbridge/core/text.hpp"
#include "mdbridge/core/value.hpp"

namespace mdbridge::core {

// Comment of an entry: `path` holds the keys from the root (array indices left out), `stack` the
// enclosing tables, outermost (the document) first. std::nullopt or "" means no comment.
using TomlDocFunction =
    std::function<std::optional<std::string>(const std::vector<std::string>& path, const std::vector<const Value*>& stack)>;

// tomlio.key / quote / value: a bare key when possible, a basic string, a TOML value.
std::string tomlKey(std::string_view name);
std::string tomlQuote(std::string_view text);
std::string tomlValue(const Value& value);

// tomlio.Writer: tables, arrays of tables and scalars with their comments.
class TomlWriter {
public:
    explicit TomlWriter(TomlDocFunction doc = {});

    // "# ..." lines wrapped to 100 columns (nothing for an empty text).
    void comment(const std::optional<std::string>& text, std::string_view indent = "");
    // An empty line, unless the last line is empty already.
    void blank();
    // The whole document: header comment lines, then the data.
    std::string document(const Value& data, const std::vector<std::string>& header = {});
    // "[a.b]" or "[[a.b]]" and its entries; the comment is written for the first array element only.
    void table(const std::vector<std::string>& path, const Value& data, std::vector<const Value*> stack,
               bool array = false, bool first = true);

    const std::vector<std::string>& lines() const { return m_lines; }

private:
    void body(const std::vector<std::string>& path, const Value& data, const std::vector<const Value*>& stack);
    std::optional<std::string> doc(const std::vector<std::string>& path, const std::vector<const Value*>& stack) const;

    TomlDocFunction m_doc;
    std::vector<std::string> m_lines;
};

// tomlio.dumps
std::string dumpToml(const Value& data, TomlDocFunction doc = {}, const std::vector<std::string>& header = {});

}  // namespace mdbridge::core
