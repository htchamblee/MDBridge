// The determinism lint (lint.hpp): a tokenizer for Lua 5.4 and checks over its tokens.

#include "mdbridge/workflow/lua/lint.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <map>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mdbridge::workflow::lua {

namespace {

enum class Kind { Name, Number, String, Op };

struct Token {
    Kind kind = Kind::Op;
    std::string text;  // a String's content between its delimiters, escapes as written
    int line = 0;
    int column = 0;
    bool isFloat = false;  // Number
};

constexpr std::array kRules = {std::string_view("library"), std::string_view("pairs"), std::string_view("tostring"),
                               std::string_view("format"),  std::string_view("concat"), std::string_view("collector"),
                               std::string_view("waiver"),  std::string_view("syntax")};

constexpr std::array kMissing = {
    std::string_view("os"),       std::string_view("io"),       std::string_view("debug"),
    std::string_view("package"),  std::string_view("coroutine"), std::string_view("require"),
    std::string_view("dofile"),   std::string_view("loadfile"), std::string_view("print"),
    std::string_view("warn"),     std::string_view("collectgarbage"), std::string_view("rawset"),
    std::string_view("loadstring"), std::string_view("setfenv"), std::string_view("getfenv"),
    std::string_view("module"),   std::string_view("unpack")};

// math members whose results are floats
constexpr std::array kFloatMath = {std::string_view("pi"),  std::string_view("huge"), std::string_view("sqrt"),
                                   std::string_view("exp"), std::string_view("log"),  std::string_view("sin"),
                                   std::string_view("cos"), std::string_view("tan"),  std::string_view("asin"),
                                   std::string_view("acos"), std::string_view("atan")};

template <std::size_t N>
bool contains(const std::array<std::string_view, N>& list, std::string_view name) {
    return std::find(list.begin(), list.end(), name) != list.end();
}

bool isNameStart(char c) {
    return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || c == '_';
}

bool isDigit(char c) {
    return c >= '0' && c <= '9';
}

bool isNameChar(char c) {
    return isNameStart(c) || isDigit(c);
}

class Lexer {
public:
    explicit Lexer(std::string_view source) : m_source(source) {}

    // Tokens of the whole file; comments become waivers; an unfinished string or comment ends the
    // scan with a "syntax" finding.
    void run() {
        while (m_pos < m_source.size()) {
            const char c = m_source[m_pos];
            if (c == '\n' || c == '\r') {
                newline();
            } else if (c == ' ' || c == '\t' || c == '\f' || c == '\v') {
                ++m_pos;
            } else if (c == '-' && peek(1) == '-') {
                if (!comment()) {
                    return;
                }
            } else if (isNameStart(c)) {
                const std::size_t start = m_pos;
                while (m_pos < m_source.size() && isNameChar(m_source[m_pos])) {
                    ++m_pos;
                }
                add(Kind::Name, std::string(m_source.substr(start, m_pos - start)), start);
            } else if (isDigit(c) || (c == '.' && isDigit(peek(1)))) {
                number();
            } else if (c == '"' || c == '\'') {
                if (!quoted(c)) {
                    return;
                }
            } else if (c == '[' && longBracketLevel() >= 0) {
                const std::size_t start = m_pos;
                const int line = m_line;
                const int column = columnOf(start);
                std::string content;
                if (!longBracket(content)) {
                    fail(line, column, "unfinished long string");
                    return;
                }
                Token token{Kind::String, std::move(content), line, column, false};
                m_tokens.push_back(std::move(token));
            } else {
                op();
            }
        }
    }

    std::map<int, std::vector<std::string>> waivers;  // line -> rules
    std::vector<LintFinding> findings;

    std::vector<Token>& result() { return m_tokens; }

private:
    char peek(std::size_t ahead) const {
        return m_pos + ahead < m_source.size() ? m_source[m_pos + ahead] : '\0';
    }

    int columnOf(std::size_t pos) const { return static_cast<int>(pos - m_lineStart) + 1; }

    void newline() {
        const char c = m_source[m_pos++];
        if (m_pos < m_source.size() && (m_source[m_pos] == '\n' || m_source[m_pos] == '\r') && m_source[m_pos] != c) {
            ++m_pos;  // \r\n or \n\r
        }
        ++m_line;
        m_lineStart = m_pos;
    }

    void add(Kind kind, std::string text, std::size_t start, bool isFloat = false) {
        m_tokens.push_back(Token{kind, std::move(text), m_line, columnOf(start), isFloat});
    }

    void fail(int line, int column, std::string message) {
        findings.push_back(LintFinding{line, column, "syntax", std::move(message)});
    }

    // At '[': the level of a long bracket ("[==[" is 2), or -1.
    int longBracketLevel() const {
        std::size_t at = m_pos + 1;
        int level = 0;
        while (at < m_source.size() && m_source[at] == '=') {
            ++level;
            ++at;
        }
        return at < m_source.size() && m_source[at] == '[' ? level : -1;
    }

    // From '[' of a long bracket to after its closing bracket; false if it does not close.
    bool longBracket(std::string& content) {
        const int level = longBracketLevel();
        m_pos += static_cast<std::size_t>(level) + 2;
        const std::string close = "]" + std::string(static_cast<std::size_t>(level), '=') + "]";
        while (m_pos < m_source.size()) {
            if (m_source.compare(m_pos, close.size(), close) == 0) {
                m_pos += close.size();
                return true;
            }
            if (m_source[m_pos] == '\n' || m_source[m_pos] == '\r') {
                const std::size_t before = m_pos;
                newline();
                content.append(m_source.substr(before, m_pos - before));
            } else {
                content += m_source[m_pos++];
            }
        }
        return false;
    }

    bool comment() {
        const int line = m_line;
        const int column = columnOf(m_pos);
        m_pos += 2;
        std::string text;
        if (m_pos < m_source.size() && m_source[m_pos] == '[' && longBracketLevel() >= 0) {
            if (!longBracket(text)) {
                fail(line, column, "unfinished long comment");
                return false;
            }
        } else {
            const std::size_t start = m_pos;
            while (m_pos < m_source.size() && m_source[m_pos] != '\n' && m_source[m_pos] != '\r') {
                ++m_pos;
            }
            text = std::string(m_source.substr(start, m_pos - start));
        }
        waiver(text, line, column);
        return true;
    }

    // "lint-allow: rule, rule ..." at the start of a comment
    void waiver(std::string_view text, int line, int column) {
        std::size_t at = text.find_first_not_of(" \t-");
        constexpr std::string_view kTag = "lint-allow:";
        if (at == std::string_view::npos || text.substr(at, kTag.size()) != kTag) {
            return;
        }
        at += kTag.size();
        std::vector<std::string>& rules = waivers[line];
        while (at < text.size()) {
            const char c = text[at];
            if (c == ' ' || c == '\t' || c == ',') {
                ++at;
                continue;
            }
            if (c < 'a' || c > 'z') {
                break;  // the rest is a reason
            }
            const std::size_t start = at;
            while (at < text.size() && text[at] >= 'a' && text[at] <= 'z') {
                ++at;
            }
            const std::string rule(text.substr(start, at - start));
            if (!contains(kRules, rule)) {
                findings.push_back(LintFinding{line, column, "waiver", "the waiver names an unknown rule '" + rule + "'"});
            }
            rules.push_back(rule);
        }
    }

    void number() {
        const std::size_t start = m_pos;
        bool isFloat = false;
        const bool hex = m_source[m_pos] == '0' && (peek(1) == 'x' || peek(1) == 'X');
        if (hex) {
            m_pos += 2;
        }
        const char exponent1 = hex ? 'p' : 'e';
        const char exponent2 = hex ? 'P' : 'E';
        while (m_pos < m_source.size()) {
            const char c = m_source[m_pos];
            if (c == exponent1 || c == exponent2) {
                isFloat = true;
                ++m_pos;
                if (m_pos < m_source.size() && (m_source[m_pos] == '+' || m_source[m_pos] == '-')) {
                    ++m_pos;
                }
            } else if (c == '.') {
                if (peek(1) == '.') {
                    break;  // "1..x" is a malformed number in Lua; stop before the operator here
                }
                isFloat = true;
                ++m_pos;
            } else if (isNameChar(c)) {
                ++m_pos;
            } else {
                break;
            }
        }
        add(Kind::Number, std::string(m_source.substr(start, m_pos - start)), start, isFloat);
    }

    bool quoted(char quote) {
        const std::size_t start = m_pos;
        const int line = m_line;
        const int column = columnOf(start);
        ++m_pos;
        std::string content;
        while (m_pos < m_source.size()) {
            const char c = m_source[m_pos];
            if (c == quote) {
                ++m_pos;
                m_tokens.push_back(Token{Kind::String, std::move(content), line, column, false});
                return true;
            }
            if (c == '\n' || c == '\r') {
                break;
            }
            if (c == '\\' && m_pos + 1 < m_source.size()) {
                const char next = m_source[m_pos + 1];
                content += c;
                if (next == '\n' || next == '\r') {
                    ++m_pos;
                    content += next;
                    newline();
                    continue;
                }
                if (next == 'z') {  // skips the following white space, line breaks included
                    m_pos += 2;
                    while (m_pos < m_source.size()) {
                        const char w = m_source[m_pos];
                        if (w == '\n' || w == '\r') {
                            newline();
                        } else if (w == ' ' || w == '\t' || w == '\f' || w == '\v') {
                            ++m_pos;
                        } else {
                            break;
                        }
                    }
                    continue;
                }
                content += next;
                m_pos += 2;
                continue;
            }
            content += c;
            ++m_pos;
        }
        fail(line, column, "unfinished string");
        return false;
    }

    void op() {
        static constexpr std::array kLong = {std::string_view("..."), std::string_view(".."), std::string_view("::"),
                                             std::string_view("=="),  std::string_view("~="), std::string_view("<="),
                                             std::string_view(">="),  std::string_view("<<"), std::string_view(">>"),
                                             std::string_view("//")};
        const std::size_t start = m_pos;
        for (const std::string_view candidate : kLong) {
            if (m_source.compare(m_pos, candidate.size(), candidate) == 0) {
                m_pos += candidate.size();
                add(Kind::Op, std::string(candidate), start);
                return;
            }
        }
        ++m_pos;
        add(Kind::Op, std::string(1, m_source[start]), start);
    }

    std::string_view m_source;
    std::size_t m_pos = 0;
    std::size_t m_lineStart = 0;
    int m_line = 1;
    std::vector<Token> m_tokens;
};

bool isOp(const std::vector<Token>& tokens, std::size_t i, std::string_view text) {
    return i < tokens.size() && tokens[i].kind == Kind::Op && tokens[i].text == text;
}

bool isName(const std::vector<Token>& tokens, std::size_t i, std::string_view text) {
    return i < tokens.size() && tokens[i].kind == Kind::Name && tokens[i].text == text;
}

// The first float conversion of a format string ("%.3f"), or empty.
std::string floatConversion(std::string_view format) {
    for (std::size_t at = 0; at < format.size(); ++at) {
        if (format[at] != '%') {
            continue;
        }
        const std::size_t start = at++;
        if (at < format.size() && format[at] == '%') {
            continue;
        }
        while (at < format.size() && std::string_view("-+ #0").find(format[at]) != std::string_view::npos) {
            ++at;
        }
        while (at < format.size() && (isDigit(format[at]) || format[at] == '.')) {
            ++at;
        }
        if (at < format.size() && std::string_view("aAeEfFgG").find(format[at]) != std::string_view::npos) {
            return std::string(format.substr(start, at - start + 1));
        }
    }
    return {};
}

class Checker {
public:
    explicit Checker(const std::vector<Token>& tokens) : m_tokens(tokens) {}

    std::vector<LintFinding> run() {
        for (std::size_t i = 0; i < m_tokens.size(); ++i) {
            const Token& token = m_tokens[i];
            if (token.kind == Kind::Op) {
                bracket(token.text);
                if (token.text == "..") {
                    concat(i);
                } else if (token.text == ":" && isName(m_tokens, i + 1, "format")) {
                    methodFormat(i);
                }
                continue;
            }
            if (token.kind == Kind::String) {
                if (token.text == "__gc" || token.text == "__mode") {
                    add(token, "collector", "'" + token.text + "' depends on when the garbage collector runs");
                }
                continue;
            }
            if (token.kind != Kind::Name) {
                continue;
            }
            if (token.text == "__gc" || token.text == "__mode") {
                add(token, "collector", "'" + token.text + "' depends on when the garbage collector runs");
            }
            if (!isFree(i)) {
                continue;
            }
            if (contains(kMissing, token.text)) {
                add(token, "library",
                    "'" + token.text + "' does not exist in the sandbox (no os, io, debug, package or coroutine "
                                       "library, no files, printing or modules)");
            } else if (token.text == "pairs" || token.text == "next") {
                add(token, "pairs",
                    "'" + token.text + "' visits keys in hash order: iterate arrays with ipairs and maps with "
                                       "mdb.sorted_pairs (or waive it with '-- lint-allow: pairs' where the order "
                                       "cannot reach the output)");
            } else if (token.text == "tostring") {
                add(token, "tostring", "tostring writes floats unlike Python (%.14g): format numbers with mdb.text");
            } else if (token.text == "string" && isOp(m_tokens, i + 1, ".") && isName(m_tokens, i + 2, "format")) {
                functionFormat(i);
            }
        }
        return std::move(m_findings);
    }

private:
    void add(const Token& token, std::string rule, std::string message) {
        m_findings.push_back(LintFinding{token.line, token.column, std::move(rule), std::move(message)});
    }

    void bracket(const std::string& text) {
        if (text == "{" || text == "(" || text == "[") {
            m_brackets.push_back(text[0]);
        } else if ((text == "}" || text == ")" || text == "]") && !m_brackets.empty()) {
            m_brackets.pop_back();
        }
    }

    // A name that stands for a variable: not a field (`x.name`, `x:name`), a label or goto target,
    // or a key of a table constructor (`{ name = ... }`).
    bool isFree(std::size_t i) const {
        if (i > 0) {
            const Token& previous = m_tokens[i - 1];
            if (previous.kind == Kind::Op && (previous.text == "." || previous.text == ":" || previous.text == "::")) {
                return false;
            }
            if (previous.kind == Kind::Name && previous.text == "goto") {
                return false;
            }
        }
        if (!m_brackets.empty() && m_brackets.back() == '{' && isOp(m_tokens, i + 1, "=") &&
            (i == 0 || isOp(m_tokens, i - 1, "{") || isOp(m_tokens, i - 1, ",") || isOp(m_tokens, i - 1, ";"))) {
            return false;
        }
        return true;
    }

    void checkLiteral(const Token& where, const Token& literal) {
        const std::string conversion = floatConversion(literal.text);
        if (!conversion.empty()) {
            add(where, "format",
                "'" + conversion + "' formats a float with C's printf, not Python's rules: use mdb.text.format");
        }
    }

    // string.format at i ("string")
    void functionFormat(std::size_t i) {
        const Token& where = m_tokens[i];
        if (i + 3 < m_tokens.size() && m_tokens[i + 3].kind == Kind::String) {
            checkLiteral(where, m_tokens[i + 3]);
        } else if (isOp(m_tokens, i + 3, "(") && i + 4 < m_tokens.size() && m_tokens[i + 4].kind == Kind::String) {
            checkLiteral(where, m_tokens[i + 4]);
        } else {
            add(where, "format", "the format of string.format is not a literal string, so its conversions cannot be checked");
        }
    }

    // ("..."):format at i (":")
    void methodFormat(std::size_t i) {
        const Token& where = m_tokens[i + 1];
        if (i >= 3 && isOp(m_tokens, i - 1, ")") && m_tokens[i - 2].kind == Kind::String && isOp(m_tokens, i - 3, "(")) {
            checkLiteral(where, m_tokens[i - 2]);
        } else if (i >= 1 && m_tokens[i - 1].kind == Kind::String) {
            checkLiteral(where, m_tokens[i - 1]);
        } else {
            add(where, "format", "the format of :format is not a literal string, so its conversions cannot be checked");
        }
    }

    // `..` at i
    void concat(std::size_t i) {
        const Token& where = m_tokens[i];
        bool isFloat = i > 0 && m_tokens[i - 1].kind == Kind::Number && m_tokens[i - 1].isFloat;
        std::size_t next = i + 1;
        if (isOp(m_tokens, next, "-")) {
            ++next;
        }
        if (next < m_tokens.size() && m_tokens[next].kind == Kind::Number && m_tokens[next].isFloat) {
            isFloat = true;
        }
        if (isName(m_tokens, next, "math") && isOp(m_tokens, next + 1, ".") && next + 2 < m_tokens.size() &&
            m_tokens[next + 2].kind == Kind::Name && contains(kFloatMath, m_tokens[next + 2].text)) {
            isFloat = true;
        }
        if (isFloat) {
            add(where, "concat", "'..' writes the float with Lua's %.14g, unlike Python: format it with mdb.text");
        }
    }

    const std::vector<Token>& m_tokens;
    std::vector<char> m_brackets;
    std::vector<LintFinding> m_findings;
};

}  // namespace

std::vector<LintFinding> lint(std::string_view source) {
    Lexer lexer(source);
    lexer.run();
    std::vector<LintFinding> findings = Checker(lexer.result()).run();
    findings.insert(findings.end(), lexer.findings.begin(), lexer.findings.end());
    std::vector<LintFinding> kept;
    for (LintFinding& finding : findings) {
        const auto waived = lexer.waivers.find(finding.line);
        if (finding.rule != "waiver" && finding.rule != "syntax" && waived != lexer.waivers.end() &&
            std::find(waived->second.begin(), waived->second.end(), finding.rule) != waived->second.end()) {
            continue;
        }
        kept.push_back(std::move(finding));
    }
    std::stable_sort(kept.begin(), kept.end(), [](const LintFinding& a, const LintFinding& b) {
        return a.line != b.line ? a.line < b.line : a.column < b.column;
    });
    return kept;
}

std::string formatFinding(std::string_view file, const LintFinding& finding) {
    return std::string(file) + ":" + std::to_string(finding.line) + ":" + std::to_string(finding.column) + ": " +
           finding.rule + ": " + finding.message;
}

}  // namespace mdbridge::workflow::lua
