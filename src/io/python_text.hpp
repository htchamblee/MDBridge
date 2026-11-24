#pragma once
// Python str operations on code points for the fixed-column PDB / GRO code (internal to io):
// slicing, strip / split with Unicode whitespace, str.upper(), isalpha / isalnum, int() and
// float(), ljust / rjust, and %-style fixed formatting of doubles.

#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "mdbridge/core/errors.hpp"

namespace mdbridge::io::detail {

using U32 = std::u32string;
using U32View = std::u32string_view;

U32 toU32(std::string_view utf8);
std::string toUtf8(U32View text);

// text[first:last] with Python's clamping (non-negative bounds; last = npos for text[first:]).
U32 slice(U32View text, std::size_t first, std::size_t last = U32::npos);

bool isSpace(char32_t c);
U32 strip(U32View text);
U32 lstrip(U32View text);
U32 rstrip(U32View text);
// str.split() without arguments
std::vector<U32> split(U32View text);
U32 ljust(U32View text, std::size_t width);
U32 rjust(U32View text, std::size_t width);

// str.upper() (full case mapping, Unicode 15.1)
U32 upper(U32View text);
bool isAlpha(char32_t c);
bool isAlnum(char32_t c);
bool isDigit(char32_t c);       // re's \d
bool isAllAlpha(U32View text);  // str.isalpha()

// int(text) / float(text) of a str; std::nullopt where Python raises ValueError (and for ints
// beyond 64 bits)
std::optional<std::int64_t> toInt(U32View text);
std::optional<double> toFloat(U32View text);

// Python's UnicodeDecodeError of a strict UTF-8 decode
class UnicodeDecodeError : public core::ValueError {
public:
    using core::ValueError::ValueError;
    const char* pythonName() const noexcept override { return "UnicodeDecodeError"; }
};

// The message of the UnicodeDecodeError of bytes.decode() ("'utf-8' codec can't decode byte 0xe9
// in position 10: invalid continuation byte"), std::nullopt for valid UTF-8.
std::optional<std::string> utf8Error(std::string_view bytes);
// bytes.decode(); throws UnicodeDecodeError
U32 decodeStrict(std::string_view bytes);
// \r\n and \r as \n (a text file's universal newlines)
U32 universalNewlines(U32View text);
// The lines of iterating a text file: each with its final \n
std::vector<U32> fileLines(U32View text);

// f"{x:{width}.{precision}f}" of a double (Python's correctly rounded fixed notation).
std::string fixed(double x, int width, int precision);
// f"{n:{width}d}"
std::string padInt(std::int64_t n, int width);

}  // namespace mdbridge::io::detail
