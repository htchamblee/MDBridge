#pragma once
// Text as Python handles it: MDATK's files are written with Python's string functions, so the
// C++ side counts, splits and wraps characters (code points) exactly as Python does.

#include <string>
#include <string_view>
#include <vector>

namespace mdbridge::core {

// UTF-8 text as code points. A byte that is not part of a valid UTF-8 sequence becomes U+FFFD,
// as an undecodable byte (surrogate-escaped by Python) becomes U+FFFD when MDATK writes it.
std::u32string decodeUtf8(std::string_view text);
// UTF-8 text as Python's surrogateescape error handler decodes it (os.fsdecode, sys.argv): a byte
// that is not part of a valid sequence becomes the lone surrogate U+DC00 + byte, which repr() and
// json.dumps show as \udcXX.
std::u32string decodeUtf8Surrogateescape(std::string_view text);
std::string encodeUtf8(std::u32string_view text);
void appendUtf8(std::string& out, char32_t code);

// Number of characters (code points), as Python's len() counts them.
std::size_t characterCount(std::string_view text);

// Python's str.isspace() / \w / \d / str.isprintable() for one character (Unicode data of
// tools/gen_unicode_tables.py). Printable: everything except the categories Cc, Cf, Cs, Co, Cn,
// Zl, Zp and Zs other than the space, as repr() decides what to escape.
bool isPythonSpace(char32_t code);
bool isPythonWordCharacter(char32_t code);
bool isPythonDigit(char32_t code);
bool isPythonPrintable(char32_t code);

// Python's str.splitlines(): splits at \n, \r\n, \r, \v, \f, \x1c-\x1e, \x85, U+2028 and U+2029.
std::vector<std::string> pythonSplitlines(std::string_view text);

// Python's str.strip() / rstrip() (Unicode whitespace).
std::string pythonStrip(std::string_view text);
std::string pythonRstrip(std::string_view text);

// Python's textwrap.wrap(text, width) with the default options: tabs expanded, whitespace
// characters turned into spaces, words broken after hyphens and when longer than a line.
std::vector<std::string> pythonTextwrap(std::string_view text, int width);

}  // namespace mdbridge::core
