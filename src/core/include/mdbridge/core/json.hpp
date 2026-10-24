#pragma once
// JSON and Python-style text of plain data, formatted like MDATK's Python side so that outputs
// can be compared byte for byte (parity checks against MDATK v0.0.4.1).

#include <string>
#include <string_view>

#include "mdbridge/core/project_file.hpp"
#include "mdbridge/core/value.hpp"

namespace mdbridge::core {

// The text of Python's json.dumps(value, indent=indent, ensure_ascii=False): dates and times as
// their isoformat() strings, floats as repr(), NaN / Infinity as Python writes them.
// indent < 0 gives the compact form of json.dumps(value, ensure_ascii=False).
std::string toJson(const Value& value, int indent = 2);

// Python's repr() of a float: the shortest text that reads back as the same number, in fixed
// notation for exponents -4 to 15 and with ".0" for integral values (1.0, 0.0001, 1e-05, 1e+16).
std::string pythonFloatRepr(double value);

// Python's repr() of a str: single quotes unless the text holds ' and no "; \\, the quote, \t, \n
// and \r escaped, other characters that are not str.isprintable() as \xHH, \uHHHH or \UHHHHHHHH.
// Bytes that are not valid UTF-8 show as the \udcXX of Python's surrogateescape.
std::string pythonStringRepr(const std::string& text);

// Python's json.loads(text): objects as tables in the order their keys first appear, holding the
// last value of a repeated key; integers as Int, other numbers as Float (also NaN, Infinity and
// -Infinity), read exactly. Throws ProjectError("<source>: not valid JSON: <json's message>"),
// e.g. "Expecting value: line 1 column 1 (char 0)". Where a Value cannot hold what Python reads,
// MDBridge deviates: a \u escape of a lone surrogate gives U+FFFD; an integer outside the 64-bit
// range that the result would hold is an error ("Integer outside the 64-bit range: ...", only
// once the text is otherwise valid), and so is nesting deeper than 1000 arrays / objects
// (Python's C scanner allows about 10000).
Value parseJson(std::string_view text, std::string_view source);

}  // namespace mdbridge::core
