#pragma once
// JSON and Python-style text of plain data, formatted like MDATK's Python side so that outputs
// can be compared byte for byte (parity checks against MDATK v0.0.4.1).

#include <string>

#include "mdbridge/core/value.hpp"

namespace mdbridge::core {

// The text of Python's json.dumps(value, indent=indent, ensure_ascii=False): dates and times as
// their isoformat() strings, floats as repr(), NaN / Infinity as Python writes them.
// indent < 0 gives the compact form of json.dumps(value, ensure_ascii=False).
std::string toJson(const Value& value, int indent = 2);

// Python's repr() of a float: the shortest text that reads back as the same number, in fixed
// notation for exponents -4 to 15 and with ".0" for integral values (1.0, 0.0001, 1e-05, 1e+16).
std::string pythonFloatRepr(double value);

// Python's repr() of a str: single quotes unless the text holds ' and no ".
std::string pythonStringRepr(const std::string& text);

}  // namespace mdbridge::core
