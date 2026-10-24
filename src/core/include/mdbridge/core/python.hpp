#pragma once
// What MDATK's Python does with plain data and paths, reproduced where its output reaches users or
// files: repr() / str() of values in messages, round(), json.dumps for hashes, and POSIX
// os.path functions for resolving project paths (parity is defined on Linux).

#include <cstdint>
#include <optional>
#include <string>
#include <string_view>

#include "mdbridge/core/value.hpp"

namespace mdbridge::core {

// repr() of plain data: 'text' (Python's quoting and escapes, incl. non-printable Unicode),
// True / False / None, floats as repr, [a, b], {'k': v}, datetime.datetime(2026, 10, 15, 14, 32, 5) ...
std::string pythonRepr(const Value& value);

// str() of plain data: strings unchanged, True / False / None, datetimes with a space
// ("2026-10-15 14:32:05"), lists and tables as their repr.
std::string pythonStr(const Value& value);

// type(value).__name__: str, int, float, bool, list, dict, datetime, date, time, NoneType.
std::string_view pythonTypeName(const Value& value);

// round(x, digits): correctly rounded on the exact binary value, ties to even. Throws
// std::overflow_error("rounded value too large to represent") as Python does (digits < 0 only).
double pythonRound(double x, int digits);

// Python's json.dumps(value, sort_keys=True, separators=(",", ":"), default=str) with
// ensure_ascii: compact, keys sorted by code point, non-ASCII as \uXXXX (surrogate pairs beyond
// U+FFFF), NaN / Infinity, dates as str(). Bytes that are not UTF-8 count as surrogateescape's
// lone surrogates (\udcXX), as in a Python str of an undecodable file name.
std::string pythonJsonCompact(const Value& value);

// int(text) and float(text) of a str, as Python parses them: surrounding Unicode whitespace, a
// sign, single underscores between digits, Unicode decimal digits; float also nan / inf /
// infinity in any case, overflow to inf and underflow to 0. Throws ValueError with Python's
// message: "invalid literal for int() with base 10: '<repr>'", "could not convert string to float:
// '<repr>'"; an int beyond 64 bits throws std::overflow_error("Python int too large to convert to
// C long").
std::int64_t pythonInt(std::string_view text);
double pythonFloat(std::string_view text);
// The same without exceptions.
std::optional<std::int64_t> tryPythonInt(std::string_view text);
std::optional<double> tryPythonFloat(std::string_view text);

// POSIX os.path (posixpath) on UTF-8 strings.
namespace posixpath {
// "~" and "~user" expansion; $HOME first for "~" (set, even empty), else the passwd entry
// (unchanged path when no home is known). Windows, where MDATK has no reference behaviour:
// "~" is %USERPROFILE% and "~user" stays as it is.
std::string expanduser(std::string_view path);
bool isabs(std::string_view path);
std::string join(std::string_view base, std::string_view path);
std::string normpath(std::string_view path);
// normpath(join(cwd, path)) for a relative path, normpath(path) for an absolute one
std::string abspath(std::string_view path, std::string_view cwd);
std::string dirname(std::string_view path);
std::string basename(std::string_view path);
// (root, ext): "a.b.toml" -> ("a.b", ".toml"), ".hidden" -> (".hidden", "")
std::pair<std::string, std::string> splitext(std::string_view path);
}  // namespace posixpath

// The process's current directory as a UTF-8 string.
std::string currentDirectory();

}  // namespace mdbridge::core
