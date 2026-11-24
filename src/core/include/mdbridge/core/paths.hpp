#pragma once
// UTF-8 text <-> file system paths. MDBridge keeps text in UTF-8; on Windows the native path
// encoding is UTF-16, so conversions must not go through the narrow (ANSI) code page.

#include <filesystem>
#include <string>
#include <string_view>

namespace mdbridge::core {

std::filesystem::path pathFromUtf8(std::string_view text);
std::string pathToUtf8(const std::filesystem::path& path);

// os.path.basename of a UTF-8 path: posixpath's rules on POSIX; on Windows '\\' separates too and
// a drive ("C:") is not part of the name, as ntpath.
std::string pathBasename(std::string_view path);

}  // namespace mdbridge::core
