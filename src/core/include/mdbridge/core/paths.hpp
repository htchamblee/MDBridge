#pragma once
// UTF-8 text <-> file system paths. MDBridge keeps text in UTF-8; on Windows the native path
// encoding is UTF-16, so conversions must not go through the narrow (ANSI) code page.

#include <filesystem>
#include <string>
#include <string_view>

namespace mdbridge::core {

std::filesystem::path pathFromUtf8(std::string_view text);
std::string pathToUtf8(const std::filesystem::path& path);

}  // namespace mdbridge::core
