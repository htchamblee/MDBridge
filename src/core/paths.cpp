#include "mdbridge/core/paths.hpp"

#include "mdbridge/core/python.hpp"

namespace mdbridge::core {

std::filesystem::path pathFromUtf8(std::string_view text) {
    return std::filesystem::path(std::u8string(text.begin(), text.end()));
}

std::string pathToUtf8(const std::filesystem::path& path) {
    const std::u8string text = path.u8string();
    return std::string(text.begin(), text.end());
}

std::string pathBasename(std::string_view path) {
#ifdef _WIN32
    return pathToUtf8(pathFromUtf8(path).filename());
#else
    return posixpath::basename(path);
#endif
}

}  // namespace mdbridge::core
