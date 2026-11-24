#pragma once
// Text files as Python reads them with open(path, "r", encoding="utf-8", errors="replace"): every
// maximal invalid UTF-8 subpart becomes one U+FFFD, and iteration yields lines ending at \n, \r\n
// or \r, each translated to a final "\n" (universal newlines). Lines are code points, so that
// slicing counts characters as Python does.

#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <string_view>
#include <system_error>
#include <vector>

namespace mdbridge::io::detail {

using Lines = std::vector<std::u32string>;

inline std::u32string decodeUtf8Replace(std::string_view bytes) {
    std::u32string out;
    out.reserve(bytes.size());
    const auto byte = [&](std::size_t i) {
        return static_cast<unsigned char>(bytes[i]);
    };
    std::size_t i = 0;
    while (i < bytes.size()) {
        const unsigned char b = byte(i);
        if (b < 0x80) {
            out += static_cast<char32_t>(b);
            ++i;
            continue;
        }
        std::size_t length = 0;
        unsigned char low = 0x80;
        unsigned char high = 0xBF;
        char32_t code = 0;
        if (b >= 0xC2 && b <= 0xDF) {
            length = 2;
            code = b & 0x1F;
        } else if (b >= 0xE0 && b <= 0xEF) {
            length = 3;
            code = b & 0x0F;
            low = b == 0xE0 ? 0xA0 : 0x80;
            high = b == 0xED ? 0x9F : 0xBF;
        } else if (b >= 0xF0 && b <= 0xF4) {
            length = 4;
            code = b & 0x07;
            low = b == 0xF0 ? 0x90 : 0x80;
            high = b == 0xF4 ? 0x8F : 0xBF;
        } else {
            out += U'�';
            ++i;
            continue;
        }
        std::size_t k = 1;
        for (; k < length && i + k < bytes.size(); ++k) {
            const unsigned char c = byte(i + k);
            if (c < (k == 1 ? low : 0x80) || c > (k == 1 ? high : 0xBF)) {
                break;
            }
            code = (code << 6) | (c & 0x3F);
        }
        if (k == length) {
            out += code;
        } else {
            out += U'�';
        }
        i += k;
    }
    return out;
}

inline Lines splitPythonLines(std::u32string_view text) {
    Lines lines;
    std::u32string line;
    for (std::size_t i = 0; i < text.size(); ++i) {
        const char32_t c = text[i];
        if (c == U'\r' || c == U'\n') {
            if (c == U'\r' && i + 1 < text.size() && text[i + 1] == U'\n') {
                ++i;
            }
            line += U'\n';
            lines.push_back(std::move(line));
            line.clear();
        } else {
            line += c;
        }
    }
    if (!line.empty()) {
        lines.push_back(std::move(line));
    }
    return lines;
}

// The path as MDATK's messages show it (UTF-8), and a path from MDATK's UTF-8 text.
inline std::string pathText(const std::filesystem::path& path) {
    const std::u8string text = path.u8string();
    return std::string(text.begin(), text.end());
}

inline std::filesystem::path pathFromUtf8(std::string_view text) {
    return std::filesystem::path(std::u8string(text.begin(), text.end()));
}

inline std::string readBytes(const std::filesystem::path& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw std::filesystem::filesystem_error(
            "cannot read file", path, std::make_error_code(std::errc::no_such_file_or_directory));
    }
    std::ostringstream text;
    text << in.rdbuf();
    return std::move(text).str();
}

inline Lines readTextLines(const std::filesystem::path& path) {
    return splitPythonLines(decodeUtf8Replace(readBytes(path)));
}

}  // namespace mdbridge::io::detail
