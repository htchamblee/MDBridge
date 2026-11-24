#pragma once
// Random-access reading of binary trajectory files and explicit-endian decoding of their values
// (independent of the host's byte order).

#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <string>
#include <system_error>

namespace mdbridge::io::detail {

class BinaryFile {
public:
    explicit BinaryFile(const std::filesystem::path& path) : m_in(path, std::ios::binary) {
        if (!m_in) {
            throw std::filesystem::filesystem_error(
                "cannot read file", path,
                std::make_error_code(std::errc::no_such_file_or_directory));
        }
        m_in.seekg(0, std::ios::end);
        m_size = static_cast<std::uint64_t>(m_in.tellg());
        m_in.seekg(0);
    }

    std::uint64_t size() const { return m_size; }

    // Up to n bytes at `offset`; the count read (less at the end of the file).
    std::size_t read(std::uint64_t offset, void* out, std::size_t n) {
        if (offset >= m_size || n == 0) {
            return 0;
        }
        m_in.clear();
        m_in.seekg(static_cast<std::streamoff>(offset));
        m_in.read(static_cast<char*>(out), static_cast<std::streamsize>(n));
        return static_cast<std::size_t>(m_in.gcount());
    }

    std::string read(std::uint64_t offset, std::size_t n) {
        std::string out(n, '\0');
        out.resize(read(offset, out.data(), n));
        return out;
    }

private:
    std::ifstream m_in;
    std::uint64_t m_size = 0;
};

inline std::uint16_t loadU16(const void* p, bool bigEndian) {
    unsigned char b[2];
    std::memcpy(b, p, 2);
    return bigEndian ? static_cast<std::uint16_t>(b[0] << 8 | b[1])
                     : static_cast<std::uint16_t>(b[1] << 8 | b[0]);
}

inline std::uint32_t loadU32(const void* p, bool bigEndian) {
    unsigned char b[4];
    std::memcpy(b, p, 4);
    if (bigEndian) {
        return std::uint32_t{b[0]} << 24 | std::uint32_t{b[1]} << 16 | std::uint32_t{b[2]} << 8 |
               b[3];
    }
    return std::uint32_t{b[3]} << 24 | std::uint32_t{b[2]} << 16 | std::uint32_t{b[1]} << 8 | b[0];
}

inline std::uint64_t loadU64(const void* p, bool bigEndian) {
    const auto* b = static_cast<const unsigned char*>(p);
    const std::uint64_t first = loadU32(b, bigEndian);
    const std::uint64_t second = loadU32(b + 4, bigEndian);
    return bigEndian ? (first << 32 | second) : (second << 32 | first);
}

inline std::int32_t loadI32(const void* p, bool bigEndian) {
    const std::uint32_t bits = loadU32(p, bigEndian);
    std::int32_t value;
    std::memcpy(&value, &bits, 4);
    return value;
}

inline std::int64_t loadI64(const void* p, bool bigEndian) {
    const std::uint64_t bits = loadU64(p, bigEndian);
    std::int64_t value;
    std::memcpy(&value, &bits, 8);
    return value;
}

inline float loadF32(const void* p, bool bigEndian) {
    const std::uint32_t bits = loadU32(p, bigEndian);
    float value;
    std::memcpy(&value, &bits, 4);
    return value;
}

inline double loadF64(const void* p, bool bigEndian) {
    const std::uint64_t bits = loadU64(p, bigEndian);
    double value;
    std::memcpy(&value, &bits, 8);
    return value;
}

inline void storeU32BE(std::string& out, std::uint32_t value) {
    out += static_cast<char>(value >> 24 & 0xFF);
    out += static_cast<char>(value >> 16 & 0xFF);
    out += static_cast<char>(value >> 8 & 0xFF);
    out += static_cast<char>(value & 0xFF);
}

inline void storeU64BE(std::string& out, std::uint64_t value) {
    storeU32BE(out, static_cast<std::uint32_t>(value >> 32));
    storeU32BE(out, static_cast<std::uint32_t>(value & 0xFFFFFFFFu));
}

inline void storeF32BE(std::string& out, float value) {
    std::uint32_t bits;
    std::memcpy(&bits, &value, 4);
    storeU32BE(out, bits);
}

inline void storeF64BE(std::string& out, double value) {
    std::uint64_t bits;
    std::memcpy(&bits, &value, 8);
    storeU64BE(out, bits);
}

}  // namespace mdbridge::io::detail
