#include "mdbridge/core/hash.hpp"

#include <array>
#include <cstdint>
#include <system_error>

#if defined(_WIN32)
#include <ratio>
#include <type_traits>
#else
#include <cerrno>
#include <sys/stat.h>
#endif

#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/python.hpp"

namespace mdbridge::core {

namespace {

// SHA-256 (FIPS 180-4)
class Sha256 {
public:
    void update(std::string_view data) {
        for (const char c : data) {
            m_block[m_used++] = static_cast<std::uint8_t>(c);
            if (m_used == 64) {
                compress();
                m_used = 0;
            }
        }
        m_bits += 8 * static_cast<std::uint64_t>(data.size());
    }

    std::string hexDigest() {
        const std::uint64_t bits = m_bits;
        m_block[m_used++] = 0x80;
        if (m_used > 56) {
            while (m_used < 64) {
                m_block[m_used++] = 0;
            }
            compress();
            m_used = 0;
        }
        while (m_used < 56) {
            m_block[m_used++] = 0;
        }
        for (int k = 7; k >= 0; --k) {
            m_block[m_used++] = static_cast<std::uint8_t>(bits >> (8 * k));
        }
        compress();
        static constexpr char kHex[] = "0123456789abcdef";
        std::string out;
        for (const std::uint32_t word : m_state) {
            for (int shift = 28; shift >= 0; shift -= 4) {
                out += kHex[(word >> shift) & 0xF];
            }
        }
        return out;
    }

private:
    static std::uint32_t rotr(std::uint32_t x, int n) { return (x >> n) | (x << (32 - n)); }

    void compress() {
        static constexpr std::array<std::uint32_t, 64> kRound = {
            0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4,
            0xab1c5ed5, 0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe,
            0x9bdc06a7, 0xc19bf174, 0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f,
            0x4a7484aa, 0x5cb0a9dc, 0x76f988da, 0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
            0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967, 0x27b70a85, 0x2e1b2138, 0x4d2c6dfc,
            0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85, 0xa2bfe8a1, 0xa81a664b,
            0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070, 0x19a4c116,
            0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
            0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7,
            0xc67178f2};
        std::array<std::uint32_t, 64> w{};
        for (std::size_t t = 0; t < 16; ++t) {
            const auto byte = [&](std::size_t k) { return std::uint32_t{m_block[4 * t + k]}; };
            w[t] = (byte(0) << 24) | (byte(1) << 16) | (byte(2) << 8) | byte(3);
        }
        for (std::size_t t = 16; t < 64; ++t) {
            const std::uint32_t s0 = rotr(w[t - 15], 7) ^ rotr(w[t - 15], 18) ^ (w[t - 15] >> 3);
            const std::uint32_t s1 = rotr(w[t - 2], 17) ^ rotr(w[t - 2], 19) ^ (w[t - 2] >> 10);
            w[t] = w[t - 16] + s0 + w[t - 7] + s1;
        }
        std::array<std::uint32_t, 8> v = m_state;
        for (std::size_t t = 0; t < 64; ++t) {
            const std::uint32_t e = v[4];
            const std::uint32_t a = v[0];
            const std::uint32_t t1 = v[7] + (rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25)) +
                                     ((e & v[5]) ^ (~e & v[6])) + kRound[t] + w[t];
            const std::uint32_t t2 = (rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22)) +
                                     ((a & v[1]) ^ (a & v[2]) ^ (v[1] & v[2]));
            v = {t1 + t2, a, v[1], v[2], v[3] + t1, e, v[5], v[6]};
        }
        for (std::size_t k = 0; k < 8; ++k) {
            m_state[k] += v[k];
        }
    }

    std::array<std::uint32_t, 8> m_state = {0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
                                            0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19};
    std::array<std::uint8_t, 64> m_block{};
    std::size_t m_used = 0;
    std::uint64_t m_bits = 0;
};

// _to_plain of a dict: None values of tables are left out at every depth, array items stay
Value withoutNulls(const Value& value) {
    if (const Table* table = value.get<Table>()) {
        Table out;
        for (const auto& [key, item] : *table) {
            if (!item.isNull()) {
                out.emplace_back(key, withoutNulls(item));
            }
        }
        return out;
    }
    if (const Array* array = value.get<Array>()) {
        Array out;
        out.reserve(array->size());
        for (const Value& item : *array) {
            out.push_back(withoutNulls(item));
        }
        return out;
    }
    return value;
}

}  // namespace

std::string sha256Hex(std::string_view data) {
    Sha256 hash;
    hash.update(data);
    return hash.hexDigest();
}

std::string stableHash(const Value& plain, std::size_t length) {
    return sha256Hex(pythonJsonCompact(withoutNulls(plain))).substr(0, length);
}

Value fileFingerprint(const std::filesystem::path& path) {
    std::string absolute;
    std::int64_t size = 0;
    std::int64_t mtimeNs = 0;
#if defined(_WIN32)
    // what Python's os.stat reports: FILETIME (100 ns since 1601) converted to the Unix epoch
    using Ticks = std::filesystem::file_time_type::period;
    static_assert(std::is_same_v<Ticks, std::ratio<1, 10'000'000>>, "MSVC's file_clock ticks");
    size = static_cast<std::int64_t>(std::filesystem::file_size(path));
    const auto ticks = std::filesystem::last_write_time(path).time_since_epoch().count();
    mtimeNs = (static_cast<std::int64_t>(ticks) - 116'444'736'000'000'000) * 100;
    absolute = pathToUtf8(std::filesystem::absolute(path));
#else
    struct stat st {};
    if (::stat(path.c_str(), &st) != 0) {
        throw std::filesystem::filesystem_error("cannot stat the file", path,
                                                std::error_code(errno, std::generic_category()));
    }
    size = static_cast<std::int64_t>(st.st_size);
#if defined(__APPLE__)
    const timespec& modified = st.st_mtimespec;
#else
    const timespec& modified = st.st_mtim;
#endif
    mtimeNs = static_cast<std::int64_t>(modified.tv_sec) * 1'000'000'000 + modified.tv_nsec;
    absolute = posixpath::abspath(pathToUtf8(path), currentDirectory());
#endif
    return Table{{"path", absolute}, {"size", size}, {"mtime_ns", mtimeNs}};
}

}  // namespace mdbridge::core
