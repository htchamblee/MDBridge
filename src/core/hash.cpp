#include "mdbridge/core/hash.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
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

constexpr std::uint32_t kRound[64] = {
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2};

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

// Plain scalars and arrays, no calls (the rotation is a macro): fast enough unoptimised (Debug
// builds hash gigabytes in the periodic search tests).
#define MDBRIDGE_ROTR(x, n) (((x) >> (n)) | ((x) << (32 - (n))))
void Sha256::compress(const std::uint8_t* block) {
    std::uint32_t w[64];
    for (int t = 0; t < 16; ++t) {
        const std::uint8_t* p = block + 4 * t;
        w[t] = (std::uint32_t{p[0]} << 24) | (std::uint32_t{p[1]} << 16) | (std::uint32_t{p[2]} << 8) |
               std::uint32_t{p[3]};
    }
    for (int t = 16; t < 64; ++t) {
        const std::uint32_t x = w[t - 15];
        const std::uint32_t y = w[t - 2];
        w[t] = w[t - 16] + (MDBRIDGE_ROTR(x, 7) ^ MDBRIDGE_ROTR(x, 18) ^ (x >> 3)) + w[t - 7] +
               (MDBRIDGE_ROTR(y, 17) ^ MDBRIDGE_ROTR(y, 19) ^ (y >> 10));
    }
    std::uint32_t* state = m_state.data();
    std::uint32_t a = state[0], b = state[1], c = state[2], d = state[3];
    std::uint32_t e = state[4], f = state[5], g = state[6], h = state[7];
    for (int t = 0; t < 64; ++t) {
        const std::uint32_t t1 = h + (MDBRIDGE_ROTR(e, 6) ^ MDBRIDGE_ROTR(e, 11) ^ MDBRIDGE_ROTR(e, 25)) +
                                 ((e & f) ^ (~e & g)) + kRound[t] + w[t];
        const std::uint32_t t2 =
            (MDBRIDGE_ROTR(a, 2) ^ MDBRIDGE_ROTR(a, 13) ^ MDBRIDGE_ROTR(a, 22)) + ((a & b) ^ (a & c) ^ (b & c));
        h = g;
        g = f;
        f = e;
        e = d + t1;
        d = c;
        c = b;
        b = a;
        a = t1 + t2;
    }
    state[0] += a;
    state[1] += b;
    state[2] += c;
    state[3] += d;
    state[4] += e;
    state[5] += f;
    state[6] += g;
    state[7] += h;
}
#undef MDBRIDGE_ROTR

void Sha256::update(std::string_view data) {
    m_bits += 8 * static_cast<std::uint64_t>(data.size());
    const auto* p = reinterpret_cast<const std::uint8_t*>(data.data());
    std::size_t left = data.size();
    if (m_used > 0) {  // complete the pending block first
        const std::size_t take = std::min(left, std::size_t{64} - m_used);
        std::memcpy(m_block.data() + m_used, p, take);
        m_used += take;
        p += take;
        left -= take;
        if (m_used < 64) {
            return;
        }
        compress(m_block.data());
        m_used = 0;
    }
    for (; left >= 64; p += 64, left -= 64) {
        compress(p);
    }
    std::memcpy(m_block.data(), p, left);
    m_used = left;
}

std::string Sha256::hexDigest() {
    const std::uint64_t bits = m_bits;
    m_block[m_used++] = 0x80;
    if (m_used > 56) {
        while (m_used < 64) {
            m_block[m_used++] = 0;
        }
        compress(m_block.data());
        m_used = 0;
    }
    while (m_used < 56) {
        m_block[m_used++] = 0;
    }
    for (int k = 7; k >= 0; --k) {
        m_block[m_used++] = static_cast<std::uint8_t>(bits >> (8 * k));
    }
    compress(m_block.data());
    static constexpr char kHex[] = "0123456789abcdef";
    std::string out;
    for (const std::uint32_t word : m_state) {
        for (int shift = 28; shift >= 0; shift -= 4) {
            out += kHex[(word >> shift) & 0xF];
        }
    }
    return out;
}

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
