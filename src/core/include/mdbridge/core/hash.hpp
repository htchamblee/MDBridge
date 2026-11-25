#pragma once
// Stable hashes of settings (cache keys) and cheap file identities, as MDATK's config.py computes
// them (stable_hash, file_fingerprint).

#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>
#include <string_view>

#include "mdbridge/core/value.hpp"

namespace mdbridge::core {

// SHA-256 of the bytes, as 64 lowercase hex digits.
std::string sha256Hex(std::string_view data);

// SHA-256 (FIPS 180-4) fed in pieces: update() with successive pieces, then hexDigest() gives
// sha256Hex of their concatenation (call it once; the object is spent afterwards). For digests of
// data too large to hold at once (the periodic search tests hash 1e8 pairs as they come).
class Sha256 {
public:
    void update(std::string_view data);
    std::string hexDigest();

private:
    void compress(const std::uint8_t* block);

    std::array<std::uint32_t, 8> m_state = {0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
                                            0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19};
    std::array<std::uint8_t, 64> m_block{};
    std::size_t m_used = 0;
    std::uint64_t m_bits = 0;
};

// stable_hash(obj, length): the first `length` hex digits of the SHA-256 of
// pythonJsonCompact(plain), after dropping null entries of tables at every depth (MDATK's
// _to_plain of a dict; every MDATK caller passes dicts). Arrays keep their null items.
std::string stableHash(const Value& plain, std::size_t length = 16);

// file_fingerprint(path): {"path": absolute path, "size": bytes, "mtime_ns": modification time in
// nanoseconds since the Unix epoch}. Throws std::filesystem::filesystem_error for a missing file.
Value fileFingerprint(const std::filesystem::path& path);

}  // namespace mdbridge::core
