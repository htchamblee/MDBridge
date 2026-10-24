#pragma once
// Stable hashes of settings (cache keys) and cheap file identities, as MDATK's config.py computes
// them (stable_hash, file_fingerprint).

#include <cstddef>
#include <filesystem>
#include <string>
#include <string_view>

#include "mdbridge/core/value.hpp"

namespace mdbridge::core {

// SHA-256 of the bytes, as 64 lowercase hex digits.
std::string sha256Hex(std::string_view data);

// stable_hash(obj, length): the first `length` hex digits of the SHA-256 of
// pythonJsonCompact(plain), after dropping null entries of tables at every depth (MDATK's
// _to_plain of a dict; every MDATK caller passes dicts). Arrays keep their null items.
std::string stableHash(const Value& plain, std::size_t length = 16);

// file_fingerprint(path): {"path": absolute path, "size": bytes, "mtime_ns": modification time in
// nanoseconds since the Unix epoch}. Throws std::filesystem::filesystem_error for a missing file.
Value fileFingerprint(const std::filesystem::path& path);

}  // namespace mdbridge::core
