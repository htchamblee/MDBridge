#pragma once
// The analysis cache (mdatk/io/cache.py CacheStore): one stored result per analysis instance and
// cache key, and a manifest describing the latest run.
//
// The cache is MDBridge's own (ARCHITECTURE "Compatibility": MDATK's cache is not read, results are
// recomputed); its keys are MDATK's (registry.hpp cacheKey). Files, in the cache folder
// (core::cacheDir, by default <out_dir>/MDATK_cache, which MDATK's files may share):
//
//   MDBridge_cache.json              the manifest (writeManifest), MDATK's schema with the version
//                                    field "mdbridge_version" and project "MDBridge"
//   <name>-<key>.mdbridge.npz        the result's arrays: a numpy .npz archive (ZIP, every member
//                                    deflated with zlib, ZIP64 where sizes need it) of .npy files
//                                    (format 1.0, or 2.0 for headers over 64 KiB; little-endian;
//                                    numpy's dtypes |b1 |i1 <i4 <i8 <f4 <f8 <U<n>), in the result's
//                                    order, then "_frames" (int64, the analysed frame numbers) and
//                                    "_time_ps" (float64, their times) - np.load reads it
//   <name>-<key>.mdbridge.json       {"project", "mdbridge_version", "created", "name", "key",
//                                    "params", "meta", "arrays" (sorted names), "n_frames"}
//
// MDATK's own files (MDATK_cache.json, <name>-<key>.npz / .json) have other names and are never
// read or written. JSON is written as Python's json.dumps(indent=2, ensure_ascii=False) writes it
// (core::toJson), atomically (<file>.tmp, then rename).

#include <cstdint>
#include <filesystem>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "mdbridge/core/dataset.hpp"
#include "mdbridge/core/value.hpp"

namespace mdbridge::analysis {

inline constexpr std::string_view kManifestName = "MDBridge_cache.json";
// Inserted before ".npz" / ".json" so that MDBridge's entries never collide with MDATK's.
inline constexpr std::string_view kEntryTag = ".mdbridge";

// A numpy .npz archive of named arrays (written deflated, as np.savez_compressed). Throws
// std::runtime_error for files it cannot write or read ("<path>: ..."); readNpz reads archives
// written by writeNpz and by numpy (stored or deflated members, ZIP64) whose arrays have the
// dtypes above (C order).
void writeNpz(const std::filesystem::path& path, const core::NamedArrays& arrays);
core::NamedArrays readNpz(const std::filesystem::path& path);

struct CachedResult {
    core::AnalysisResult result;
    std::vector<std::int64_t> frames;
    std::vector<double> timePs;
};

// What a manifest records about one stored result: entries maps instance -> key.
using CacheEntries = std::vector<std::pair<std::string, std::string>>;

class CacheStore {
public:
    // `cacheDir`: UTF-8; made absolute against the current directory (generic form on Windows).
    explicit CacheStore(std::string cacheDir);

    const std::string& cacheDir() const { return m_dir; }
    std::string manifestPath() const;
    // <cache>/<name>-<key>.mdbridge (add ".npz" / ".json")
    std::string entryBase(std::string_view name, std::string_view key) const;
    // The manifest's "file" of an entry: "<name>-<key>.mdbridge.npz".
    static std::string entryFile(std::string_view name, std::string_view key);
    void ensure() const;

    // Both files of the entry exist.
    bool has(std::string_view name, std::string_view key) const;
    // save_result: the arrays plus _frames / _time_ps, and the description; returns the .npz path.
    std::string saveResult(const core::AnalysisResult& result, std::span<const std::int64_t> frames,
                           std::span<const double> timePs) const;
    // load_result: std::nullopt when the entry is missing or unreadable (logged as a warning
    // "ignoring unreadable cache entry <base>: <reason>").
    std::optional<CachedResult> loadResult(std::string_view name, std::string_view key) const;

    // The manifest, or std::nullopt when there is none. Throws core::ProjectError for invalid JSON.
    std::optional<core::Value> readManifest() const;
    // write_manifest: {project, mdbridge_version, updated (core::clockNow), system, frames
    // (frameSelection + n_frames, first_frame, last_frame), time_axis {unit "ps", n, start, end,
    // dt = median of the steps}, groups, available (dataset.describe()), entries [{name, key, file,
    // params, outputs (per kind, sorted names)}] of `entries`, all_entries} then `extra`.
    // all_entries: the previous manifest's all_entries plus these, by name and key (an unreadable
    // previous manifest counts as none; MDATK fails there). Returns the manifest's path.
    std::string writeManifest(const core::Value& system, const core::Table& frameSelection,
                              const core::Dataset& dataset, const core::Array& groups, const CacheEntries& entries,
                              const core::Table& extra = {}) const;
    // load_dataset_from_manifest: the results the manifest's entries name, frames and times of the
    // first. Throws core::NamedError("FileNotFoundError", ...): "no MDBridge_cache.json in <dir>",
    // "cache entry <file> is missing", "the cache manifest lists no analyses".
    core::Dataset loadDatasetFromManifest(const std::string& timeUnit = "ns") const;

private:
    std::string m_dir;
};

}  // namespace mdbridge::analysis
