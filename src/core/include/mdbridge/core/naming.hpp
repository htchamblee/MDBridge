#pragma once
// File and folder names of generated workflow files: MDATK's mdatk/workflows/naming.py.

#include <cstdint>
#include <optional>
#include <string>
#include <string_view>

#include "mdbridge/core/value.hpp"

namespace mdbridge::core {

// strftime format of workflow time stamps, e.g. 20261015-1432.
inline constexpr std::string_view kStampFormat = "%Y%m%d-%H%M";

// The stamp of `now` (a local date-time), or of the current local time.
std::string newStamp(const std::optional<DateTime>& now = std::nullopt);

// Whitespace stripped, every run of characters other than A-Z a-z 0-9 . _ - replaced by "-", then
// "-", "." and "_" stripped from both ends; "project" when nothing is left.
std::string sanitizeProject(std::string_view name);

// Like sanitizeProject with the characters A-Z a-z 0-9 - only (so "_" becomes "-"); "step" when
// nothing is left.
std::string sanitizeTag(std::string_view tag);

// At least two digits: 3 -> "03", 123 -> "123", -1 -> "-1".
std::string stepLabel(std::int64_t step);

// <project>_<step>_<tag>[_<stamp>][.<ext>], e.g. apt1_03_amber-build_20261015-1432.sh; leading
// dots of `ext` are dropped.
std::string fileName(std::string_view project, std::int64_t step, std::string_view tag,
                     std::string_view stamp = {}, std::string_view ext = {});

// <step>_<tag>, e.g. 03_amber-build.
std::string stepFolder(std::int64_t step, std::string_view tag);

struct ParsedName {
    std::string project, step, tag, stamp, ext;
    bool operator==(const ParsedName&) const = default;
};

// The parts of a name made by fileName, or std::nullopt.
std::optional<ParsedName> parseName(std::string_view name);

}  // namespace mdbridge::core
