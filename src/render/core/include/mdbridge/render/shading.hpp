#pragma once
// Residue colouring and proximity glow of the structure view, shared with the sequence and
// secondary-structure panels so that all three show the same colours and glow (MDATK's
// mdatk/render/coloring.py and residues.py): glow intensities per mode, static colours per colour
// mode, per-frame proximity colours, one-letter codes, trace rows of residue ranges and the rows
// that get residue-number tags. Qt-free.
//
// Differences from MDATK (docs/DIFFERENCES.md): tagRows honours a chain in visual.label_residues
// ("B:12" tags residue 12 of chain B only; RN-24) and warns once about a spec it cannot read
// (MDATK: silently no tags from it; RN-23).

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "mdbridge/core/dataset.hpp"
#include "mdbridge/core/project.hpp"
#include "mdbridge/render/color.hpp"

namespace mdbridge::render {

inline constexpr std::array<std::string_view, 4> kGlowModes{"none", "ever", "fade", "current"};
inline constexpr std::array<std::string_view, 3> kColorModes{"chain", "groups", "proximity"};

// rows x frames of float32 (numpy's (rows, frames) float32 array), row-major.
struct FrameMatrix {
    std::size_t rows = 0;
    std::size_t frames = 0;
    std::vector<float> values;

    float at(std::size_t row, std::size_t frame) const { return values[row * frames + frame]; }
    // The column of one frame (every row).
    std::vector<float> column(std::size_t frame) const;
};

// glow_intensity(min_distance, member, cutoff, mode, fade_frames): per row and frame 0..1.
//   ever    1 for the whole movie if the residue ever comes within the cutoff
//   fade    1 while within the cutoff, then a linear fade over fadeFrames (fadeFrames <= 0: current)
//   current 1 only while within the cutoff
//   none    always 0
// Members of the proximity set never glow. minDistance: rows x frames, row-major; member: one per
// row (std::invalid_argument for other sizes). Throws core::ValueError("unknown glow mode 'x'
// (choose from none, ever, fade, current)").
FrameMatrix glowIntensity(std::span<const double> minDistance, std::size_t rows, std::size_t frames,
                          const std::vector<bool>& member, double cutoff, std::string_view mode,
                          std::int64_t fadeFrames = 0);

// static_colors: RGBA per residue that do not change over time. "chain": the chain palette's
// colours by sorted distinct chain index; "groups" / "proximity": the colour of groupIndex's group;
// residues without one (and every residue in other modes) get baseColor. Colour strings are
// Matplotlib's (toRgba; core::ValueError for invalid ones).
std::vector<Rgba> staticColors(std::string_view mode, std::span<const std::int64_t> chain,
                               std::span<const std::int64_t> groupIndex, const std::vector<std::string>& groupColors,
                               std::string_view baseColor, std::string_view chainPalette);

// frame_colors: in proximity mode non-members blend baseColor -> proximityColor by `intensity`
// (one per residue), members keep their static colour; other modes: the static colours.
std::vector<Rgba> frameColors(const std::vector<Rgba>& colors, std::string_view mode, std::string_view baseColor,
                              std::string_view proximityColor, std::span<const float> intensity,
                              const std::vector<bool>& member);

// warn_once: a warning (core::logWarning) once per process for each message: the GUI rebuilds its
// panels on every preview refresh. resetWarnOnce() forgets the messages (tests).
void warnOnce(const std::string& message);
void resetWarnOnce();

// one_letter: the base letter of a nucleotide in any naming scheme (io::parseResidueName: G5, RA,
// LG, 0A, DT ...), the amino-acid code (incl. HID / HIE / HIP / CYX), else "X".
std::string oneLetter(std::string_view resname);

// chain_names: the chain letters of the trace rows ("chain_name" when any is non-empty; else
// letters A..Z from the chain index, modulo 26).
std::vector<std::string> chainNames(const core::AnalysisResult& trace);
// resnames: "resname" of the trace, else names from meta["residue_labels"] ("A:G12" -> "G"; "X"
// for an empty one) when there is one per row, else "X" per row.
std::vector<std::string> residueNames(const core::AnalysisResult& trace);
// residue_rows: trace rows of residue ranges such as "12", "B:12" or "A:1-20, B:5-30" (row order);
// all rows for an empty spec. Throws core::SelectionError for a malformed spec and
// "no residue of the structure view matches '<spec>'" for one matching nothing.
std::vector<std::int64_t> residueRows(const core::AnalysisResult& trace, std::string_view spec);

// The colours and glow of the trace rows, exactly as the structure view computes them.
struct Shading {
    std::vector<Rgba> staticColors;  // colour mode
    std::vector<bool> member;        // residues of the proximity set (never glow)
    FrameMatrix glow;                // rows x frames
    FrameMatrix colorIntensity;      // rows x frames proximity-colour intensity
    std::array<double, 3> glowRgb{};
};

// residue_shading(ctx, trace, need_proximity): needProximity defaults to color_mode "proximity"
// or a glow mode other than "none". Proximity rows are matched to the trace by resindex (missing
// rows: infinitely far, not members). Throws std::out_of_range("proximity") when the proximity
// result is needed but missing.
Shading residueShading(const core::Project& project, const core::Dataset& dataset, const core::AnalysisResult& trace,
                       std::optional<bool> needProximity = std::nullopt);

// tag_rows: rows that get a residue-number tag: resid % label_every == 0 (label_every > 0) and the
// residues of label_residues (with a chain only on that chain, RN-24); a spec that cannot be read
// adds none and warns once "visual.label_residues: <reason>; no residue-number tags from it" (RN-23).
std::vector<std::int64_t> tagRows(const core::AnalysisResult& trace, const core::VisualParams& visual);

}  // namespace mdbridge::render
