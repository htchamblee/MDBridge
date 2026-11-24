// Amber ASCII coordinates as MDAnalysis 2.11.0 reads them: TRJReader (mdcrd: 10F8.3 lines, a
// box line after each frame when the line after the first frame has exactly three F8.3 fields)
// and INPReader (rst7 / inpcrd: one frame, time from the second line, never a box).

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

#include "binary_file.hpp"
#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/core/text.hpp"
#include "readers.hpp"
#include "text_file.hpp"

namespace mdbridge::io::detail {

namespace {

// CLEAN ROOM: reimplement from scratchpad/cleanroom/amber-ascii.md

}  // namespace

std::unique_ptr<core::Trajectory> openAmberTrajectory(const std::filesystem::path& path,
                                                      std::optional<double> dt, std::size_t atoms) {
    (void)path; (void)dt; (void)atoms;
    throw std::logic_error("clean room: openAmberTrajectory not implemented");
}

std::unique_ptr<core::Trajectory> openAmberRestart(const std::filesystem::path& path,
                                                   std::optional<double> dt, std::size_t) {
    (void)path; (void)dt;
    throw std::logic_error("clean room: openAmberRestart not implemented");
}

}  // namespace mdbridge::io::detail
