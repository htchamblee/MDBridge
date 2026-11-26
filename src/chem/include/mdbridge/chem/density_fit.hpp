#pragma once
// The rigid-body density fit of MDATK's mdatk/density.py (the owner's code, ported; unchanged in
// MDATK 0.1.1): near-uniform rotations (super-Fibonacci spirals, Alexa 2022), the fit problem (the
// model's simulated envelope scored against the map by the correlation about the mean, ChimeraX
// fitmap's convention), its local refinement (Powell's method on six rigid-body parameters), and the
// global search (every rotation x FFT translations on a coarse grid scored by Roseman's locally
// normalised correlation, the best distinct poses refined against an 8 A low-pass copy of the map
// first, then the best three and the 180-degree flips of the best two at the map's resolution).
// Owner decision D2 (docs/specs/m6-plan.md): MDATK's algorithm and defaults (30 degrees, 12 poses,
// 60 000 envelope points, 8 A coarse stage), results within D10's tolerance of MDATK's (0.5 degrees,
// 0.2 A, CC 1e-3).
//
// Differences from MDATK v0.0.4.1 (docs/DIFFERENCES.md, "Density fit"):
// - CH-10: a non-finite resolution is refused ("the resolution must be a positive number of Å")
//   where MDATK failed in numpy ("negative dimensions are not allowed").
// - CH-11: non-finite map voxels carry no density (they are read as 0, so the smoothing of the
//   coarse search and of the coarse-to-fine copy cannot spread them), and fitModel's message ends
//   "; N non-finite map voxel(s) read as no density".
// - The map values sampled for a score are float32 (scipy's map_coordinates of a float32 map), but
//   every sum, mean and FFT is in double precision (MDATK: float32 means and BLAS sdot of the
//   centred samples, float32 FFTs of the coarse map), and sums run in a fixed order, so the result
//   does not depend on the thread count or the CPU's BLAS kernels.
// - The global search scores rotations and refines poses on several threads; the candidates, their
//   order and every number are those of one thread.

#include <array>
#include <atomic>
#include <complex>
#include <cstddef>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "mdbridge/core/progress.hpp"
#include "mdbridge/core/rigid_transform.hpp"
#include "mdbridge/io/mrc.hpp"  // DensityMap, FitError, kSigmaFactor

namespace mdbridge::chem {

using io::FitError;  // density.FitError

inline constexpr double kEnvelopeFraction = 0.05;  // simulated-map points scored: value >= fraction * max
inline constexpr double kInsideFull = 0.5;         // envelope fraction inside the map box: full score
inline constexpr double kInsideMin = 0.25;         // ... no score
inline constexpr double kCoarseResolution = 8.0;   // global-search candidates are refined at no finer first
inline constexpr double kMaxAngularStep = 45.0;    // coarser global searches miss the right orientation
inline constexpr std::size_t kMaxEnvelopePoints = 60000;
inline constexpr std::size_t kTopPoses = 12;
inline constexpr double kDefaultAngularStep = 30.0;

// rotation_samples(step): max(24, ceil(8 pi^2 / step^3)) rotations (step in radians) covering SO(3)
// about `stepDeg` apart. FitError "the angular step must be at least 3 degrees (got 2.0)" for a step
// below 3 or NaN (the value as Python's repr writes it).
std::vector<core::Mat3> rotationSamples(double stepDeg);

// _fft_size(n): the smallest 2^a 3^b 5^c >= n (n itself for n < 1).
std::size_t fftSize(std::size_t n);

// centre_transform: the placement without rotation that puts the model centre on the centre of the
// density above `level` (numerics::densityCentroid; non-finite voxels carry none).
core::RigidTransform centreTransform(const io::DensityMap& map, double level, const core::Vec3& modelCentre);

struct FitResult {
    core::RigidTransform transform;  // world -> map
    double cc = 0.0;
    std::string start;  // "manual", "global" or "identity"
    std::size_t candidates = 1;
    std::size_t evaluations = 0;
    std::string message;
};

// The work of a fit, shared by its steps: the reporter's progress ("global search: rotation 1/551",
// "refining pose 1/12", "final refinement", "global search done", "local fit done") and its cancel
// flag (checked before each rotation batch and each refinement: core::Cancelled), and the threads
// the global search may use (0: the hardware's, at most 16; 1: the calling thread only).
struct FitControl {
    const core::Reporter* reporter = nullptr;
    std::size_t threads = 0;
};

// FitProblem: the correlation of a model (atoms in world coordinates) with an experimental map.
// The problem keeps a pointer to `map` unless the map has non-finite voxels (then it keeps its own
// copy with them as 0, CH-11): the map must outlive the problem. Const member functions may be
// called from several threads at once.
class FitProblem {
public:
    // FitError "the resolution must be a positive number of Å" (NaN, infinite; CH-10), "enter the
    // map resolution (Å) for the simulated map" (zero, negative), "at least three atoms are needed
    // for a fit"; std::invalid_argument without one weight per atom.
    FitProblem(const io::DensityMap& map, std::span<const core::Vec3> xyz, std::span<const double> weights,
               double resolution, double sigmaFactor = io::kSigmaFactor, double envelope = kEnvelopeFraction,
               std::size_t maxPoints = kMaxEnvelopePoints);
    FitProblem(FitProblem&& other) noexcept;
    FitProblem(const FitProblem&) = delete;
    FitProblem& operator=(const FitProblem&) = delete;
    FitProblem& operator=(FitProblem&&) = delete;
    ~FitProblem();

    // cc: the correlation about the mean over the envelope points (or the given subset) inside the
    // map box after `transform`; poses with less than kInsideFull of them inside are scaled down, to
    // 0 at kInsideMin. Counts one evaluation.
    double cc(const core::RigidTransform& transform) const;
    double cc(const core::RigidTransform& transform, std::span<const core::Vec3> points,
              std::span<const double> values) const;
    double insideFraction(const core::RigidTransform& transform) const;
    // at_resolution: the same problem against a copy of the map low-pass filtered to `resolution`
    // (default envelope and point count, as MDATK builds it).
    FitProblem atResolution(double resolution) const;
    // 180-degree rotations about the model's principal axes through its centre (eigenvalues
    // ascending).
    std::vector<core::RigidTransform> flips() const;
    // refine: maximise the correlation from `start` (Powell, xtol, ftol 1e-7, maxfev 4000, initial
    // steps max(0.5, resolution / 4) for the rotation vector times rg and the shift); with maxPoints
    // the envelope thinned to that many points. Returns the pose and its cc over all points.
    std::pair<core::RigidTransform, double> refine(const core::RigidTransform& start,
                                                   std::optional<std::size_t> maxPoints = std::nullopt,
                                                   double xtol = 0.02,
                                                   std::optional<double> maxShift = std::nullopt) const;
    // global_search: FitError for a step below 3 degrees; "no placement of the model fits inside the
    // map at this resolution (the map is smaller than the coarse search grid)"; core::Cancelled.
    FitResult globalSearch(double angularStep = kDefaultAngularStep, std::size_t top = kTopPoses,
                           const FitControl& control = {}) const;
    // "; 60% of the model inside the map" below 99.5 %, else ""
    std::string insideNote(const core::RigidTransform& transform) const;

    const io::DensityMap& map() const { return *m_map; }
    double resolution() const { return m_resolution; }
    double sigmaFactor() const { return m_sigmaFactor; }
    std::span<const core::Vec3> xyz() const { return m_xyz; }
    std::span<const double> weights() const { return m_weights; }
    const core::Vec3& centre() const { return m_centre; }
    double radius() const { return m_radius; }
    double rg() const { return m_rg; }
    const io::DensityMap& simulated() const { return m_sim; }
    std::span<const core::Vec3> points() const { return m_points; }  // world coordinates
    std::span<const double> values() const { return m_values; }
    std::size_t evaluations() const { return m_evaluations.load(); }
    // the map's non-finite voxels, read as 0 (CH-11)
    std::size_t nonFiniteVoxels() const { return m_nonFinite; }

private:
    // _subset(n): the envelope thinned evenly to n points (all of them when there are no more)
    void subset(std::size_t n, std::vector<core::Vec3>& points, std::vector<double>& values) const;

    std::unique_ptr<io::DensityMap> m_owned;  // the map with non-finite voxels as 0, when it had any
    const io::DensityMap* m_map = nullptr;
    std::size_t m_nonFinite = 0;
    double m_resolution = 0.0;
    double m_sigmaFactor = io::kSigmaFactor;
    std::vector<core::Vec3> m_xyz;
    std::vector<double> m_weights;
    core::Vec3 m_centre{};
    double m_radius = 0.0;
    double m_rg = 1.0;
    io::DensityMap m_sim;
    std::vector<core::Vec3> m_points;
    std::vector<double> m_values;
    mutable std::atomic<std::size_t> m_evaluations{0};
};

// _CoarseSearch: the FFT translation search for one rotation at a time on a coarse, padded grid
// (spacing about resolution / 2, at least the map voxel, at most ~96 points per axis, 2^a 3^b 5^c
// points per axis).
class CoarseSearch {
public:
    CoarseSearch(const FitProblem& problem, double angularStep);

    struct Peak {
        double score = 0.0;
        core::Vec3 centre{};  // the model centre's position in map coordinates
    };
    // peaks(r, count): the best `count` translations for rotation r by the locally normalised
    // correlation, each peak's neighbourhood (5^3 points) suppressed before the next one.
    std::vector<Peak> peaks(const core::Mat3& rotation, std::size_t count = 1) const;

    double spacing() const { return m_h; }
    const core::Vec3& lowCorner() const { return m_lo; }
    const std::array<std::size_t, 3>& gridXyz() const { return m_n; }  // nx, ny, nz
    double varFloor() const { return m_varFloor; }

private:
    double m_h = 0.0;
    core::Vec3 m_lo{};
    std::array<std::size_t, 3> m_n{};  // x, y, z
    std::vector<std::complex<double>> m_eFt;
    std::vector<std::complex<double>> m_e2Ft;
    double m_varFloor = 0.0;
    std::vector<unsigned char> m_valid;  // [z, y, x]
    std::vector<double> m_gauss;         // [z, y, x / 2 + 1]
    std::vector<core::Vec3> m_rel;
    std::vector<double> m_w;
};

enum class FitSearch { Global, Manual };

// fit_model: fit atoms (world coordinates) into `map`. Global runs the exhaustive search; Manual
// refines locally from `start` (identity when not given: start "identity"). Messages: "551
// rotations (30° step), 12 poses at 4 sites refined", "start CC 0.320", each with insideNote and the
// CH-11 note.
FitResult fitModel(const io::DensityMap& map, std::span<const core::Vec3> xyz, std::span<const std::string> elements,
                   double resolution, const std::optional<core::RigidTransform>& start = std::nullopt,
                   FitSearch search = FitSearch::Global, double angularStep = kDefaultAngularStep,
                   const FitControl& control = {});

}  // namespace mdbridge::chem
