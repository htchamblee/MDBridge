// The rigid-body density fit of mdatk/density.py (rotation_samples, FitProblem, _CoarseSearch,
// _fft_size, centre_transform, fit_model), ported from MDATK's own code (the owner's; plan D2).
// SciPy's pieces come from src/numerics (map_coordinates, gaussian_filter, rfftn / irfftn, Powell).

#include "mdbridge/chem/density_fit.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <thread>

#include "mdbridge/core/detail/numpy_sum.hpp"
#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/json.hpp"  // pythonFloatRepr
#include "mdbridge/core/linalg.hpp"
#include "mdbridge/numerics/density.hpp"
#include "mdbridge/numerics/fft.hpp"
#include "mdbridge/numerics/ndimage.hpp"
#include "mdbridge/numerics/powell.hpp"

namespace mdbridge::chem {

using core::Mat3;
using core::RigidTransform;
using core::Vec3;
using Complex = std::complex<double>;

namespace {

constexpr std::size_t kMaxDefaultThreads = 16;
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

std::size_t resolveThreads(std::size_t requested) {
    if (requested > 0) {
        return requested;
    }
    const unsigned hardware = std::thread::hardware_concurrency();
    return hardware == 0 ? 1 : std::min<std::size_t>(hardware, kMaxDefaultThreads);
}

// fn(i) for every i in [begin, end) on up to `threads` threads (each index exactly once, so results
// written per index do not depend on the thread count); the first exception is rethrown.
template <class F>
void parallelRange(std::size_t begin, std::size_t end, std::size_t threads, const F& fn) {
    const std::size_t count = end > begin ? end - begin : 0;
    threads = std::min(threads, count);
    if (threads <= 1) {
        for (std::size_t i = begin; i < end; ++i) {
            fn(i);
        }
        return;
    }
    std::atomic<std::size_t> next{begin};
    std::vector<std::exception_ptr> errors(threads);
    std::vector<std::thread> pool;
    pool.reserve(threads);
    for (std::size_t t = 0; t < threads; ++t) {
        pool.emplace_back([&, t]() {
            try {
                for (;;) {
                    const std::size_t i = next.fetch_add(1);
                    if (i >= end) {
                        break;
                    }
                    fn(i);
                }
            } catch (...) {
                errors[t] = std::current_exception();
                next.store(end);
            }
        });
    }
    for (std::thread& th : pool) {
        th.join();
    }
    for (const std::exception_ptr& e : errors) {
        if (e) {
            std::rethrow_exception(e);
        }
    }
}

// work(i) for i in [0, count) in batches of `batch` (parallel within a batch), cancellation checked
// before each batch, then after(i) in index order on the calling thread (progress)
template <class Work, class After>
void runBatches(std::size_t count, std::size_t batch, std::size_t threads, const core::Reporter* reporter,
                const Work& work, const After& after) {
    batch = std::max<std::size_t>(1, batch);
    for (std::size_t first = 0; first < count; first += batch) {
        if (reporter) {
            reporter->checkCancel();
        }
        const std::size_t last = std::min(count, first + batch);
        parallelRange(first, last, threads, work);
        for (std::size_t i = first; i < last; ++i) {
            after(i);
        }
    }
}

void progress(const FitControl& control, double fraction, const std::string& message) {
    if (control.reporter) {
        control.reporter->progress(fraction, message);
    }
}

// numpy's sum of n values (pairwise)
template <class F>
double sum(std::size_t n, const F& value) {
    return core::detail::numpySum<double>(value, n);
}

std::string format(const char* spec, double value) {
    char buffer[64];
    std::snprintf(buffer, sizeof buffer, spec, value);
    return buffer;
}

// np.linspace(0, n - 1, count).astype(np.int64)
std::vector<std::size_t> linspaceIndices(std::size_t n, std::size_t count) {
    std::vector<std::size_t> out(count);
    if (count == 0) {
        return out;
    }
    const double stop = static_cast<double>(n - 1);
    if (count == 1) {
        out[0] = 0;
        return out;
    }
    const double step = stop / static_cast<double>(count - 1);
    for (std::size_t i = 0; i < count; ++i) {
        out[i] = static_cast<std::size_t>(static_cast<double>(i) * step);
    }
    out[count - 1] = n - 1;
    return out;
}

numerics::Shape3 shapeZyx(const io::DensityMap& map) {
    return {static_cast<std::size_t>(map.shapeXyz[2]), static_cast<std::size_t>(map.shapeXyz[1]),
            static_cast<std::size_t>(map.shapeXyz[0])};
}

// sample_map(map, [p], cval=NaN, order=1): the trilinear value (float32, as scipy returns it for a
// float32 map), NaN outside the grid
double sample(const io::DensityMap& map, const numerics::Shape3& shape, const Vec3& p) {
    const Vec3 ijk = map.xyzToIndex(p);
    const double v = numerics::mapCoordinatesLinearConstant(map.data, shape, ijk[2], ijk[1], ijk[0], kNaN);
    return static_cast<double>(static_cast<float>(v));
}

bool poseClose(const RigidTransform& a, const RigidTransform& b, const Vec3& centre, double shift, double angle) {
    const Vec3 pa = a.apply(centre), pb = b.apply(centre);
    const double d = std::sqrt((pa[0] - pb[0]) * (pa[0] - pb[0]) + (pa[1] - pb[1]) * (pa[1] - pb[1]) +
                               (pa[2] - pb[2]) * (pa[2] - pb[2]));
    if (!(d < shift)) {
        return false;
    }
    Mat3 rel{};
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            rel[i][j] = a.rotation[i][0] * b.rotation[j][0] + a.rotation[i][1] * b.rotation[j][1] +
                        a.rotation[i][2] * b.rotation[j][2];
        }
    }
    return core::rotationAngle(rel) < angle;
}

// np.fft.fftfreq(n, d) (rfft: the first n / 2 + 1)
double fftFrequency(std::size_t i, std::size_t n, double d) {
    const double val = 1.0 / (static_cast<double>(n) * d);
    const std::size_t half = (n - 1) / 2 + 1;
    const double index = i < half ? static_cast<double>(i) : -static_cast<double>(n - i);
    return index * val;
}

// conj(a) * b, as numpy multiplies complex numbers
Complex conjTimes(const Complex& a, const Complex& b) {
    return {a.real() * b.real() + a.imag() * b.imag(), a.real() * b.imag() - a.imag() * b.real()};
}

}  // namespace

// ---- rotations, grid sizes, centring --------------------------------------------------------

std::vector<Mat3> rotationSamples(double stepDeg) {
    if (!(stepDeg >= 3.0)) {
        throw FitError("the angular step must be at least 3 degrees (got " + core::pythonFloatRepr(stepDeg) + ")");
    }
    const double pi = std::numbers::pi;
    const double step = stepDeg * (pi / 180.0);  // math.radians
    const double estimate = std::ceil(8 * std::pow(pi, 2.0) / std::pow(step, 3.0));
    const std::size_t n = std::max<std::size_t>(24, static_cast<std::size_t>(estimate));
    const double phi = std::sqrt(2.0);
    const double psi = 1.533751168755204288118041;
    const double twoPi = 2 * pi;
    std::vector<Mat3> out;
    out.reserve(n);
    const double nd = static_cast<double>(n);
    for (std::size_t k = 0; k < n; ++k) {
        const double s = static_cast<double>(k) + 0.5;
        const double r = std::sqrt(s / nd);
        const double big = std::sqrt(1.0 - s / nd);
        const double alpha = twoPi * s / phi;
        const double beta = twoPi * s / psi;
        out.push_back(core::quaternionMatrix(
            {r * std::sin(alpha), r * std::cos(alpha), big * std::sin(beta), big * std::cos(beta)}));
    }
    return out;
}

std::size_t fftSize(std::size_t n) {
    std::size_t best = 0;
    for (std::size_t p2 = 1; p2 < 2 * n; p2 *= 2) {
        for (std::size_t p3 = p2; p3 < 2 * n; p3 *= 3) {
            for (std::size_t p5 = p3; p5 < 2 * n; p5 *= 5) {
                if (p5 >= n && (best == 0 || p5 < best)) {
                    best = p5;
                }
            }
        }
    }
    return best != 0 ? best : n;
}

RigidTransform centreTransform(const io::DensityMap& map, double level, const Vec3& modelCentre) {
    const Vec3 c = numerics::densityCentroid(map, level);
    RigidTransform out;
    out.translation = {c[0] - modelCentre[0], c[1] - modelCentre[1], c[2] - modelCentre[2]};
    return out;
}

// ---- FitProblem --------------------------------------------------------------------------------

FitProblem::FitProblem(const io::DensityMap& densityMap, std::span<const Vec3> atomXyz,
                       std::span<const double> atomWeights, double res, double sigma, double envelope,
                       std::size_t maxPoints)
    : m_resolution(res), m_sigmaFactor(sigma) {
    if (std::isnan(res) || std::isinf(res)) {  // CH-10
        throw FitError("the resolution must be a positive number of Å");
    }
    if (res <= 0) {
        throw FitError("enter the map resolution (Å) for the simulated map");
    }
    if (atomXyz.size() < 3) {
        throw FitError("at least three atoms are needed for a fit");
    }
    if (atomWeights.size() != atomXyz.size()) {
        throw std::invalid_argument("FitProblem: one weight per atom is needed");
    }
    // CH-11: non-finite voxels carry no density
    for (const float v : densityMap.data) {
        if (!std::isfinite(v)) {
            ++m_nonFinite;
        }
    }
    if (m_nonFinite > 0) {
        m_owned = std::make_unique<io::DensityMap>(densityMap);
        for (float& v : m_owned->data) {
            v = std::isfinite(v) ? v : 0.0f;
        }
        m_map = m_owned.get();
    } else {
        m_map = &densityMap;
    }
    m_xyz.assign(atomXyz.begin(), atomXyz.end());
    m_weights.assign(atomWeights.begin(), atomWeights.end());
    const std::size_t n = m_xyz.size();
    // xyz.mean(axis=0): the rows added in order
    Vec3 total{0, 0, 0};
    for (const Vec3& p : m_xyz) {
        for (std::size_t k = 0; k < 3; ++k) {
            total[k] += p[k];
        }
    }
    for (std::size_t k = 0; k < 3; ++k) {
        m_centre[k] = total[k] / static_cast<double>(n);
    }
    std::vector<double> d2(n);
    for (std::size_t i = 0; i < n; ++i) {
        const double dx = m_xyz[i][0] - m_centre[0], dy = m_xyz[i][1] - m_centre[1], dz = m_xyz[i][2] - m_centre[2];
        d2[i] = dx * dx + dy * dy + dz * dz;
        const double r = std::sqrt(d2[i]);
        if (i == 0 || r > m_radius) {  // the largest distance (np.max)
            m_radius = r;
        }
    }
    m_rg = std::max(std::sqrt(sum(n, [&](std::size_t i) { return d2[i]; }) / static_cast<double>(n)), 1.0);
    const double minVoxel = std::min({m_map->voxel[0], m_map->voxel[1], m_map->voxel[2]});
    const double spacing = std::max(res / 3.0, minVoxel * 0.5);
    m_sim = io::simulateMap(m_xyz, m_weights, res, spacing, sigma);
    double vmax = -std::numeric_limits<double>::infinity();
    for (const float v : m_sim.data) {
        vmax = std::max(vmax, static_cast<double>(v));
    }
    const double threshold = envelope * vmax;
    std::vector<std::size_t> keep;
    for (std::size_t i = 0; i < m_sim.data.size(); ++i) {
        if (static_cast<double>(m_sim.data[i]) >= threshold) {
            keep.push_back(i);
        }
    }
    if (keep.size() > maxPoints) {  // thin evenly; the envelope is smooth
        std::vector<std::size_t> thinned;
        thinned.reserve(maxPoints);
        for (const std::size_t j : linspaceIndices(keep.size(), maxPoints)) {
            thinned.push_back(keep[j]);
        }
        keep.swap(thinned);
    }
    const auto ny = static_cast<std::size_t>(m_sim.shapeXyz[1]), nx = static_cast<std::size_t>(m_sim.shapeXyz[0]);
    m_points.reserve(keep.size());
    m_values.reserve(keep.size());
    for (const std::size_t i : keep) {
        const std::size_t z = i / (ny * nx), rem = i % (ny * nx);
        const std::size_t y = rem / nx, x = rem % nx;
        m_points.push_back(
            m_sim.indexToXyz({static_cast<double>(x), static_cast<double>(y), static_cast<double>(z)}));
        m_values.push_back(static_cast<double>(m_sim.data[i]));
    }
}

FitProblem::FitProblem(FitProblem&& other) noexcept
    : m_owned(std::move(other.m_owned)),
      m_map(other.m_map),
      m_nonFinite(other.m_nonFinite),
      m_resolution(other.m_resolution),
      m_sigmaFactor(other.m_sigmaFactor),
      m_xyz(std::move(other.m_xyz)),
      m_weights(std::move(other.m_weights)),
      m_centre(other.m_centre),
      m_radius(other.m_radius),
      m_rg(other.m_rg),
      m_sim(std::move(other.m_sim)),
      m_points(std::move(other.m_points)),
      m_values(std::move(other.m_values)),
      m_evaluations(other.m_evaluations.load()) {}

FitProblem::~FitProblem() = default;

double FitProblem::cc(const RigidTransform& transform) const {
    return cc(transform, m_points, m_values);
}

double FitProblem::cc(const RigidTransform& transform, std::span<const Vec3> pts,
                      std::span<const double> vals) const {
    ++m_evaluations;
    const numerics::Shape3 shape = shapeZyx(*m_map);
    const std::size_t n = pts.size();
    std::vector<double> e, s;
    e.reserve(n);
    s.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        const double v = sample(*m_map, shape, transform.apply(pts[i]));
        if (std::isfinite(v)) {
            e.push_back(v);
            s.push_back(vals[i]);
        }
    }
    const double frac = n > 0 ? static_cast<double>(e.size()) / static_cast<double>(n) : 0.0;
    if (frac <= kInsideMin) {
        return 0.0;
    }
    const std::size_t m = e.size();
    const double md = static_cast<double>(m);
    const double sMean = sum(m, [&](std::size_t i) { return s[i]; }) / md;
    const double eMean = sum(m, [&](std::size_t i) { return e[i]; }) / md;
    for (std::size_t i = 0; i < m; ++i) {
        s[i] -= sMean;
        e[i] -= eMean;
    }
    const double bNorm = std::sqrt(sum(m, [&](std::size_t i) { return s[i] * s[i]; }));
    const double aNorm = std::sqrt(sum(m, [&](std::size_t i) { return e[i] * e[i]; }));
    double out = 0.0;
    if (!(aNorm <= 0 || bNorm <= 0)) {
        out = sum(m, [&](std::size_t i) { return e[i] * s[i]; }) / (aNorm * bNorm);
    }
    if (frac < kInsideFull) {
        out *= (frac - kInsideMin) / (kInsideFull - kInsideMin);
    }
    return out;
}

double FitProblem::insideFraction(const RigidTransform& transform) const {
    if (m_points.empty()) {
        return 0.0;
    }
    const numerics::Shape3 shape = shapeZyx(*m_map);
    std::size_t inside = 0;
    for (const Vec3& p : m_points) {
        if (std::isfinite(sample(*m_map, shape, transform.apply(p)))) {
            ++inside;
        }
    }
    return static_cast<double>(inside) / static_cast<double>(m_points.size());
}

std::string FitProblem::insideNote(const RigidTransform& transform) const {
    const double frac = insideFraction(transform);
    return frac < 0.995 ? "; " + format("%.0f%%", frac * 100.0) + " of the model inside the map" : "";
}

FitProblem FitProblem::atResolution(double target) const {
    const io::DensityMap& source = *m_map;
    const double extra =
        m_sigmaFactor * std::sqrt(std::max(target * target - m_resolution * m_resolution, 0.0));
    auto filtered = std::make_unique<io::DensityMap>();
    filtered->shapeXyz = source.shapeXyz;
    filtered->origin = source.origin;
    filtered->voxel = source.voxel;
    filtered->path = source.path;
    filtered->data = source.data;
    const std::array<double, 3> sigma{extra / source.voxel[2], extra / source.voxel[1], extra / source.voxel[0]};
    numerics::gaussianFilter(filtered->data, shapeZyx(source), sigma, numerics::ExtendMode::Reflect);
    FitProblem out(*filtered, m_xyz, m_weights, target, m_sigmaFactor);
    out.m_owned = std::move(filtered);  // the filtered map is finite: the problem pointed at it
    out.m_map = out.m_owned.get();
    return out;
}

std::vector<RigidTransform> FitProblem::flips() const {
    Mat3 inertia{};
    for (const Vec3& p : m_xyz) {
        const Vec3 rel{p[0] - m_centre[0], p[1] - m_centre[1], p[2] - m_centre[2]};
        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t j = 0; j < 3; ++j) {
                inertia[i][j] += rel[i] * rel[j];
            }
        }
    }
    // np.linalg.eigh of the symmetric positive semi-definite matrix: its singular vectors, the
    // eigenvalues ascending (the sign of an axis does not change its flip)
    const core::Svd3 svd = core::svd3x3(inertia);
    std::vector<RigidTransform> out;
    for (std::size_t k = 0; k < 3; ++k) {
        const Vec3& a = svd.vt[2 - k];
        RigidTransform t;
        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t j = 0; j < 3; ++j) {
                t.rotation[i][j] = 2.0 * a[i] * a[j] - (i == j ? 1.0 : 0.0);
            }
        }
        for (std::size_t i = 0; i < 3; ++i) {
            const double rc = t.rotation[i][0] * m_centre[0] + t.rotation[i][1] * m_centre[1] +
                              t.rotation[i][2] * m_centre[2];
            t.translation[i] = m_centre[i] - rc;
        }
        out.push_back(t);
    }
    return out;
}

void FitProblem::subset(std::size_t n, std::vector<Vec3>& points, std::vector<double>& values) const {
    if (m_points.size() <= n) {
        points = m_points;
        values = m_values;
        return;
    }
    points.clear();
    values.clear();
    for (const std::size_t i : linspaceIndices(m_points.size(), n)) {
        points.push_back(m_points[i]);
        values.push_back(m_values[i]);
    }
}

std::pair<RigidTransform, double> FitProblem::refine(const RigidTransform& start, std::optional<std::size_t> maxPoints,
                                                     double xtol, std::optional<double> maxShift) const {
    std::vector<Vec3> pts;
    std::vector<double> sim;
    if (maxPoints && *maxPoints > 0) {
        subset(*maxPoints, pts, sim);
    } else {
        pts = m_points;
        sim = m_values;
    }
    const Vec3 cMap = start.apply(m_centre);
    const double rg = m_rg;
    const auto make = [&](std::span<const double> p) {
        RigidTransform rot;
        rot.rotation = core::rotvecMatrix({p[0] / rg, p[1] / rg, p[2] / rg});
        rot.translation = {p[3], p[4], p[5]};
        return rot.about(cMap).compose(start);
    };
    const auto cost = [&](std::span<const double> p) {
        if (maxShift && std::sqrt(p[3] * p[3] + p[4] * p[4] + p[5] * p[5]) > *maxShift) {
            return 1.0;
        }
        return -cc(make(p), pts, sim);
    };
    const double step = std::max(0.5, 0.25 * m_resolution);
    numerics::PowellOptions options;
    options.xtol = xtol;
    options.ftol = 1e-7;
    options.maxfev = 4000;
    options.direc.assign(6, std::vector<double>(6, 0.0));
    for (std::size_t i = 0; i < 6; ++i) {
        options.direc[i][i] = step;
    }
    const std::vector<double> x0(6, 0.0);
    const numerics::PowellResult res = numerics::minimizePowell(cost, x0, options);
    const RigidTransform best = make(res.x);
    return {best, cc(best)};
}

FitResult FitProblem::globalSearch(double angularStep, std::size_t top, const FitControl& control) const {
    const std::vector<Mat3> rots = rotationSamples(angularStep);
    const CoarseSearch coarse(*this, angularStep);
    const std::size_t threads = resolveThreads(control.threads);
    const std::size_t n = rots.size();
    // every rotation's two best translations; MDATK's progress after rotations 1, 17, 33, ... and the last
    std::vector<std::vector<CoarseSearch::Peak>> peaks(n);
    const std::size_t rotationBatch = 16 * ((threads + 15) / 16);
    runBatches(
        n, rotationBatch, threads, control.reporter, [&](std::size_t k) { peaks[k] = coarse.peaks(rots[k], 2); },
        [&](std::size_t k) {
            if (k % 16 == 0 || k == n - 1) {
                progress(control, 0.7 * static_cast<double>(k + 1) / static_cast<double>(n),
                         "global search: rotation " + std::to_string(k + 1) + "/" + std::to_string(n));
            }
        });
    struct Candidate {
        double score;
        std::size_t rotation;
        Vec3 centre;
    };
    std::vector<Candidate> cands;
    for (std::size_t k = 0; k < n; ++k) {
        for (const CoarseSearch::Peak& p : peaks[k]) {
            cands.push_back({p.score, k, p.centre});
        }
    }
    if (cands.empty()) {
        throw FitError(
            "no placement of the model fits inside the map at this resolution "
            "(the map is smaller than the coarse search grid)");
    }
    std::stable_sort(cands.begin(), cands.end(),
                     [](const Candidate& a, const Candidate& b) { return -a.score < -b.score; });
    std::vector<RigidTransform> chosen;
    const double minShift = std::max(2.0, 0.5 * m_resolution);
    const double siteRadius = std::max(4.0, m_resolution);
    struct Site {
        Vec3 centre;
        int poses;
    };
    std::vector<Site> sites;
    for (const Candidate& c : cands) {
        const Mat3& r = rots[c.rotation];
        RigidTransform t;
        t.rotation = r;
        for (std::size_t i = 0; i < 3; ++i) {
            t.translation[i] = c.centre[i] - (r[i][0] * m_centre[0] + r[i][1] * m_centre[1] + r[i][2] * m_centre[2]);
        }
        bool close = false;
        for (const RigidTransform& other : chosen) {
            if (poseClose(t, other, m_centre, minShift, angularStep * 0.75)) {
                close = true;
                break;
            }
        }
        if (close) {
            continue;
        }
        // at most 3 orientations per site, so the refined poses cover several sites
        Site* near = nullptr;
        for (Site& site : sites) {
            const double dx = site.centre[0] - c.centre[0], dy = site.centre[1] - c.centre[1],
                         dz = site.centre[2] - c.centre[2];
            if (std::sqrt(dx * dx + dy * dy + dz * dz) < siteRadius) {
                near = &site;
                break;
            }
        }
        if (near && near->poses >= 3) {
            continue;
        }
        if (near) {
            ++near->poses;
        } else {
            sites.push_back({c.centre, 1});
        }
        chosen.push_back(t);
        if (chosen.size() >= top) {
            break;
        }
    }
    // coarse-to-fine: candidates first meet a low-pass copy of the map (wider basins), then the best
    // ones and their 180-degree flips are refined at the real resolution
    std::optional<FitProblem> lowPass;
    if (m_resolution < kCoarseResolution) {
        lowPass.emplace(atResolution(kCoarseResolution));
    }
    const FitProblem& stage = lowPass ? *lowPass : *this;
    std::vector<std::pair<double, RigidTransform>> results(chosen.size());
    runBatches(
        chosen.size(), threads, threads, control.reporter,
        [&](std::size_t k) {
            auto [fit, value] = stage.refine(chosen[k], 4000, 0.1);
            results[k] = {value, fit};
        },
        [&](std::size_t k) {
            progress(control, 0.7 + 0.2 * static_cast<double>(k + 1) / static_cast<double>(chosen.size()),
                     "refining pose " + std::to_string(k + 1) + "/" + std::to_string(chosen.size()));
        });
    const auto byCc = [](const std::pair<double, RigidTransform>& a, const std::pair<double, RigidTransform>& b) {
        return -a.first < -b.first;
    };
    std::stable_sort(results.begin(), results.end(), byCc);
    std::vector<RigidTransform> tries;
    for (std::size_t k = 0; k < std::min<std::size_t>(3, results.size()); ++k) {
        tries.push_back(results[k].second);
    }
    const std::vector<RigidTransform> flipped = flips();
    for (std::size_t k = 0; k < std::min<std::size_t>(2, results.size()); ++k) {
        for (const RigidTransform& f : flipped) {
            tries.push_back(results[k].second.compose(f));
        }
    }
    std::vector<std::pair<double, RigidTransform>> finals(tries.size());
    runBatches(
        tries.size(), threads, threads, control.reporter,
        [&](std::size_t k) {
            RigidTransform fit = tries[k];
            if (lowPass) {
                fit = stage.refine(fit, 4000, 0.1).first;
            }
            auto [fit2, cc2] = refine(fit);
            finals[k] = {cc2, fit2};
        },
        [&](std::size_t k) {
            progress(control, 0.9 + 0.1 * static_cast<double>(k + 1) / static_cast<double>(tries.size()),
                     "final refinement");
        });
    std::stable_sort(finals.begin(), finals.end(), byCc);
    progress(control, 1.0, "global search done");
    FitResult out;
    out.cc = finals[0].first;
    out.transform = finals[0].second;
    out.start = "global";
    out.candidates = chosen.size();
    out.evaluations = evaluations();
    out.message = std::to_string(n) + " rotations (" + format("%g", angularStep) + "° step), " +
                  std::to_string(chosen.size()) + " poses at " + std::to_string(sites.size()) + " sites refined" +
                  insideNote(out.transform);
    return out;
}

// ---- CoarseSearch -------------------------------------------------------------------------------

CoarseSearch::CoarseSearch(const FitProblem& problem, double /*angularStep*/) {
    const io::DensityMap& map = problem.map();
    const double res = problem.resolution();
    const double sigma = problem.sigmaFactor() * res;
    Vec3 extent{};
    for (std::size_t k = 0; k < 3; ++k) {
        extent[k] = static_cast<double>(map.shapeXyz[k]) * map.voxel[k];
    }
    const double radius = problem.radius();
    // coarse spacing: about resolution / 2, at least the map voxel, and <= ~96 points per axis
    const double maxVoxel = std::max({map.voxel[0], map.voxel[1], map.voxel[2]});
    const double maxExtent = std::max({extent[0] + 2 * radius, extent[1] + 2 * radius, extent[2] + 2 * radius});
    m_h = std::max({res / 2.0, maxVoxel, maxExtent / 96.0});
    const double pad = radius + 3 * sigma;
    for (std::size_t k = 0; k < 3; ++k) {
        m_lo[k] = map.origin[k] - pad;
        m_n[k] = fftSize(static_cast<std::size_t>(std::ceil((extent[k] + 2 * pad) / m_h)) + 1);
    }
    const std::size_t nx = m_n[0], ny = m_n[1], nz = m_n[2];
    const numerics::FftShape3 shape{nz, ny, nx};
    // the experimental map resampled onto the coarse grid after anti-alias smoothing
    io::DensityMap smooth;
    smooth.shapeXyz = map.shapeXyz;
    smooth.origin = map.origin;
    smooth.voxel = map.voxel;
    smooth.data = map.data;
    const numerics::Shape3 mapShape = shapeZyx(map);
    numerics::gaussianFilter(smooth.data, mapShape,
                             std::array<double, 3>{0.5 * m_h / map.voxel[2], 0.5 * m_h / map.voxel[1],
                                                   0.5 * m_h / map.voxel[0]},
                             numerics::ExtendMode::Reflect);
    const std::size_t total = nx * ny * nz;
    std::vector<double> vals(total);
    std::size_t inside = 0;
    for (std::size_t z = 0; z < nz; ++z) {
        for (std::size_t y = 0; y < ny; ++y) {
            for (std::size_t x = 0; x < nx; ++x) {
                const Vec3 p{m_lo[0] + m_h * static_cast<double>(x), m_lo[1] + m_h * static_cast<double>(y),
                             m_lo[2] + m_h * static_cast<double>(z)};
                const double v = sample(smooth, mapShape, p);
                vals[(z * ny + y) * nx + x] = v;
                if (std::isfinite(v)) {
                    ++inside;
                }
            }
        }
    }
    double mean = 0.0;
    if (inside > 0) {
        mean = sum(total, [&](std::size_t i) { return std::isfinite(vals[i]) ? vals[i] : 0.0; }) /
               static_cast<double>(inside);
    }
    std::vector<double> e(total), e2(total);
    for (std::size_t i = 0; i < total; ++i) {
        e[i] = std::isfinite(vals[i]) ? vals[i] - mean : 0.0;
        e2[i] = e[i] * e[i];
    }
    m_eFt = numerics::rfftn(e, shape);
    m_e2Ft = numerics::rfftn(e2, shape);
    if (inside > 0) {  // e[inside].var()
        std::vector<double> in;
        in.reserve(inside);
        for (std::size_t i = 0; i < total; ++i) {
            if (std::isfinite(vals[i])) {
                in.push_back(e[i]);
            }
        }
        const double m = static_cast<double>(in.size());
        const double inMean = sum(in.size(), [&](std::size_t i) { return in[i]; }) / m;
        const double var = sum(in.size(), [&](std::size_t i) { return (in[i] - inMean) * (in[i] - inMean); }) / m;
        m_varFloor = 1e-3 * var;
    } else {
        m_varFloor = 1e-3 * 1.0;
    }
    // translations whose model centre lies inside the map box or up to radius / 2 outside it
    const Vec3 hi = map.extent();
    const double slack = 0.5 * radius;
    std::array<std::vector<unsigned char>, 3> ok;
    for (std::size_t k = 0; k < 3; ++k) {
        ok[k].resize(m_n[k]);
        for (std::size_t i = 0; i < m_n[k]; ++i) {
            const double c = m_lo[k] + m_h * static_cast<double>(i);
            ok[k][i] = c >= map.origin[k] - slack && c <= hi[k] + slack ? 1 : 0;
        }
    }
    m_valid.resize(total);
    for (std::size_t z = 0; z < nz; ++z) {
        for (std::size_t y = 0; y < ny; ++y) {
            for (std::size_t x = 0; x < nx; ++x) {
                m_valid[(z * ny + y) * nx + x] = static_cast<unsigned char>(ok[2][z] & ok[1][y] & ok[0][x]);
            }
        }
    }
    // the Gaussian transfer function of the model blur (sigma) on the coarse grid
    const std::size_t hx = nx / 2 + 1;
    const double factor = -2.0 * std::pow(std::numbers::pi, 2.0) * std::pow(sigma, 2.0);
    m_gauss.resize(nz * ny * hx);
    for (std::size_t z = 0; z < nz; ++z) {
        const double kz = fftFrequency(z, nz, m_h);
        for (std::size_t y = 0; y < ny; ++y) {
            const double ky = fftFrequency(y, ny, m_h);
            for (std::size_t x = 0; x < hx; ++x) {
                const double kx = static_cast<double>(x) * (1.0 / (static_cast<double>(nx) * m_h));
                const double k2 = kz * kz + ky * ky + kx * kx;
                m_gauss[(z * ny + y) * hx + x] = std::exp(factor * k2);
            }
        }
    }
    const Vec3& c = problem.centre();
    for (const Vec3& p : problem.xyz()) {
        m_rel.push_back({p[0] - c[0], p[1] - c[1], p[2] - c[2]});
    }
    m_w.assign(problem.weights().begin(), problem.weights().end());
}

std::vector<CoarseSearch::Peak> CoarseSearch::peaks(const Mat3& r, std::size_t count) const {
    const std::size_t nx = m_n[0], ny = m_n[1], nz = m_n[2];
    const numerics::FftShape3 shape{nz, ny, nx};
    const std::size_t total = nx * ny * nz;
    // the work arrays of this thread, reused from rotation to rotation
    thread_local std::vector<double> tmpl, mask, inMask, tz, corr;
    thread_local std::vector<Complex> product;
    tmpl.assign(total, 0.0);
    const std::size_t atoms = m_rel.size();
    std::vector<std::array<std::int64_t, 3>> base(atoms);
    std::vector<Vec3> f(atoms);
    for (std::size_t i = 0; i < atoms; ++i) {
        const Vec3& p = m_rel[i];
        for (std::size_t j = 0; j < 3; ++j) {
            const double g = (p[0] * r[j][0] + p[1] * r[j][1] + p[2] * r[j][2]) / m_h;
            const double b = std::floor(g);
            base[i][j] = static_cast<std::int64_t>(b);
            f[i][j] = g - b;
        }
    }
    const auto wrap = [](std::int64_t v, std::size_t n) {
        const auto m = static_cast<std::int64_t>(n);
        return static_cast<std::size_t>(((v % m) + m) % m);
    };
    for (int dx = 0; dx < 2; ++dx) {
        for (int dy = 0; dy < 2; ++dy) {
            for (int dz = 0; dz < 2; ++dz) {
                for (std::size_t i = 0; i < atoms; ++i) {
                    const double wx = dx ? f[i][0] : 1 - f[i][0];
                    const double wy = dy ? f[i][1] : 1 - f[i][1];
                    const double wz = dz ? f[i][2] : 1 - f[i][2];
                    const std::size_t z = wrap(base[i][2] + dz, nz), y = wrap(base[i][1] + dy, ny),
                                      x = wrap(base[i][0] + dx, nx);
                    tmpl[(z * ny + y) * nx + x] += m_w[i] * wx * wy * wz;
                }
            }
        }
    }
    std::vector<Complex> spectrum = numerics::rfftn(tmpl, shape);
    for (std::size_t i = 0; i < spectrum.size(); ++i) {
        spectrum[i] = {spectrum[i].real() * m_gauss[i], spectrum[i].imag() * m_gauss[i]};
    }
    const std::vector<double> t = numerics::irfftn(spectrum, shape);
    double tmax = t[0];
    for (const double v : t) {
        tmax = v > tmax || std::isnan(v) ? v : tmax;
    }
    const double level = kEnvelopeFraction * tmax;
    mask.resize(total);
    inMask.clear();
    for (std::size_t i = 0; i < total; ++i) {
        mask[i] = t[i] > level ? 1.0 : 0.0;
        if (t[i] > level) {
            inMask.push_back(t[i]);
        }
    }
    const double npts = static_cast<double>(inMask.size());
    const double tMean = sum(inMask.size(), [&](std::size_t i) { return inMask[i]; }) / npts;
    tz.resize(total);
    for (std::size_t i = 0; i < total; ++i) {
        tz[i] = mask[i] != 0.0 ? t[i] - tMean : 0.0;
    }
    const double tnorm = std::sqrt(sum(total, [&](std::size_t i) { return tz[i] * tz[i]; }));
    const std::vector<Complex> tzFt = numerics::rfftn(tz, shape);
    const std::vector<Complex> mFt = numerics::rfftn(mask, shape);
    product.resize(tzFt.size());
    for (std::size_t i = 0; i < product.size(); ++i) {
        product[i] = conjTimes(tzFt[i], m_eFt[i]);
    }
    const std::vector<double> num = numerics::irfftn(product, shape);
    for (std::size_t i = 0; i < product.size(); ++i) {
        product[i] = conjTimes(mFt[i], m_eFt[i]);
    }
    const std::vector<double> s1 = numerics::irfftn(product, shape);
    for (std::size_t i = 0; i < product.size(); ++i) {
        product[i] = conjTimes(mFt[i], m_e2Ft[i]);
    }
    const std::vector<double> s2 = numerics::irfftn(product, shape);
    const double nDiv = std::max(npts, 1.0);
    const double varFloor = m_varFloor * npts;
    const double tDiv = std::max(tnorm, 1e-12);
    corr.resize(total);
    for (std::size_t i = 0; i < total; ++i) {
        if (!m_valid[i]) {
            corr[i] = -std::numeric_limits<double>::infinity();
            continue;
        }
        const double v = s2[i] - s1[i] * s1[i] / nDiv;
        const double var = std::isnan(v) || v > varFloor ? v : varFloor;  // np.maximum
        corr[i] = num[i] / (std::sqrt(var) * tDiv);
    }
    std::vector<Peak> out;
    for (std::size_t c = 0; c < count; ++c) {
        std::size_t flat = 0;  // np.argmax: the first maximum, or the first NaN
        for (std::size_t i = 1; i < total && !std::isnan(corr[flat]); ++i) {
            if (corr[i] > corr[flat] || std::isnan(corr[i])) {
                flat = i;
            }
        }
        const double score = corr[flat];
        if (!std::isfinite(score)) {
            break;
        }
        const std::size_t z = flat / (ny * nx), y = (flat / nx) % ny, x = flat % nx;
        out.push_back({score,
                       {m_lo[0] + m_h * static_cast<double>(x), m_lo[1] + m_h * static_cast<double>(y),
                        m_lo[2] + m_h * static_cast<double>(z)}});
        // suppress the neighbourhood of this peak before the next one
        const std::size_t rad = 2;
        for (std::size_t zz = z >= rad ? z - rad : 0; zz < std::min(nz, z + rad + 1); ++zz) {
            for (std::size_t yy = y >= rad ? y - rad : 0; yy < std::min(ny, y + rad + 1); ++yy) {
                for (std::size_t xx = x >= rad ? x - rad : 0; xx < std::min(nx, x + rad + 1); ++xx) {
                    corr[(zz * ny + yy) * nx + xx] = -std::numeric_limits<double>::infinity();
                }
            }
        }
    }
    return out;
}

// ---- fit_model ----------------------------------------------------------------------------------

FitResult fitModel(const io::DensityMap& map, std::span<const Vec3> xyz, std::span<const std::string> elements,
                   double resolution, const std::optional<RigidTransform>& start, FitSearch search,
                   double angularStep, const FitControl& control) {
    const std::vector<double> weights = io::atomWeights(elements);
    const FitProblem prob(map, xyz, weights, resolution);
    FitResult out;
    if (search == FitSearch::Global) {
        out = prob.globalSearch(angularStep, kTopPoses, control);
    } else {
        const RigidTransform t0 = start ? *start : RigidTransform::identity();
        if (control.reporter) {
            control.reporter->checkCancel();
        }
        const double cc0 = prob.cc(t0);
        auto [fit, value] = prob.refine(t0);
        progress(control, 1.0, "local fit done");
        out.transform = fit;
        out.cc = value;
        out.start = start ? "manual" : "identity";
        out.evaluations = prob.evaluations();
        out.message = "start CC " + format("%.3f", cc0) + prob.insideNote(fit);
    }
    if (prob.nonFiniteVoxels() > 0) {  // CH-11
        out.message += (out.message.empty() ? "" : "; ") + std::to_string(prob.nonFiniteVoxels()) +
                       " non-finite map voxel(s) read as no density";
    }
    return out;
}

}  // namespace mdbridge::chem
