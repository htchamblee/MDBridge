// The density shadow (src/render/molecule/density_shadow.cpp, the affine resampling of ndimage.cpp)
// against MDATK's mdatk/density.py and scipy (data/render/density_shadow,
// tests/parity/make_density_shadow_cases.py): images and RGBA within 1e-6, placements within 1e-9,
// shapes, extents (1e-9 Å), None / errors exact; MDBridge's improvements RN-30 ... RN-32.

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <map>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/value.hpp"
#include "mdbridge/io/mrc.hpp"
#include "mdbridge/render/molecule/density_shadow.hpp"
#include "mdbridge/render/molecule/ndimage.hpp"

namespace mdbridge::render::molecule {
namespace {

using core::Array;
using core::Value;

constexpr double kImageTol = 1e-6;
constexpr double kPlacementTol = 1e-9;
const std::string kDemoMap = "render/inputs/demo_dry_map.mrc";
const std::string kSmallMap = "io/maps/files/axes_321le.mrc";

std::filesystem::path dataRoot() { return std::filesystem::path(MDBRIDGE_TEST_DATA); }

const Value& loadJson(std::string_view file) {
    static std::map<std::string, Value> cache;
    const std::string name(file);
    auto it = cache.find(name);
    if (it == cache.end()) {
        const auto path = dataRoot() / "render" / "density_shadow" / name;
        std::ifstream in(path, std::ios::binary);
        if (!in) {
            throw std::runtime_error("cannot read " + path.string());
        }
        std::stringstream text;
        text << in.rdbuf();
        it = cache.emplace(name, core::parseJson(text.str(), name)).first;
    }
    return it->second;
}

const Value& at(const Value& v, std::string_view key) {
    const Value* p = v.find(key);
    if (!p) {
        throw std::runtime_error("missing key " + std::string(key));
    }
    return *p;
}
bool has(const Value& v, std::string_view key) { return v.find(key) != nullptr; }
const Array& list(const Value& v) { return *v.get<Array>(); }
double number(const Value& v) {
    if (const auto* d = v.get<double>()) {
        return *d;
    }
    return static_cast<double>(*v.get<std::int64_t>());
}
std::size_t size(const Value& v) { return static_cast<std::size_t>(*v.get<std::int64_t>()); }
const std::string& str(const Value& v) { return *v.get<std::string>(); }
std::vector<double> numbers(const Value& v) {
    std::vector<double> out;
    for (const Value& x : list(v)) {
        out.push_back(number(x));
    }
    return out;
}
core::Mat3 mat3(const Value& v) {
    const auto n = numbers(v);
    return {{{n[0], n[1], n[2]}, {n[3], n[4], n[5]}, {n[6], n[7], n[8]}}};
}
core::Vec3 vec3(const Value& v) {
    const auto n = numbers(v);
    return {n[0], n[1], n[2]};
}

const io::DensityMap& mapOf(const std::string& rel) {
    static std::map<std::string, io::DensityMap> maps;
    auto it = maps.find(rel);
    if (it == maps.end()) {
        it = maps.emplace(rel, io::readMap(dataRoot() / std::filesystem::path(rel))).first;
    }
    return it->second;
}

template <class T>
void expectValues(const std::vector<T>& got, const std::vector<double>& want, double tol, const std::string& what) {
    ASSERT_EQ(got.size(), want.size()) << what;
    int reported = 0;
    for (std::size_t i = 0; i < got.size(); ++i) {
        const double g = static_cast<double>(got[i]);
        const double w = want[i];
        const bool ok = std::isnan(w) ? std::isnan(g)
                                      : (std::isinf(w) ? g == w : std::fabs(g - w) <= tol * std::max(1.0, std::fabs(w)));
        if (!ok && reported++ < 5) {
            ADD_FAILURE() << what << " [" << i << "]: " << g << " != " << w;
        }
    }
}

// ---- scipy.ndimage.affine_transform ------------------------------------------------------------------

TEST(DensityShadow, AffineTransformMatchesScipy) {
    const Array& cases = list(at(loadJson("affine.json"), "cases"));
    ASSERT_GE(cases.size(), 10u);
    for (const Value& c : cases) {
        const auto s = numbers(at(c, "shape"));
        const auto o = numbers(at(c, "output_shape"));
        const Shape3 shape = {static_cast<std::size_t>(s[0]), static_cast<std::size_t>(s[1]),
                              static_cast<std::size_t>(s[2])};
        const Shape3 outShape = {static_cast<std::size_t>(o[0]), static_cast<std::size_t>(o[1]),
                                 static_cast<std::size_t>(o[2])};
        std::vector<float> data;
        for (double x : numbers(at(c, "data"))) {
            data.push_back(static_cast<float>(x));
        }
        const auto mv = numbers(at(c, "matrix"));
        const Matrix3 m = {{{mv[0], mv[1], mv[2]}, {mv[3], mv[4], mv[5]}, {mv[6], mv[7], mv[8]}}};
        const auto ov = numbers(at(c, "offset"));
        const auto out = affineTransformLinearConstant(data, shape, m, {ov[0], ov[1], ov[2]}, outShape,
                                                       number(at(c, "cval")));
        // float32 results: scipy's arithmetic order is followed, so they agree to the last bit or so
        expectValues(out, numbers(at(c, "output")), 1e-6, str(at(c, "name")));
    }
}

TEST(DensityShadow, AffineTransformChecksTheShape) {
    const std::vector<float> data(5, 1.0f);
    EXPECT_THROW(affineTransformLinearConstant(data, {1, 2, 3}, Matrix3{}, {0, 0, 0}, {1, 1, 1}),
                 std::invalid_argument);
}

// ---- _above_level, project_density, shadow_rgba ---------------------------------------------------------

TEST(DensityShadow, AboveLevelMatchesMdatk) {
    int n = 0;
    for (const Value& c : list(at(loadJson("shadows.json"), "cases"))) {
        if (str(at(c, "kind")) != "above_level") {
            continue;
        }
        ++n;
        const auto& map = mapOf(str(at(c, "map")));
        const double level = number(at(c, "level"));
        const auto got = aboveLevel(map, level);
        // exact float32 values
        expectValues(got, numbers(at(c, "values")), 0.0, str(at(c, "map")) + " level " + std::to_string(level));
    }
    EXPECT_GE(n, 10);
}

void checkProject(const Value& c, const std::optional<Shadow>& sh, const std::string& name) {
    const Value& want = at(c, "shadow");
    if (want.find("rows") == nullptr) {
        EXPECT_FALSE(sh.has_value()) << name << ": MDATK gives None";
        return;
    }
    ASSERT_TRUE(sh.has_value()) << name;
    ASSERT_EQ(sh->rows, size(at(want, "rows"))) << name;
    ASSERT_EQ(sh->cols, size(at(want, "cols"))) << name;
    expectValues(sh->image, numbers(at(want, "image")), kImageTol, name + " image");
    const auto extent = numbers(at(want, "extent"));
    for (int k = 0; k < 4; ++k) {
        EXPECT_NEAR(sh->extent[k], extent[k], 1e-9 * std::max(1.0, std::fabs(extent[k]))) << name << " extent " << k;
    }
    for (const Value& r : list(at(c, "rgba"))) {
        const auto color = numbers(at(r, "color"));
        const auto got = shadowRgba(*sh, {color[0], color[1], color[2]}, number(at(r, "opacity")),
                                    number(at(r, "gamma")));
        expectValues(got, numbers(at(r, "values")), kImageTol, name + " rgba");
    }
}

TEST(DensityShadow, ProjectDensityMatchesMdatk) {
    int n = 0, images = 0, errors = 0;
    for (const Value& c : list(at(loadJson("shadows.json"), "cases"))) {
        if (str(at(c, "kind")) != "project") {
            continue;
        }
        ++n;
        const std::string name = str(at(c, "name"));
        SCOPED_TRACE(name);
        const auto& map = mapOf(str(at(c, "map")));
        core::RigidTransform t;
        t.rotation = mat3(at(c, "rotation"));
        t.translation = vec3(at(c, "translation"));
        std::optional<ShadowWindow> window;
        if (const Array* w = at(c, "window").get<Array>()) {
            window = ShadowWindow{number((*w)[0]), number((*w)[1]), number((*w)[2]), number((*w)[3])};
        }
        const int maxPixels = static_cast<int>(size(at(c, "max_pixels")));
        const int maxDepth = static_cast<int>(size(at(c, "max_depth")));
        const auto call = [&] {
            return projectDensity(map, t, mat3(at(c, "view")), vec3(at(c, "centre")), number(at(c, "level")), window,
                                  maxPixels, maxDepth);
        };
        if (has(c, "error")) {
            ++errors;
            try {
                call();
                ADD_FAILURE() << name << ": no error";
            } catch (const core::ValueError& e) {
                EXPECT_EQ(e.what(), str(at(at(c, "error"), "message")));
                EXPECT_STREQ(e.pythonName(), str(at(at(c, "error"), "type")).c_str());
            }
            continue;
        }
        const auto sh = call();
        images += sh.has_value() ? 1 : 0;
        checkProject(c, sh, name);
        // the cached overload gives the same image
        const auto shared = std::make_shared<const io::DensityMap>(map);
        const auto again = projectDensity(shared, t, mat3(at(c, "view")), vec3(at(c, "centre")), number(at(c, "level")),
                                          window, maxPixels, maxDepth);
        ASSERT_EQ(again.has_value(), sh.has_value());
        if (sh) {
            EXPECT_EQ(again->image, sh->image);
        }
    }
    EXPECT_GE(n, 35);
    EXPECT_GE(images, 25);
    EXPECT_GE(errors, 4);
}

// MDATK's test_projected_shadow_lies_under_the_model, on the demo map (simulated from demo_dry's first
// frame): the shadow's centroid is where the map's own density centroid projects, for every view
TEST(DensityShadow, ShadowCentroidFollowsTheView) {
    const auto& map = mapOf(kDemoMap);
    const core::Vec3 centre = {(map.origin[0] + map.extent()[0]) / 2, (map.origin[1] + map.extent()[1]) / 2,
                               (map.origin[2] + map.extent()[2]) / 2};
    // density-weighted centroid of the voxels above the level, in world (= map) coordinates
    const double level = io::mapLevel(map, 1.0);
    const auto above = aboveLevel(map, level);
    core::Vec3 mu{};
    double total = 0;
    std::size_t i = 0;
    for (std::int64_t z = 0; z < map.shapeXyz[2]; ++z) {
        for (std::int64_t y = 0; y < map.shapeXyz[1]; ++y) {
            for (std::int64_t x = 0; x < map.shapeXyz[0]; ++x, ++i) {
                const core::Vec3 p = map.indexToXyz({double(x), double(y), double(z)});
                for (int k = 0; k < 3; ++k) {
                    mu[k] += above[i] * (p[k] - centre[k]);
                }
                total += above[i];
            }
        }
    }
    for (int k = 0; k < 3; ++k) {
        mu[k] /= total;
    }
    for (const core::Vec3 rv : {core::Vec3{0, 0, 0}, core::Vec3{0.7, -0.5, 0.2}, core::Vec3{1.5, 0.1, 0.8}}) {
        const core::Mat3 v = core::rotvecMatrix(rv);
        const auto sh = projectDensity(map, core::RigidTransform::identity(), v, centre, level);
        ASSERT_TRUE(sh);
        double su = 0, sv = 0, st = 0;
        const double du = (sh->extent[1] - sh->extent[0]) / double(sh->cols);
        const double dv = (sh->extent[3] - sh->extent[2]) / double(sh->rows);
        for (std::size_t r = 0; r < sh->rows; ++r) {
            for (std::size_t col = 0; col < sh->cols; ++col) {
                const double w = sh->image[r * sh->cols + col];
                su += w * (sh->extent[0] + (double(col) + 0.5) * du);
                sv += w * (sh->extent[2] + (double(r) + 0.5) * dv);
                st += w;
            }
        }
        const double pu = v[0][0] * mu[0] + v[0][1] * mu[1] + v[0][2] * mu[2];
        const double pv = v[1][0] * mu[0] + v[1][1] * mu[1] + v[1][2] * mu[2];
        EXPECT_LT(std::hypot(su / st - pu, sv / st - pv), 1.0);
        const auto rgba = shadowRgba(*sh, {0.1, 0.2, 0.9}, 0.5);
        float amax = 0;
        for (std::size_t k = 3; k < rgba.size(); k += 4) {
            amax = std::max(amax, rgba[k]);
        }
        EXPECT_GT(amax, 0.0f);
        EXPECT_LE(amax, 0.5f);
    }
}

TEST(DensityShadow, RejectsBadLimits) {
    const auto& map = mapOf(kDemoMap);
    const core::Mat3 eye{{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
    EXPECT_THROW(projectDensity(map, {}, eye, {0, 0, 0}, 0.0, std::nullopt, 0, 10), std::invalid_argument);
    EXPECT_THROW(projectDensity(std::shared_ptr<const io::DensityMap>(), {}, eye, {0, 0, 0}, 0.0),
                 std::invalid_argument);
}

TEST(DensityShadow, AboveLevelCacheKeepsTwoEntries) {
    const auto a = std::make_shared<const io::DensityMap>(mapOf(kSmallMap));
    const auto b = std::make_shared<const io::DensityMap>(mapOf(kSmallMap));
    const auto first = aboveLevelCached(a, 1.0);
    EXPECT_EQ(aboveLevelCached(a, 1.0), first);   // same map object and level: shared
    EXPECT_NE(aboveLevelCached(b, 1.0), first);   // another map object (same contents): its own entry
    EXPECT_EQ(*aboveLevelCached(b, 1.0), *first);
    aboveLevelCached(a, 2.0);                     // a third entry evicts the oldest
    EXPECT_NE(aboveLevelCached(a, 1.0), first);
}

// ---- anchored ---------------------------------------------------------------------------------------------

TEST(DensityShadow, AnchoredMatchesMdatk) {
    int n = 0;
    for (const Value& c : list(at(loadJson("anchored.json"), "cases"))) {
        ++n;
        const std::string name = str(at(c, "name"));
        SCOPED_TRACE(name);
        const auto values = numbers(at(c, "values"));
        std::vector<std::vector<double>> anchor;
        for (const Value& r : list(at(c, "anchor"))) {
            anchor.push_back(numbers(r));
        }
        std::vector<core::Vec3> reference;
        for (const Value& r : list(at(c, "reference"))) {
            reference.push_back(vec3(r));
        }
        if (has(c, "error")) {
            try {
                anchored(values, anchor, reference);
                ADD_FAILURE() << "no error";
            } catch (const core::ValueError& e) {
                EXPECT_EQ(e.what(), str(at(at(c, "error"), "message")));
            }
            continue;
        }
        const AnchoredPlacement got = anchored(values, anchor, reference);
        const auto list12 = got.transform.toList(12);
        const auto want = numbers(at(c, "transform"));
        ASSERT_EQ(list12.size(), want.size());
        for (std::size_t k = 0; k < want.size(); ++k) {
            EXPECT_NEAR(list12[k], want[k], kPlacementTol) << k;
        }
        EXPECT_EQ(got.note, has(c, "note") ? str(at(c, "note")) : std::string());
    }
    EXPECT_GE(n, 15);
}

}  // namespace
}  // namespace mdbridge::render::molecule
