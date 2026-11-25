#pragma once
// The molecule render cases of tests/data/render/molecule (tests/parity/make_molecule_cases.py):
// molecule_cases.json (scenes, renders) with the items and MDATK's images in scenes.npz /
// renders.npz. The same reader as test_molecule_render.cpp's, for the Vulkan renderer's tests.

#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "mdbridge/core/json.hpp"
#include "mdbridge/core/value.hpp"
#include "mdbridge/render/molecule/geometry.hpp"
#include "mdbridge/render/molecule/scene.hpp"
#include "mdbridge/render/molecule/software.hpp"

namespace mdbridge::render::molecule::molecule_cases {

using core::Array;
using core::Value;

inline std::filesystem::path dataDir() { return std::filesystem::path(MDBRIDGE_TEST_DATA) / "render" / "molecule"; }

inline const Value& cases() {
    static const Value doc = [] {
        std::ifstream in(dataDir() / "molecule_cases.json", std::ios::binary);
        if (!in) {
            throw std::runtime_error("cannot read molecule_cases.json");
        }
        std::stringstream text;
        text << in.rdbuf();
        return core::parseJson(text.str(), "molecule_cases.json");
    }();
    return doc;
}

// An array of a .npz archive written by render_common.write_npz (ZIP, members stored), widened to
// double: '<f4', '<f8', '<i4', '<i8', '|i1', '|u1', '|b1' in C order.
struct Arr {
    std::vector<std::size_t> shape;
    std::vector<double> values;
    std::vector<double> toDouble() const { return values; }
};

inline std::uint32_t le32(const unsigned char* p) {
    return std::uint32_t{p[0]} | (std::uint32_t{p[1]} << 8) | (std::uint32_t{p[2]} << 16) | (std::uint32_t{p[3]} << 24);
}
inline std::uint16_t le16(const unsigned char* p) { return static_cast<std::uint16_t>(p[0] | (p[1] << 8)); }

inline Arr parseNpy(const unsigned char* d, std::size_t n, const std::string& name) {
    if (n < 10 || d[0] != 0x93 || std::string(reinterpret_cast<const char*>(d) + 1, 5) != "NUMPY") {
        throw std::runtime_error(name + ": not a .npy");
    }
    const std::size_t hl = d[6] == 1 ? le16(d + 8) : le32(d + 8);
    const std::size_t off = d[6] == 1 ? 10 : 12;
    const std::string header(reinterpret_cast<const char*>(d) + off, hl);
    auto field = [&](const std::string& key) {
        const std::size_t k = header.find("'" + key + "'");
        if (k == std::string::npos) {
            throw std::runtime_error(name + ": no " + key);
        }
        return header.substr(header.find(':', k) + 1);
    };
    std::string descr = field("descr");
    descr = descr.substr(descr.find('\'') + 1);
    descr = descr.substr(0, descr.find('\''));
    if (field("fortran_order").find("False") != 1 && field("fortran_order").find("False") != 0) {
        throw std::runtime_error(name + ": Fortran order");
    }
    std::string sh = field("shape");
    sh = sh.substr(sh.find('(') + 1, sh.find(')') - sh.find('(') - 1);
    Arr a;
    std::size_t count = 1;
    std::stringstream ss(sh);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        if (tok.find_first_not_of(' ') != std::string::npos) {
            a.shape.push_back(static_cast<std::size_t>(std::stoull(tok)));
            count *= a.shape.back();
        }
    }
    const unsigned char* p = d + off + hl;
    a.values.resize(count);
    for (std::size_t i = 0; i < count; ++i) {
        if (descr == "<f4") {
            float f;
            std::memcpy(&f, p + 4 * i, 4);
            a.values[i] = f;
        } else if (descr == "<f8") {
            std::memcpy(&a.values[i], p + 8 * i, 8);
        } else if (descr == "<i4") {
            std::int32_t v;
            std::memcpy(&v, p + 4 * i, 4);
            a.values[i] = v;
        } else if (descr == "<i8") {
            std::int64_t v;
            std::memcpy(&v, p + 8 * i, 8);
            a.values[i] = static_cast<double>(v);
        } else if (descr == "|i1") {
            a.values[i] = static_cast<std::int8_t>(p[i]);
        } else if (descr == "|u1" || descr == "|b1") {
            a.values[i] = p[i];
        } else {
            throw std::runtime_error(name + ": dtype " + descr);
        }
    }
    return a;
}

inline const std::map<std::string, Arr>& npzFile(const std::string& file) {
    static std::map<std::string, std::map<std::string, Arr>> cache;
    auto it = cache.find(file);
    if (it == cache.end()) {
        std::ifstream in(dataDir() / file, std::ios::binary);
        const std::vector<unsigned char> bytes((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        std::map<std::string, Arr> arrays;
        std::size_t pos = 0;
        while (pos + 30 <= bytes.size() && le32(bytes.data() + pos) == 0x04034b50) {
            const unsigned char* h = bytes.data() + pos;
            if (le16(h + 8) != 0) {
                throw std::runtime_error(file + ": compressed member");
            }
            const std::size_t size = le32(h + 18), nl = le16(h + 26), el = le16(h + 28);
            std::string name(reinterpret_cast<const char*>(h) + 30, nl);
            const unsigned char* data = h + 30 + nl + el;
            if (name.size() > 4 && name.substr(name.size() - 4) == ".npy") {
                name.resize(name.size() - 4);
            }
            arrays.emplace(name, parseNpy(data, size, file + "/" + name));
            pos += 30 + nl + el + size;
        }
        it = cache.emplace(file, std::move(arrays)).first;
    }
    return it->second;
}

inline const Arr& npz(const std::string& file, const std::string& key) {
    const auto& all = npzFile(file);
    const auto it = all.find(key);
    if (it == all.end()) {
        throw std::runtime_error(file + " has no " + key);
    }
    return it->second;
}

inline const Value& at(const Value& v, std::string_view key) {
    const Value* p = v.find(key);
    if (!p) {
        throw std::runtime_error("missing key " + std::string(key));
    }
    return *p;
}
inline double number(const Value& v) {
    if (const auto* d = v.get<double>()) {
        return *d;
    }
    return static_cast<double>(*v.get<std::int64_t>());
}
inline std::vector<double> doubles(const Value& v) {
    std::vector<double> out;
    for (const Value& x : *v.get<Array>()) {
        out.push_back(number(x));
    }
    return out;
}
inline Vec3 vec3(const Value& v) {
    const std::vector<double> d = doubles(v);
    return {d.at(0), d.at(1), d.at(2)};
}

inline std::vector<Vec3> vec3s(const Arr& a) {
    const std::vector<double> d = a.toDouble();
    std::vector<Vec3> out(d.size() / 3);
    for (std::size_t i = 0; i < out.size(); ++i) {
        out[i] = {d[3 * i], d[3 * i + 1], d[3 * i + 2]};
    }
    return out;
}
inline std::vector<std::int32_t> ints(const Arr& a) {
    std::vector<std::int32_t> out;
    for (double v : a.values) {
        out.push_back(static_cast<std::int32_t>(v));
    }
    return out;
}

// The items as make_molecule_cases.py stored them ("<prefix>/<k>/<field>", float32 widened to double).
inline std::vector<Colored> storedItems(const std::string& file, const std::string& prefix, const Value& records) {
    std::vector<Colored> items;
    std::size_t k = 0;
    for (const Value& rec : *records.get<Array>()) {
        const std::string p = prefix + "/" + std::to_string(k++) + "/";
        const std::string& kind = *at(rec, "kind").get<std::string>();
        Colored it;
        it.opacity = number(at(rec, "opacity"));
        it.rgb = vec3s(npz(file, p + "rgb"));
        if (kind == "mesh") {
            Mesh m;
            m.verts = vec3s(npz(file, p + "verts"));
            m.normals = vec3s(npz(file, p + "normals"));
            const std::vector<std::int32_t> t = ints(npz(file, p + "tris"));
            for (std::size_t i = 0; i + 2 < t.size(); i += 3) {
                m.tris.push_back({t[i], t[i + 1], t[i + 2]});
            }
            m.owner = ints(npz(file, p + "owner"));
            m.role.assign(m.verts.size(), 0);
            it.prim = std::move(m);
        } else if (kind == "spheres") {
            Spheres s;
            s.xyz = vec3s(npz(file, p + "xyz"));
            s.radius = npz(file, p + "radius").toDouble();
            s.owner = ints(npz(file, p + "owner"));
            s.atom = ints(npz(file, p + "atom"));
            it.prim = std::move(s);
        } else {
            Sticks s;
            s.a = vec3s(npz(file, p + "a"));
            s.b = vec3s(npz(file, p + "b"));
            s.radius = npz(file, p + "radius").toDouble();
            s.owner = ints(npz(file, p + "owner"));
            s.atom = ints(npz(file, p + "atom"));
            it.prim = std::move(s);
        }
        items.push_back(std::move(it));
    }
    return items;
}

struct RenderCase {
    std::string name;
    std::vector<Colored> items;
    Viewport vp;
    Shading shading;
    RenderQuality quality = RenderQuality::High;
    std::vector<std::uint8_t> expected;
    std::size_t width = 0, height = 0;
};

inline RenderCase renderCase(const std::string& name) {
    const Value& rec = at(at(cases(), "renders"), name);
    RenderCase c;
    c.name = name;
    const Value& src = at(rec, "source");
    if (const Value* scene = src.find("scene")) {
        const std::string& s = *scene->get<std::string>();
        c.items = storedItems("scenes.npz", s, at(at(at(cases(), "scenes"), s), "items"));
    } else {
        c.items = storedItems("renders.npz", name, at(rec, "items"));
    }
    const Value& vp = at(rec, "viewport");
    const Array* rot = at(vp, "rotation").get<Array>();
    for (std::size_t i = 0; i < 3; ++i) {
        c.vp.rotation[i] = vec3((*rot)[i]);
    }
    const std::vector<double> xl = doubles(at(vp, "xlim")), yl = doubles(at(vp, "ylim")), zr = doubles(at(vp, "zrange"));
    c.vp.xlim = {xl[0], xl[1]};
    c.vp.ylim = {yl[0], yl[1]};
    c.vp.zrange = {zr[0], zr[1]};
    c.vp.width = static_cast<int>(number(at(vp, "width")));
    c.vp.height = static_cast<int>(number(at(vp, "height")));
    const Value& sh = at(rec, "shading");
    c.shading.fog_rgb = vec3(at(sh, "fog_rgb"));
    c.shading.fog_strength = number(at(sh, "fog_strength"));
    if (!at(sh, "glow").isNull()) {
        c.shading.glow = doubles(at(sh, "glow"));
    }
    c.shading.glow_rgb = vec3(at(sh, "glow_rgb"));
    c.quality = *at(rec, "quality").get<std::string>() == "fast" ? RenderQuality::Fast : RenderQuality::High;
    const std::string key = name + "/rgba";
    const Arr img = npz("renders.npz", key);
    c.height = img.shape.at(0);
    c.width = img.shape.at(1);
    for (double v : img.values) {
        c.expected.push_back(static_cast<std::uint8_t>(v));
    }
    return c;
}

inline std::vector<std::string> renderNames() {
    std::vector<std::string> out;
    for (const auto& [k, v] : *at(cases(), "renders").get<core::Table>()) {
        out.push_back(k);
    }
    return out;
}


}  // namespace mdbridge::render::molecule::molecule_cases
