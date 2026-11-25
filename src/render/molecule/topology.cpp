// Port of MDATK's mdatk/render/molecule/topology.py (MDATK v0.0.4.1, the owner's code).

#include "mdbridge/render/molecule/topology.hpp"

#include <algorithm>
#include <cmath>
#include <string>
#include <utility>

#include "mdbridge/core/text.hpp"
#include "mdbridge/io/pdb.hpp"

namespace mdbridge::render::molecule {
namespace {

std::string upperAscii(std::string_view text) {
    std::string out(text);
    for (char& c : out) {
        if (c >= 'a' && c <= 'z') {
            c = static_cast<char>(c - 'a' + 'A');
        }
    }
    return out;
}

bool startsWith(std::string_view text, std::string_view prefix) {
    return text.substr(0, prefix.size()) == prefix;
}

struct ElementEntry {
    std::string_view element;
    double value;
};

// Bondi van der Waals radii (Å)
constexpr std::array<ElementEntry, 13> kVdw = {{{"H", 1.10}, {"C", 1.70}, {"N", 1.55}, {"O", 1.52},
                                                 {"P", 1.80}, {"S", 1.80}, {"F", 1.47}, {"MG", 1.73},
                                                 {"NA", 2.27}, {"K", 2.75}, {"CL", 1.75}, {"CA", 2.31},
                                                 {"ZN", 1.39}}};

struct RgbEntry {
    std::string_view key;
    Vec3 rgb;
};

constexpr std::array<RgbEntry, 7> kElementRgb = {{{"C", {0.56, 0.56, 0.56}},
                                                   {"N", {0.19, 0.31, 0.97}},
                                                   {"O", {0.90, 0.13, 0.13}},
                                                   {"P", {1.00, 0.55, 0.00}},
                                                   {"H", {0.92, 0.92, 0.92}},
                                                   {"S", {0.90, 0.78, 0.20}},
                                                   {"MG", {0.0, 0.51, 0.0}}}};

constexpr std::array<RgbEntry, 6> kNucleotideRgb = {{{"A", {0.27, 0.53, 0.80}},
                                                      {"G", {0.33, 0.66, 0.41}},
                                                      {"C", {0.85, 0.65, 0.13}},
                                                      {"U", {0.80, 0.33, 0.33}},
                                                      {"T", {0.80, 0.33, 0.33}},
                                                      {"?", {0.6, 0.6, 0.6}}}};

using Table = std::map<std::string, std::int32_t, std::less<>>;

template <std::size_t N>
bool hasAll(const Table& t, const std::array<std::string_view, N>& names) {
    return std::all_of(names.begin(), names.end(), [&](std::string_view n) { return t.find(n) != t.end(); });
}

bool has(const Table& t, std::string_view name) { return t.find(name) != t.end(); }

char letterFromAtoms(const Table& t) {
    if (hasAll(t, kPurineRing)) {
        return has(t, "O6") ? 'G' : has(t, "N6") ? 'A' : '?';
    }
    if (hasAll(t, kPyrimidineRing)) {
        if (has(t, "N4")) {
            return 'C';
        }
        if (has(t, "C7") || has(t, "C5M")) {
            return 'T';
        }
        return has(t, "O4") ? 'U' : '?';
    }
    return '?';
}

}  // namespace

std::string normaliseName(std::string_view name) {
    std::string n = core::pythonStrip(name);
    std::replace(n.begin(), n.end(), '*', '\'');
    if (n == "O1P") {
        return "OP1";
    }
    if (n == "O2P") {
        return "OP2";
    }
    if (n == "O3P") {
        return "OP3";
    }
    return n;
}

char baseLetter(std::string_view resname) {
    const std::string r = upperAscii(core::pythonStrip(resname));
    if (r == "PSU" || r == "PU") {
        return 'U';
    }
    if (r.size() <= 4) {
        const std::string_view head = std::string_view(r).substr(0, 3);
        static constexpr std::array<std::pair<std::string_view, char>, 6> charmm = {
            {{"ADE", 'A'}, {"GUA", 'G'}, {"CYT", 'C'}, {"URA", 'U'}, {"URI", 'U'}, {"THY", 'T'}}};
        for (const auto& [name, letter] : charmm) {
            if (head == name) {
                return letter;
            }
        }
        for (std::string_view prefix : {"RA", "RG", "RC", "RU", "DA", "DG", "DC", "DT", "LA", "LG", "LC", "LU"}) {
            if (startsWith(r, prefix)) {
                return prefix[1];
            }
        }
    }
    for (char k : {'A', 'G', 'C', 'U', 'T'}) {
        if (r.size() == 1 && r[0] == k) {
            return k;
        }
        if (!r.empty() && r[0] == k) {
            const std::string_view rest = std::string_view(r).substr(1);
            if (rest == "3" || rest == "5" || rest == "N") {
                return k;
            }
        }
    }
    return '?';
}

double vdwRadius(std::string_view element) {
    for (const auto& e : kVdw) {
        if (e.element == element) {
            return e.value;
        }
    }
    return 1.7;
}

Vec3 elementRgb(std::string_view element) {
    for (const auto& e : kElementRgb) {
        if (e.key == element) {
            return e.rgb;
        }
    }
    return {0.6, 0.4, 0.7};
}

Vec3 nucleotideRgb(char letter) {
    for (const auto& e : kNucleotideRgb) {
        if (e.key[0] == letter) {
            return e.rgb;
        }
    }
    return kNucleotideRgb.back().rgb;
}

Topology::Topology(const TopologyInput& input)
    : n_res(input.n_res), row(input.row), bonds(input.bonds), chain(input.chain), resindex(input.resindex) {
    const std::size_t n_atoms = row.size();
    names.reserve(n_atoms);
    // a CCD component's own atom names (0DA's N9A ...) read as the standard ones (MDATK 0.1.x)
    for (std::size_t i = 0; i < input.names.size(); ++i) {
        const std::string& n = input.names[i];
        names.push_back(normaliseName(i < input.resnames.size() ? io::standardAtomName(input.resnames[i], n) : n));
    }
    elements.reserve(n_atoms);
    for (const std::string& e : input.elements) {
        std::string el = upperAscii(core::pythonStrip(e));
        elements.push_back(el.empty() ? std::string("C") : std::move(el));
    }
    const auto nr = static_cast<std::size_t>(std::max(n_res, 0));
    table.assign(nr, {});
    std::map<std::int64_t, std::int32_t> first;
    for (std::size_t i = 0; i < n_atoms && i < names.size(); ++i) {
        const std::int32_t r = row[i];
        if (r >= 0 && r < n_res) {
            table[static_cast<std::size_t>(r)].emplace(names[i], static_cast<std::int32_t>(i));
        }
        first.emplace(r, static_cast<std::int32_t>(i));
    }
    resname.resize(nr);
    letter.resize(nr);
    for (std::size_t r = 0; r < nr; ++r) {
        const auto it = first.find(static_cast<std::int64_t>(r));
        resname[r] = it != first.end() ? input.resnames[static_cast<std::size_t>(it->second)] : std::string("?");
        letter[r] = baseLetter(resname[r]);
        if (letter[r] == '?') {  // unknown names: tell the base by its exocyclic atoms
            letter[r] = letterFromAtoms(table[r]);
        }
    }
    radius.reserve(n_atoms);
    element_rgb.reserve(n_atoms);
    is_carbon.reserve(n_atoms);
    is_hydrogen.reserve(n_atoms);
    for (const std::string& e : elements) {
        radius.push_back(vdwRadius(e));
        element_rgb.push_back(elementRgb(e));
        is_carbon.push_back(e == "C");
        is_hydrogen.push_back(e == "H");
    }

    auto pick = [&](std::size_t r, std::initializer_list<std::string_view> candidates) -> std::int32_t {
        for (std::string_view n : candidates) {
            const auto it = table[r].find(n);
            if (it != table[r].end()) {
                return it->second;
            }
        }
        return -1;
    };
    const bool c4_backbone = input.backbone_atom == "C4'";
    control.resize(nr);
    anchor_atom.resize(nr);
    side.resize(nr);
    c1.resize(nr);
    c4.resize(nr);
    for (std::size_t r = 0; r < nr; ++r) {
        control[r] = c4_backbone ? pick(r, {"C4'", "C3'", "C5'", "C1'", "CA"})
                                 : pick(r, {"P", "O5'", "C5'", "C4'", "CA"});
        if (control[r] >= 0) {
            anchor_atom[r] = control[r];
        } else {
            const auto it = first.find(static_cast<std::int64_t>(r));
            anchor_atom[r] = it != first.end() ? it->second : -1;
        }
        side[r] = pick(r, {"C1'", "C4'", "CB", "O"});
        c1[r] = pick(r, {"C1'"});
        c4[r] = pick(r, {"C4'"});
    }
    for (std::size_t r = 0; r < nr; ++r) {
        const Table& t = table[r];
        const auto rr = static_cast<std::int32_t>(r);
        if (hasAll(t, kPurineRing)) {
            std::array<std::int32_t, 9> ring{};
            for (std::size_t k = 0; k < ring.size(); ++k) {
                ring[k] = t.find(kPurineRing[k])->second;
            }
            purine_res.push_back(rr);
            purine_ring.push_back(ring);
            purine_gly.push_back(t.find("N9")->second);
        } else if (hasAll(t, kPyrimidineRing)) {
            std::array<std::int32_t, 6> ring{};
            for (std::size_t k = 0; k < ring.size(); ++k) {
                ring[k] = t.find(kPyrimidineRing[k])->second;
            }
            // pseudouridine is C-glycosidic (C1'-C5)
            const std::string up = upperAscii(resname[r]);
            const bool psu = up == "PSU" || up == "PU";
            pyrimidine_res.push_back(rr);
            pyrimidine_ring.push_back(ring);
            pyrimidine_gly.push_back(t.find(psu ? "C5" : "N1")->second);
        }
        if (hasAll(t, kSugarRing)) {
            std::array<std::int32_t, 5> ring{};
            for (std::size_t k = 0; k < ring.size(); ++k) {
                ring[k] = t.find(kSugarRing[k])->second;
            }
            sugar_res.push_back(rr);
            sugar_ring.push_back(ring);
        }
    }
    segments = computeSegments(input.ref_xyz);
}

std::vector<std::vector<std::int32_t>> Topology::computeSegments(std::span<const Vec3> xyz) const {
    std::vector<std::vector<std::int32_t>> segs;
    std::vector<std::int32_t> cur;
    auto close = [&] {
        if (cur.size() > 1) {
            segs.push_back(cur);
        }
        cur.clear();
    };
    for (std::int32_t r = 0; r < n_res; ++r) {
        const std::int32_t c = control[static_cast<std::size_t>(r)];
        if (c < 0) {
            close();
            continue;
        }
        if (!cur.empty()) {
            const std::int32_t p = cur.back();
            const auto ur = static_cast<std::size_t>(r);
            const auto up = static_cast<std::size_t>(p);
            const Vec3& a = xyz[static_cast<std::size_t>(c)];
            const Vec3& b = xyz[static_cast<std::size_t>(control[up])];
            const double dx = a[0] - b[0];
            const double dy = a[1] - b[1];
            const double dz = a[2] - b[2];
            const bool broken = chain[ur] != chain[up] || resindex[ur] - resindex[up] != 1 ||
                                std::sqrt(dx * dx + dy * dy + dz * dz) > kSegmentBreak;
            if (broken) {
                close();
            }
        }
        cur.push_back(r);
    }
    close();
    return segs;
}

std::vector<std::int32_t> Topology::nucleic() const {
    std::vector<std::int32_t> out(purine_res);
    out.insert(out.end(), pyrimidine_res.begin(), pyrimidine_res.end());
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

std::int32_t Topology::atomOf(std::int32_t residue, std::string_view name) const {
    if (residue < 0 || residue >= n_res) {
        return -1;
    }
    const auto& t = table[static_cast<std::size_t>(residue)];
    const auto it = t.find(name);
    return it != t.end() ? it->second : -1;
}

}  // namespace mdbridge::render::molecule
