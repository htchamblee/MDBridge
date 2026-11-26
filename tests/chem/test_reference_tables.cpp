// The generated reference tables (src/chem/helix_tables.inc, component_tables.inc from
// tools/gen_chem_tables.py; their equality with MDATK's embedded tables is checked by the generator)
// against chemistry they must satisfy: bond lengths of the fibre monomer and of the 4JRD nucleotide,
// the O3'-P link the screws make, the L components' composition and charges.

#include <array>
#include <cmath>
#include <cstddef>
#include <map>
#include <span>
#include <string>
#include <string_view>

#include <gtest/gtest.h>

#include "mdbridge/chem/reference_tables.hpp"

using namespace mdbridge::chem::reference;

namespace {

constexpr double kPi = 3.14159265358979323846;

using Vec = std::array<double, 3>;

Vec cartesian(double r, double phiDegrees, double z) {
    const double p = phiDegrees * kPi / 180.0;
    return {r * std::cos(p), r * std::sin(p), z};
}

// One level on: a turn by -twist about z and a move of -rise along z (MDATK's _screw / _Screw).
Vec screwed(const Vec& v, double twist, double rise) {
    const double t = -twist * kPi / 180.0;
    return {std::cos(t) * v[0] - std::sin(t) * v[1], std::sin(t) * v[0] + std::cos(t) * v[1], v[2] - rise};
}

double distance(const Vec& a, const Vec& b) {
    return std::hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}

std::map<std::string, Vec, std::less<>> atoms(std::span<const CylindricalAtom> rows) {
    std::map<std::string, Vec, std::less<>> out;
    for (const CylindricalAtom& a : rows) {
        out[std::string(a.name)] = cartesian(a.r, a.phi, a.z);
    }
    return out;
}

}  // namespace

TEST(ReferenceTables, ArnottMonomerHasA_RnaBonds) {
    auto backbone = atoms(arnottBackbone());
    ASSERT_EQ(backbone.size(), 12U);
    EXPECT_NEAR(distance(backbone["P"], backbone["O5'"]), 1.6, 0.05);
    EXPECT_NEAR(distance(backbone["C1'"], backbone["O4'"]), 1.42, 0.05);
    EXPECT_NEAR(distance(backbone["C2'"], backbone["O2'"]), 1.42, 0.05);
    for (const char base : {'A', 'C', 'G', 'U'}) {
        const auto rows = arnottBase(base);
        ASSERT_FALSE(rows.empty()) << base;
        EXPECT_EQ(rows.front().name, (base == 'A' || base == 'G') ? "N9" : "N1") << base;
        const Vec n = cartesian(rows.front().r, rows.front().phi, rows.front().z);
        EXPECT_NEAR(distance(backbone["C1'"], n), 1.47, 0.03) << base;
    }
    EXPECT_EQ(arnottBase('A').size(), 10U);
    EXPECT_EQ(arnottBase('G').size(), 11U);
    EXPECT_EQ(arnottBase('C').size(), 8U);
    EXPECT_EQ(arnottBase('U').size(), 8U);
    EXPECT_TRUE(arnottBase('T').empty());
    // the 11-fold A-RNA screw of the builder (32.7 degrees, 2.81 A) links O3' to the next P at 1.59 A
    EXPECT_NEAR(distance(backbone["O3'"], screwed(backbone["P"], 32.7, 2.81)), 1.59, 0.01);
}

TEST(ReferenceTables, JrdNucleotideAndScrewCloseTheBackbone) {
    auto jrd = atoms(jrdNucleotide());
    ASSERT_EQ(jrd.size(), 22U);
    EXPECT_EQ(jrdNucleotide().front().name, "P");
    // frame: x through C1'
    EXPECT_EQ(jrd["C1'"][1], 0.0);
    EXPECT_EQ(jrd["C1'"][2], 0.0);
    EXPECT_NEAR(jrdTwist(), 41.75, 0.01);  // 8.6 residues per turn
    EXPECT_NEAR(jrdRise(), 3.79, 0.01);
    EXPECT_NEAR(distance(jrd["C1'"], jrd["N9"]), 1.47, 0.03);
    EXPECT_NEAR(distance(jrd["P"], jrd["OP1"]), 1.49, 0.03);
    // the rigid screw of the duplex links this O3' to the next residue's P at 1.48 A (crystal, A10 O3' -
    // A11 P: 1.60 A; MDATK's frame gives 1.77 A, CHD-02); the builder closes the links by torsions
    EXPECT_NEAR(distance(jrd["O3'"], screwed(jrd["P"], jrdTwist(), jrdRise())), 1.485, 0.001);
}

TEST(ReferenceTables, LComponentsInDictionaryForm) {
    ASSERT_EQ(components().size(), 4U);
    const std::map<std::string_view, std::array<int, 6>> formula = {  // C H N O P, net charge
        {"0A", {10, 12, 5, 7, 1}}, {"0C", {9, 12, 3, 8, 1}}, {"0G", {10, 12, 5, 8, 1}}, {"0U", {9, 11, 2, 9, 1}}};
    for (const ComponentInfo& c : components()) {
        std::map<std::string_view, int> count;
        int charge = 0;
        for (const ComponentAtom& a : componentAtoms(c.code)) {
            ++count[a.element];
            charge += a.charge;
            if (a.name == "OP2" || a.name == "OP3") {
                EXPECT_EQ(a.charge, -1) << c.code << " " << a.name;
            }
        }
        count["H"] = static_cast<int>(componentHydrogens(c.code).size());
        const auto& f = formula.at(c.code);
        EXPECT_EQ(count["C"], f[0]) << c.code;
        EXPECT_EQ(count["H"], f[1]) << c.code;
        EXPECT_EQ(count["N"], f[2]) << c.code;
        EXPECT_EQ(count["O"], f[3]) << c.code;
        EXPECT_EQ(count["P"], f[4]) << c.code;
        EXPECT_EQ(charge, -2) << c.code;
        EXPECT_EQ(componentAtoms(c.code).front().name, "OP3") << c.code;
        // a ring: atoms - 1 bonds would be a tree; ribose + bases add one ring per cycle
        const std::size_t heavy = componentAtoms(c.code).size();
        const std::size_t rings = (c.base == 'A' || c.base == 'G') ? 3 : 2;
        EXPECT_EQ(componentBonds(c.code).size(), heavy - 1 + rings) << c.code;
    }
    EXPECT_TRUE(componentAtoms("0T").empty());
    // the CCD's 0A names HN6 / HN6A / H5'A come under the standard names
    std::map<std::string_view, std::string_view> parent;
    for (const ComponentHydrogen& h : componentHydrogens("0A")) {
        parent[h.hydrogen] = h.heavy;
    }
    EXPECT_EQ(parent["H61"], "N6");
    EXPECT_EQ(parent["H62"], "N6");
    EXPECT_EQ(parent["H5''"], "C5'");
    EXPECT_EQ(parent.count("HOP2"), 0U);
}

TEST(ReferenceTables, HydrogenPositionsAndChirality) {
    std::map<std::string_view, std::string_view> g;
    for (const HydrogenPosition& p : hydrogenPositions("0G")) {
        g[p.hydrogen] = p.kind;
    }
    EXPECT_EQ(g["H21"], "trans");  // CHD-01: the CCD's G, H21 on the N3 side
    EXPECT_EQ(g["H22"], "cis");
    EXPECT_EQ(g["H5'"], "-");  // the mirror image of D's pro-S H5'
    for (const HydrogenPosition& p : hydrogenPositions("A")) {
        if (p.hydrogen == "H5'") {
            EXPECT_EQ(p.kind, "+");
        }
        if (p.hydrogen == "H61") {
            EXPECT_EQ(p.kind, "cis");
        }
    }
    EXPECT_EQ(hydrogenPositions("0U").size(), 2U);
    ASSERT_EQ(lRiboseChirality().size(), 4U);
    EXPECT_EQ(lRiboseChirality()[0].centre, "C4'");
    EXPECT_EQ(lRiboseChirality()[0].sign, 1);
    EXPECT_EQ(lRiboseChirality()[1].sign, -1);
}
