// Restraint dictionaries and SMILES (src/chem/dictionaries.cpp, smiles.cpp) against what MDATK's
// mdatk/dictionaries.py gives (tests/parity/make_dict_cases.py -> data/chem/dict_cases.json, with
// improvements_chem.py's CHD-01), and the counterparts of MDATK's tests/test_dictionaries.py
// (named after them; the mmCIF and find_dictionary tests are in tests/io, the PubChem and command
// line tests with the `dict` commands, the planner-step tests with milestone 7). Everything exact:
// written dictionaries byte for byte, messages equal.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "mdbridge/chem/dictionaries.hpp"
#include "mdbridge/chem/smiles.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/io/cif.hpp"
#include "test_support.hpp"

using namespace mdbridge;
using namespace mdbridge::chem;
using core::Array;
using core::Table;
using core::Value;
using test::kData;

namespace {

const Value& cases() {
    static const Value data =
        core::parseJson(test::readFile(kData / "chem" / "dict_cases.json"), "dict_cases.json");
    return data;
}

const Value& at(const Value& value, std::string_view key) {
    const Value* v = value.find(key);
    if (!v) {
        throw std::runtime_error("missing key " + std::string(key));
    }
    return *v;
}
const Array& list(const Value& value) { return *value.get<Array>(); }
std::string str(const Value& value) { return *value.get<std::string>(); }
bool flag(const Value& value) { return *value.get<bool>(); }

Value strings(const std::vector<std::string>& items) {
    Array out;
    for (const auto& s : items) {
        out.emplace_back(s);
    }
    return out;
}

Value errorValue(const DictionaryError& e) {
    return Table{{"type", Value("DictionaryError")}, {"message", Value(std::string(e.what()))}};
}

Value graphValue(const MolGraph& g) {
    Array charges, hcount, adj, stereo;
    for (const auto& q : g.charges) {
        charges.emplace_back(q ? Value(*q) : Value());
    }
    for (const int h : g.hcount) {
        hcount.emplace_back(h);
    }
    for (const auto& a : g.adj) {
        std::vector<std::pair<int, int>> sorted(a.begin(), a.end());
        std::sort(sorted.begin(), sorted.end());
        Array row;
        for (const auto& [j, o] : sorted) {
            row.emplace_back(Array{Value(j), Value(o)});
        }
        adj.emplace_back(std::move(row));
    }
    for (const auto& [k, s] : g.stereo) {
        Array order;
        for (const int x : s.order) {
            order.emplace_back(x == MolGraph::kHydrogen ? Value("H") : Value(x));
        }
        stereo.emplace_back(Array{Value(k), Value(std::move(order)), Value(s.tag)});
    }
    return Table{{"elements", strings(g.elements)}, {"charges", charges}, {"hcount", hcount},
                 {"adj", adj},  {"stereo", stereo},  {"formula", Value(g.formula())}};
}

Value componentValue(const Component& c) {
    Array atoms, bonds, hydrogens, hpos, chirality;
    for (const auto& a : c.atoms) {
        atoms.emplace_back(Array{Value(a.name), Value(a.element), Value(a.charge)});
    }
    for (const auto& b : c.bonds) {
        bonds.emplace_back(Array{Value(b.atom1), Value(b.atom2), Value(b.order)});
    }
    for (const auto& [x, hs] : c.hydrogens) {
        hydrogens.emplace_back(Array{Value(x), strings(hs)});
    }
    for (const auto& r : c.hpos) {
        hpos.emplace_back(Array{Value(r.hydrogen), Value(r.x), Value(r.a), Value(r.b), Value(r.kind)});
    }
    for (const auto& r : c.chirality()) {
        chirality.emplace_back(Array{Value(r.centre), strings({r.order.begin(), r.order.end()}), Value(r.sign)});
    }
    return Table{{"code", Value(c.code)},
                 {"parent", Value(c.parent)},
                 {"base", Value(c.base)},
                 {"terminal", Value(c.terminal)},
                 {"name", Value(c.name)},
                 {"smiles", Value(c.smiles)},
                 {"atoms", atoms},
                 {"bonds", bonds},
                 {"hydrogens", hydrogens},
                 {"hpos", hpos},
                 {"leaving", strings(c.leaving)},
                 {"description", Value(c.description)},
                 {"formula", Value(c.formula())},
                 {"residue_atoms", strings(c.residueAtoms())},
                 {"chirality", chirality}};
}

FinishOptions optionsOf(const Value& options) {
    FinishOptions o;
    if (const Value* v = options.find("group")) {
        o.group = str(*v);
    }
    if (const Value* v = options.find("name")) {
        o.name = str(*v);
    }
    if (const Value* v = options.find("source_id")) {
        o.sourceId = str(*v);
    }
    if (const Value* v = options.find("allow_protonation")) {
        o.allowProtonation = flag(*v);
    }
    if (const Value* v = options.find("keep_names")) {
        o.keepNames = flag(*v);
    }
    return o;
}

const Value& finishCase(std::string_view label) {
    for (const Value& c : list(at(cases(), "finish"))) {
        if (str(at(c, "label")) == label) {
            return c;
        }
    }
    throw std::runtime_error("no finish case " + std::string(label));
}

FinishResult finish(const Value& c) {
    return finishDictionary(str(at(c, "input")), str(at(c, "code")), optionsOf(at(c, "options")));
}

std::string mappingOf(const FinishResult& res, const std::string& old) {
    for (const auto& [o, n] : res.mapping) {
        if (o == old) {
            return n;
        }
    }
    throw std::runtime_error("no mapping of " + old);
}

const io::CifLoop& loopOf(const io::CifDocument& doc, std::string_view block, std::string_view category) {
    const io::CifBlock* b = doc.block(block);
    if (!b || !b->get(category)) {
        throw std::runtime_error("no " + std::string(category) + " in " + std::string(block));
    }
    return *b->get(category);
}

// The CCD's own isomeric SMILES of 0A 0C 0G 0U (wwPDB CCD, CC0: _pdbx_chem_comp_descriptor
// SMILES_CANONICAL of CACTVS and OpenEye OEToolkits; components.cif of 2026-11-14).
const std::map<std::string, std::array<std::string, 2>> kCcdSmiles = {
    {"0A",
     {"Nc1ncnc2n(cnc12)[C@H]3O[C@@H](CO[P](O)(O)=O)[C@H](O)[C@@H]3O",
      "c1nc(c2c(n1)n(cn2)[C@@H]3[C@H]([C@H]([C@@H](O3)COP(=O)(O)O)O)O)N"}},
    {"0C",
     {"NC1=NC(=O)N(C=C1)[C@H]2O[C@@H](CO[P](O)(O)=O)[C@H](O)[C@@H]2O",
      "C1=CN(C(=O)N=C1N)[C@@H]2[C@H]([C@H]([C@@H](O2)COP(=O)(O)O)O)O"}},
    {"0G",
     {"NC1=Nc2n(cnc2C(=O)N1)[C@H]3O[C@@H](CO[P](O)(O)=O)[C@H](O)[C@@H]3O",
      "c1nc2c(n1[C@@H]3[C@H]([C@H]([C@@H](O3)COP(=O)(O)O)O)O)N=C(NC2=O)N"}},
    {"0U",
     {"O[C@@H]1[C@H](O)[C@H](O[C@H]1CO[P](O)(O)=O)N2C=CC(=O)NC2=O",
      "C1=CN(C(=O)NC1=O)[C@@H]2[C@H]([C@H]([C@@H](O2)COP(=O)(O)O)O)O"}},
};

// MDATK's test strings (tests/test_dictionaries.py)
const std::string kDAdenosine = "C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)CO)O)O)N";
const std::string kLAdenosine = "C1=NC(=C2C(=N1)N(C=N2)[C@@H]3[C@H]([C@H]([C@@H](O3)CO)O)O)N";
const std::string kMonomerLibrary0A = "Nc1ncnc2n(cnc12)[C@@H]3O[C@H](CO[P](O)(O)=O)[C@@H](O)[C@H]3O";

using Vec = std::array<double, 3>;

double volume(const Vec& c, const Vec& a, const Vec& b, const Vec& d) {
    const Vec va{a[0] - c[0], a[1] - c[1], a[2] - c[2]}, vb{b[0] - c[0], b[1] - c[1], b[2] - c[2]},
        vd{d[0] - c[0], d[1] - c[1], d[2] - c[2]};
    return va[0] * (vb[1] * vd[2] - vb[2] * vd[1]) - va[1] * (vb[0] * vd[2] - vb[2] * vd[0]) +
           va[2] * (vb[0] * vd[1] - vb[1] * vd[0]);
}

double torsion(const Vec& p0, const Vec& p1, const Vec& p2, const Vec& p3) {
    const auto sub = [](const Vec& a, const Vec& b) { return Vec{a[0] - b[0], a[1] - b[1], a[2] - b[2]}; };
    const auto dot = [](const Vec& a, const Vec& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; };
    const Vec b0 = sub(p0, p1), b1 = sub(p2, p1), b2 = sub(p3, p2);
    const double n = std::sqrt(dot(b1, b1));
    const Vec u{b1[0] / n, b1[1] / n, b1[2] / n};
    const Vec v = sub(b0, Vec{dot(b0, u) * u[0], dot(b0, u) * u[1], dot(b0, u) * u[2]});
    const Vec w = sub(b2, Vec{dot(b2, u) * u[0], dot(b2, u) * u[1], dot(b2, u) * u[2]});
    const Vec c{u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]};
    return std::atan2(dot(c, w), dot(v, w)) * 180.0 / 3.141592653589793;
}

}  // namespace

// --- the generated cases ------------------------------------------------------------------------

TEST(DictionaryCases, ComponentsEqualMdatk) {
    const Array& want = list(at(cases(), "components"));
    const auto& comps = dictionaryComponents();
    ASSERT_EQ(comps.size(), want.size());
    for (std::size_t k = 0; k < comps.size(); ++k) {
        const Value got = componentValue(comps[k]);
        EXPECT_EQ(got, want[k]) << core::toJson(got, -1) << "\n" << core::toJson(want[k], -1);
    }
}

TEST(DictionaryCases, SmilesParser) {
    for (const Value& c : list(at(cases(), "smiles"))) {
        const std::string smiles = str(at(c, "smiles"));
        Value got;
        try {
            got = graphValue(parseSmiles(smiles));
            EXPECT_EQ(got, at(c, "graph")) << smiles << "\n" << core::toJson(got, -1);
        } catch (const DictionaryError& e) {
            ASSERT_NE(c.find("error"), nullptr) << smiles << ": " << e.what();
            EXPECT_EQ(errorValue(e), at(c, "error")) << smiles;
        }
    }
}

TEST(DictionaryCases, GraphMatching) {
    for (const Value& c : list(at(cases(), "match"))) {
        const Component& comp = *findComponent(str(at(c, "code")));
        const MolGraph ref = comp.graph();
        std::vector<int> perm;
        for (const Value& p : list(at(c, "perm"))) {
            perm.push_back(static_cast<int>(*p.get<std::int64_t>()));
        }
        std::vector<int> inverse(perm.size());
        for (std::size_t k = 0; k < perm.size(); ++k) {
            inverse[static_cast<std::size_t>(perm[k])] = static_cast<int>(k);
        }
        MolGraph tgt;
        for (const int old : perm) {
            const auto o = static_cast<std::size_t>(old);
            tgt.addAtom(ref.elements[o], ref.charges[o], ref.hcount[o]);
        }
        for (std::size_t i = 0; i < ref.size(); ++i) {
            for (const auto& [j, order] : ref.adj[i]) {
                if (static_cast<int>(i) < j) {
                    tgt.addBond(inverse[i], inverse[static_cast<std::size_t>(j)], order);
                }
            }
        }
        if (flag(at(c, "removed"))) {
            for (std::size_t i = 0; i < tgt.size(); ++i) {
                if (!tgt.adj[i].empty()) {
                    tgt.removeBond(static_cast<int>(i), tgt.adj[i].front().first);
                    break;
                }
            }
        }
        const auto m = matchGraph(ref, tgt);
        Value got;
        if (m) {
            Array a;
            for (const int x : *m) {
                a.emplace_back(x);
            }
            got = a;
        }
        EXPECT_EQ(got, at(c, "match")) << str(at(c, "code"));
    }
}

TEST(DictionaryCases, CheckAndIdentifySmiles) {
    const Value& data = at(cases(), "check");
    for (const Value& row : list(at(data, "check"))) {
        const Array& r = list(row);
        Value got;
        try {
            got = strings(checkSmiles(str(r[0]), str(r[1])));
        } catch (const DictionaryError& e) {
            got = errorValue(e);
        }
        EXPECT_EQ(got, r[2]) << str(r[0]) << " " << str(r[1]) << ": " << core::toJson(got, -1);
    }
    for (const Value& row : list(at(data, "identify"))) {
        const Array& r = list(row);
        Value got;
        try {
            got = strings(identifySmiles(str(r[0])));
        } catch (const DictionaryError& e) {
            got = errorValue(e);
        }
        EXPECT_EQ(got, r[1]) << str(r[0]);
    }
}

TEST(DictionaryCases, FinishDictionary) {
    std::size_t n = 0;
    for (const Value& c : list(at(cases(), "finish"))) {
        const std::string label = str(at(c, "label"));
        try {
            const FinishResult res = finish(c);
            ASSERT_EQ(c.find("error"), nullptr) << label << " gave no error";
            EXPECT_EQ(res.text, str(at(c, "text"))) << label;
            EXPECT_EQ(strings(res.notes), at(c, "notes")) << label << ": " << core::toJson(strings(res.notes), -1);
            Table mapping;
            for (const auto& [o, nn] : res.mapping) {
                mapping.emplace_back(o, Value(nn));
            }
            EXPECT_EQ(Value(mapping), at(c, "mapping")) << label;
            EXPECT_EQ(res.sourceId, str(at(c, "source_id"))) << label;
            EXPECT_EQ(res.program, str(at(c, "program"))) << label;
        } catch (const DictionaryError& e) {
            ASSERT_NE(c.find("error"), nullptr) << label << ": " << e.what();
            EXPECT_EQ(errorValue(e), at(c, "error")) << label;
        }
        ++n;
    }
    EXPECT_EQ(n, 103U);
}

// --- MDATK's tests/test_dictionaries.py --------------------------------------------------------

TEST(Dictionaries, ComponentTableAndVariantConventions) {
    std::vector<std::string> codes;
    for (const Component& c : dictionaryComponents()) {
        codes.push_back(c.code);
    }
    EXPECT_EQ(codes, (std::vector<std::string>{"0A5", "0A3", "0C5", "0C3", "0G5", "0G3", "0U5", "0U3", "0A", "0C",
                                               "0G", "0U"}));
    const std::map<std::string, std::pair<std::size_t, std::size_t>> ccdCounts = {
        {"0A", {23, 35}}, {"0C", {21, 33}}, {"0G", {24, 36}}, {"0U", {21, 32}}};
    for (const auto& [code, counts] : ccdCounts) {
        const Component& comp = *findComponent(code);
        const auto atoms = comp.residueAtoms();
        const std::set<std::string> all(atoms.begin(), atoms.end());
        EXPECT_EQ(atoms.size(), counts.second);
        EXPECT_EQ(comp.atoms.size(), counts.first);
        EXPECT_TRUE(comp.leaving.empty());
        for (const char* a : {"OP3", "P", "OP1", "OP2", "HO3'"}) {
            EXPECT_EQ(all.count(a), 1U) << code << " " << a;
        }
        EXPECT_EQ(all.count("HO5'"), 0U);
        const Component& five = *findComponent(code + "5");
        const Component& three = *findComponent(code + "3");
        const auto a5v = five.residueAtoms(), a3v = three.residueAtoms();
        std::set<std::string> a5(a5v.begin(), a5v.end()), a3(a3v.begin(), a3v.end());
        std::set<std::string> want5 = all;
        for (const char* a : {"P", "OP1", "OP2", "OP3"}) {
            want5.erase(a);
        }
        want5.insert("HO5'");
        EXPECT_EQ(a5, want5) << code;
        std::set<std::string> want3 = all;
        want3.erase("OP3");
        EXPECT_EQ(a3, want3) << code;
        EXPECT_EQ(three.leaving, std::vector<std::string>{"OP3"});
        EXPECT_EQ(three.smiles, comp.smiles);  // generated from the parent, OP3 removed afterwards
        EXPECT_EQ(five.parent, code);
        EXPECT_EQ(three.parent, code);
        EXPECT_EQ(five.terminal, "5");
        EXPECT_EQ(three.terminal, "3");
    }
    EXPECT_EQ(findComponent("0A5")->formula(), "C10H13N5O4");
    EXPECT_EQ(findComponent("0G3")->formula(), "C10H12N5O8P-2");
    try {
        getComponent("XYZ");
        ADD_FAILURE() << "no error";
    } catch (const DictionaryError& e) {
        EXPECT_NE(std::string(e.what()).find("unknown residue code"), std::string::npos);
    }
}

TEST(Dictionaries, BuiltinSmilesMatchReferenceGraphAndStereo) {
    for (const Component& comp : dictionaryComponents()) {
        EXPECT_TRUE(checkSmiles(comp.code, comp.smiles).empty()) << comp.code;
        EXPECT_EQ(parseSmiles(comp.smiles).formula(), comp.formula()) << comp.code;
        const auto ids = identifySmiles(comp.smiles);
        EXPECT_NE(std::find(ids.begin(), ids.end(), comp.code), ids.end()) << comp.code;
    }
}

TEST(Dictionaries, SmilesChecksRejectWrongStereoAndCompounds) {
    const auto message = [](auto&& f) {
        try {
            f();
        } catch (const DictionaryError& e) {
            return std::string(e.what());
        }
        return std::string("no error");
    };
    EXPECT_TRUE(checkSmiles("0A5", kLAdenosine).empty());
    EXPECT_NE(message([&] { checkSmiles("0A5", kDAdenosine); }).find("opposite configuration"), std::string::npos);
    // MDATK's copy of the monomer library's 0A descriptor (D); the CCD's own one is L (below)
    EXPECT_NE(message([&] { checkSmiles("0A", kMonomerLibrary0A); }).find("opposite configuration"),
              std::string::npos);
    EXPECT_NE(message([&] { checkSmiles("0A5", "Nc1ncnc2n(cnc12)C3OC(CO)C(O)C3O"); }).find("does not specify"),
              std::string::npos);
    EXPECT_NE(message([&] { checkSmiles("0G5", kLAdenosine); }).find("does not match 0G5"), std::string::npos);
    EXPECT_NE(message([&] { checkSmiles("0G5", "OC[C@@H]1O[C@H](N2C=[NH+]C3=C2N=C(N)NC3=O)[C@@H](O)[C@H]1O"); })
                  .find("tautomer"),
              std::string::npos);
    std::string neutral = findComponent("0U")->smiles;
    neutral.replace(neutral.find("[O-]P(=O)([O-])"), 15, "OP(=O)(O)");
    const auto notes = checkSmiles("0U3", neutral);
    EXPECT_TRUE(std::any_of(notes.begin(), notes.end(),
                            [](const std::string& n) { return n.find("protonation") != std::string::npos; }));
    EXPECT_EQ(identifySmiles(kLAdenosine), std::vector<std::string>{"0A5"});
    EXPECT_TRUE(identifySmiles(kDAdenosine).empty());
}

TEST(Dictionaries, SmilesParserDetails) {
    MolGraph g = parseSmiles("c1ccncc1");  // aromatic: implicit H from the lowest valence
    EXPECT_EQ(g.formula(), "C5H5N");
    EXPECT_EQ(g.hcount, (std::vector<int>{1, 1, 1, 0, 1, 1}));
    g = parseSmiles("[NH4+].[Cl-]");
    EXPECT_EQ(g.formula(), "ClH4N");
    EXPECT_EQ(g.charges, (std::vector<std::optional<int>>{1, -1}));
    g = parseSmiles("[H][C@](F)(Cl)Br");  // explicit H atom folded into its neighbour
    EXPECT_EQ(g.formula(), "CHBrClF");
    EXPECT_EQ(g.stereo.at(0).order.at(0), MolGraph::kHydrogen);
    g = parseSmiles("C1CC%10CC1.C%10");
    EXPECT_EQ(g.size(), 6U);
    EXPECT_TRUE(g.bond(2, 5).has_value());
    for (const char* bad : {"C1CC", "C(C", "C)", "[Xx", "", "Q"}) {
        EXPECT_THROW(parseSmiles(bad), DictionaryError) << bad;
    }
}

TEST(Dictionaries, GraphMatchingRecoversNamesDeterministically) {
    const Component& comp = *findComponent("0G3");
    const MolGraph ref = comp.graph();
    // the permutation of MDATK's test (random.Random(3).shuffle of range(24)): the first match case
    const Value& c = list(at(cases(), "match")).front();
    ASSERT_EQ(str(at(c, "code")), "0G3");
    std::vector<int> perm;
    for (const Value& p : list(at(c, "perm"))) {
        perm.push_back(static_cast<int>(*p.get<std::int64_t>()));
    }
    std::vector<int> inverse(perm.size());
    for (std::size_t k = 0; k < perm.size(); ++k) {
        inverse[static_cast<std::size_t>(perm[k])] = static_cast<int>(k);
    }
    MolGraph tgt;
    for (const int old : perm) {
        tgt.addAtom(ref.elements[static_cast<std::size_t>(old)], ref.charges[static_cast<std::size_t>(old)],
                    ref.hcount[static_cast<std::size_t>(old)]);
    }
    for (std::size_t i = 0; i < ref.size(); ++i) {
        for (const auto& [j, order] : ref.adj[i]) {
            if (static_cast<int>(i) < j) {
                tgt.addBond(inverse[i], inverse[static_cast<std::size_t>(j)], order);
            }
        }
    }
    const auto m = matchGraph(ref, tgt);
    ASSERT_TRUE(m.has_value());
    EXPECT_EQ(m, matchGraph(ref, tgt));
    const std::set<std::string> symmetric = {"OP2", "OP3"};  // either way round
    std::set<std::string> swapped;
    for (std::size_t i = 0; i < ref.size(); ++i) {
        const int back = perm[static_cast<std::size_t>((*m)[i])];
        if (symmetric.count(ref.names[i]) != 0) {
            swapped.insert(ref.names[static_cast<std::size_t>(back)]);
        } else {
            EXPECT_EQ(back, static_cast<int>(i)) << ref.names[i];
        }
    }
    EXPECT_EQ(swapped, symmetric);
    for (std::size_t i = 0; i < tgt.size(); ++i) {  // a bond removed: not isomorphic any more
        if (!tgt.adj[i].empty()) {
            tgt.removeBond(static_cast<int>(i), tgt.adj[i].front().first);
            break;
        }
    }
    EXPECT_FALSE(matchGraph(ref, tgt).has_value());
}

// MDATK checks its SMILES with RDKit (formula, graph, CIP labels S S R S, an embedded conformer's
// chiral volumes); MDBridge checks them against the wwPDB CCD itself: the reference graph is the
// CCD's (reference_tables.hpp), the CCD's own isomeric SMILES are accepted as the same L molecule,
// and MDATK's SMILES agree with the CCD's bond orders except the Kekule structure of adenine.
TEST(Dictionaries, SmilesAgreeWithTheCcd) {
    for (const Component& comp : dictionaryComponents()) {
        const MolGraph ref = comp.graph();
        const MolGraph tgt = parseSmiles(comp.smiles);
        const auto m = matchGraph(ref, tgt);
        ASSERT_TRUE(m.has_value()) << comp.code;
        EXPECT_EQ(allMatches(ref, tgt).size(), comp.terminal == "5" ? 1U : 6U) << comp.code;
        // adenine: C5=C6 N1=C2 N3=C4 in MDATK's SMILES, C5=C4 C6=N1 C2=N3 in the CCD (same molecule)
        EXPECT_EQ(matchScore(ref, tgt, *m), comp.base == "A" ? 6 : 0) << comp.code;
    }
    for (const auto& [parent, smiles] : kCcdSmiles) {
        for (const std::string& s : smiles) {
            for (const std::string& code : {parent, parent + "3"}) {
                EXPECT_EQ(checkSmiles(code, s),
                          (std::vector<std::string>{"OP3 has 1 H (reference 0): different phosphate protonation",
                                                    "OP2 has 1 H (reference 0): different phosphate protonation"}))
                    << code << " " << s;
            }
            EXPECT_EQ(identifySmiles(s), (std::vector<std::string>{parent + "3", parent})) << s;
        }
    }
}

TEST(Dictionaries, RenameRoundTrip) {
    for (const char* code : {"0A5", "0A3", "0C5", "0C3", "0G5", "0G3", "0U5", "0U3", "0G"}) {
        const Component& comp = *findComponent(code);
        for (const std::string& style : {std::string("acedrg"), std::string("elbow")}) {
            const Value& c = finishCase(std::string(code) + " " + style);
            const FinishResult res = finish(c);
            const io::CifDocument out = io::parseCif(res.text);
            ASSERT_GE(out.blocks.size(), 2U);
            EXPECT_EQ(out.blocks[0].name, "comp_list");
            EXPECT_EQ(out.blocks[1].name, "comp_" + std::string(code));
            const io::CifLoop& clist = loopOf(out, "comp_list", "_chem_comp");
            const auto cell = [&](std::string_view item) { return clist.column(item).at(0); };
            EXPECT_EQ(cell("id"), code);
            EXPECT_EQ(cell("three_letter_code"), code);
            EXPECT_EQ(cell("group"), "RNA");
            EXPECT_EQ(cell("name"), comp.name);
            EXPECT_EQ(cell("number_atoms_all"), std::to_string(comp.residueAtoms().size()));
            const io::CifLoop& atoms = loopOf(out, "comp_" + std::string(code), "_chem_comp_atom");
            auto names = atoms.column("atom_id"), want = comp.residueAtoms();
            std::sort(names.begin(), names.end());
            std::sort(want.begin(), want.end());
            EXPECT_EQ(names, want) << code << " " << style;
            const auto ids = atoms.column("comp_id");
            EXPECT_TRUE(std::all_of(ids.begin(), ids.end(), [&](const std::string& id) { return id == code; }));
            const io::CifLoop& angles = loopOf(out, "comp_" + std::string(code), "_chem_comp_angle");
            for (const auto& row : angles.rows) {
                EXPECT_TRUE(comp.leaving.empty() || std::find(row.begin(), row.end(), "OP3") == row.end());
            }
            const auto centres = loopOf(out, "comp_" + std::string(code), "_chem_comp_chir").column("atom_id_centre");
            for (const char* centre : {"C1'", "C2'", "C3'", "C4'"}) {
                EXPECT_NE(std::find(centres.begin(), centres.end(), centre), centres.end());
            }
            const bool keepsOp3 = std::find(want.begin(), want.end(), "OP3") != want.end();
            EXPECT_EQ(std::find(centres.begin(), centres.end(), "P") != centres.end(), keepsOp3) << code;
            const auto plane = loopOf(out, "comp_" + std::string(code), "_chem_comp_plane_atom").column("atom_id");
            EXPECT_EQ(std::find(plane.begin(), plane.end(), "C1'"), plane.end());
            if (style == "acedrg") {
                EXPECT_EQ(loopOf(out, "extra_block", "_extra").rows,
                          (std::vector<std::vector<std::string>>{{"42"}}));
                EXPECT_EQ(loopOf(out, "comp_" + std::string(code), "_pdbx_chem_comp_descriptor").column("descriptor"),
                          std::vector<std::string>{comp.smiles});
                EXPECT_EQ(res.program, "acedrg");
            } else {
                EXPECT_EQ(res.program, "phenix.elbow");
                EXPECT_EQ(res.sourceId, "LIG");
            }
            EXPECT_EQ(res.text.rfind("# " + std::string(code) + ": restraint dictionary from " + res.program, 0), 0U);
            // finishing again changes nothing
            EXPECT_EQ(finishDictionary(res.text, code).text, res.text) << code << " " << style;
            EXPECT_EQ(renameDictionary(str(at(c, "input")), code), res.text);
        }
    }
}

TEST(Dictionaries, RenameGroupNameSourceIdAndKeepNames) {
    const FinishResult named = finish(finishCase("group and name"));
    const auto row = loopOf(io::parseCif(named.text), "comp_list", "_chem_comp").rows.at(0);
    EXPECT_NE(std::find(row.begin(), row.end(), "non-polymer"), row.end());
    EXPECT_NE(std::find(row.begin(), row.end(), "my uridine"), row.end());
    EXPECT_NE(finish(finishCase("two components")).text.find("0U5"), std::string::npos);
    EXPECT_FALSE(finish(finishCase("source id")).mapping.empty());
    EXPECT_THROW(finish(finishCase("missing source id")), DictionaryError);
    const FinishResult kept = finish(finishCase("no reference"));
    EXPECT_NE(kept.notes.at(0).find("XY1 has no reference atom names"), std::string::npos);
    const io::CifDocument atomsDoc = io::parseCif(kept.text);
    const io::CifLoop& atoms = loopOf(atomsDoc, "comp_XY1", "_chem_comp_atom");
    const auto names = atoms.column("atom_id");
    EXPECT_NE(std::find(names.begin(), names.end(), "C1"), names.end());
    for (const std::string& id : atoms.column("comp_id")) {
        EXPECT_EQ(id, "XY1");
    }
    EXPECT_THROW(finish(finishCase("bad code")), DictionaryError);
}

TEST(Dictionaries, RenameRefusesWrongCompoundAndDStereo) {
    const auto message = [](const Value& c) {
        try {
            finish(c);
        } catch (const DictionaryError& e) {
            return std::string(e.what());
        }
        return std::string("no error");
    };
    EXPECT_NE(message(finishCase("wrong compound")).find("does not match 0G5"), std::string::npos);
    EXPECT_NE(message(finishCase("flipped chirality")).find("opposite configuration"), std::string::npos);
    EXPECT_NE(message(finishCase("not a dictionary")).find("restraint dictionary"), std::string::npos);
    const auto notes = finish(finishCase("no chirality rows")).notes;
    EXPECT_TRUE(std::any_of(notes.begin(), notes.end(),
                            [](const std::string& n) { return n.find("not checked") != std::string::npos; }));
}

TEST(Dictionaries, RenamePhosphateProtonationAndOtherHydrogens) {
    const FinishResult res = finish(finishCase("phosphate protonation"));
    const auto names = loopOf(io::parseCif(res.text), "comp_0C3", "_chem_comp_atom").column("atom_id");
    const auto has = [&](const char* n) { return std::find(names.begin(), names.end(), n) != names.end(); };
    EXPECT_FALSE(has("OP3"));
    EXPECT_FALSE(has("HOP3"));
    EXPECT_NE(has("HOP2"), has("HOP1"));  // the protonated non-bridging oxygen keeps its H
    EXPECT_TRUE(std::any_of(res.notes.begin(), res.notes.end(),
                            [](const std::string& n) { return n.find("carries 1 H") != std::string::npos; }));
    try {
        finish(finishCase("base protonation"));
        ADD_FAILURE() << "no error";
    } catch (const DictionaryError& e) {
        EXPECT_NE(std::string(e.what()).find("protonation"), std::string::npos);
    }
    const FinishResult allowed = finish(finishCase("base protonation allowed"));
    const auto n7 = loopOf(io::parseCif(allowed.text), "comp_0G5", "_chem_comp_atom").column("atom_id");
    EXPECT_NE(std::find(n7.begin(), n7.end(), "HN7"), n7.end());
}

// MDATK embeds the components with RDKit; the generator places them with the structure builder
// (L helices, CCD names), mirrored for the D case. CHD-01: guanine's H22 is the amino hydrogen cis
// to N1 (MDATK: H21).
TEST(Dictionaries, HydrogenNamesFollowTheLConventions) {
    const std::map<std::string, std::array<const char*, 4>> amino = {
        {"A", {"H61", "N6", "C6", "N1"}}, {"C", {"H41", "N4", "C4", "N3"}}, {"G", {"H22", "N2", "C2", "N1"}}};
    for (const Component& comp : dictionaryComponents()) {
        for (const bool mirror : {false, true}) {
            for (const char* seed : {"1", "2"}) {
                const std::string label =
                    comp.code + " coords mirror=" + (mirror ? "True" : "False") + " seed=" + seed;
                const FinishResult res = finish(finishCase(label));
                const io::CifDocument atomsDoc = io::parseCif(res.text);
                const io::CifLoop& atoms = loopOf(atomsDoc, "comp_" + comp.code, "_chem_comp_atom");
                std::map<std::string, Vec> got;
                const auto names = atoms.column("atom_id"), xs = atoms.column("x"), ys = atoms.column("y"),
                           zs = atoms.column("z");
                for (std::size_t k = 0; k < names.size(); ++k) {
                    got[names[k]] = {std::stod(xs[k]), std::stod(ys[k]), std::stod(zs[k])};
                }
                const double v1 = volume(got["C5'"], got["O5'"], got["C4'"], got["H5'"]);
                const double v2 = volume(got["C5'"], got["O5'"], got["C4'"], got["H5''"]);
                EXPECT_NE(v1 < 0 && 0 < v2, mirror) << label;  // L: H5' negative, H5'' positive
                if (const auto it = amino.find(comp.base); it != amino.end()) {
                    const auto& [h, n, c, r] = it->second;
                    EXPECT_LT(std::fabs(torsion(got[h], got[n], got[c], got[r])), 90.0) << label;
                }
                const bool warned = std::any_of(res.notes.begin(), res.notes.end(), [](const std::string& x) {
                    return x.find("D configuration") != std::string::npos;
                });
                EXPECT_EQ(warned, mirror) << label;
            }
        }
        // without chirality rows, D coordinates are refused
        try {
            finish(finishCase(comp.code + " mirrored coords without chirality rows"));
            ADD_FAILURE() << comp.code;
        } catch (const DictionaryError& e) {
            EXPECT_NE(std::string(e.what()).find("opposite configuration"), std::string::npos);
        }
    }
}

TEST(Dictionaries, CcdStyleDictionaryIsLeftUnchanged) {
    const FinishResult again = finish(finishCase("CCD-style input"));
    for (const auto& [o, n] : again.mapping) {
        EXPECT_TRUE(n.empty() || o == n) << o << " -> " << n;
    }
    EXPECT_EQ(mappingOf(again, "OP3"), "");
}

// DIC-01: a charge that is no finite number reads as no formal charge (MDATK: OverflowError).
TEST(Dictionaries, NonFiniteChargesAreNoFormalCharges) {
    const Value& c = finishCase("0A5 acedrg");
    const FinishResult base = finish(c);
    std::string text = str(at(c, "input"));
    // the acedrg dictionary's charge column holds 0 for every atom: make one infinite, one NaN
    const std::size_t first = text.find(" OT 0 ");
    ASSERT_NE(first, std::string::npos);
    text.replace(first, 6, " OT inf ");
    const std::size_t second = text.find(" NT 0 ");
    ASSERT_NE(second, std::string::npos);
    text.replace(second, 6, " NT nan ");
    const FinishResult res = finishDictionary(text, "0A5");
    EXPECT_EQ(res.notes, base.notes);
    EXPECT_EQ(res.mapping, base.mapping);
}

TEST(Dictionaries, ResidueCodesAreCheckedAfterStrippingAndUpperCasing) {
    try {
        finishDictionary(str(at(finishCase("0A5 acedrg"), "input")), " 0a5x77 ");
        ADD_FAILURE() << "no error";
    } catch (const DictionaryError& e) {
        EXPECT_EQ(std::string(e.what()), "residue code '0A5X77': use 1-5 letters / digits");
    }
    EXPECT_EQ(&getComponent(" 0a5\t"), findComponent("0A5"));
    // DIC-02: ASCII upper case only; bracket counts and charges up to 10^6
    try {
        finishDictionary(str(at(finishCase("0A5 acedrg"), "input")), "0\xc3\xa9" "5");
        ADD_FAILURE() << "no error";
    } catch (const DictionaryError& e) {
        EXPECT_EQ(std::string(e.what()), "residue code '0\xc3\xa9" "5': use 1-5 letters / digits");
    }
    EXPECT_EQ(parseSmiles("[CH4]").hcount.at(0), 4);
    EXPECT_EQ(parseSmiles("[Fe+1000000]").charges.at(0), 1000000);
    for (const char* big : {"[CH1000001]", "[CH99999999999999999999]", "[Fe+1000001]"}) {
        try {
            parseSmiles(big);
            ADD_FAILURE() << big;
        } catch (const DictionaryError& e) {
            EXPECT_EQ(std::string(e.what()), std::string("SMILES: cannot read the atom ") + big);
        }
    }
}
