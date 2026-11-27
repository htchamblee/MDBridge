// The other pairings and triplexes of the structure builder (src/chem/helix_fit.cpp): the fitted
// builds against the patched MDATK's (tests/parity/make_builder_cases.py -> data/chem/
// builder_fitted_cases.json "rebuilt", improvements_chem.py CH-16) within docs/specs/m6-plan.md D10
// (heavy-atom RMSD <= 0.01 Å, every atom <= 0.05 Å, names / records / header and notes equal but for
// their numbers within one unit of the printed last digit), and the counterparts of MDATK's
// tests/test_builder.py for them: test_hh_duplex_matches_4jrd, test_pairings_build_the_intended_pairs,
// test_triplex_third_strand, test_l_pairings_are_mirror_images,
// test_pairing_mismatches_messages_and_errors, test_pairing_cli_and_planner_step (the build part; the
// planner step is milestone 7's), test_rebuilt_helices_have_no_close_contacts,
// test_reverse_wc_bases_clear_their_own_sugars, test_sugar_steps_keep_c2_away_from_the_next_c5,
// test_closed_strands_keep_the_pair_geometry and test_contact_report_includes_hydrogens_and_close_heavy_atoms.
// Barnaba's Leontis-Westhof codes (MDATK skips them without barnaba) are not checked; tleap's and
// the window's tests have other counterparts (tests/gui/test_tools.cpp).

#include <algorithm>
#include <array>
#include <atomic>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <optional>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <tuple>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "../../src/chem/helix_internal.hpp"
#include "mdbridge/chem/helix.hpp"
#include "mdbridge/chem/pdbtools.hpp"
#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/linalg.hpp"
#include "mdbridge/core/progress.hpp"
#include "mdbridge/io/pdb.hpp"
#include "test_support.hpp"

using namespace mdbridge;
using namespace mdbridge::chem;
using core::Array;
using core::Value;
using test::kData;

namespace {

using Vec = std::array<double, 3>;

const Value& fittedCases() {
    static const Value data = core::parseJson(test::readFile(kData / "chem" / "builder_fitted_cases.json"),
                                              "builder_fitted_cases.json");
    return data;
}

const Array& list(const Value& value) {
    return *value.get<Array>();
}
const Value& at(const Value& value, std::string_view key) {
    const Value* v = value.find(key);
    if (!v) {
        throw std::runtime_error("missing key " + std::string(key));
    }
    return *v;
}
std::string str(const Value& value) {
    return *value.get<std::string>();
}

std::vector<std::string> lines(const std::string& text) {
    std::vector<std::string> out;
    std::istringstream in(text);
    for (std::string line; std::getline(in, line);) {
        out.push_back(line);
    }
    return out;
}

bool isAtomLine(const std::string& line) {
    return line.rfind("ATOM  ", 0) == 0 || line.rfind("HETATM", 0) == 0;
}

// D10 (regenerate_all.py _d10_text_close): equal but for numbers (decimals and integers) within one
// unit of their last printed digit (the spaces before a number may differ with its sign)
bool stringsClose(const std::string& x, const std::string& y) {
    static const std::regex decimal(R"(-?\d+(?:\.\d+)?)");
    const auto split = [](const std::string& s, std::vector<std::string>& texts, std::vector<std::string>& numbers) {
        std::size_t last = 0;
        for (auto it = std::sregex_iterator(s.begin(), s.end(), decimal); it != std::sregex_iterator(); ++it) {
            texts.push_back(s.substr(last, static_cast<std::size_t>(it->position()) - last));
            numbers.push_back(it->str());
            last = static_cast<std::size_t>(it->position() + it->length());
        }
        texts.push_back(s.substr(last));
    };
    std::vector<std::string> ta, tb, na, nb;
    split(x, ta, na);
    split(y, tb, nb);
    if (ta.size() != tb.size()) {
        return false;
    }
    for (std::size_t k = 0; k < ta.size(); ++k) {
        std::string p = ta[k], q = tb[k];
        if (k + 1 < ta.size()) {
            p.erase(p.find_last_not_of(' ') + 1);
            q.erase(q.find_last_not_of(' ') + 1);
        }
        if (p != q) {
            return false;
        }
    }
    const auto digitsOf = [](const std::string& v) {
        const std::size_t dot = v.find('.');
        return dot == std::string::npos ? std::size_t{0} : v.size() - dot - 1;
    };
    for (std::size_t k = 0; k < na.size(); ++k) {
        const std::size_t digits = digitsOf(na[k]);
        if (digitsOf(nb[k]) != digits ||
            std::fabs(std::stod(na[k]) - std::stod(nb[k])) > 1.000001 * std::pow(10.0, -static_cast<double>(digits))) {
            return false;
        }
    }
    return true;
}

HelixOptions optionsOf(const Value& c) {
    HelixOptions o;
    o.pairing = str(at(c, "pairing"));
    if (!at(c, "triplex").isNull()) {
        o.triplex = str(at(c, "triplex"));
    }
    for (const auto& [key, value] : *at(c, "options").get<core::Table>()) {
        if (key == "partner") o.partner = str(value);
        else if (key == "chirality") o.chirality = str(value);
        else if (key == "chains") o.chains = str(value);
        else if (key == "start") o.start = str(value);
        else if (key == "scheme") o.scheme = str(value);
        else if (key == "third") o.third = str(value);
        else if (key == "phosphate5") o.phosphate5 = *value.get<bool>();
        else if (key == "hydrogens") o.hydrogens = *value.get<bool>();
        else throw std::runtime_error("unknown helix option " + key);
    }
    return o;
}

struct Built {
    std::string text;
    std::vector<std::string> notes;
};

struct Request {
    std::string seq;
    HelixOptions options;
};

HelixOptions fitted(std::string_view pairing, std::optional<std::string> triplex = std::nullopt, bool hydrogens = false) {
    HelixOptions o;
    o.pairing = std::string(pairing);
    o.triplex = std::move(triplex);
    o.hydrogens = hydrogens;
    return o;
}

// the builds of a test, at the same time (each closure on one thread; the result does not depend on it)
std::vector<Built> buildAll(std::vector<Request> requests) {
    std::vector<Built> out(requests.size());
    std::atomic<std::size_t> next{0};
    std::vector<std::exception_ptr> errors(requests.size());
    const std::size_t workers = std::max<std::size_t>(1, std::min<std::size_t>(std::thread::hardware_concurrency(),
                                                                                requests.size()));
    std::vector<std::thread> pool;
    for (std::size_t w = 0; w < workers; ++w) {
        pool.emplace_back([&] {
            for (std::size_t k = next++; k < requests.size(); k = next++) {
                try {
                    requests[k].options.threads = 1;
                    out[k].text = buildHelix(requests[k].seq, requests[k].options, &out[k].notes);
                } catch (...) {
                    errors[k] = std::current_exception();
                }
            }
        });
    }
    for (std::thread& t : pool) {
        t.join();
    }
    for (const std::exception_ptr& e : errors) {
        if (e) {
            std::rethrow_exception(e);
        }
    }
    return out;
}

Built buildOne(const std::string& seq, HelixOptions options) {
    return buildAll({{seq, std::move(options)}})[0];
}

// --- model helpers (as tests/test_builder.py's) ------------------------------------------------------

Vec sub(const Vec& a, const Vec& b) {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}
double dot(const Vec& a, const Vec& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vec crossv(const Vec& a, const Vec& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
double dist(const Vec& a, const Vec& b) {
    return std::sqrt(dot(sub(a, b), sub(a, b)));
}
double dihedral(const Vec& a, const Vec& b, const Vec& c, const Vec& d) {
    const Vec b0 = sub(a, b);
    Vec b1 = sub(c, b);
    const Vec b2 = sub(d, c);
    const double n = std::sqrt(dot(b1, b1));
    b1 = {b1[0] / n, b1[1] / n, b1[2] / n};
    const double s0 = dot(b0, b1), s2 = dot(b2, b1);
    const Vec v = {b0[0] - s0 * b1[0], b0[1] - s0 * b1[1], b0[2] - s0 * b1[2]};
    const Vec w = {b2[0] - s2 * b1[0], b2[1] - s2 * b1[1], b2[2] - s2 * b1[2]};
    return std::atan2(dot(crossv(b1, v), w), dot(v, w)) * 180.0 / 3.14159265358979323846;
}

using Coords = std::map<std::string, Vec>;

Coords coordsOf(const io::PdbResidue& r) {
    Coords out;
    for (const auto& [name, xyz] : r.coords()) {
        out[name] = xyz;
    }
    return out;
}

char baseOf(const io::PdbResidue& r) {
    return io::parseResidueName(r.name)->base;
}

// [(residue of strand 1, partner residue)] at each level
std::vector<std::pair<const io::PdbResidue*, const io::PdbResidue*>> strandPairs(const io::PdbStructure& s, bool parallel,
                                                                                std::size_t strand = 1) {
    const auto& one = s.segments[0];
    const auto& other = s.segments[strand];
    const std::size_t n = one.size();
    std::vector<std::pair<const io::PdbResidue*, const io::PdbResidue*>> out;
    for (std::size_t k = 0; k < n; ++k) {
        out.emplace_back(&one[k], parallel ? &other[k] : &other[n - 1 - k]);
    }
    return out;
}

using HbondTable = std::map<std::string, std::vector<std::pair<std::string, std::string>>>;

// donor-acceptor pairs of the intended pairs, written out independently of the builder
const HbondTable& hbonds(std::string_view rule) {
    static const std::map<std::string, HbondTable, std::less<>> tables = {
        {"reverse-wc",
         {{"AU", {{"N1", "N3"}, {"N6", "O2"}}}, {"UA", {{"N3", "N1"}, {"O2", "N6"}}},
          {"GC", {{"N1", "O2"}, {"N2", "N3"}}}, {"CG", {{"O2", "N1"}, {"N3", "N2"}}}}},
        {"hoogsteen", {{"AU", {{"N7", "N3"}, {"N6", "O4"}}}, {"GC", {{"N7", "N3"}, {"O6", "N4"}}}}},
        {"reverse-hoogsteen",
         {{"AU", {{"N7", "N3"}, {"N6", "O2"}}}, {"AA", {{"N7", "N6"}, {"N6", "N1"}}},
          {"GG", {{"N7", "N1"}, {"O6", "N2"}}}}},
        {"hh", {{"AA", {{"N6", "N7"}, {"N7", "N6"}}}}},
        {"wc", {{"AU", {{"N1", "N3"}, {"N6", "O4"}}}, {"GC", {{"N1", "N3"}, {"N2", "O2"}, {"O6", "N4"}}}}},
    };
    return tables.find(rule)->second;
}

bool parallelPairing(std::string_view pairing) {
    return pairing == "reverse-wc" || pairing == "hoogsteen" || pairing == "hh";
}

std::vector<double> o3pDistances(const io::PdbStructure& s) {
    std::vector<double> out;
    for (const auto& seg : s.segments) {
        for (std::size_t k = 0; k + 1 < seg.size(); ++k) {
            out.push_back(dist(coordsOf(seg[k]).at("O3'"), coordsOf(seg[k + 1]).at("P")));
        }
    }
    return out;
}

struct ModelAtom {
    std::size_t residue;
    const io::PdbResidue* res;
    std::string name;
    Vec xyz;
};

std::vector<ModelAtom> allAtoms(const io::PdbStructure& s, bool heavyOnly) {
    std::vector<ModelAtom> out;
    std::size_t k = 0;
    for (const io::PdbResidue* r : s.residues()) {
        for (const io::PdbAtom& a : r->atoms) {
            if (!heavyOnly || a.name[0] != 'H') {
                out.push_back({k, r, a.name, a.xyz()});
            }
        }
        ++k;
    }
    return out;
}

std::vector<std::pair<std::size_t, std::size_t>> backboneLinks(const io::PdbStructure& s, const std::vector<ModelAtom>& atoms) {
    std::map<std::pair<std::size_t, std::string>, std::size_t> pos;
    for (std::size_t i = 0; i < atoms.size(); ++i) {
        pos[{atoms[i].residue, atoms[i].name}] = i;
    }
    std::vector<std::pair<std::size_t, std::size_t>> links;
    std::size_t k = 0;
    for (const auto& seg : s.segments) {
        for (std::size_t a = k; a + 1 < k + seg.size(); ++a) {
            const auto o3 = pos.find({a, "O3'"});
            const auto p = pos.find({a + 1, "P"});
            if (o3 != pos.end() && p != pos.end()) {
                links.emplace_back(o3->second, p->second);
            }
        }
        k += seg.size();
    }
    return links;
}

// bond neighbours: within residues by distance (H: < 1.25 Å, no H-H), and the O3'-P links
std::vector<std::set<std::size_t>> neighbours(const io::PdbStructure& s, const std::vector<ModelAtom>& atoms) {
    std::vector<std::set<std::size_t>> nb(atoms.size());
    for (std::size_t i = 0; i < atoms.size(); ++i) {
        for (std::size_t j = i + 1; j < atoms.size() && atoms[j].residue == atoms[i].residue; ++j) {
            const bool hi = atoms[i].name[0] == 'H', hj = atoms[j].name[0] == 'H';
            const double d = dist(atoms[i].xyz, atoms[j].xyz);
            if (d < 1.75 && !(hi && hj) && (d < 1.25 || !(hi || hj))) {
                nb[i].insert(j);
                nb[j].insert(i);
            }
        }
    }
    for (const auto& [i, j] : backboneLinks(s, atoms)) {
        nb[i].insert(j);
        nb[j].insert(i);
    }
    return nb;
}

// the closest distance of heavy atoms three or more bonds apart (topology_distances' topo >= 3)
double closestHeavy(const io::PdbStructure& s) {
    const std::vector<ModelAtom> atoms = allAtoms(s, true);
    const std::vector<std::set<std::size_t>> nb = neighbours(s, atoms);
    double best = 1e9;
    for (std::size_t i = 0; i < atoms.size(); ++i) {
        std::set<std::size_t> near = {i};
        std::set<std::size_t> frontier = {i};
        for (int step = 0; step < 2; ++step) {
            std::set<std::size_t> next;
            for (const std::size_t a : frontier) {
                for (const std::size_t b : nb[a]) {
                    if (!near.count(b)) next.insert(b);
                }
            }
            near.insert(next.begin(), next.end());
            frontier = next;
        }
        for (std::size_t j = i + 1; j < atoms.size(); ++j) {
            if (!near.count(j)) {
                best = std::min(best, dist(atoms[i].xyz, atoms[j].xyz));
            }
        }
    }
    return best;
}

const std::set<std::string>& polarH() {
    static const std::set<std::string> names = {"H1", "H3", "H21", "H22", "H41", "H42", "H61", "H62", "HO2'", "HO5'", "HO3'"};
    return names;
}

// close_contacts of tests/test_builder.py (independent of the builder): atoms four or more bonds apart
// H-H < 2.0, polar H - N/O < 1.6, other H-heavy < 1.9, heavy < 2.8 (N/O of two residues < 2.6); heavy
// atoms three bonds apart < 2.4
std::vector<std::string> closeContactsOf(const io::PdbStructure& s) {
    const std::vector<ModelAtom> atoms = allAtoms(s, false);
    const std::vector<std::set<std::size_t>> nb = neighbours(s, atoms);
    std::vector<std::string> out;
    for (std::size_t i = 0; i < atoms.size(); ++i) {
        std::set<std::size_t> two;
        for (const std::size_t c : nb[i]) two.insert(nb[c].begin(), nb[c].end());
        std::set<std::size_t> three;
        for (const std::size_t c : two) three.insert(nb[c].begin(), nb[c].end());
        for (std::size_t j = i + 1; j < atoms.size(); ++j) {
            const double d = dist(atoms[i].xyz, atoms[j].xyz);
            if (d > 2.8 || nb[i].count(j) || two.count(j)) continue;
            const bool hi = atoms[i].name[0] == 'H', hj = atoms[j].name[0] == 'H';
            double limit = 0.0;
            if (three.count(j)) {
                limit = hi || hj ? 0.0 : 2.4;
            } else if (hi && hj) {
                limit = 2.0;
            } else if (hi || hj) {
                const std::string& h = hi ? atoms[i].name : atoms[j].name;
                const std::string& x = hi ? atoms[j].name : atoms[i].name;
                limit = polarH().count(h) && (x[0] == 'N' || x[0] == 'O') ? 1.6 : 1.9;
            } else {
                const bool no = (atoms[i].name[0] == 'N' || atoms[i].name[0] == 'O') &&
                                (atoms[j].name[0] == 'N' || atoms[j].name[0] == 'O');
                limit = no && atoms[i].residue != atoms[j].residue ? 2.6 : 2.8;
            }
            if (d < limit) {
                out.push_back(atoms[i].res->label() + " " + atoms[i].name + " - " + atoms[j].res->label() + " " +
                              atoms[j].name + " " + std::to_string(d));
            }
        }
    }
    std::sort(out.begin(), out.end());
    return out;
}

bool startsWith(const std::string& s, std::string_view p) {
    return s.rfind(p, 0) == 0;
}
bool anyNote(const std::vector<std::string>& notes, const std::function<bool(const std::string&)>& f) {
    return std::any_of(notes.begin(), notes.end(), f);
}

// check_model: H-bonds 2.6-3.3 Å, O3'-P links 1.5-1.7 Å, no heavy-atom contacts < 2.4 Å three or more
// bonds apart, the report with hydrogens
void checkModel(const Built& b, const std::vector<std::pair<const io::PdbResidue*, const io::PdbResidue*>>& pairs,
                const HbondTable& table, const std::string& where) {
    const io::PdbStructure s = io::readPdbText(b.text);
    for (const auto& [x, y] : pairs) {
        const Coords ca = coordsOf(*x), cb = coordsOf(*y);
        const std::string key = std::string(1, baseOf(*x)) + baseOf(*y);
        const auto it = table.find(key);
        ASSERT_NE(it, table.end()) << where << " " << key;
        for (const auto& [p, q] : it->second) {
            const double d = dist(ca.at(p), cb.at(q));
            EXPECT_TRUE(2.6 <= d && d <= 3.3) << where << " " << x->label() << " " << y->label() << " " << p << " " << q << " " << d;
        }
    }
    const std::vector<double> d = o3pDistances(s);
    EXPECT_GE(*std::min_element(d.begin(), d.end()), 1.5) << where;
    EXPECT_LE(*std::max_element(d.begin(), d.end()), 1.7) << where;
    EXPECT_TRUE(anyNote(b.notes, [](const std::string& n) { return startsWith(n, "O3'-P links:"); })) << where;
    EXPECT_GE(closestHeavy(s), 2.4) << where;
    EXPECT_TRUE(anyNote(b.notes, [](const std::string& n) {
        return startsWith(n, "No close contacts (") && n.find("hydrogens included") != std::string::npos;
    })) << where;
}

// --- the fitted builds against the patched MDATK's (D10) ----------------------------------------------

class FittedCase : public ::testing::TestWithParam<std::size_t> {};

TEST_P(FittedCase, WithinD10OfMdatk) {
    const Value& c = list(at(fittedCases(), "rebuilt")).at(GetParam());
    const std::string seq = str(at(c, "sequence"));
    const std::string where = "rebuilt[" + std::to_string(GetParam()) + "] " + str(at(c, "pairing")) + " " + seq;
    HelixOptions options = optionsOf(c);
    std::vector<std::string> notes;
    if (const Value* error = c.find("error")) {
        try {
            buildHelix(seq, options, &notes);
            FAIL() << where << ": expected " << str(at(*error, "message"));
        } catch (const core::ValueError& e) {
            EXPECT_EQ(std::string(e.pythonName()), str(at(*error, "type"))) << where;
            EXPECT_EQ(std::string(e.what()), str(at(*error, "message"))) << where;
        }
        return;
    }
    const std::string text = buildHelix(seq, options, &notes);
    const std::vector<std::string> got = lines(text), want = lines(str(at(c, "text")));
    ASSERT_EQ(got.size(), want.size()) << where;
    double sum = 0.0, worst = 0.0, worstAll = 0.0;
    std::size_t heavy = 0;
    for (std::size_t k = 0; k < got.size(); ++k) {
        if (!isAtomLine(want[k])) {
            EXPECT_TRUE(stringsClose(got[k], want[k])) << where << "\n  " << got[k] << "\n  " << want[k];
            continue;
        }
        ASSERT_TRUE(isAtomLine(got[k])) << where << ": " << got[k];
        EXPECT_EQ(got[k].substr(0, 30), want[k].substr(0, 30)) << where;
        EXPECT_EQ(got[k].substr(54), want[k].substr(54)) << where;
        Vec a{}, b{};
        for (std::size_t i = 0; i < 3; ++i) {
            a[i] = std::stod(got[k].substr(30 + 8 * i, 8));
            b[i] = std::stod(want[k].substr(30 + 8 * i, 8));
        }
        const double d = dist(a, b);
        worstAll = std::max(worstAll, d);
        const std::string name = got[k].substr(12, 4);
        if (name[name.find_first_not_of(' ')] != 'H') {
            sum += d * d;
            worst = std::max(worst, d);
            ++heavy;
        }
    }
    const double rmsd = std::sqrt(sum / static_cast<double>(heavy));
    EXPECT_LE(rmsd, 0.01) << where << ": heavy-atom RMSD";
    EXPECT_LE(worstAll, 0.05) << where << ": largest atom deviation";
    std::cout << where << ": heavy-atom RMSD " << rmsd << " Å, largest heavy " << worst << " Å, any atom " << worstAll
              << " Å\n";
    const Array& wantNotes = list(at(c, "notes"));
    ASSERT_EQ(notes.size(), wantNotes.size()) << where;
    for (std::size_t k = 0; k < notes.size(); ++k) {
        EXPECT_TRUE(stringsClose(notes[k], str(wantNotes[k]))) << where << "\n  " << notes[k] << "\n  " << str(wantNotes[k]);
    }
}

INSTANTIATE_TEST_SUITE_P(Rebuilt, FittedCase, ::testing::Range<std::size_t>(0, 24));

TEST(FittedCaseData, AllRebuiltCasesAreTested) {
    EXPECT_EQ(list(at(fittedCases(), "rebuilt")).size(), 24u);
}

// --- MDATK's tests -------------------------------------------------------------------------------

// test_hh_duplex_matches_4jrd: the hh A10 duplex against PDB 4JRD (data/sources/pdb/4jrd_atoms.tsv)
TEST(BuilderMdatk, HhDuplexMatches4jrd) {
    HelixOptions o = fitted("hh");
    o.phosphate5 = true;
    const Built b = buildOne(std::string(10, 'A'), o);
    const io::PdbStructure s = io::readPdbText(b.text);
    checkModel(b, strandPairs(s, true), hbonds("hh"), "hh A10");
    // 4JRD atoms (chain, residue number) -> name -> xyz
    std::map<std::pair<std::string, int>, Coords> ref;
    const std::vector<std::string> rows = lines(test::readFile(kData / ".." / ".." / "data" / "sources" / "pdb" / "4jrd_atoms.tsv"));
    std::vector<std::string> header;
    for (std::size_t k = 0; k < rows.size(); ++k) {
        std::vector<std::string> f;
        std::istringstream in(rows[k]);
        for (std::string cell; std::getline(in, cell, '\t');) f.push_back(cell);
        if (k == 0) {
            header = f;
            continue;
        }
        const auto col = [&](std::string_view name) {
            return f[static_cast<std::size_t>(std::find(header.begin(), header.end(), name) - header.begin())];
        };
        if (col("label_comp_id") != "A" || col("pdbx_PDB_model_num") != "1" || (col("label_alt_id") != "." && col("label_alt_id") != "A")) {
            continue;
        }
        ref[{col("auth_asym_id"), std::stoi(col("auth_seq_id"))}][col("label_atom_id")] =
            Vec{std::stod(col("Cartn_x")), std::stod(col("Cartn_y")), std::stod(col("Cartn_z"))};
    }
    std::map<std::pair<std::string, int>, Coords> mine;
    for (const io::PdbResidue* r : s.residues()) {
        mine[{r->chain, static_cast<int>(r->resseq)}] = coordsOf(*r);
    }
    std::vector<std::pair<std::pair<std::string, int>, std::pair<std::string, int>>> keys;
    for (int k = 1; k <= 10; ++k) keys.push_back({{"A", k}, {"A", k + 1}});
    for (int k = 1; k <= 10; ++k) keys.push_back({{"B", k}, {"B", k}});
    for (const std::vector<std::string>& names :
         {std::vector<std::string>{"P"}, std::vector<std::string>{"P", "C4'", "C1'", "N9", "N1"}}) {
        std::vector<Vec> x, y;
        for (const auto& [m, r] : keys) {
            for (const std::string& a : names) {
                if (ref[r].count(a) && mine[m].count(a)) {
                    x.push_back(mine[m][a]);
                    y.push_back(ref[r][a]);
                }
            }
        }
        ASSERT_GE(x.size(), 15u);
        // Kabsch (unweighted) with core::svd3x3
        Vec mx{}, my{};
        for (std::size_t k = 0; k < x.size(); ++k)
            for (std::size_t c = 0; c < 3; ++c) {
                mx[c] += x[k][c] / static_cast<double>(x.size());
                my[c] += y[k][c] / static_cast<double>(x.size());
            }
        core::Mat3 h{};
        for (std::size_t k = 0; k < x.size(); ++k)
            for (std::size_t p = 0; p < 3; ++p)
                for (std::size_t q = 0; q < 3; ++q) h[p][q] += (x[k][p] - mx[p]) * (y[k][q] - my[q]);
        const core::Svd3 svd = core::svd3x3(h);
        core::Mat3 r{};
        for (std::size_t i = 0; i < 3; ++i)
            for (std::size_t j = 0; j < 3; ++j)
                for (std::size_t k = 0; k < 3; ++k) r[i][j] += svd.vt[k][i] * svd.u[j][k];
        const double det = core::det3(r);
        if (det < 0) {
            for (std::size_t i = 0; i < 3; ++i)
                for (std::size_t j = 0; j < 3; ++j) r[i][j] -= 2.0 * svd.vt[2][i] * svd.u[j][2];
        }
        double sum = 0.0;
        for (std::size_t k = 0; k < x.size(); ++k) {
            const Vec d = sub(x[k], mx);
            Vec moved{};
            for (std::size_t i = 0; i < 3; ++i) moved[i] = r[i][0] * d[0] + r[i][1] * d[1] + r[i][2] * d[2] + my[i];
            sum += dot(sub(moved, y[k]), sub(moved, y[k]));
        }
        const double rmsd = std::sqrt(sum / static_cast<double>(x.size()));
        EXPECT_LT(rmsd, 1.5) << names.size() << " atom names";
    }
    EXPECT_NE(b.text.find("REMARK 951 NUCLEOTIDE OF PDB 4JRD"), std::string::npos);
}

// test_pairings_build_the_intended_pairs
TEST(BuilderMdatk, PairingsBuildTheIntendedPairs) {
    const std::vector<std::pair<std::string, std::string>> cases = {
        {"reverse-wc", "GAUCAGUAGC"}, {"hoogsteen", "AAGAGAAGGA"}, {"reverse-hoogsteen", "AAGAAGGAAG"}, {"hh", "AAAAAAAA"}};
    std::vector<Request> requests;
    for (const auto& [pairing, seq] : cases) requests.push_back({seq, fitted(pairing)});
    const std::vector<Built> built = buildAll(requests);
    for (std::size_t k = 0; k < cases.size(); ++k) {
        const auto& [pairing, seq] = cases[k];
        const io::PdbStructure s = io::readPdbText(built[k].text);
        checkModel(built[k], strandPairs(s, parallelPairing(pairing)), hbonds(pairing), pairing);
        std::string partner;
        for (const io::PdbResidue& r : s.segments[1]) partner += r.name;
        EXPECT_EQ(partner, autoPartner(seq, pairing)) << pairing;
        EXPECT_TRUE(anyNote(built[k].notes, [&](const std::string& n) { return startsWith(n, "Pairing " + pairing + ":"); }));
        const std::string lw = pairing == "reverse-wc" ? "tWW" : pairing == "hoogsteen" ? "cHW" : pairing == "hh" ? "tHH" : "tHW";
        EXPECT_NE(built[k].text.find("PAIRING " + lw), std::string::npos) << pairing;
        EXPECT_FALSE(anyNote(built[k].notes, [](const std::string& n) { return n.find("mismatch") != std::string::npos; }));
    }
}

// test_triplex_third_strand
TEST(BuilderMdatk, TriplexThirdStrand) {
    const std::string seq = "AAGAGAAGGA";
    const std::vector<Built> built = buildAll({{seq, fitted("wc", "hoogsteen")}, {seq, fitted("wc", "reverse-hoogsteen")}});
    for (const std::string motif : {"hoogsteen", "reverse-hoogsteen"}) {
        const Built& b = built[motif == "hoogsteen" ? 0 : 1];
        const io::PdbStructure s = io::readPdbText(b.text);
        EXPECT_EQ(s.chainIds(), (std::vector<std::string>{"A", "B", "C"}));
        ASSERT_EQ(s.segments.size(), 3u);
        for (const auto& seg : s.segments) EXPECT_EQ(seg.size(), 10u);
        std::string third;
        for (const io::PdbResidue& r : s.segments[2]) third += r.name;
        EXPECT_EQ(third, resolveThird(seq, "auto", motif).sequence);
        EXPECT_EQ(third, motif == "hoogsteen" ? "UUCUCUUCCU" : "AGGAAGAGAA");
        checkModel(b, strandPairs(s, false, 1), hbonds("wc"), motif + " duplex");
        checkModel(b, strandPairs(s, motif == "hoogsteen", 2), hbonds(motif), motif + " triples");
        if (motif == "hoogsteen") {
            EXPECT_TRUE(anyNote(b.notes, [](const std::string& n) {
                return n.find("G·C+ Hoogsteen pair(s): the cytosine is written as plain C") != std::string::npos;
            }));
        }
        EXPECT_TRUE(anyNote(b.notes, [](const std::string& n) { return startsWith(n, "Triplex, "); }));
        EXPECT_NE(b.text.find("THIRD STRAND"), std::string::npos);
    }
}

// test_l_pairings_are_mirror_images
TEST(BuilderMdatk, LPairingsAreMirrorImages) {
    const std::string seq = "AAGAGA";
    const std::vector<HelixOptions> kinds = {fitted("hoogsteen"), fitted("reverse-wc"), fitted("wc", "reverse-hoogsteen")};
    std::vector<Request> requests;
    for (const HelixOptions& o : kinds) {
        requests.push_back({seq, o});
        HelixOptions l = o;
        l.chirality = "L";
        requests.push_back({seq, l});
    }
    const std::vector<Built> built = buildAll(requests);
    for (std::size_t k = 0; k < kinds.size(); ++k) {
        const io::PdbStructure d = io::readPdbText(built[2 * k].text);
        const io::PdbStructure lh = io::readPdbText(built[2 * k + 1].text);
        const auto rd = d.residues(), rl = lh.residues();
        ASSERT_EQ(rd.size(), rl.size());
        for (std::size_t r = 0; r < rd.size(); ++r) {
            const Coords cd = coordsOf(*rd[r]), cl = coordsOf(*rl[r]);
            ASSERT_EQ(cd.size(), cl.size());
            for (const auto& [name, x] : cd) {
                ASSERT_TRUE(cl.count(name)) << name;
                EXPECT_NEAR(-x[0], cl.at(name)[0], 2e-3);
                EXPECT_NEAR(x[1], cl.at(name)[1], 2e-3);
                EXPECT_NEAR(x[2], cl.at(name)[2], 2e-3);
            }
        }
        const std::string n = std::to_string(rl.size());
        io::PdbStructure lcopy = lh, dcopy = d;
        EXPECT_EQ(assignChirality(lcopy).at(0), "Chirality (sugar geometry): " + n + " L and 0 D nucleotides");
        EXPECT_EQ(assignChirality(dcopy).at(0), "Chirality (sugar geometry): 0 L and " + n + " D nucleotides");
        for (const io::PdbResidue* r : rl) {
            EXPECT_TRUE(r->name == "0A" || r->name == "0C" || r->name == "0G" || r->name == "0U") << r->name;
        }
    }
}

// test_pairing_mismatches_messages_and_errors (the build parts; the messages: test_helix.cpp's
// check_helix cases)
TEST(BuilderMdatk, PairingMismatchesMessagesAndErrors) {
    HelixOptions hoog = fitted("hoogsteen");
    hoog.partner = "UUUU";
    HelixOptions pq = fitted("wc", "hoogsteen");
    pq.chains = "P,Q";
    pq.start = "5,10";
    const std::vector<Built> built = buildAll({{"AAGA", hoog}, {"AAUAG", fitted("wc", "hoogsteen")}, {"AAGAG", pq}});
    EXPECT_TRUE(anyNote(built[0].notes, [](const std::string& n) { return startsWith(n, "WARNING: mismatch A:3 · B:3: no cHW pair"); }));
    EXPECT_EQ(io::readPdbText(built[0].text).segments.size(), 2u);
    EXPECT_TRUE(anyNote(built[1].notes, [](const std::string& n) {
        return n.find("a strand-1 pyrimidine has no purine Hoogsteen edge") != std::string::npos;
    }));
    const io::PdbStructure s = io::readPdbText(built[2].text);
    EXPECT_EQ(s.chainIds(), (std::vector<std::string>{"P", "Q", "C"}));
    EXPECT_EQ(s.segments[0][0].resseq, 5);
    EXPECT_EQ(s.segments[1][0].resseq, 10);
    EXPECT_EQ(s.segments[2][0].resseq, 1);
    const auto fails = [](HelixOptions o, const std::string& match) {
        try {
            buildHelix("AAGAG", o);
        } catch (const BuildError& e) {
            return std::string(e.what()).find(match) != std::string::npos;
        }
        return false;
    };
    HelixOptions o = fitted("hh", "hoogsteen");
    EXPECT_TRUE(fails(o, "Watson-Crick duplex"));
    o = fitted("hoogsteen");
    o.partner = "none";
    EXPECT_TRUE(fails(o, "single strand"));
    o = fitted("wc", "hoogsteen");
    o.third = "UU";
    EXPECT_TRUE(fails(o, "same length"));
    EXPECT_TRUE(fails(fitted("zipper"), "unknown pairing"));
    EXPECT_TRUE(fails(fitted("wc", "tetraplex"), "unknown triplex"));
    o = fitted("wc", "hoogsteen");
    o.chains = "A,B,A";
    EXPECT_TRUE(fails(o, "three strands"));
    EXPECT_EQ(helixKind(), "A-form");
    EXPECT_EQ(helixKind("hh"), "parallel tHH");
    EXPECT_EQ(helixKind("wc", "reverse-hoogsteen"), "reverse-hoogsteen triplex");
}

// test_pairing_cli_and_planner_step: the builds of its command lines (the command line itself:
// tests/data/cli cases; the planner step: milestone 7)
TEST(BuilderMdatk, PairingCliBuilds) {
    HelixOptions l = fitted("hoogsteen");
    l.chirality = "L";
    HelixOptions tri = fitted("wc", "reverse-hoogsteen");
    tri.third = "GAGAA";
    tri.chains = "A,B,T";
    tri.start = "1,1,7";
    const std::vector<Built> built = buildAll({{"AAGAG", l}, {"AAGAG", tri}});
    std::vector<std::string> names;
    const io::PdbStructure hoog = io::readPdbText(built[0].text);
    for (const io::PdbResidue& r : hoog.segments[1]) names.push_back(r.name);
    EXPECT_EQ(names, (std::vector<std::string>{"0U", "0U", "0C", "0U", "0C"}));
    const io::PdbStructure s = io::readPdbText(built[1].text);
    EXPECT_EQ(s.chainIds(), (std::vector<std::string>{"A", "B", "T"}));
    EXPECT_EQ(s.segments[2][0].resseq, 7);
    EXPECT_THROW(buildHelix("AAGAG", fitted("hh", "hoogsteen")), BuildError);
}

struct Mixed {
    std::string pairing;
    std::optional<std::string> motif;
    std::string seq;
};
// MIXED: other pairs than the reference pair, strand-1 pyrimidines, mismatches
const std::vector<Mixed>& mixed() {
    static const std::vector<Mixed> cases = {
        {"reverse-wc", {}, "GAGGAUCA"}, {"reverse-wc", {}, "AUGCAU"}, {"hoogsteen", {}, "GAGGAUCA"},
        {"reverse-hoogsteen", {}, "GAAGAG"}, {"reverse-hoogsteen", {}, "AUGCAU"}, {"hh", {}, "GAGGAUCA"},
        {"wc", "hoogsteen", "GAGGAUCA"}, {"wc", "reverse-hoogsteen", "AUGCAU"}};
    return cases;
}

std::vector<Built> buildWithHydrogens(const std::vector<Mixed>& cases) {
    std::vector<Request> requests;
    for (const Mixed& m : cases) requests.push_back({m.seq, fitted(m.pairing, m.motif, true)});
    return buildAll(requests);
}

// test_rebuilt_helices_have_no_close_contacts
TEST(BuilderMdatk, RebuiltHelicesHaveNoCloseContacts) {
    const std::vector<Built> built = buildWithHydrogens(mixed());
    for (std::size_t k = 0; k < built.size(); ++k) {
        const std::string where = mixed()[k].pairing + " " + mixed()[k].seq;
        const io::PdbStructure s = io::readPdbText(built[k].text);
        EXPECT_EQ(closeContactsOf(s), std::vector<std::string>{}) << where;
        const std::vector<double> d = o3pDistances(s);
        EXPECT_GE(*std::min_element(d.begin(), d.end()), 1.5) << where;
        EXPECT_LE(*std::max_element(d.begin(), d.end()), 1.7) << where;
        std::vector<std::string> report;
        for (const std::string& n : built[k].notes) {
            std::string lower = n;
            std::transform(lower.begin(), lower.end(), lower.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
            if (lower.find("close contact") != std::string::npos) report.push_back(n);
        }
        ASSERT_EQ(report.size(), 1u) << where;
        EXPECT_NE(report[0].find("hydrogens included"), std::string::npos) << where;
        for (const std::string_view k2 : {"(H-H <", "(H-heavy <", "(H-bond H <"}) {
            EXPECT_EQ(report[0].find(k2), std::string::npos) << where << ": " << report[0];
        }
    }
}

// test_reverse_wc_bases_clear_their_own_sugars
TEST(BuilderMdatk, ReverseWcBasesClearTheirOwnSugars) {
    const std::vector<std::string> seqs = {"UCUCAG", "GAGGAUCA", "AUGCAU"};
    std::vector<Mixed> cases;
    for (const std::string& seq : seqs) cases.push_back({"reverse-wc", {}, seq});
    const std::vector<Built> built = buildWithHydrogens(cases);
    for (std::size_t k = 0; k < seqs.size(); ++k) {
        const io::PdbStructure model = io::readPdbText(built[k].text);
        for (const io::PdbResidue* r : model.residues()) {
            const Coords c = coordsOf(*r);
            const bool purine = baseOf(*r) == 'A' || baseOf(*r) == 'G';
            const std::string ring = purine ? "C8" : "C6", h = purine ? "H8" : "H6", n = purine ? "N9" : "N1",
                              x = purine ? "C4" : "C2";
            const double chi = dihedral(c.at("O4'"), c.at("C1'"), c.at(n), c.at(x));
            EXPECT_TRUE(chi <= -135.0 || chi >= 170.0) << r->label() << " " << chi;
            EXPECT_GE(dist(c.at(ring), c.at("C3'")), 3.1) << r->label();
            EXPECT_GE(dist(c.at(h), c.at("H3'")), 2.0) << r->label();
            EXPECT_GE(dist(c.at(h), c.at("H2'")), 2.0) << r->label();
        }
        const auto note = std::find_if(built[k].notes.begin(), built[k].notes.end(),
                                       [](const std::string& s) { return startsWith(s, "Backbone torsions"); });
        ASSERT_NE(note, built[k].notes.end());
        std::vector<double> chis;
        static const std::regex chiRe(R"(chi (-?\d+))");
        for (auto it = std::sregex_iterator(note->begin(), note->end(), chiRe); it != std::sregex_iterator(); ++it) {
            chis.push_back(std::stod((*it)[1].str()));
        }
        for (std::size_t i = 1; i < chis.size(); ++i) {  // the first one is the 4JRD reference
            EXPECT_TRUE(chis[i] <= -135 || chis[i] >= 170) << seqs[k] << " " << chis[i];
        }
    }
}

// test_sugar_steps_keep_c2_away_from_the_next_c5
TEST(BuilderMdatk, SugarStepsKeepC2AwayFromTheNextC5) {
    const std::vector<Mixed> cases = {{"reverse-hoogsteen", {}, "AAGAGG"}, {"reverse-hoogsteen", {}, "GAAGAG"},
                                      {"wc", "hoogsteen", "GAGGAG"}};
    const std::vector<Built> built = buildWithHydrogens(cases);
    for (std::size_t k = 0; k < cases.size(); ++k) {
        const io::PdbStructure model = io::readPdbText(built[k].text);
        for (const auto& seg : model.segments) {
            for (std::size_t i = 0; i + 1 < seg.size(); ++i) {
                const Coords ca = coordsOf(seg[i]), cb = coordsOf(seg[i + 1]);
                const std::string where = cases[k].seq + " " + seg[i].label() + " " + seg[i + 1].label();
                EXPECT_GE(dist(ca.at("C2'"), cb.at("C5'")), 3.1) << where;
                EXPECT_GE(dist(ca.at("O2'"), cb.at("C5'")), 3.0) << where;
                for (const std::string hb : {"H5'", "H5''"}) {
                    EXPECT_GE(dist(ca.at("H2'"), cb.at(hb)), 2.0) << where << " " << hb;
                    EXPECT_GE(dist(ca.at("C2'"), cb.at(hb)), 2.4) << where << " " << hb;
                }
            }
        }
    }
}

// the base-plane normal of a residue's base heavy atoms (np.linalg.svd(ring - mean)[2][2])
Vec baseNormal(const Coords& c, char base) {
    std::vector<Vec> ring;
    for (const std::string_view n : baseAtomNames(base)) ring.push_back(c.at(std::string(n)));
    Vec m{};
    for (const Vec& v : ring)
        for (std::size_t i = 0; i < 3; ++i) m[i] += v[i] / static_cast<double>(ring.size());
    for (Vec& v : ring) v = sub(v, m);
    return core::svdRows3(ring).vt[2];
}

// test_closed_strands_keep_the_pair_geometry
TEST(BuilderMdatk, ClosedStrandsKeepThePairGeometry) {
    const std::vector<Mixed> cases = {{"reverse-wc", {}, "AUGCAU"}, {"reverse-wc", {}, "GAGGAUCA"},
                                      {"hoogsteen", {}, "GAGGAUCA"}, {"reverse-hoogsteen", {}, "GAGGAUCA"},
                                      {"wc", "hoogsteen", "GAGGAUCA"}, {"wc", "reverse-hoogsteen", "GAAGAG"}};
    const std::vector<Built> built = buildWithHydrogens(cases);
    for (std::size_t k = 0; k < cases.size(); ++k) {
        const Mixed& m = cases[k];
        const io::PdbStructure s = io::readPdbText(built[k].text);
        std::vector<std::pair<std::vector<std::pair<const io::PdbResidue*, const io::PdbResidue*>>, const HbondTable*>> groups;
        if (!m.motif) {
            groups.push_back({strandPairs(s, parallelPairing(m.pairing)), &hbonds(m.pairing)});
        } else {
            groups.push_back({strandPairs(s, false, 1), &hbonds("wc")});
            groups.push_back({strandPairs(s, *m.motif == "hoogsteen", 2), &hbonds(*m.motif)});
        }
        int checked = 0;
        for (const auto& [pairs, table] : groups) {
            for (const auto& [a, b] : pairs) {
                const std::string key = std::string(1, baseOf(*a)) + baseOf(*b);
                const auto it = table->find(key);
                if (it == table->end()) continue;
                const Coords ca = coordsOf(*a), cb = coordsOf(*b);
                const std::string where = m.pairing + " " + m.seq + " " + a->label() + " " + b->label();
                for (const auto& [x, y] : it->second) {
                    const double d = dist(ca.at(x), cb.at(y));
                    EXPECT_TRUE(2.8 <= d && d <= 3.0) << where << " " << x << " " << y << " " << d;
                    for (const auto& [don, acc, cd, cc] :
                         {std::make_tuple(x, y, &ca, &cb), std::make_tuple(y, x, &cb, &ca)}) {
                        std::vector<std::string> hs;
                        for (const auto& [h, v] : *cd) {
                            if (h[0] == 'H' && dist(v, cd->at(don)) < 1.1) hs.push_back(h);
                        }
                        if (!hs.empty() && (acc[0] == 'N' || acc[0] == 'O') && don[0] == 'N') {
                            const std::string h = *std::min_element(hs.begin(), hs.end(), [&](const std::string& p, const std::string& q) {
                                return dist(cd->at(p), cc->at(acc)) < dist(cd->at(q), cc->at(acc));
                            });
                            if (dist(cd->at(h), cc->at(acc)) < 2.6) {
                                const Vec u = sub(cd->at(don), cd->at(h)), v = sub(cc->at(acc), cd->at(h));
                                const double angle = std::acos(dot(u, v) / std::sqrt(dot(u, u)) / std::sqrt(dot(v, v))) *
                                                     180.0 / 3.14159265358979323846;
                                EXPECT_LE(dist(cd->at(h), cc->at(acc)), 2.05) << where << " " << h << " " << acc;
                                EXPECT_GE(angle, 150.0) << where << " " << h << " " << acc;
                            }
                        }
                    }
                }
                const Vec n1 = baseNormal(ca, key[0]), n2 = baseNormal(cb, key[1]);
                EXPECT_LE(std::acos(std::min(1.0, std::fabs(dot(n1, n2)))) * 180.0 / 3.14159265358979323846, 10.0) << where;
                ++checked;
            }
        }
        EXPECT_GE(checked, 5) << m.pairing << " " << m.seq;
        EXPECT_EQ(closeContactsOf(s), std::vector<std::string>{}) << m.pairing << " " << m.seq;
    }
}

// test_contact_report_includes_hydrogens_and_close_heavy_atoms
TEST(BuilderMdatk, ContactReportIncludesHydrogensAndCloseHeavyAtoms) {
    HelixOptions wc;
    wc.partner = "none";
    wc.hydrogens = true;
    const io::PdbStructure s = io::readPdbText(buildHelix("GGCAUC", wc));
    const auto residues = s.residues();
    std::vector<std::string> labels;
    std::vector<char> bases;
    std::vector<helix_detail::ResAtoms> resAtoms;
    for (const io::PdbResidue* r : residues) {
        labels.push_back(r->label());
        bases.push_back(baseOf(*r));
        helix_detail::ResAtoms atoms;
        for (const io::PdbAtom& a : r->atoms) atoms.emplace_back(a.name, a.xyz());
        resAtoms.push_back(std::move(atoms));
    }
    std::vector<std::vector<std::size_t>> segments(1);
    for (std::size_t k = 0; k < residues.size(); ++k) segments[0].push_back(k);
    std::vector<helix_detail::Contact> found = helix_detail::closeContacts(resAtoms, segments, bases, {});
    const auto pairOf = [](const helix_detail::Contact& c) { return std::set<std::string>{c.a.second, c.b.second}; };
    EXPECT_TRUE(std::any_of(found.begin(), found.end(), [&](const helix_detail::Contact& c) {
        const std::set<std::string> p = pairOf(c);
        return (p == std::set<std::string>{"H2'", "H6"} || p == std::set<std::string>{"H2'", "H8"}) &&
               startsWith(c.what, "H-H") && c.d < 2.0;
    }));
    const std::string note = helix_detail::geometryNotes(labels, resAtoms, segments, {}, bases).back();
    EXPECT_TRUE(startsWith(note, "WARNING:")) << note;
    EXPECT_NE(note.find("close contact(s)"), std::string::npos);
    EXPECT_NE(note.find("(H-H < 2.0)"), std::string::npos);
    // residue 2 moved (rigidly) so that its C5' is 2.7 Å from C2' of residue 1: reported (heavy < 3.0)
    Vec c2{}, c5{};
    for (const auto& [n, x] : resAtoms[0]) if (n == "C2'") c2 = x;
    for (const auto& [n, x] : resAtoms[1]) if (n == "C5'") c5 = x;
    const Vec u = sub(c2, c5);
    const double f = 1.0 - 2.7 / std::sqrt(dot(u, u));
    for (auto& [n, x] : resAtoms[1]) x = {x[0] + u[0] * f, x[1] + u[1] * f, x[2] + u[2] * f};
    found = helix_detail::closeContacts(resAtoms, segments, bases, {});
    EXPECT_TRUE(std::any_of(found.begin(), found.end(), [](const helix_detail::Contact& c) {
        return c.a.second == "C2'" && c.b.second == "C5'" && startsWith(c.what, "heavy") && std::fabs(c.d - 2.7) < 1e-6;
    }));
    // a heavy-atom model gets hydrogens for the check: the same report as with hydrogens
    const std::vector<Mixed> two(mixed().begin(), mixed().begin() + 2);
    std::vector<Request> requests;
    for (const Mixed& m : two) {
        requests.push_back({m.seq, fitted(m.pairing, m.motif, true)});
        requests.push_back({m.seq, fitted(m.pairing, m.motif, false)});
    }
    const std::vector<Built> built = buildAll(requests);
    for (std::size_t k = 0; k < two.size(); ++k) {
        const auto report = [](const std::vector<std::string>& notes) {
            std::vector<std::string> out;
            for (const std::string& n : notes) {
                if (n.find("close contact") != std::string::npos || n.find("Close contact") != std::string::npos) out.push_back(n);
            }
            return out;
        };
        EXPECT_EQ(report(built[2 * k + 1].notes), report(built[2 * k].notes)) << two[k].seq;
    }
}

// cancellation through the reporter: a fitted build stops with core::Cancelled
TEST(Builder, FittedBuildIsCancellable) {
    auto flag = std::make_shared<std::atomic<bool>>(false);
    int calls = 0;
    const core::Reporter reporter({}, [&](double, const std::string&) {
        if (++calls == 3) flag->store(true);
    }, flag);
    HelixOptions o = fitted("hoogsteen");
    o.reporter = &reporter;
    o.threads = 1;
    EXPECT_THROW(buildHelix("GAGGAUCA", o), core::Cancelled);
}

// the same model for any number of finite-difference threads
TEST(Builder, FittedBuildDoesNotDependOnTheThreadCount) {
    HelixOptions one = fitted("reverse-wc", std::nullopt, true);
    one.threads = 1;
    HelixOptions four = one;
    four.threads = 4;
    EXPECT_EQ(buildHelix("AUGCAU", one), buildHelix("AUGCAU", four));
}

}  // namespace
