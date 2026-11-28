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
// CH-18 (docs/DIFFERENCES.md): MDATK's tests of reverse-wc, reverse-hoogsteen and the purine-motif
// triplex check here that RNA builds of them are refused, and run on the kinds still built; the
// Hoogsteen duplex and the pyrimidine-motif triplex have experimental (not 2.9 Å, coplanar) pairs; and
// uniform helices over a turn of every kind are checked against the experimental geometry of
// docs/reports/2026-11-27-non-wc-geometry.md (contacts, stacking, twist, rise).
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
#include <mutex>
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
#include "mdbridge/chem/reference_tables.hpp"
#include "mdbridge/chem/annotate.hpp"
#include "mdbridge/chem/helix.hpp"
#include "mdbridge/chem/mutate.hpp"
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
    static const Value data =
        core::parseJson(test::readFile(kData / "chem" / "builder_fitted_cases.json"), "builder_fitted_cases.json");
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
        if (key == "partner")
            o.partner = str(value);
        else if (key == "chirality")
            o.chirality = str(value);
        else if (key == "chains")
            o.chains = str(value);
        else if (key == "start")
            o.start = str(value);
        else if (key == "scheme")
            o.scheme = str(value);
        else if (key == "third")
            o.third = str(value);
        else if (key == "phosphate5")
            o.phosphate5 = *value.get<bool>();
        else if (key == "hydrogens")
            o.hydrogens = *value.get<bool>();
        else
            throw std::runtime_error("unknown helix option " + key);
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

HelixOptions fitted(std::string_view pairing, std::optional<std::string> triplex = std::nullopt,
                    bool hydrogens = false) {
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
    const std::size_t workers =
        std::max<std::size_t>(1, std::min<std::size_t>(std::thread::hardware_concurrency(), requests.size()));
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
    if (r.name == mdbridge::chem::reference::protonatedCytidine().code) return 'C';  // D19 / CH-21: the C+ (CCD CH)
    return io::parseResidueName(r.name)->base;
}

// the one-letter bases of a strand (CH, the C+, as C)
std::string strandBases(const std::vector<io::PdbResidue>& seg) {
    std::string out;
    for (const io::PdbResidue& r : seg) out += baseOf(r);
    return out;
}

// [(residue of strand 1, partner residue)] at each level
std::vector<std::pair<const io::PdbResidue*, const io::PdbResidue*>> strandPairs(const io::PdbStructure& s,
                                                                                 bool parallel,
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
         {{"AU", {{"N1", "N3"}, {"N6", "O2"}}},
          {"UA", {{"N3", "N1"}, {"O2", "N6"}}},
          {"GC", {{"N1", "O2"}, {"N2", "N3"}}},
          {"CG", {{"O2", "N1"}, {"N3", "N2"}}}}},
        {"hoogsteen", {{"AU", {{"N7", "N3"}, {"N6", "O4"}}}, {"GC", {{"N7", "N3"}, {"O6", "N4"}}}}},
        {"reverse-hoogsteen",
         {{"AU", {{"N7", "N3"}, {"N6", "O2"}}},
          {"AA", {{"N7", "N6"}, {"N6", "N1"}}},
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

std::vector<std::pair<std::size_t, std::size_t>> backboneLinks(const io::PdbStructure& s,
                                                               const std::vector<ModelAtom>& atoms) {
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
    static const std::set<std::string> names = {"H1",  "H3",  "H21",  "H22",  "H41", "H42",
                                                "H61", "H62", "HO2'", "HO5'", "HO3'", "HN3"};  // HN3: CH's (D19)
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
            EXPECT_TRUE(2.6 <= d && d <= 3.3)
                << where << " " << x->label() << " " << y->label() << " " << p << " " << q << " " << d;
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
        EXPECT_TRUE(stringsClose(notes[k], str(wantNotes[k])))
            << where << "\n  " << notes[k] << "\n  " << str(wantNotes[k]);
    }
}

INSTANTIATE_TEST_SUITE_P(Rebuilt, FittedCase, ::testing::Range<std::size_t>(0, 26));

TEST(FittedCaseData, AllRebuiltCasesAreTested) {
    EXPECT_EQ(list(at(fittedCases(), "rebuilt")).size(), 26u);
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
    const std::vector<std::string> rows =
        lines(test::readFile(kData / ".." / ".." / "data" / "sources" / "pdb" / "4jrd_atoms.tsv"));
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
        if (col("label_comp_id") != "A" || col("pdbx_PDB_model_num") != "1" ||
            (col("label_alt_id") != "." && col("label_alt_id") != "A")) {
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
    const std::vector<std::pair<std::string, std::string>> cases = {{"hoogsteen", "AAGAGAAGGA"}, {"hh", "AAAAAAAA"}};
    std::vector<Request> requests;
    for (const auto& [pairing, seq] : cases) requests.push_back({seq, fitted(pairing)});
    const std::vector<Built> built = buildAll(requests);
    for (std::size_t k = 0; k < cases.size(); ++k) {
        const auto& [pairing, seq] = cases[k];
        const io::PdbStructure s = io::readPdbText(built[k].text);
        checkModel(built[k], strandPairs(s, parallelPairing(pairing)), hbonds(pairing), pairing);
        EXPECT_EQ(strandBases(s.segments[1]), autoPartner(seq, pairing)) << pairing;
        for (const io::PdbResidue& r : s.segments[1]) {  // D19: the C+ of the G·C+ Hoogsteen pairs is CH
            if (baseOf(r) == 'C') {
                EXPECT_EQ(r.name, pairing == "hoogsteen" ? "CH" : "C") << r.label();
            }
        }
        EXPECT_TRUE(
            anyNote(built[k].notes, [&](const std::string& n) { return startsWith(n, "Pairing " + pairing + ":"); }));
        const std::string lw = pairing == "reverse-wc"  ? "tWW"
                               : pairing == "hoogsteen" ? "cHW"
                               : pairing == "hh"        ? "tHH"
                                                        : "tHW";
        EXPECT_NE(built[k].text.find("PAIRING " + lw), std::string::npos) << pairing;
        EXPECT_FALSE(
            anyNote(built[k].notes, [](const std::string& n) { return n.find("mismatch") != std::string::npos; }));
    }
    // CH-18: RNA is not known to form the parallel reverse-WC and the reverse-Hoogsteen duplexes
    for (const auto& [pairing, what] :
         {std::pair<std::string, std::string>{"reverse-wc", "the parallel reverse Watson-Crick (tWW) duplex"},
          {"reverse-hoogsteen", "the reverse-Hoogsteen (tHW) duplex"}}) {
        const std::string message = "pairing " + pairing + ": RNA is not known to form " + what +
                                    "; DNA and DNA-RNA hybrid builds of it will come with DNA structure building";
        ASSERT_EQ(rnaRefusal(pairing), message);
        try {
            buildHelix("AAGAGAAGGA", fitted(pairing));
            ADD_FAILURE() << pairing << " was built";
        } catch (const BuildError& e) {
            EXPECT_EQ(std::string(e.what()), message);
        }
        EXPECT_EQ(rnaRefusal(pairing), message);
        EXPECT_EQ(checkHelix("AAGAGAAGGA", "auto", pairing), (std::vector<HelixMessage>{{"error", message}}));
    }
    EXPECT_EQ(rnaRefusal("hoogsteen"), std::nullopt);
    EXPECT_EQ(rnaRefusal("hh"), std::nullopt);
    EXPECT_EQ(rnaRefusal("wc"), std::nullopt);
}

// test_triplex_third_strand
TEST(BuilderMdatk, TriplexThirdStrand) {
    const std::string seq = "AAGAGAAGGA";
    const std::vector<Built> built = buildAll({{seq, fitted("wc", "hoogsteen")}});
    for (const std::string motif : {"hoogsteen"}) {
        const Built& b = built[0];
        const io::PdbStructure s = io::readPdbText(b.text);
        EXPECT_EQ(s.chainIds(), (std::vector<std::string>{"A", "B", "C"}));
        ASSERT_EQ(s.segments.size(), 3u);
        for (const auto& seg : s.segments) EXPECT_EQ(seg.size(), 10u);
        const std::string third = strandBases(s.segments[2]);
        EXPECT_EQ(third, resolveThird(seq, "auto", motif).sequence);
        for (const io::PdbResidue& r : s.segments[2]) {  // D19: the third strand's C+ are CH
            if (baseOf(r) == 'C') {
                EXPECT_EQ(r.name, "CH") << r.label();
            }
        }
        for (const io::PdbResidue& r : s.segments[1]) EXPECT_NE(r.name, "CH") << r.label();  // Watson-Crick C
        EXPECT_EQ(third, motif == "hoogsteen" ? "UUCUCUUCCU" : "AGGAAGAGAA");
        checkModel(b, strandPairs(s, false, 1), hbonds("wc"), motif + " duplex");
        checkModel(b, strandPairs(s, motif == "hoogsteen", 2), hbonds(motif), motif + " triples");
        if (motif == "hoogsteen") {
            EXPECT_TRUE(anyNote(b.notes, [](const std::string& n) {
                return n.find("G·C+ Hoogsteen pair(s): the cytosine is N3-protonated (C+), residue CH (the CCD's "
                              "N3-protonated cytidine; HN3 comes with the hydrogens)") != std::string::npos;
            }));
        }
        EXPECT_TRUE(anyNote(b.notes, [](const std::string& n) { return startsWith(n, "Triplex, "); }));
        EXPECT_NE(b.text.find("THIRD STRAND"), std::string::npos);
    }
    // CH-18: RNA is not known to form the purine-motif triplex
    const std::string message =
        "triplex reverse-hoogsteen: RNA is not known to form the purine-motif "
        "(reverse-Hoogsteen) triplex; DNA builds of it will come with DNA structure building";
    try {
        buildHelix(seq, fitted("wc", "reverse-hoogsteen"));
        ADD_FAILURE() << "the purine-motif triplex was built";
    } catch (const BuildError& e) {
        EXPECT_EQ(std::string(e.what()), message);
    }
    EXPECT_EQ(rnaRefusal("wc", "reverse-hoogsteen"), message);
    EXPECT_EQ(rnaRefusal("wc", "hoogsteen"), std::nullopt);
    EXPECT_EQ(checkHelix(seq, "auto", "wc", "reverse-hoogsteen"), (std::vector<HelixMessage>{{"error", message}}));
}

// test_l_pairings_are_mirror_images
TEST(BuilderMdatk, LPairingsAreMirrorImages) {
    const std::vector<HelixOptions> kinds = {fitted("hoogsteen"), fitted("hh"), fitted("wc", "hoogsteen")};
    std::vector<Request> requests;
    for (const HelixOptions& o : kinds) {
        const std::string seq = o.pairing == std::optional<std::string>("hh") ? "AAAAAA" : "AAGAGA";  // D19
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
        // the D build's C+ are CH (CH-21), which the PDB tools read as C nucleotides (CH-22)
        const auto ch = std::count_if(rd.begin(), rd.end(), [](const io::PdbResidue* r) { return r->name == "CH"; });
        EXPECT_EQ(ch, k == 1 ? 0 : 2);
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
    const std::vector<Built> built = buildAll({{"AAGA", hoog}, {"AAGAG", pq}});
    EXPECT_TRUE(anyNote(built[0].notes, [](const std::string& n) {
        return startsWith(n, "WARNING: mismatch A:3 · B:3: no cHW pair");
    }));
    EXPECT_EQ(io::readPdbText(built[0].text).segments.size(), 2u);
    // D19: a strand-1 pyrimidine has no Hoogsteen edge: refused (MDATK: placed and reported)
    try {
        buildHelix("AAUAG", fitted("wc", "hoogsteen"));
        ADD_FAILURE() << "AAUAG was built";
    } catch (const BuildError& e) {
        EXPECT_EQ(std::string(e.what()), "triplex hoogsteen: the third strand pairs with the Hoogsteen edge of strand "
                                         "1, which must be a purine (A, G); strand 1 has U or C at 3:U");
    }
    const io::PdbStructure s = io::readPdbText(built[1].text);
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
    HelixOptions tri = fitted("wc", "hoogsteen");  // MDATK's: the purine motif (third GAGAA), refused (CH-18)
    tri.third = "UUCUC";
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
// MIXED: other pairs than the reference pair, G·C+ and C+·G-C levels (MDATK's GAGGAUCA and AUGCAU of
// every kind: reverse-wc, reverse-hoogsteen and the purine-motif triplex refused, CH-18; their strand-1
// pyrimidines refused in the Hoogsteen kinds, hh's non-A levels refused, D19)
const std::vector<Mixed>& mixed() {
    static const std::vector<Mixed> cases = {{"hoogsteen", {}, "AAGAGAAGGA"},
                                             {"wc", "hoogsteen", "GAGGAG"},
                                             {"hoogsteen", {}, "GGAGGA"},
                                             {"wc", "hoogsteen", "AAGAGAAGGA"}};
    return cases;
}

// the builds of mixed sequences, each once in the test program (with or without hydrogens)
std::vector<Built> buildCached(const std::vector<Mixed>& cases, bool hydrogens) {
    static std::mutex mutex;
    static std::map<std::string, Built> cache;
    const auto key = [&](const Mixed& m) {
        return m.pairing + "/" + m.motif.value_or("") + "/" + m.seq + (hydrogens ? "/H" : "");
    };
    std::vector<Request> requests;
    std::vector<std::string> keys;
    {
        const std::lock_guard<std::mutex> lock(mutex);
        for (const Mixed& m : cases) {
            if (!cache.count(key(m)) && std::find(keys.begin(), keys.end(), key(m)) == keys.end()) {
                requests.push_back({m.seq, fitted(m.pairing, m.motif, hydrogens)});
                keys.push_back(key(m));
            }
        }
    }
    const std::vector<Built> built = buildAll(requests);
    const std::lock_guard<std::mutex> lock(mutex);
    for (std::size_t k = 0; k < keys.size(); ++k) cache[keys[k]] = built[k];
    std::vector<Built> out;
    for (const Mixed& m : cases) out.push_back(cache.at(key(m)));
    return out;
}

std::vector<Built> buildWithHydrogens(const std::vector<Mixed>& cases) {
    return buildCached(cases, true);
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
            std::transform(lower.begin(), lower.end(), lower.begin(),
                           [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
            if (lower.find("close contact") != std::string::npos) report.push_back(n);
        }
        ASSERT_EQ(report.size(), 1u) << where;
        EXPECT_NE(report[0].find("hydrogens included"), std::string::npos) << where;
        for (const std::string_view k2 : {"(H-H <", "(H-heavy <", "(H-bond H <"}) {
            EXPECT_EQ(report[0].find(k2), std::string::npos) << where << ": " << report[0];
        }
    }
}

// test_reverse_wc_bases_clear_their_own_sugars (reverse-wc is refused, CH-18: the kinds still built)
TEST(BuilderMdatk, ReverseWcBasesClearTheirOwnSugars) {
    const std::vector<Mixed>& cases = mixed();
    std::vector<std::string> seqs;
    for (const Mixed& m : cases) seqs.push_back(m.pairing + " " + m.motif.value_or("") + " " + m.seq);
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
    // MDATK's: reverse-hoogsteen AAGAGG, GAAGAG (refused, CH-18) and the triplex GAGGAG; hh A only (D19)
    const std::vector<Mixed> cases = {{"hoogsteen", {}, "AAGAGG"}, {"hh", {}, "AAAAAA"}, {"wc", "hoogsteen", "GAGGAG"}};
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
    // MDATK's also reverse-wc AUGCAU, GAGGAUCA, reverse-hoogsteen GAGGAUCA and the purine-motif triplex
    // GAAGAG (refused, CH-18), and GAGGAUCA of the Hoogsteen kinds (strand-1 pyrimidines, refused: D19);
    // hh added. hh (the 4JRD A·A pair, 2.74 Å), the Hoogsteen duplex and the pyrimidine-motif triplex keep
    // their experimental pairs: H-bonds 2.6-3.3 Å (MDATK's 2.8-3.0: its pairs are fitted to 2.9 Å),
    // H...acceptor within 2.4 Å at 133 degrees or more (the RNA triple helices' smallest, hydrogens placed
    // by the builder's rules; D19's closure moves the partners within their levels: 134° at a helix end),
    // and the propeller of the RNA triple helices (bases within 20 degrees; MDATK's: 10, coplanar pairs)
    const std::vector<Mixed> cases = {{"hoogsteen", {}, "GGAGGA"},
                                      {"hoogsteen", {}, "AAGAGAAGGA"},
                                      {"hh", {}, "AAAAAAAA"},
                                      {"wc", "hoogsteen", "GAGGAG"}};
    const std::vector<Built> built = buildWithHydrogens(cases);
    for (std::size_t k = 0; k < cases.size(); ++k) {
        const Mixed& m = cases[k];
        const io::PdbStructure s = io::readPdbText(built[k].text);
        std::vector<std::pair<std::vector<std::pair<const io::PdbResidue*, const io::PdbResidue*>>, const HbondTable*>>
            groups;
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
                const bool experimental =
                    m.pairing == "hoogsteen" || m.pairing == "hh" || m.motif == std::optional<std::string>("hoogsteen");
                for (const auto& [x, y] : it->second) {
                    const double d = dist(ca.at(x), cb.at(y));
                    if (experimental) {
                        EXPECT_TRUE(2.6 <= d && d <= 3.3) << where << " " << x << " " << y << " " << d;
                    } else {
                        EXPECT_TRUE(2.8 <= d && d <= 3.0) << where << " " << x << " " << y << " " << d;
                    }
                    for (const auto& [don, acc, cd, cc] :
                         {std::make_tuple(x, y, &ca, &cb), std::make_tuple(y, x, &cb, &ca)}) {
                        std::vector<std::string> hs;
                        for (const auto& [h, v] : *cd) {
                            if (h[0] == 'H' && dist(v, cd->at(don)) < 1.1) hs.push_back(h);
                        }
                        if (!hs.empty() && (acc[0] == 'N' || acc[0] == 'O') && don[0] == 'N') {
                            const std::string h = *std::min_element(
                                hs.begin(), hs.end(), [&](const std::string& p, const std::string& q) {
                                    return dist(cd->at(p), cc->at(acc)) < dist(cd->at(q), cc->at(acc));
                                });
                            if (dist(cd->at(h), cc->at(acc)) < 2.6) {
                                const Vec u = sub(cd->at(don), cd->at(h)), v = sub(cc->at(acc), cd->at(h));
                                const double angle =
                                    std::acos(dot(u, v) / std::sqrt(dot(u, u)) / std::sqrt(dot(v, v))) * 180.0 /
                                    3.14159265358979323846;
                                EXPECT_LE(dist(cd->at(h), cc->at(acc)), experimental ? 2.4 : 2.05)
                                    << where << " " << h << " " << acc;
                                EXPECT_GE(angle, experimental ? 133.0 : 150.0) << where << " " << h << " " << acc;
                            }
                        }
                    }
                }
                const Vec n1 = baseNormal(ca, key[0]), n2 = baseNormal(cb, key[1]);
                EXPECT_LE(std::acos(std::min(1.0, std::fabs(dot(n1, n2)))) * 180.0 / 3.14159265358979323846,
                          experimental ? 20.0 : 10.0)
                    << where;
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
    for (const auto& [n, x] : resAtoms[0])
        if (n == "C2'") c2 = x;
    for (const auto& [n, x] : resAtoms[1])
        if (n == "C5'") c5 = x;
    const Vec u = sub(c2, c5);
    const double f = 1.0 - 2.7 / std::sqrt(dot(u, u));
    for (auto& [n, x] : resAtoms[1]) x = {x[0] + u[0] * f, x[1] + u[1] * f, x[2] + u[2] * f};
    found = helix_detail::closeContacts(resAtoms, segments, bases, {});
    EXPECT_TRUE(std::any_of(found.begin(), found.end(), [](const helix_detail::Contact& c) {
        return c.a.second == "C2'" && c.b.second == "C5'" && startsWith(c.what, "heavy") && std::fabs(c.d - 2.7) < 1e-6;
    }));
    // a heavy-atom model gets hydrogens for the check: the same report as with hydrogens
    const std::vector<Mixed> two(mixed().begin(), mixed().begin() + 2);
    const std::vector<Built> withH = buildCached(two, true), heavyOnly = buildCached(two, false);
    std::vector<Built> built;
    for (std::size_t k = 0; k < two.size(); ++k) {
        built.push_back(withH[k]);
        built.push_back(heavyOnly[k]);
    }
    for (std::size_t k = 0; k < two.size(); ++k) {
        const auto report = [](const std::vector<std::string>& notes) {
            std::vector<std::string> out;
            for (const std::string& n : notes) {
                if (n.find("close contact") != std::string::npos || n.find("Close contact") != std::string::npos)
                    out.push_back(n);
            }
            return out;
        };
        EXPECT_EQ(report(built[2 * k + 1].notes), report(built[2 * k].notes)) << two[k].seq;
    }
}

// --- CH-18: the geometry of uniform helices over a turn, against experiment ----------------------------

// the atoms of level k (strand 1 residue k, each partner's residue at that level): C1' and the base heavy
// atoms, as docs/reports/2026-11-27-non-wc-geometry.md measures them
std::vector<Vec> levelAtoms(const io::PdbStructure& s, std::size_t k, const std::vector<bool>& parallel, bool bases) {
    std::vector<Vec> out;
    const std::size_t n = s.segments[0].size();
    for (std::size_t strand = 0; strand < s.segments.size(); ++strand) {
        const io::PdbResidue& r = s.segments[strand][parallel[strand] ? k : n - 1 - k];
        const Coords c = coordsOf(r);
        if (!bases) out.push_back(c.at("C1'"));
        for (const std::string_view a : baseAtomNames(baseOf(r))) out.push_back(c.at(std::string(a)));
    }
    return out;
}

// the screw (twist degrees, rise Å) of the rigid fit of level k onto level k + 1
struct StepScrew {
    double twist = 0.0, rise = 0.0;
    Vec axis{};
};
StepScrew stepScrewOf(const std::vector<Vec>& x, const std::vector<Vec>& y);
std::pair<double, double> stepScrew(const std::vector<Vec>& x, const std::vector<Vec>& y) {
    const StepScrew s = stepScrewOf(x, y);
    return {s.twist, s.rise};
}
StepScrew stepScrewOf(const std::vector<Vec>& x, const std::vector<Vec>& y) {
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
    if (core::det3(r) < 0) {
        for (std::size_t i = 0; i < 3; ++i)
            for (std::size_t j = 0; j < 3; ++j) r[i][j] -= 2.0 * svd.vt[2][i] * svd.u[j][2];
    }
    Vec t{};  // y = r x + t
    for (std::size_t i = 0; i < 3; ++i) t[i] = my[i] - (r[i][0] * mx[0] + r[i][1] * mx[1] + r[i][2] * mx[2]);
    const double angle = std::acos(std::clamp((r[0][0] + r[1][1] + r[2][2] - 1.0) / 2.0, -1.0, 1.0));
    Vec axis = {r[2][1] - r[1][2], r[0][2] - r[2][0], r[1][0] - r[0][1]};
    const double len = std::sqrt(dot(axis, axis));
    axis = {axis[0] / len, axis[1] / len, axis[2] / len};
    return {angle * 180.0 / 3.14159265358979323846, std::fabs(dot(axis, t)), axis};
}

// the separation of the planes of the base heavy atoms of two levels along their mean normal
double planeSeparation(const std::vector<Vec>& a, const std::vector<Vec>& b) {
    const auto plane = [](std::vector<Vec> p) {
        Vec m{};
        for (const Vec& v : p)
            for (std::size_t i = 0; i < 3; ++i) m[i] += v[i] / static_cast<double>(p.size());
        for (Vec& v : p) v = sub(v, m);
        return std::make_pair(m, core::svdRows3(p).vt[2]);
    };
    auto [ca, na] = plane(a);
    auto [cb, nb] = plane(b);
    if (dot(na, nb) < 0) nb = {-nb[0], -nb[1], -nb[2]};
    Vec n = {na[0] + nb[0], na[1] + nb[1], na[2] + nb[2]};
    const double len = std::sqrt(dot(n, n));
    n = {n[0] / len, n[1] / len, n[2] / len};
    return std::fabs(dot(sub(cb, ca), n));
}

// the stacking of consecutive bases of one strand: the perpendicular distance of each base's centroid from
// the other's plane, averaged (the report's "stack perp")
double baseStacking(const io::PdbResidue& a, const io::PdbResidue& b) {
    const auto plane = [](const io::PdbResidue& r) {
        const Coords c = coordsOf(r);
        std::vector<Vec> p;
        for (const std::string_view n : baseAtomNames(baseOf(r))) p.push_back(c.at(std::string(n)));
        Vec m{};
        for (const Vec& v : p)
            for (std::size_t i = 0; i < 3; ++i) m[i] += v[i] / static_cast<double>(p.size());
        for (Vec& v : p) v = sub(v, m);
        return std::make_pair(m, core::svdRows3(p).vt[2]);
    };
    const auto [ca, na] = plane(a);
    const auto [cb, nb] = plane(b);
    const Vec v = sub(cb, ca);
    return 0.5 * (std::fabs(dot(v, na)) + std::fabs(dot(v, nb)));
}

// C1' and the base atoms of levels k and k + 1 that both residues of a strand have (an A level and a G level
// share C1', N9, C8, N7, C5, C6, N1, C2, N3, C4), as the report pairs them
std::pair<std::vector<Vec>, std::vector<Vec>> matchedLevelAtoms(const io::PdbStructure& s, std::size_t k,
                                                                const std::vector<bool>& parallel) {
    std::vector<Vec> x, y;
    const std::size_t n = s.segments[0].size();
    for (std::size_t strand = 0; strand < s.segments.size(); ++strand) {
        const io::PdbResidue& a = s.segments[strand][parallel[strand] ? k : n - 1 - k];
        const io::PdbResidue& b = s.segments[strand][parallel[strand] ? k + 1 : n - 2 - k];
        const Coords ca = coordsOf(a), cb = coordsOf(b);
        std::vector<std::string> names = {"C1'"};
        for (const std::string_view nm : baseAtomNames(baseOf(a))) names.emplace_back(nm);
        for (const std::string& nm : names) {
            if (cb.count(nm)) {
                x.push_back(ca.at(nm));
                y.push_back(cb.at(nm));
            }
        }
    }
    return {x, y};
}

// D19: the geometry of the experimental kinds against the RNA triple helices they come from (the
// steps of data/sources/pdb/triplex_levels.tsv measured as docs/reports/2026-11-27-non-wc-geometry.md
// measures them; tests/chem/test_builder.cpp's helpers). Ranges: the experimental steps' extremes of twist,
// rise, inclination and each strand's stacking (the generated tables, tools/gen_chem_tables.py, which the
// closure's walls hold the steps within), the largest base-plane separation, the means of the tables (the
// uniform helices).
struct ExperimentRange {
    double twistLow, twistHigh, riseLow, riseHigh, planesHigh, inclinationHigh;
    std::vector<double> stackHigh;  // per strand (strand 1, partners in the model's chain order)
    double twistMean, riseMean, planesMean;
};
const ExperimentRange& experimentRange(bool triplex) {
    static const auto of = [](std::string_view kind) {
        const reference::TriplexScrew& sc = *reference::triplexScrew(kind);
        const std::size_t strands = kind == "hoogsteen" ? 2 : 3;
        return ExperimentRange{sc.twistRange[0], sc.twistRange[1], sc.riseRange[0], sc.riseRange[1],
                               3.49,  // purine + Hoogsteen strands and the triples: tests/data/README.md
                               sc.inclinationRange[1], std::vector<double>(sc.stackMax, sc.stackMax + strands),
                               sc.twist, sc.rise, sc.stack};
    };
    static const ExperimentRange duplex = of("hoogsteen");
    static const ExperimentRange triple = of("triplex-hoogsteen");
    return triplex ? triple : duplex;
}

// a model's steps are measured on its PDB text (coordinates to 0.001 Å) and the closure's walls are stiff,
// not hard: within these of the experimental extremes
constexpr double kAngleTolerance = 0.05;   // degrees
constexpr double kLengthTolerance = 0.005;  // Å

double angleDeg(const Vec& a, const Vec& b, const Vec& c) {
    const Vec u = sub(a, b), v = sub(c, b);
    return std::acos(std::clamp(dot(u, v) / std::sqrt(dot(u, u) * dot(v, v)), -1.0, 1.0)) * 180.0 /
           3.14159265358979323846;
}

// the checks of one model of an experimental kind (uniform: also its screw against the tables' means)
void checkAgainstExperiment(const Built& built, bool triplex, const std::string& where, bool uniform) {
    const ExperimentRange& e = experimentRange(triplex);
    const io::PdbStructure s = io::readPdbText(built.text);
    EXPECT_EQ(closeContactsOf(s), std::vector<std::string>{}) << where;
    EXPECT_TRUE(anyNote(built.notes, [](const std::string& n) { return startsWith(n, "No close contacts ("); }))
        << where;
    EXPECT_FALSE(anyNote(built.notes, [](const std::string& n) { return startsWith(n, "WARNING"); })) << where;
    const std::vector<bool> parallel = triplex ? std::vector<bool>{true, false, true} : std::vector<bool>{true, true};
    ASSERT_EQ(s.segments.size(), parallel.size()) << where;
    const std::size_t n = s.segments[0].size();
    for (std::size_t level = 0; level + 1 < n; ++level) {
        const std::string at = where + " step " + std::to_string(level + 1);
        const auto [from, to] = matchedLevelAtoms(s, level, parallel);
        const StepScrew st = stepScrewOf(from, to);
        const std::vector<Vec> here = levelAtoms(s, level, parallel, true);
        const double layers = planeSeparation(here, levelAtoms(s, level + 1, parallel, true));
        EXPECT_TRUE(e.twistLow - kAngleTolerance <= st.twist && st.twist <= e.twistHigh + kAngleTolerance)
            << at << ": twist " << st.twist;
        EXPECT_TRUE(e.riseLow - kLengthTolerance <= st.rise && st.rise <= e.riseHigh + kLengthTolerance)
            << at << ": rise " << st.rise;
        EXPECT_TRUE(3.0 <= layers && layers <= e.planesHigh) << at << ": level planes " << layers;
        // inclination: the level's base plane against the step's axis
        Vec m{};
        std::vector<Vec> p = here;
        for (const Vec& v : p)
            for (std::size_t i = 0; i < 3; ++i) m[i] += v[i] / static_cast<double>(p.size());
        for (Vec& v : p) v = sub(v, m);
        const Vec normal = core::svdRows3(p).vt[2];
        const double incl =
            std::acos(std::min(1.0, std::fabs(dot(normal, st.axis)))) * 180.0 / 3.14159265358979323846;
        EXPECT_LE(incl, e.inclinationHigh + kAngleTolerance) << at << ": inclination";
        for (std::size_t strand = 0; strand < s.segments.size(); ++strand) {
            const auto& seg = s.segments[strand];
            const double stack = baseStacking(seg[level], seg[level + 1]);
            EXPECT_TRUE(3.0 <= stack && stack <= e.stackHigh[strand] + kLengthTolerance)
                << at << " " << seg[level].label() << ": stacking " << stack;
        }
        // stacked bases do not H-bond: their N/O atoms 3.05 Å apart or more (experiment: 2.82 Å, the U O4 pairs
        // of consecutive Hoogsteen pyrimidines; A N1 - U N3 of the next level 3.15 Å)
        for (std::size_t a = 0; a < s.segments.size(); ++a) {
            for (std::size_t b = 0; b < s.segments.size(); ++b) {
                const auto& sa = s.segments[a];
                const auto& sb = s.segments[b];
                const io::PdbResidue& ra = sa[parallel[a] ? level : n - 1 - level];
                const io::PdbResidue& rb = sb[parallel[b] ? level + 1 : n - 2 - level];
                const Coords ca = coordsOf(ra), cb = coordsOf(rb);
                for (const std::string_view x : baseAtomNames(baseOf(ra))) {
                    for (const std::string_view y : baseAtomNames(baseOf(rb))) {
                        if ((x[0] == 'N' || x[0] == 'O') && (y[0] == 'N' || y[0] == 'O')) {
                            EXPECT_GE(dist(ca.at(std::string(x)), cb.at(std::string(y))), 3.05)
                                << at << " " << ra.label() << " " << x << " - " << rb.label() << " " << y;
                        }
                    }
                }
            }
        }
        if (uniform && level == 0) {
            std::cout << where << ": twist " << st.twist << " rise " << st.rise << " planes " << layers
                      << " inclination " << incl << "\n";
            // the uniform helix on the experimental screw (the backbone torsion terms keep it up to 1.5° off)
            EXPECT_NEAR(st.twist, e.twistMean, 1.5) << where;
            EXPECT_NEAR(st.rise, e.riseMean, 0.05) << where;
            EXPECT_NEAR(layers, e.planesMean, 0.03) << where;
        }
    }
    // the links: lengths, bond angles (RNA crystal structures C3'-O3'-P 119.7 ± 1.1°, 117.1-124.6°;
    // O3'-P-O5' 103.1 ± 1.1°, 100.2-106.4°) and epsilon (the experimental strands: -188 ... -101°)
    for (const auto& seg : s.segments) {
        for (std::size_t k = 0; k + 1 < seg.size(); ++k) {
            const Coords a = coordsOf(seg[k]), b = coordsOf(seg[k + 1]);
            const std::string at = where + " " + seg[k].label() + "-" + seg[k + 1].label();
            const double d = dist(a.at("O3'"), b.at("P"));
            EXPECT_TRUE(1.58 <= d && d <= 1.64) << at << ": O3'-P " << d;
            const double c3 = angleDeg(a.at("C3'"), a.at("O3'"), b.at("P"));
            EXPECT_TRUE(117.0 <= c3 && c3 <= 123.0) << at << ": C3'-O3'-P " << c3;
            const double o5 = angleDeg(a.at("O3'"), b.at("P"), b.at("O5'"));
            EXPECT_TRUE(100.2 <= o5 && o5 <= 106.4) << at << ": O3'-P-O5' " << o5;
            double eps = dihedral(a.at("C4'"), a.at("C3'"), a.at("O3'"), b.at("P"));
            if (eps > 0) eps -= 360.0;
            EXPECT_TRUE(-200.0 <= eps && eps <= -125.0) << at << ": epsilon " << eps;
        }
    }
    const auto linkNote = std::find_if(built.notes.begin(), built.notes.end(),
                                       [](const std::string& x) { return startsWith(x, "O3'-P links: "); });
    ASSERT_NE(linkNote, built.notes.end()) << where;
    EXPECT_NE(linkNote->find("; angles C3'-O3'-P "), std::string::npos) << *linkNote;
}

// Uniform helices of 14 levels (more than a turn) and mixed sequences of the experimental kinds, against the
// experiment; hh (PDB 4JRD, MDATK's reference fit): twist, rise and stacking against 4JRD.
TEST(BuilderGeometry, ExperimentalKindsMatchTheRnaTripleHelices) {
    struct Case {
        std::optional<std::string> motif;
        std::string seq;
    };
    const std::vector<Case> cases = {{{}, std::string(14, 'A')},
                                     {"hoogsteen", std::string(14, 'A')},
                                     {{}, "AAGAGAAGGA"},
                                     {"hoogsteen", "AAGAGAAGGA"},
                                     {{}, "GGAGGA"},
                                     {"hoogsteen", "GAGGAG"}};
    std::vector<Mixed> requests;
    for (const Case& c : cases) requests.push_back({c.motif ? "wc" : "hoogsteen", c.motif, c.seq});
    const std::vector<Built> built = buildWithHydrogens(requests);
    for (std::size_t k = 0; k < cases.size(); ++k) {
        const std::string where = std::string(cases[k].motif ? "triplex " : "hoogsteen ") + cases[k].seq;
        checkAgainstExperiment(built[k], cases[k].motif.has_value(), where,
                               cases[k].seq.find_first_not_of('A') == std::string::npos);
    }
    // hh: the 4JRD screw (42.3° / 3.76 Å; MDATK's fit 41.7° / 3.79 Å), strands stacked 3.4-3.7 Å
    const Built hh = buildOne(std::string(14, 'A'), fitted("hh", std::nullopt, true));
    const io::PdbStructure s = io::readPdbText(hh.text);
    EXPECT_EQ(closeContactsOf(s), std::vector<std::string>{});
    for (std::size_t level = 0; level + 1 < 14; ++level) {
        const auto [twist, rise] =
            stepScrew(levelAtoms(s, level, {true, true}, false), levelAtoms(s, level + 1, {true, true}, false));
        EXPECT_NEAR(twist, 42.3, 1.0) << level;
        EXPECT_NEAR(rise, 3.76, 0.05) << level;
        for (const auto& seg : s.segments) {
            const double stack = baseStacking(seg[level], seg[level + 1]);
            EXPECT_TRUE(3.4 <= stack && stack <= 3.7) << level << ": stacking " << stack;
        }
    }
}

// D19: a triplex of 14 C+·G-C triples, whose steps went past the experimental extremes before the closure's
// walls (twist 35.6°, inclination 17.6°: docs/reports/2026-11-27-non-wc-geometry.md); its own test, as it
// takes the longest of the builds here (5-7 minutes on one thread)
TEST(BuilderGeometry, TriplexGTractMatchesTheRnaTripleHelices) {
    const std::vector<Built> built = buildWithHydrogens({{"wc", "hoogsteen", std::string(14, 'G')}});
    checkAgainstExperiment(built[0], true, "triplex GGGGGGGGGGGGGG", false);
}

// D19: the kinds that are built refuse the sequences they cannot pair (the positions named)
TEST(BuilderGeometry, SequencesAKindCannotPairAreRefused) {
    const auto refused = [](const std::string& seq, HelixOptions o) {
        try {
            buildHelix(seq, o);
        } catch (const BuildError& e) {
            return std::string(e.what());
        }
        return std::string("built");
    };
    EXPECT_EQ(refused("AAGAGAUGGA", fitted("hoogsteen")),
              "pairing hoogsteen: strand 2 pairs with the Hoogsteen edge of strand 1, which must be a purine (A, G); "
              "strand 1 has U or C at 7:U");
    EXPECT_EQ(refused("UCUCUCUCUC", fitted("wc", "hoogsteen")),
              "triplex hoogsteen: the third strand pairs with the Hoogsteen edge of strand 1, which must be a purine "
              "(A, G); strand 1 has U or C at 1:U, 2:C, 3:U, 4:C, 5:U, 6:C, 7:U, 8:C ...");
    EXPECT_EQ(refused("AAGAG", fitted("hh")),
              "pairing hh: the parallel duplex of PDB 4JRD has A·A pairs only; other bases at 3:G·G, 5:G·G");
    HelixOptions hh = fitted("hh");
    hh.partner = "AACA";
    EXPECT_EQ(refused("AAAA", hh), "pairing hh: the parallel duplex of PDB 4JRD has A·A pairs only; other bases at 3:A·C");
    // the check reports the same error; the requests MDATK refuses keep MDATK's errors
    EXPECT_EQ(checkHelix("AAGAG", "auto", "hh"),
              (std::vector<HelixMessage>{
                  {"error", "pairing hh: the parallel duplex of PDB 4JRD has A·A pairs only; other bases at 3:G·G, 5:G·G"}}));
    EXPECT_EQ(sequenceRefusal("hoogsteen", std::nullopt, "AAGAG", "UUCUC"), std::nullopt);
    EXPECT_EQ(sequenceRefusal("wc", std::nullopt, "UCUC", "GAGA"), std::nullopt);
    HelixOptions single = fitted("hoogsteen");
    single.partner = "none";
    EXPECT_NE(refused("UUUU", single).find("a single strand has no pairing"), std::string::npos);
}

// D19 / CH-21: the C+ of the G·C+ Hoogsteen pairs is the CCD's N3-protonated cytidine CH, with its HN3 when
// hydrogens are written (the H-bond of the pair: N7 ... HN3 under 2.1 Å); an L build writes plain C
TEST(BuilderGeometry, ProtonatedCytidinesAreTheCcdsCh) {
    HelixOptions h = fitted("hoogsteen", std::nullopt, true);
    HelixOptions heavy = fitted("hoogsteen");
    HelixOptions mirror = fitted("hoogsteen");
    mirror.chirality = "L";
    const std::vector<Built> built = buildAll({{"AGGA", h}, {"AGGA", heavy}, {"AGGA", mirror}});
    const io::PdbStructure s = io::readPdbText(built[0].text);
    for (std::size_t k = 0; k < 4; ++k) {
        const io::PdbResidue& r = s.segments[1][k];
        const bool plus = k == 1 || k == 2;
        EXPECT_EQ(r.name, plus ? "CH" : "U") << k;
        const Coords c = coordsOf(r);
        EXPECT_EQ(c.count("HN3"), plus ? 1u : 0u) << k;
        if (plus) {
            EXPECT_NEAR(dist(c.at("HN3"), c.at("N3")), 1.01, 0.02);
            EXPECT_LT(dist(c.at("HN3"), coordsOf(s.segments[0][k]).at("N7")), 2.1);
        }
    }
    EXPECT_NE(built[0].text.find("ATOM  "), std::string::npos);
    EXPECT_TRUE(anyNote(built[0].notes, [](const std::string& n) {
        return n == "2 G·C+ Hoogsteen pair(s): the cytosine is N3-protonated (C+), residue CH (the CCD's "
                    "N3-protonated cytidine, with HN3); at neutral pH the pair is weak: restrain it for the simulation";
    }));
    const io::PdbStructure noH = io::readPdbText(built[1].text);
    EXPECT_EQ(noH.segments[1][1].name, "CH");
    EXPECT_EQ(coordsOf(noH.segments[1][1]).count("HN3"), 0u);
    const io::PdbStructure l = io::readPdbText(built[2].text);
    EXPECT_EQ(l.segments[1][1].name, "0C");
    EXPECT_TRUE(anyNote(built[2].notes, [](const std::string& n) {
        return n.find("the CCD has no N3-protonated L-cytidine") != std::string::npos;
    }));
    HelixOptions pdb = fitted("hoogsteen");
    pdb.scheme = "pdb";
    const std::string text = buildHelix("AGGA", pdb);
    EXPECT_NE(text.find("HETATM"), std::string::npos);  // a modified residue: HETATM records (PDB scheme)
    EXPECT_EQ(checkHelix("AGGA", "auto", "hoogsteen").back(),
              (HelixMessage{"info", "G·C+ Hoogsteen pairs need N3-protonated cytosine (C+): written as CH, the CCD's "
                                    "N3-protonated cytidine (HN3 with the hydrogens; L: plain C; restrain the pair for "
                                    "the simulation)"}));
}

// CH-22: what the builder writes, the PDB tools read: a G·C+ Hoogsteen helix with hydrogens, its CH residues
// C nucleotides for pdb check (sequence, chirality count, not "other residues"), annotate (the sequence),
// clean's naming schemes (CH kept) and mutate (CH -> G drops HN3)
TEST(BuilderGeometry, ProtonatedCytidinesReadBackAsC) {
    const Built built = buildOne("AGGA", fitted("hoogsteen", std::nullopt, true));
    const io::PdbStructure s = io::readPdbText(built.text);
    ASSERT_EQ(s.segments[1][1].name, "CH");
    std::string partner;
    for (const io::PdbResidue& r : s.segments[1]) {
        partner += r.name == "CH" ? 'C' : io::parseResidueName(r.name)->base;
    }
    const std::string info = core::toJsonAscii(check(PdbSource::text(built.text)));
    EXPECT_NE(info.find("\"sequence\": \"" + partner + "\""), std::string::npos) << info;
    EXPECT_NE(info.find("\"other_residues\": {}"), std::string::npos) << info;
    EXPECT_NE(info.find("Chirality (sugar geometry): 0 L and 8 D nucleotides"), std::string::npos) << info;
    const AnnotationResult a = annotate(s);
    EXPECT_NE(std::find(a.records.begin(), a.records.end(), "REMARK 950 SEQ   " + partner), a.records.end());
    for (const std::string& n : a.notes) {
        EXPECT_EQ(n.find("CH"), std::string::npos) << n;  // no note of a residue it cannot read
    }
    io::PdbStructure named = s;
    convertNames(named, "amber");
    EXPECT_EQ(named.segments[1][1].name, "CH");
    EXPECT_EQ(named.segments[1][0].name, "U5");
    const io::PdbResidue& c = s.segments[1][1];
    const MutateResult m = mutate(PdbSource::text(built.text), c.chain + ":" + std::to_string(c.resseq) + "G");
    const io::PdbStructure mutated = io::readPdbText(m.text);
    EXPECT_EQ(mutated.segments[1][1].name, "G");
    EXPECT_EQ(coordsOf(mutated.segments[1][1]).count("HN3"), 0u);
}

// CH-18: a built pair whose bases have no N/O donor-acceptor pair within 3.5 Å is reported (MDATK
// reports only the mismatch, not that it has no H-bond); hh's U·U and C·C levels are refused (D19), the
// case is a Hoogsteen G·A mismatch
TEST(BuilderGeometry, PairsWithoutABaseBaseHbondAreReported) {
    HelixOptions o = fitted("hoogsteen");
    o.partner = "UUAU";
    const Built built = buildOne("AAGA", o);
    const auto note = [](const std::vector<std::string>& notes) {
        const auto it = std::find_if(notes.begin(), notes.end(), [](const std::string& n) {
            return startsWith(n, "WARNING: ") && n.find("without a base-base H-bond") != std::string::npos;
        });
        return it == notes.end() ? std::string() : *it;
    };
    const std::string hoog = note(built.notes);
    EXPECT_EQ(hoog, "WARNING: 1 base pair(s) without a base-base H-bond (no N/O donor-acceptor pair within 3.5 Å): "
                    "A:3 G · B:3 A");
    // the note comes before the closing "Built ..." note
    const auto pos = std::find(built.notes.begin(), built.notes.end(), hoog);
    ASSERT_NE(pos, built.notes.end());
    ASSERT_NE(pos + 1, built.notes.end());
    EXPECT_TRUE(startsWith(*(pos + 1), "Built "));
    // canonical pairs only: no note; a Watson-Crick mismatch without an H-bond: the note
    EXPECT_EQ(note(buildOne("AAAAAAAA", fitted("hh")).notes), "");
    EXPECT_EQ(note(buildOne("GGCAUC", HelixOptions{}).notes), "");
    HelixOptions mism;
    mism.partner = "GAUCCC";
    const std::string wc = note(buildOne("GGCAUC", mism).notes);
    EXPECT_TRUE(startsWith(wc, "WARNING: ")) << wc;
}

// cancellation through the reporter: a fitted build stops with core::Cancelled
TEST(Builder, FittedBuildIsCancellable) {
    auto flag = std::make_shared<std::atomic<bool>>(false);
    int calls = 0;
    const core::Reporter reporter(
        {},
        [&](double, const std::string&) {
            if (++calls == 3) flag->store(true);
        },
        flag);
    HelixOptions o = fitted("hoogsteen");
    o.reporter = &reporter;
    o.threads = 1;
    EXPECT_THROW(buildHelix("GAGGAGAAGA", o), core::Cancelled);
}

// the same model for any number of finite-difference threads
TEST(Builder, FittedBuildDoesNotDependOnTheThreadCount) {
    HelixOptions one = fitted("hoogsteen", std::nullopt, true);
    one.threads = 1;
    HelixOptions four = one;
    four.threads = 4;
    EXPECT_EQ(buildHelix("AGGAAG", one), buildHelix("AGGAAG", four));
    // the reference fit itself, fitted afresh (the builds above share the cached one)
    for (const std::string kind : {"hoogsteen", "triplex-hoogsteen"}) {
        EXPECT_EQ(helix_detail::referenceFitFresh(kind, 1), helix_detail::referenceFitFresh(kind, 4)) << kind;
    }
}

// rotationVector: angles up to pi, the near-pi branch included
TEST(Builder, RotationVectorUpToPi) {
    for (const double deg : {0.0, 1e-9, 30.0, 120.0, 179.0, 179.9999, 180.0}) {
        for (const mdbridge::chem::Vec3& axis0 : {mdbridge::chem::Vec3{1.0, 0.0, 0.0}, mdbridge::chem::Vec3{0.3, -0.5, 0.81},
                                                 mdbridge::chem::Vec3{-0.2, -0.7, -0.4}}) {
            const mdbridge::chem::Vec3 axis = helix_detail::unit(axis0);
            const helix_detail::Mat3 r = helix_detail::turnMatrix(axis, deg);
            const mdbridge::chem::Vec3 w = helix_detail::rotationVector(r);
            const helix_detail::Mat3 back = helix_detail::turnMatrix(
                helix_detail::norm(w) > 0 ? helix_detail::unit(w) : axis, helix_detail::norm(w) * 180.0 / 3.14159265358979323846);
            for (std::size_t i = 0; i < 3; ++i)
                for (std::size_t j = 0; j < 3; ++j) EXPECT_NEAR(back[i][j], r[i][j], 1e-9) << deg;
            EXPECT_NEAR(helix_detail::norm(w) * 180.0 / 3.14159265358979323846, deg, 1e-6) << deg;
        }
    }
}

}  // namespace
