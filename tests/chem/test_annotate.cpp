// chem::annotate against MDATK's base-pair annotation (mdatk/annotate.py with MDBridge's D13 / CH-14,
// tests/parity/improvements_chem.py): every case of tests/data/chem/annotate_cases.json
// (tests/parity/make_annotate_cases.py) byte for byte, and counterparts of MDATK v0.0.4.1's
// tests/test_annotate.py (8 tests, one TEST each, named after it; the fixtures pdb_text, duplex,
// nucleotide, partner and atom_line of MDATK's tests are ported). MDBridge has no Barnaba (CH-14):
// the Barnaba tests check that engine "barnaba" runs the builtin detector with the warning;
// test_cli_annotate_logs_step is checked through the library call `pdb annotate` makes (the command
// and its step log: W2-C).

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <initializer_list>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "mdbridge/chem/annotate.hpp"
#include "mdbridge/chem/pdbtools.hpp"
#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/core/linalg.hpp"
#include "mdbridge/io/pdb.hpp"
#include "test_support.hpp"

using namespace mdbridge;
using namespace mdbridge::test;

namespace {

// --- annotate_cases.json -----------------------------------------------------------------------

const core::Value& cases() {
    static const core::Value data = readJson(kData / "chem" / "annotate_cases.json");
    return data;
}

const core::Array& list(const core::Value& value) {
    return *value.get<core::Array>();
}

const std::string& str(const core::Value& value) {
    return *value.get<std::string>();
}

std::vector<std::string> strings(const core::Value& value) {
    std::vector<std::string> out;
    for (const core::Value& v : list(value)) {
        out.push_back(str(v));
    }
    return out;
}

const std::string& input(std::size_t k) {
    return str(list(*cases().find("inputs"))[k]);
}

chem::AnnotateOptions options(const core::Value& o) {
    chem::AnnotateOptions out;
    if (const core::Value* v = o.find("pairs")) {
        out.pairs = *v->get<bool>();
    }
    if (const core::Value* v = o.find("stacking")) {
        out.stacking = *v->get<bool>();
    }
    if (const core::Value* v = o.find("dotbracket")) {
        out.dotbracket = *v->get<bool>();
    }
    if (const core::Value* v = o.find("engine")) {
        out.engine = str(*v);
    }
    if (const core::Value* v = o.find("wobble")) {
        out.wobble = *v->get<bool>();
    }
    if (const core::Value* v = o.find("cutoff")) {
        if (const auto* d = v->get<double>()) {
            out.cutoff = *d;
        } else {
            out.cutoff = static_cast<double>(*v->get<std::int64_t>());
        }
    }
    return out;
}

TEST(Annotate, EveryCaseOfAnnotateCases) {
    const core::Array& all = list(*cases().find("cases"));
    ASSERT_EQ(all.size(), 216u);
    for (std::size_t n = 0; n < all.size(); ++n) {
        const core::Value& c = all[n];
        const auto k = static_cast<std::size_t>(*c.find("input")->get<std::int64_t>());
        const std::string context = "case " + std::to_string(n) + " (input " + std::to_string(k) + ", " +
                                    core::pythonJsonCompact(*c.find("options")) + ")";
        const chem::AnnotateOptions o = options(*c.find("options"));
        if (const core::Value* error = c.find("error")) {
            try {
                chem::annotate(chem::PdbSource::text(input(k)), o);
                ADD_FAILURE() << context << ": no error";
            } catch (const core::ValueError& e) {
                EXPECT_EQ(std::string(e.pythonName()), str(*error->find("type"))) << context;
                EXPECT_EQ(std::string(e.what()), str(*error->find("message"))) << context;
            }
            continue;
        }
        const chem::AnnotationResult res = chem::annotate(chem::PdbSource::text(input(k)), o);
        EXPECT_EQ(res.records, strings(*c.find("records"))) << context;
        EXPECT_EQ(res.notes, strings(*c.find("notes"))) << context;
        EXPECT_EQ(res.engine, str(*c.find("engine"))) << context;
        if (c.find("insert")) {
            EXPECT_EQ(res.text, io::insertHeaderRecords(input(k), res.records)) << context;
        } else {
            EXPECT_EQ(res.text, str(*c.find("text"))) << context;
        }
    }
}

TEST(Annotate, StructureCasesOfAnnotateCases) {
    for (const core::Value& c : list(*cases().find("structure"))) {
        const auto k = static_cast<std::size_t>(*c.find("input")->get<std::int64_t>());
        const io::PdbStructure s = io::readPdbText(input(k));
        const io::PdbStructure before = s;
        const chem::AnnotationResult res = chem::annotate(s);
        EXPECT_EQ(res.text, str(*c.find("text"))) << "input " << k;
        EXPECT_EQ(res.notes, strings(*c.find("notes"))) << "input " << k;
        EXPECT_EQ(res.records, strings(*c.find("records"))) << "input " << k;
        EXPECT_EQ(s, before);  // read_pdb of a Structure is a deep copy
    }
}

TEST(Annotate, RecordHelpers) {
    io::PdbResidue r;
    r.name = "G5";
    r.chain = "";
    r.resseq = -12;
    EXPECT_EQ(chem::annotationResidueField(r), "  -12 G5 ");
    r.chain = "AB";
    r.icode = "x";
    r.name = "LG5X";
    r.resseq = 12345;
    EXPECT_EQ(chem::annotationResidueField(r), "A12345xLG5X");
    EXPECT_EQ(chem::annotationWrapped("SEQ", ""), std::vector<std::string>{"REMARK 950 SEQ   "});
    EXPECT_EQ(chem::annotationWrapped("CHAINS", "A"), std::vector<std::string>{"REMARK 950 CHAINS A"});
    const std::vector<std::vector<std::size_t>> strands = {{0, 1}, {2}, {3, 4}, {5}};
    using Groups = std::vector<std::vector<std::size_t>>;  // MDATK's _strand_groups of the same input
    EXPECT_EQ(chem::annotationStrandGroups(strands, {{0, 5}, {3, 4}}), (Groups{{0, 3}, {1}, {2}}));
    EXPECT_EQ(chem::annotationStrandGroups(strands, {{2, 5}, {0, 3}}), (Groups{{0, 2}, {1, 3}}));
    chem::AnnotateOptions unknown;
    unknown.engine = "rnaview";  // refused before the text is read
    EXPECT_THROW(chem::annotate(chem::PdbSource::text("REMARK only\n"), unknown), core::ValueError);
}

TEST(Annotate, ReadsFiles) {
    const TempDir dir;
    const std::filesystem::path path = dir.path() / "in.pdb";
    std::string crlf;
    for (const char ch : input(0)) {
        crlf += ch == '\n' ? std::string("\r\n") : std::string(1, ch);
    }
    writeFile(path, crlf);
    const chem::AnnotationResult file = chem::annotate(chem::PdbSource::file(path));
    const chem::AnnotationResult text = chem::annotate(chem::PdbSource::text(input(0)));
    EXPECT_EQ(file.text, text.text);
    EXPECT_EQ(file.notes, text.notes);
    EXPECT_THROW(chem::annotate(chem::PdbSource::file(dir.path() / "missing.pdb")),
                 std::filesystem::filesystem_error);
}

// --- fixtures of MDATK's tests (test_pdbtools.py BACKBONE / BASES / atom_line, test_basepairs.py
// nucleotide / partner / moved, test_annotate.py pdb_text / duplex) ------------------------------

using Vec3 = std::array<double, 3>;
using Atoms = std::vector<std::pair<std::string, Vec3>>;  // a dict in insertion order

struct AtomXyz {
    const char* name;
    double x, y, z;
};

const std::vector<AtomXyz> kBackbone = {
    {"P", 2.314, 1.777, -0.478},  {"OP1", 1.933, 3.073, -1.082}, {"OP2", 3.157, 0.872, -1.292},
    {"O5'", 2.994, 2.142, 0.923}, {"C5'", 2.156, 2.552, 2.020},  {"H5''", 1.605, 3.450, 1.741},
    {"C4'", 3.008, 2.848, 3.243}, {"O4'", 3.228, 1.605, 3.977},  {"C1'", 4.559, 1.151, 3.778},
    {"C3'", 4.416, 3.378, 2.961}, {"C2'", 5.230, 2.096, 2.779},  {"O2'", 6.273, 2.050, 3.757},
    {"HO2'", 6.759, 1.235, 3.611}, {"O3'", 5.030, 4.129, 3.998},
};

const std::vector<AtomXyz>& baseAtoms(char base) {
    static const std::vector<AtomXyz> a = {
        {"N9", 4.539, -0.297, 3.316}, {"C8", 4.337, -0.846, 2.077}, {"N7", 4.453, -2.141, 2.049},
        {"C5", 4.754, -2.480, 3.363}, {"C6", 4.998, -3.708, 3.989}, {"N6", 4.977, -4.883, 3.346},
        {"N1", 5.265, -3.689, 5.302}, {"C2", 5.285, -2.520, 5.935}, {"N3", 5.073, -1.309, 5.456},
        {"C4", 4.806, -1.356, 4.133}};
    static const std::vector<AtomXyz> c = {
        {"N1", 4.501, -0.268, 3.327}, {"C6", 4.258, -0.582, 2.026}, {"C5", 4.206, -1.878, 1.616},
        {"C4", 4.416, -2.881, 2.614}, {"N4", 4.379, -4.159, 2.286}, {"N3", 4.653, -2.564, 3.884},
        {"C2", 4.701, -1.262, 4.273}, {"O2", 4.919, -0.928, 5.447}};
    static const std::vector<AtomXyz> g = {
        {"N9", 4.500, -0.268, 3.327}, {"C8", 4.268, -0.773, 2.063}, {"N7", 4.279, -2.081, 2.002},
        {"C5", 4.535, -2.475, 3.310}, {"C6", 4.665, -3.772, 3.873}, {"O6", 4.578, -4.865, 3.316},
        {"N1", 4.926, -3.709, 5.247}, {"C2", 5.048, -2.547, 5.982}, {"N2", 5.299, -2.702, 7.289},
        {"N3", 4.928, -1.333, 5.453}, {"C4", 4.674, -1.377, 4.118}};
    static const std::vector<AtomXyz> u = {
        {"N1", 4.539, -0.307, 3.313}, {"C6", 4.317, -0.673, 2.011}, {"C5", 4.360, -1.958, 1.634},
        {"C4", 4.643, -2.997, 2.590}, {"O4", 4.703, -4.199, 2.337}, {"N3", 4.853, -2.538, 3.874},
        {"C2", 4.815, -1.224, 4.292}, {"O2", 5.015, -0.915, 5.455}};
    return base == 'A' ? a : base == 'C' ? c : base == 'G' ? g : u;
}

Atoms moved(const Atoms& res, const Vec3& shift) {
    Atoms out = res;
    for (auto& [name, p] : out) {
        for (std::size_t d = 0; d < 3; ++d) {
            p[d] += shift[d];
        }
    }
    return out;
}

Atoms nucleotide(char base, const Vec3& shift = {0.0, 0.0, 0.0}) {
    Atoms out;
    const auto add = [&out](const AtomXyz& a) { out.emplace_back(a.name, Vec3{a.x, a.y, a.z}); };
    for (std::size_t k = 0; k < 9; ++k) {
        add(kBackbone[k]);
    }
    for (const AtomXyz& a : baseAtoms(base)) {
        add(a);
    }
    for (std::size_t k = 9; k < kBackbone.size(); ++k) {
        add(kBackbone[k]);
    }
    return moved(out, shift);
}

const Vec3& at(const Atoms& res, const std::string& name) {
    return std::find_if(res.begin(), res.end(), [&name](const auto& a) { return a.first == name; })->second;
}

Vec3 sub(const Vec3& a, const Vec3& b) {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}
Vec3 add(const Vec3& a, const Vec3& b) {
    return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}
Vec3 scale(const Vec3& a, double s) {
    return {a[0] * s, a[1] * s, a[2] * s};
}
double dot(const Vec3& a, const Vec3& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
double norm(const Vec3& a) {
    return std::sqrt(dot(a, a));
}
Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
Vec3 mean(const std::vector<Vec3>& p) {
    Vec3 m{0.0, 0.0, 0.0};
    for (const Vec3& v : p) {
        m = add(m, v);
    }
    return scale(m, 1.0 / static_cast<double>(p.size()));
}

using AtomPairs = std::vector<std::pair<std::string, std::string>>;

AtomPairs pairAtoms(char first, char second) {
    if (first == 'G' && second == 'C') {
        return {{"N1", "N3"}, {"N2", "O2"}, {"O6", "N4"}};
    }
    if (first == 'A' && second == 'U') {
        return {{"N6", "O4"}, {"N1", "N3"}};
    }
    return {{"O6", "N3"}, {"N1", "O2"}};  // G-U
}

// _plane: ring centre and normal of the base atoms
std::pair<Vec3, Vec3> plane(const Atoms& res, char base) {
    std::vector<Vec3> p;
    for (const AtomXyz& a : baseAtoms(base)) {
        p.push_back(at(res, a.name));
    }
    const Vec3 c = mean(p);
    for (Vec3& v : p) {
        v = sub(v, c);
    }
    const core::SvdRows3 svd = core::svdRows3(p);
    return {c, svd.vt[2]};
}

// _kabsch: weighted rotation r and centroids (src -> dst)
struct Kabsch {
    core::Mat3 r{};
    Vec3 cs{};
    Vec3 cd{};
};
Kabsch kabsch(const std::vector<Vec3>& src, const std::vector<Vec3>& dst, const std::vector<double>& w) {
    double total = 0.0;
    for (const double x : w) {
        total += x;
    }
    Kabsch k;
    for (std::size_t n = 0; n < src.size(); ++n) {
        k.cs = add(k.cs, scale(src[n], w[n] / total));
        k.cd = add(k.cd, scale(dst[n], w[n] / total));
    }
    core::Mat3 h{};
    for (std::size_t n = 0; n < src.size(); ++n) {
        const Vec3 a = scale(sub(src[n], k.cs), w[n] / total);
        const Vec3 b = sub(dst[n], k.cd);
        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t j = 0; j < 3; ++j) {
                h[i][j] += a[i] * b[j];
            }
        }
    }
    const core::Svd3 svd = core::svd3x3(h);
    core::Mat3 vu{};  // vt.T @ u.T
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            for (std::size_t m = 0; m < 3; ++m) {
                vu[i][j] += svd.vt[m][i] * svd.u[j][m];
            }
        }
    }
    const double sign = core::det3(vu) < 0.0 ? -1.0 : 1.0;
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            for (std::size_t m = 0; m < 3; ++m) {
                k.r[i][j] += svd.vt[m][i] * (m == 2 ? sign : 1.0) * svd.u[j][m];
            }
        }
    }
    return k;
}

// partner: nucleotide `sbase` coplanar with `first`, its H-bond atoms about d from first's
Atoms partner(const Atoms& first, char fbase, char sbase, double d = 2.9) {
    const AtomPairs atoms = pairAtoms(fbase, sbase);
    const Atoms second = nucleotide(sbase);
    std::vector<Vec3> a;
    std::vector<Vec3> b;
    for (const auto& [x, y] : atoms) {
        a.push_back(at(first, x));
        b.push_back(at(second, y));
    }
    const auto [ring1, n1] = plane(first, fbase);
    const auto [ring2, n2] = plane(second, sbase);
    const Vec3 edgeVector = sub(a.back(), a.front());
    const Vec3 edge = scale(edgeVector, 1.0 / norm(edgeVector));
    Vec3 out = cross(n1, edge);
    out = scale(out, 1.0 / norm(out));
    if (dot(out, sub(mean(a), ring1)) < 0.0) {
        out = scale(out, -1.0);
    }
    std::vector<Vec3> targets;
    for (const Vec3& p : a) {
        targets.push_back(add(p, scale(out, d)));
    }
    const Vec3 centre = add(mean(targets), scale(out, norm(sub(mean(b), ring2))));
    std::optional<std::pair<double, Atoms>> best;
    for (const double sign : {1.0, -1.0}) {
        std::vector<Vec3> src = b;
        src.push_back(ring2);
        src.push_back(add(ring2, n2));
        std::vector<Vec3> dst = targets;
        dst.push_back(centre);
        dst.push_back(add(centre, scale(n1, sign)));
        std::vector<double> w(b.size(), 10.0);
        w.push_back(1.0);
        w.push_back(1.0);
        const Kabsch k = kabsch(src, dst, w);
        Atoms placed;
        for (const auto& [name, v] : second) {
            const Vec3 c = sub(v, k.cs);
            Vec3 p{};
            for (std::size_t i = 0; i < 3; ++i) {
                p[i] = k.r[i][0] * c[0] + k.r[i][1] * c[1] + k.r[i][2] * c[2] + k.cd[i];
            }
            placed.emplace_back(name, p);
        }
        double err = 0.0;
        for (const auto& [x, y] : atoms) {
            err += std::abs(norm(sub(at(first, x), at(placed, y))) - d);
        }
        if (!best || err < best->first) {
            best.emplace(err, std::move(placed));
        }
    }
    return best->second;
}

std::string fmt(const char* format, auto... args) {
    char buf[256];
    std::snprintf(buf, sizeof buf, format, args...);
    return buf;
}

std::string atomLine(int serial, const std::string& name, const std::string& resname, const std::string& chain,
                     int resseq, const Vec3& p) {
    const std::string element = name.substr(name.find_first_not_of("0123456789"), 1);
    const std::string field =
        name.size() == 4 ? name : " " + name + std::string(3 - std::min<std::size_t>(3, name.size()), ' ');
    return fmt("%-6s%5d %s %3s %1s%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s", "ATOM", serial, field.c_str(),
               resname.c_str(), chain.empty() ? " " : chain.c_str(), resseq, p[0], p[1], p[2], element.c_str());
}

struct Residue {
    std::string name;
    int resseq;
    Atoms atoms;
};
using Chains = std::vector<std::pair<std::string, std::vector<Residue>>>;

std::string pdbText(const Chains& chains, const std::vector<std::string>& header = {}) {
    std::vector<std::string> lines = header;
    int serial = 0;
    for (const auto& [chain, residues] : chains) {
        for (const Residue& r : residues) {
            for (const auto& [name, p] : r.atoms) {
                lines.push_back(atomLine(++serial, name, r.name, chain, r.resseq, p));
            }
        }
        ++serial;
        lines.push_back(fmt("TER   %5d      %3s %1s%4d", serial, residues.back().name.c_str(),
                            chain.empty() ? " " : chain.c_str(), residues.back().resseq));
    }
    lines.emplace_back("END");
    std::string text;
    for (const std::string& line : lines) {
        text += line + "\n";
    }
    return text;
}

Chains duplex(bool lRna = false) {
    const Atoms g = nucleotide('G');
    const Atoms c = partner(g, 'G', 'C');
    const Atoms a = nucleotide('A', {0.0, 0.0, 25.0});
    const Atoms u = partner(a, 'A', 'U');
    const Atoms lone = nucleotide('G', {0.0, 60.0, 0.0});
    const std::array<std::string, 4> names = lRna ? std::array<std::string, 4>{"LG5", "LA3", "LU5", "LC3"}
                                                  : std::array<std::string, 4>{"G5", "A3", "U5", "C3"};
    return {{"A", {{names[0], 1, g}, {names[1], 2, a}}},
            {"B", {{names[2], 3, u}, {names[3], 4, c}}},
            {"C", {{"G", 5, lone}, {"MG", 6, {{"MG", {0.0, 80.0, 0.0}}}}}}};
}

const std::vector<std::string> kHeader = {"HEADER    RNA                                     01-JAN-26   XXXX",
                                          "REMARK   2 RESOLUTION. 3.00 ANGSTROMS.",
                                          "REMARK 950 OLD ANNOTATION LINE"};

std::vector<std::string> splitLines(const std::string& text) {
    std::vector<std::string> out;
    std::size_t first = 0;
    for (std::size_t k = 0; k < text.size(); ++k) {
        if (text[k] == '\n') {
            out.push_back(text.substr(first, k - first));
            first = k + 1;
        }
    }
    if (first < text.size()) {
        out.push_back(text.substr(first));
    }
    return out;
}

std::vector<std::string> coordinateColumns(const std::string& text) {
    std::vector<std::string> out;
    for (const std::string& line : splitLines(text)) {
        if (line.starts_with("ATOM") || line.starts_with("HETATM")) {
            out.push_back(line.substr(30, 24));
        }
    }
    return out;
}

bool contains(const std::vector<std::string>& lines, const std::string& line) {
    return std::find(lines.begin(), lines.end(), line) != lines.end();
}

bool anyContains(const std::vector<std::string>& lines, const std::string& part) {
    return std::any_of(lines.begin(), lines.end(),
                       [&part](const std::string& l) { return l.find(part) != std::string::npos; });
}

chem::AnnotationResult annotateText(const std::string& text, chem::AnnotateOptions o = {}) {
    return chem::annotate(chem::PdbSource::text(text), o);
}

chem::AnnotateOptions engine(const std::string& name) {
    chem::AnnotateOptions o;
    o.engine = name;
    return o;
}

// --- counterparts of tests/test_annotate.py ------------------------------------------------------

TEST(AnnotateMdatk, test_builtin_records_and_placement) {
    std::string text = pdbText(duplex(), kHeader);
    const chem::AnnotationResult res = annotateText(text, engine("auto"));  // auto: the builtin engine
    const std::vector<std::string> lines = splitLines(res.text);
    EXPECT_EQ(std::vector<std::string>(lines.begin(), lines.begin() + 2),
              std::vector<std::string>(kHeader.begin(), kHeader.begin() + 2));
    EXPECT_EQ(res.text.find("REMARK 950 OLD ANNOTATION LINE"), std::string::npos);
    const auto firstAtom =
        std::find_if(lines.begin(), lines.end(), [](const std::string& l) { return l.starts_with("ATOM"); });
    EXPECT_EQ(std::vector<std::string>(lines.begin() + 2, firstAtom), res.records);
    EXPECT_EQ(res.records[0], "REMARK 950 MDBridge BASE-PAIR ANNOTATION");  // D13
    EXPECT_TRUE(res.records[1].starts_with("REMARK 950 ENGINE builtin"));
    EXPECT_TRUE(contains(res.records, "REMARK 950 PAIRS 2 (WATSON-CRICK 2, WOBBLE 0, OTHER 0)"));
    EXPECT_TRUE(contains(res.records, "REMARK 950 PAIR  A   1 G5  B   4 C3  WCc"));
    EXPECT_TRUE(contains(res.records, "REMARK 950 PAIR  A   2 A3  B   3 U5  WCc"));
    EXPECT_FALSE(anyContains(res.records, "STACK"));
    EXPECT_TRUE(anyContains(res.notes, "no stacking"));
    EXPECT_FALSE(anyContains(res.notes, "pip install"));  // CH-14
    const std::vector<std::string> tail(res.records.end() - 6, res.records.end());
    EXPECT_EQ(tail, (std::vector<std::string>{"REMARK 950 CHAINS A,B", "REMARK 950 SEQ   GA&UC",
                                              "REMARK 950 DBN   ((&))", "REMARK 950 CHAINS C", "REMARK 950 SEQ   G",
                                              "REMARK 950 DBN   ."}));
    EXPECT_EQ(coordinateColumns(res.text), coordinateColumns(text));
    for (const std::string& line : lines) {
        EXPECT_LE(line.size(), 80u) << line;
    }
    const chem::AnnotationResult again = annotateText(res.text, engine("builtin"));
    EXPECT_EQ(again.text, res.text);
}

TEST(AnnotateMdatk, test_options_and_l_rna) {
    chem::AnnotateOptions o = engine("builtin");
    o.pairs = false;
    o.dotbracket = false;
    chem::AnnotationResult res = annotateText(pdbText(duplex(true)), o);
    EXPECT_FALSE(anyContains(res.records, " PAIR ") || anyContains(res.records, "DBN"));
    EXPECT_TRUE(contains(res.records, "REMARK 950 PAIRS 2 (WATSON-CRICK 2, WOBBLE 0, OTHER 0)"));
    res = annotateText(pdbText(duplex(true)), engine("builtin"));
    EXPECT_TRUE(contains(res.records, "REMARK 950 PAIR  A   1 LG5 B   4 LC3 WCc"));
    try {
        annotateText(pdbText(duplex()), engine("rnaview"));
        ADD_FAILURE() << "no error";
    } catch (const core::ValueError& e) {
        EXPECT_NE(std::string(e.what()).find("unknown annotation engine"), std::string::npos);
    }
}

TEST(AnnotateMdatk, test_dotbracket_wrapped_to_60_characters) {
    const Atoms g = nucleotide('G');
    std::vector<Residue> residues;
    for (int k = 0; k < 70; ++k) {
        residues.push_back({"U", k + 1, moved(nucleotide('U'), {0.0, 0.0, 8.0 * k})});
    }
    residues[0] = {"G", 1, g};
    residues[69] = {"C", 70, partner(g, 'G', 'C')};
    const chem::AnnotationResult res = annotateText(pdbText({{"A", residues}}), engine("builtin"));
    std::vector<std::size_t> seq;
    std::vector<std::string> dbn;
    for (const std::string& r : res.records) {
        if (r.starts_with("REMARK 950 SEQ")) {
            seq.push_back(r.size() - 17);
        }
        if (r.starts_with("REMARK 950 DBN")) {
            dbn.push_back(r);
        }
    }
    EXPECT_EQ(seq, (std::vector<std::size_t>{60, 10}));
    ASSERT_EQ(dbn.size(), 2u);
    EXPECT_EQ(dbn[0].size() - 17, 60u);
    EXPECT_EQ(dbn[1].size() - 17, 10u);
    EXPECT_EQ(dbn[0][17], '(');
    EXPECT_TRUE(dbn[1].ends_with(")"));
}

// MDATK (with a stand-in Barnaba) writes Barnaba's records; MDBridge has no Barnaba (CH-14): "auto"
// and "barnaba" give the builtin engine's records, "barnaba" with the warning, never STACK records.
TEST(AnnotateMdatk, test_barnaba_records) {
    chem::AnnotationResult res = annotateText(pdbText(duplex()), engine("auto"));
    EXPECT_TRUE(res.records[1].starts_with("REMARK 950 ENGINE builtin"));
    EXPECT_FALSE(contains(res.notes, std::string(chem::kBarnabaNote)));
    res = annotateText(pdbText(duplex()), engine("barnaba"));
    EXPECT_TRUE(res.records[1].starts_with("REMARK 950 ENGINE builtin"));
    EXPECT_EQ(res.engine, "builtin");
    EXPECT_TRUE(contains(res.records, "REMARK 950 PAIR  A   1 G5  B   4 C3  WCc"));
    EXPECT_TRUE(contains(res.records, "REMARK 950 PAIRS 2 (WATSON-CRICK 2, WOBBLE 0, OTHER 0)"));
    EXPECT_TRUE(contains(res.records, "REMARK 950 CHAINS A,B") && contains(res.records, "REMARK 950 DBN   ((&))"));
    EXPECT_FALSE(anyContains(res.records, "STACK"));
    EXPECT_EQ(res.notes, (chem::Notes{std::string(chem::kBarnabaNote),
                                      "The builtin engine finds canonical base pairs only (no stacking)",
                                      "Annotated 2 base pair(s) with the builtin engine"}));
    chem::AnnotateOptions o = engine("barnaba");
    o.stacking = false;
    res = annotateText(pdbText(duplex()), o);
    EXPECT_FALSE(anyContains(res.records, "STACK"));
    EXPECT_EQ(res.notes,
              (chem::Notes{std::string(chem::kBarnabaNote), "Annotated 2 base pair(s) with the builtin engine"}));
}

// The command and its step log come with W2-C; here the library call `pdb annotate in.pdb -o out.pdb
// --no-stacking` makes, and `--engine barnaba` succeeding with the warning (MDATK: exit 1, not installed).
TEST(AnnotateMdatk, test_cli_annotate_logs_step) {
    const TempDir dir;
    const std::filesystem::path pdb = dir.path() / "in.pdb";
    writeFile(pdb, pdbText(duplex()));
    chem::AnnotateOptions o;
    o.stacking = false;
    const chem::AnnotationResult res = chem::annotate(chem::PdbSource::file(pdb), o);
    EXPECT_NE(res.text.find("REMARK 950 DBN   ((&))"), std::string::npos);
    EXPECT_EQ(io::readPdbText(res.text).residues().size(), 6u);
    EXPECT_EQ(res.engine, "builtin");
    const chem::AnnotationResult barnaba = chem::annotate(chem::PdbSource::file(pdb), engine("barnaba"));
    EXPECT_EQ(barnaba.notes.front(), chem::kBarnabaNote);
}

// MDATK skips this without Barnaba; MDBridge's barnaba engine (the builtin detector) gives the
// pairs real Barnaba finds.
TEST(AnnotateMdatk, test_real_barnaba_annotation) {
    const chem::AnnotationResult res = annotateText(pdbText(duplex(false)), engine("barnaba"));
    EXPECT_TRUE(contains(res.records, "REMARK 950 PAIR  A   1 G5  B   4 C3  WCc"));
    EXPECT_TRUE(contains(res.records, "REMARK 950 PAIR  A   2 A3  B   3 U5  WCc"));
    EXPECT_TRUE(contains(res.records, "REMARK 950 DBN   ((&))"));
}

TEST(AnnotateMdatk, test_strands_from_ter_segments_ions_and_verbatim_models) {
    Chains blank = duplex();
    for (auto& segment : blank) {  // Amber-style: no chain ids, TER between strands
        segment.first.clear();
    }
    blank[2].second = {{"G", 5, blank[2].second[0].atoms}, {"LA", 9, {{"LA", {0.0, 80.0, 0.0}}}}};
    const std::string text = pdbText(blank);
    const chem::AnnotationResult res = annotateText(text, engine("builtin"));
    const std::vector<std::string> tail(res.records.end() - 6, res.records.end());
    EXPECT_EQ(tail, (std::vector<std::string>{"REMARK 950 CHAINS _", "REMARK 950 SEQ   GA&UC", "REMARK 950 DBN   ((&))",
                                              "REMARK 950 CHAINS _", "REMARK 950 SEQ   G", "REMARK 950 DBN   ."}));
    EXPECT_FALSE(anyContains(res.notes, "numbering gaps"));  // LA is no nucleotide
    std::vector<std::string> body;
    for (const std::string& line : splitLines(text)) {
        if (line.starts_with("ATOM") || line.starts_with("HETATM") || line.starts_with("TER")) {
            body.push_back(line);
        }
    }
    std::vector<std::string> modelLines = {"MODEL        1"};
    modelLines.insert(modelLines.end(), body.begin(), body.end());
    modelLines.insert(modelLines.end(), {"ENDMDL", "MODEL        2"});
    modelLines.insert(modelLines.end(), body.begin(), body.end());
    modelLines.insert(modelLines.end(), {"ENDMDL", "END"});
    std::string models;
    for (const std::string& line : modelLines) {
        models += line + "\n";
    }
    const chem::AnnotationResult out = annotateText(models, engine("builtin"));
    std::vector<std::string> kept;
    for (const std::string& line : splitLines(out.text)) {
        if (!line.starts_with("REMARK 950")) {
            kept.push_back(line);
        }
    }
    EXPECT_EQ(kept, modelLines);
    EXPECT_EQ(splitLines(out.text)[0], "REMARK 950 MDBridge BASE-PAIR ANNOTATION");
    EXPECT_TRUE(anyContains(out.notes, "model 1 of 2"));
}

TEST(AnnotateMdatk, test_clean_drops_stale_annotation_master_and_missing_barnaba) {
    std::string text = pdbText(duplex(), {kHeader[0], kHeader[1]});
    text.replace(text.rfind("END\n"), 4,
                 "MASTER        1    0    0    0    0    0    0    0   10    3   0    0\nEND\n");
    const chem::AnnotationResult ann = annotateText(text, engine("builtin"));
    EXPECT_EQ(ann.text.find("MASTER"), std::string::npos);
    EXPECT_TRUE(anyContains(ann.notes, "MASTER"));
    chem::CleanOptions o;
    o.renumber = true;
    o.start = 100;
    const chem::CleanResult cleaned = chem::clean(chem::PdbSource::text(ann.text), o);
    EXPECT_EQ(cleaned.text.find("REMARK 950"), std::string::npos);
    EXPECT_TRUE(anyContains(cleaned.notes, "REMARK 950"));
    // MDATK: RuntimeError "pip install barnaba" even with nothing to pair; MDBridge (CH-14): the
    // builtin detector with the warning
    const std::string one = pdbText({{"A", {{"G", 1, nucleotide('G')}}}});
    const chem::AnnotationResult res = annotateText(one, engine("barnaba"));
    EXPECT_EQ(res.notes.front(), chem::kBarnabaNote);
    EXPECT_TRUE(contains(res.records, "REMARK 950 PAIRS 0 (WATSON-CRICK 0, WOBBLE 0, OTHER 0)"));
}

}  // namespace
