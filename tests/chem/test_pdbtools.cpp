// chem::pdbtools against MDATK's PDB tools (mdatk/pdbtools.py with MDBridge's CH fixes,
// tests/parity/improvements_chem.py): every case of tests/data/chem/pdbtools_cases.json
// (tests/parity/make_pdbtools_cases.py) byte for byte, and counterparts of MDATK's
// tests/test_pdbtools.py (test_pdbtools_mdatk.cpp).

#include <array>
#include <cstdlib>
#include <optional>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "mdbridge/chem/pdbtools.hpp"
#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/json.hpp"
#include "test_support.hpp"

using namespace mdbridge;
using namespace mdbridge::test;

namespace {

const core::Value& cases() {
    static const core::Value data = readJson(kData / "chem" / "pdbtools_cases.json");
    return data;
}

const core::Array& list(const core::Value& value) {
    return *value.get<core::Array>();
}

const core::Array& items(const core::Value& value, std::string_view key) {
    return list(*value.find(key));
}

const std::string& str(const core::Value& value) {
    return *value.get<std::string>();
}

std::int64_t integer(const core::Value& value) {
    return *value.get<std::int64_t>();
}

bool boolean(const core::Value& value) {
    return *value.get<bool>();
}

double hexFloat(const core::Value& v) {
    return std::strtod(str(v).c_str(), nullptr);
}

std::vector<std::string> strings(const core::Value& value) {
    std::vector<std::string> out;
    for (const core::Value& v : list(value)) {
        out.push_back(str(v));
    }
    return out;
}

const std::string& input(std::size_t k) {
    return str(items(cases(), "inputs")[k]);
}

chem::RangeSpec rangeSpec(const core::Value& v) {
    if (v.isNull()) {
        return {};
    }
    if (v.isString()) {
        return chem::RangeSpec(str(v));
    }
    std::vector<chem::RangeSpec::Item> out;
    for (const core::Value& item : list(v)) {
        if (item.isString()) {
            out.emplace_back(str(item));
            continue;
        }
        const core::Array& t = list(item);
        chem::ResidueRange r{str(t[0]), std::nullopt, std::nullopt};
        if (t.size() > 1 && !t[1].isNull()) {
            r.start = integer(t[1]);
        }
        if (t.size() > 2 && !t[2].isNull()) {
            r.end = integer(t[2]);
        }
        out.emplace_back(r);
    }
    return chem::RangeSpec(std::move(out));
}

std::optional<std::vector<std::string>> optionalStrings(const core::Value* v) {
    if (v == nullptr || v->isNull()) {
        return std::nullopt;
    }
    return strings(*v);
}

// run `fn`; an error must be the expected one (or none expected)
template <class Fn>
bool runExpecting(Fn&& fn, const core::Value& c, const std::string& context, std::string_view key = "error") {
    const core::Value* error = c.find(key);
    try {
        fn();
    } catch (const core::ValueError& e) {
        if (error == nullptr) {
            ADD_FAILURE() << context << ": unexpected error " << e.what();
        } else {
            EXPECT_EQ(e.what(), text(*error, "message")) << context;
            EXPECT_EQ(std::string(e.pythonName()), text(*error, "type")) << context;
        }
        return false;
    }
    if (error != nullptr) {
        ADD_FAILURE() << context << ": no error, expected " << text(*error, "message");
        return false;
    }
    return true;
}

core::Value residuesValue(const io::PdbStructure& s) {
    core::Array segments;
    for (const auto& seg : s.segments) {
        core::Array residues;
        for (const io::PdbResidue& r : seg) {
            residues.emplace_back(core::Array{
                core::Value(r.chain), core::Value(r.resseq), core::Value(r.icode), core::Value(r.name),
                r.chirality ? core::Value(*r.chirality) : core::Value(),
                core::Value(r.atoms.empty() ? std::string() : r.atoms.front().record)});
        }
        segments.emplace_back(std::move(residues));
    }
    return core::Value(std::move(segments));
}

core::Value notesValue(const chem::Notes& notes) {
    return core::Value(core::Array(notes.begin(), notes.end()));
}

std::string caseContext(const char* kind, std::size_t n, const core::Value& c) {
    return std::string(kind) + " case " + std::to_string(n) + ": " + core::toJson(c, -1).substr(0, 300);
}

}  // namespace

TEST(PdbTools, SchemeName) {
    for (const core::Value& c : items(cases(), "scheme_name")) {
        const core::Array& e = list(c);
        const std::string context = core::toJson(c, -1);
        const char base = str(e[0])[0];
        const char hand = str(e[1])[0];
        if (e[4].isString()) {
            EXPECT_EQ(chem::schemeName(base, hand, str(e[2]), str(e[3])), str(e[4])) << context;
        } else {
            runExpecting([&] { chem::schemeName(base, hand, str(e[2]), str(e[3])); },
                         core::Value(core::Table{{"error", e[4]}}), context);
        }
    }
}

TEST(PdbTools, ChiralVolumes) {
    const core::Value& chiral = *cases().find("chiral");
    for (const core::Value& c : items(chiral, "volume")) {
        const core::Array& pts = list(list(c)[0]);
        std::array<chem::Vec3, 4> p{};
        for (std::size_t k = 0; k < 4; ++k) {
            for (std::size_t i = 0; i < 3; ++i) {
                p[k][i] = *list(pts[k])[i].get<double>();
            }
        }
        EXPECT_EQ(chem::chiralVolume(p[0], p[1], p[2], p[3]), hexFloat(list(c)[1])) << core::toJson(c, -1);
    }
    for (const core::Value& c : items(chiral, "sugar")) {
        chem::AtomCoords coords;
        for (const auto& [name, xyz] : *c.find("coords")->get<core::Table>()) {
            coords.emplace_back(name, chem::Vec3{hexFloat(list(xyz)[0]), hexFloat(list(xyz)[1]), hexFloat(list(xyz)[2])});
        }
        const std::optional<char> base =
            boolean(*c.find("base_arg")) ? std::optional<char>(text(c, "base")[0]) : std::nullopt;
        const chem::SugarChirality got = chem::sugarChirality(coords, base);
        const std::string context = text(c, "base") + " " + text(c, "label") + (boolean(*c.find("mirror")) ? " L" : " D");
        const core::Value& call = *c.find("call");
        if (call.isNull()) {
            EXPECT_FALSE(got.call.has_value()) << context;
        } else {
            EXPECT_EQ(got.call, str(call)[0]) << context;
        }
        const core::Array& volumes = items(c, "volumes");
        ASSERT_EQ(got.volumes.size(), volumes.size()) << context;
        for (std::size_t k = 0; k < volumes.size(); ++k) {
            EXPECT_EQ(got.volumes[k].first, str(list(volumes[k])[0])) << context;
            EXPECT_EQ(got.volumes[k].second, hexFloat(list(volumes[k])[1])) << context;
        }
    }
}

TEST(PdbTools, Ranges) {
    std::size_t n = 0;
    for (const core::Value& c : items(cases(), "ranges")) {
        const std::string context = caseContext("range", n++, c);
        std::vector<chem::ResidueRange> got;
        if (!runExpecting([&] { got = chem::parseRanges(rangeSpec(*c.find("spec"))); }, c, context)) {
            continue;
        }
        core::Array value;
        for (const chem::ResidueRange& r : got) {
            value.emplace_back(core::Array{core::Value(r.chain), r.start ? core::Value(*r.start) : core::Value(),
                                           r.end ? core::Value(*r.end) : core::Value()});
        }
        EXPECT_EQ(core::Value(std::move(value)), *c.find("ranges")) << context;
    }
}

namespace {

// one function alone on a read input, as make_pdbtools_cases.op_case runs it
chem::Notes runOperation(io::PdbStructure& s, const std::string& op, const core::Value& args, core::Value& value) {
    const auto arg = [&](std::string_view key) -> const core::Value& { return *args.find(key); };
    if (op == "assign_chirality") {
        return chem::assignChirality(s, boolean(arg("detect")), rangeSpec(arg("l_ranges")));
    }
    if (op == "convert_names") {
        return chem::convertNames(s, text(args, "scheme"), rangeSpec(arg("only")), boolean(arg("terminals")));
    }
    if (op == "set_records") {
        return chem::setRecords(s, text(args, "policy"));
    }
    if (op == "strip_5prime_phosphates") {
        return chem::strip5primePhosphates(s);
    }
    if (op == "renumber_residues") {
        return chem::renumberResidues(s, integer(arg("start")), boolean(arg("per_chain")), boolean(arg("keep_gaps")));
    }
    if (op == "seqres_records") {
        const std::vector<std::string> lines = chem::seqresRecords(s);
        value = core::Value(core::Array(lines.begin(), lines.end()));
        return {};
    }
    if (op == "fasta") {
        value = core::Value(chem::fasta(s, text(args, "title")));
        return {};
    }
    if (op == "summary") {
        value = chem::summary(s);
        return {};
    }
    if (op == "split_chains") {
        return chem::splitChains(s, rangeSpec(arg("at")));
    }
    if (op == "combine_chains") {
        const core::Value& into = arg("into");
        return chem::combineChains(s, strings(arg("chains")),
                                   into.isNull() ? std::nullopt : std::optional<std::string>(str(into)));
    }
    if (op == "order_chains") {
        return chem::orderChains(s, optionalStrings(&arg("order")));
    }
    if (op == "restore_chains") {
        if (const core::Value* ref = args.find("reference")) {
            const io::PdbStructure reference = io::readPdbText(input(static_cast<std::size_t>(integer(*ref))));
            return chem::restoreChainsFromReference(s, reference, boolean(arg("copy_numbers")));
        }
        return chem::restoreChains(s, optionalStrings(args.find("chain_ids")));
    }
    ADD_FAILURE() << "unknown operation " << op;
    return {};
}

}  // namespace

TEST(PdbTools, Operations) {
    std::size_t n = 0;
    for (const core::Value& c : items(cases(), "operations")) {
        const std::string context = caseContext("operation", n++, c);
        const std::string& op = text(c, "op");
        const std::string& in = input(static_cast<std::size_t>(integer(*c.find("input"))));
        const core::Value& args = *c.find("args");
        if (op == "to_lrna") {
            chem::CleanResult r;
            if (runExpecting([&] {
                    r = chem::toLrna(chem::PdbSource::text(in), rangeSpec(*args.find("ranges")), text(args, "scheme"),
                                     true, boolean(*args.find("terminals")));
                }, c, context)) {
                EXPECT_EQ(r.text, text(c, "text")) << context;
                EXPECT_EQ(notesValue(r.notes), *c.find("notes")) << context;
            }
            continue;
        }
        io::PdbStructure s;
        chem::Notes notes;
        core::Value value;
        if (!runExpecting([&] {
                s = io::readPdbText(in);
                if (op != "assign_chirality") {
                    chem::assignChirality(s);
                }
                notes = runOperation(s, op, args, value);
            }, c, context)) {
            continue;
        }
        if (const core::Value* expected = c.find("value")) {
            EXPECT_EQ(value, *expected) << context;
        } else {
            EXPECT_EQ(notesValue(notes), *c.find("notes")) << context;
        }
        EXPECT_EQ(residuesValue(s), *c.find("residues")) << context;
        if (c.find("write") || c.find("write_error")) {
            std::string written;
            if (runExpecting([&] { written = chem::writePdb(s); }, c, context, "write_error")) {
                EXPECT_EQ(written, text(c, "write")) << context;
            }
        }
    }
}

namespace {

chem::CleanOptions cleanOptions(const core::Value& options) {
    chem::CleanOptions o;
    for (const auto& [key, v] : *options.get<core::Table>()) {
        const core::Value& entry = v;
        const auto optionalString = [&entry]() {
            return entry.isNull() ? std::nullopt : std::optional<std::string>(str(entry));
        };
        if (key == "scheme") {
            o.scheme = optionalString();
        } else if (key == "records") {
            o.records = optionalString();
        } else if (key == "strip_conect") {
            o.stripConect = boolean(v);
        } else if (key == "detect_l") {
            o.detectL = boolean(v);
        } else if (key == "l_ranges") {
            o.lRanges = rangeSpec(v);
        } else if (key == "chain_reference") {
            o.chainReference = chem::PdbSource::text(input(static_cast<std::size_t>(integer(v))));
        } else if (key == "chain_ids") {
            o.chainIds = optionalStrings(&v);
        } else if (key == "assign_chains") {
            o.assignChains = boolean(v);
        } else if (key == "copy_numbers") {
            o.copyNumbers = boolean(v);
        } else if (key == "split_segments") {
            o.splitSegments = boolean(v);
        } else if (key == "split_at") {
            o.splitAt = rangeSpec(v);
        } else if (key == "combine") {
            o.combine = v.isNull() ? std::vector<std::string>{} : strings(v);
        } else if (key == "combine_into") {
            o.combineInto = optionalString();
        } else if (key == "sort_chains") {
            if (v.isArray()) {
                o.chainOrder = strings(v);
            } else {
                o.sortChains = boolean(v);
            }
        } else if (key == "renumber") {
            o.renumber = boolean(v);
        } else if (key == "start") {
            o.start = integer(v);
        } else if (key == "per_chain") {
            o.perChain = boolean(v);
        } else if (key == "keep_gaps") {
            o.keepGaps = boolean(v);
        } else if (key == "seqres") {
            o.seqres = boolean(v);
        } else if (key == "make_fasta") {
            o.makeFasta = boolean(v);
        } else if (key == "title") {
            o.title = optionalString();
        } else if (key == "terminal_names") {
            o.terminalNames = boolean(v);
        } else if (key == "strip_5p") {
            o.strip5p = boolean(v);
        } else {
            ADD_FAILURE() << "unknown clean option " << key;
        }
    }
    return o;
}

}  // namespace

TEST(PdbTools, Clean) {
    std::size_t n = 0;
    for (const core::Value& c : items(cases(), "clean")) {
        const std::string context = caseContext("clean", n++, c);
        const std::string& in = input(static_cast<std::size_t>(integer(*c.find("input"))));
        chem::CleanResult r;
        if (!runExpecting([&] { r = chem::clean(chem::PdbSource::text(in), cleanOptions(*c.find("options"))); }, c,
                          context)) {
            continue;
        }
        EXPECT_EQ(r.text, text(c, "text")) << context;
        EXPECT_EQ(notesValue(r.notes), *c.find("notes")) << context;
        const core::Value& fasta = *c.find("fasta");
        if (fasta.isNull()) {
            EXPECT_FALSE(r.fasta.has_value()) << context;
        } else {
            EXPECT_EQ(r.fasta.value_or("<none>"), str(fasta)) << context;
        }
    }
}
