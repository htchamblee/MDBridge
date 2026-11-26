#include "mdbridge/chem/generator_api.hpp"

#include <algorithm>
#include <string>

#include "mdbridge/core/python.hpp"
#include "mdbridge/core/text.hpp"

namespace mdbridge::chem::generator {

std::vector<ResidueRange> parseRanges(const std::optional<std::string>& text) {
    return chem::parseRanges(text && !text->empty() ? RangeSpec(*text) : RangeSpec());
}

Model readModel(const std::filesystem::path& path, bool assign, const std::optional<std::string>& lRanges,
                const std::optional<std::string>& scheme) {
    io::PdbStructure s = PdbSource::file(path).read();
    Model out;
    out.header = s.header;
    out.readNotes = s.notes;
    // _nucleotide of the names as read (refine.py looks at them before the chirality is assigned)
    std::vector<std::optional<io::ResidueName>> read;
    for (const io::PdbResidue* r : s.residues()) read.push_back(nucleotide(*r));
    if (assign) {
        out.chiralityNotes =
            chem::assignChirality(s, true, lRanges && !lRanges->empty() ? RangeSpec(*lRanges) : RangeSpec());
    }
    if (scheme) {
        out.nameNotes = convertNames(s, *scheme);
    }
    std::size_t i = 0;
    for (const io::PdbResidue* r : s.residues()) {
        ModelResidue m;
        m.chain = r->chain;
        m.number = r->resseq;
        m.icode = r->icode;
        m.name = r->name;
        m.label = r->label();
        m.nucleotide = read[i++];
        m.chirality = residueChirality(*r);
        for (const io::PdbAtom& a : r->atoms) m.atoms.push_back(a.name);
        out.residues.push_back(std::move(m));
    }
    return out;
}

std::vector<std::string> codesInModel(const Model& model, bool variantsOnly) {
    std::vector<std::string> out;
    for (const Component& c : dictionaryComponents()) {
        if (variantsOnly && !isVariant(c)) continue;
        const bool used = std::any_of(model.residues.begin(), model.residues.end(),
                                      [&c](const ModelResidue& r) { return r.name == c.code; });
        if (used) out.push_back(c.code);
    }
    return out;
}

namespace {

// _remark_residue: (chain, resseq in 4 columns, icode, name in 3) -> nullopt when unreadable
std::optional<RemarkResidue> remarkResidue(std::string_view field) {
    const auto part = [&field](std::size_t start, std::size_t length) {
        return start < field.size() ? field.substr(start, length) : std::string_view();
    };
    const std::optional<std::int64_t> number = core::tryPythonInt(part(1, 4));
    if (!number) return std::nullopt;
    return RemarkResidue{core::pythonStrip(part(0, 1)), *number, core::pythonStrip(part(5, 1)),
                         core::pythonStrip(part(6, 3))};
}

}  // namespace

Annotation readAnnotation(std::span<const std::string> header) {
    Annotation out;
    for (const std::string& line : header) {
        if (line.starts_with("REMARK 950")) out.annotated = true;
        const bool pair = line.starts_with("REMARK 950 PAIR ");
        if (!pair && !line.starts_with("REMARK 950 STACK")) continue;
        const std::string_view text(line);
        const auto field = [&text](std::size_t start, std::size_t end) {
            return start < text.size() ? text.substr(start, end - start) : std::string_view();
        };
        const std::optional<RemarkResidue> a = remarkResidue(field(17, 26));
        const std::optional<RemarkResidue> b = remarkResidue(field(27, 36));
        if (!a || !b) continue;
        RemarkPair item{*a, *b, core::pythonStrip(text.size() > 37 ? text.substr(37) : std::string_view())};
        (pair ? out.pairs : out.stacks).push_back(std::move(item));
    }
    return out;
}

}  // namespace mdbridge::chem::generator
