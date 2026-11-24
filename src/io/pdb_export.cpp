// The PDB files MDATK's frame export writes through MDAnalysis 2.11.0's PDBWriter (AtomGroup.write:
// one frame, reindexed serials, CONECT records of the bonds read from the topology file).

#include "mdbridge/io/pdb_export.hpp"

#include <algorithm>
#include <array>
#include <fstream>
#include <map>
#include <system_error>

#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/core/selections.hpp"
#include "mdbridge/core/text.hpp"
#include "python_text.hpp"
#include "readers.hpp"
#include "text_file.hpp"

namespace mdbridge::io {

namespace {

using namespace detail;

#include "pdb_tables.inc"

template <std::size_t N>
bool contains(const std::string_view (&list)[N], std::string_view value) {
    return std::find(std::begin(list), std::end(list), value) != std::end(list);
}

template <std::size_t N>
bool containsPair(const std::pair<std::string_view, std::string_view> (&list)[N],
                  std::string_view first, std::string_view second) {
    return std::any_of(std::begin(list), std::end(list), [&](const auto& pair) {
        return pair.first == first && pair.second == second;
    });
}
// The atom-name field (columns 13-16). CLEAN ROOM: reimplement from scratchpad/cleanroom/pdb-atom-name.md
U32 atomNameField(const std::string& name, const std::string& resname) {
    (void)name; (void)resname;
    throw std::logic_error("clean room: atom name field not implemented");
}

// util.ltruncate_int: the last `digits` characters of str(value) as an int
std::int64_t lastDigits(std::int64_t value, std::size_t digits) {
    const std::string text = std::to_string(value);
    return core::pythonInt(text.size() > digits ? text.substr(text.size() - digits) : text);
}

std::string chainIdOf(const std::string& chain) {
    const U32 text = toU32(chain);
    if (text.size() != 1 || !isAlnum(text[0])) {
        return "X";
    }
    return chain;
}

std::string formalCharge(std::int64_t charge) {
    if (charge < -9) {
        throw core::ValueError("formal charge < -9 is not supported by PDB standard");
    }
    if (charge > 9) {
        throw core::ValueError("formal charge > 9 is not supported by PDB standard");
    }
    if (charge == 0) {
        return "";
    }
    return std::to_string(charge < 0 ? -charge : charge) + (charge < 0 ? "-" : "+");
}

std::string line(std::string_view text) {
    return std::string(text) + "\n";
}

// PDBWriter.has_valid_coordinates: numpy compares the float32 positions with the limits as float32
bool validCoordinates(const std::vector<float>& positions, const core::AtomGroup& atoms) {
    const float low = static_cast<float>(-999.9995);
    const float high = static_cast<float>(9999.9995);
    for (const std::int32_t atom : atoms.indices()) {
        for (std::size_t i = 0; i < 3; ++i) {
            const float x = positions[3 * static_cast<std::size_t>(atom) + i];
            if (!(low < x && x <= high)) {
                return false;
            }
        }
    }
    return true;
}

std::string conectRecords(const core::AtomGroup& atoms) {
    const core::Topology& top = atoms.topology();
    if (!top.bonds) {
        return {};
    }
    std::map<std::int32_t, std::int64_t> serial;
    for (std::size_t i = 0; i < atoms.size(); ++i) {
        serial[atoms[i]] = static_cast<std::int64_t>(i) + 1;
    }
    std::map<std::int32_t, std::vector<std::int64_t>> neighbours;
    for (const core::Bond& bond : *top.bonds) {
        const auto a = serial.find(bond[0]);
        const auto b = serial.find(bond[1]);
        if (a == serial.end() || b == serial.end()) {
            continue;
        }
        if (a->second >= 100000 || b->second >= 100000) {
            return {};  // no CONECT records at all
        }
        neighbours[bond[1]].push_back(a->second);
        neighbours[bond[0]].push_back(b->second);
    }
    std::vector<std::int32_t> order = atoms.indices();
    std::sort(order.begin(), order.end());
    std::string out;
    for (const std::int32_t atom : order) {
        const auto it = neighbours.find(atom);
        if (it == neighbours.end()) {
            continue;
        }
        std::vector<std::int64_t> bonded = it->second;
        std::sort(bonded.begin(), bonded.end());
        std::string record = "CONECT" + padInt(serial[atom], 5);
        for (const std::int64_t s : bonded) {
            record += padInt(s, 5);
        }
        out += line(record);
    }
    return out;
}

// What the writer has written when it fails, and the error (the file is closed with CONECT / END;
// a coordinate error removes it).
struct PartialText : public core::ValueError {
    PartialText(const core::ValueError& error, std::string written, bool removed)
        : core::ValueError(error.what()), text(std::move(written)), fileRemoved(removed) {}
    std::string text;
    bool fileRemoved;
};

std::string frameText(const core::AtomGroup& atoms, const core::Frame& frame,
                      const core::TrajectoryInfo& info) {
    if (atoms.empty()) {
        throw ReaderError("IndexError", "Cannot write empty AtomGroup");
    }
    const core::Topology& top = atoms.topology();
    std::string out;
    if (info.header) {
        out += line("HEADER    " + toUtf8(strip(toU32(*info.header))));
    }
    out += line("TITLE     MDANALYSIS FRAME 0: Created by PDBWriter");
    if (info.compounds) {
        for (const std::string& c : *info.compounds) {
            out += line("COMPND    " + c);
        }
    }
    if (info.remarks) {
        for (const std::string& remark : *info.remarks) {
            for (const std::string& wrapped : core::pythonTextwrap(remark, 66)) {
                out += line("REMARK     " + wrapped);
            }
        }
    }
    if (frame.box) {
        const auto& b = *frame.box;
        out += line("CRYST1" + fixed(b[0], 9, 3) + fixed(b[1], 9, 3) + fixed(b[2], 9, 3) +
                    fixed(b[3], 7, 2) + fixed(b[4], 7, 2) + fixed(b[5], 7, 2) + " P 1           1");
    } else {
        out += line("CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1");
        out += line("REMARK     285 UNITARY VALUES FOR THE UNIT CELL AUTOMATICALLY SET");
        out += line("REMARK     285 BY MDANALYSIS PDBWRITER BECAUSE UNIT CELL INFORMATION");
        out += line("REMARK     285 WAS MISSING.");
        out += line("REMARK     285 PROTEIN DATA BANK CONVENTIONS REQUIRE THAT");
        out += line("REMARK     285 CRYST1 RECORD IS INCLUDED, BUT THE VALUES ON");
        out += line("REMARK     285 THIS RECORD ARE MEANINGLESS.");
    }
    const auto close = [&] { return out + conectRecords(atoms) + line("END"); };
    if (!validCoordinates(frame.positions, atoms)) {
        throw PartialText(
            core::ValueError("PDB files must have coordinate values between -1000.000 and 9999.999 "
                             "Angstroem: file writing was aborted."),
            close(), true);
    }
    std::vector<std::string> charges(atoms.size());
    try {
        for (std::size_t i = 0; i < atoms.size(); ++i) {
            charges[i] = top.formalCharges ? formalCharge((*top.formalCharges)[atoms[i]]) : "";
        }
    } catch (const core::ValueError& error) {
        throw PartialText(error, close(), false);
    }
    const auto str = [&](const std::optional<std::vector<std::uint32_t>>& values, std::int32_t atom,
                         const std::string& fallback) {
        return values ? top.strings[(*values)[atom]] : fallback;
    };
    for (std::size_t i = 0; i < atoms.size(); ++i) {
        const std::int32_t atom = atoms[i];
        const std::int32_t residue = top.atomResidue[atom];
        const std::string record = str(top.recordTypes, atom, "ATOM");
        if (record != "ATOM" && record != "HETATM") {
            throw PartialText(
                core::ValueError("Found " + record +
                                 " for the record type, but only allowed types are ATOM or HETATM"),
                close(), false);
        }
        const std::string& resname = top.resname(residue);
        const std::string& name = top.name(atom);
        const U32 icode =
            top.icodes ? slice(toU32(top.strings[(*top.icodes)[residue]]), 0, 1) : U32(U" ");
        const U32 element = top.elements
                                ? upper(slice(toU32(top.strings[(*top.elements)[atom]]), 0, 2))
                                : U32(U" ");
        const double occupancy = top.occupancies ? (*top.occupancies)[atom] : 1.0;
        const double tempFactor = top.tempFactors ? (*top.tempFactors)[atom] : 0.0;
        U32 text = ljust(toU32(record), 6) +
                   toU32(padInt(lastDigits(static_cast<std::int64_t>(i) + 1, 5), 5)) + U" " +
                   atomNameField(name, resname) +
                   ljust(slice(toU32(str(top.altLocs, atom, " ")), 0, 1), 1) +
                   ljust(slice(toU32(resname), 0, 4), 4) +
                   toU32(chainIdOf(str(top.chainIDs, atom, ""))) +
                   toU32(padInt(lastDigits(top.resids[residue], 4), 4)) + ljust(icode, 1) + U"   ";
        for (std::size_t k = 0; k < 3; ++k) {
            text += toU32(fixed(frame.positions[3 * static_cast<std::size_t>(atom) + k], 8, 3));
        }
        text += toU32(fixed(occupancy, 6, 2) + fixed(tempFactor, 6, 2)) + U"      " +
                ljust(slice(toU32(top.segid(top.segmentOfAtom(atom))), 0, 4), 4) +
                rjust(element, 2) + ljust(toU32(charges[i]), 2);
        out += toUtf8(text) + "\n";
    }
    return close();
}

// f"{value:0{width}d}"
std::string zeroPadded(std::int64_t value, std::size_t width) {
    std::string digits = std::to_string(value < 0 ? -value : value);
    const std::size_t sign = value < 0 ? 1 : 0;
    if (digits.size() + sign < width) {
        digits.insert(0, width - digits.size() - sign, '0');
    }
    return (value < 0 ? "-" : "") + digits;
}

void writeText(const std::filesystem::path& path, const std::string& text) {
    std::ofstream file(path, std::ios::binary | std::ios::trunc);
    if (!file) {
        throw std::filesystem::filesystem_error("cannot write file", path,
                                                std::make_error_code(std::errc::permission_denied));
    }
    file.write(text.data(), static_cast<std::streamsize>(text.size()));
}

}  // namespace

std::string pdbFrameText(const core::AtomGroup& atoms, const core::Frame& frame,
                         const core::TrajectoryInfo& info) {
    try {
        return frameText(atoms, frame, info);
    } catch (const PartialText& error) {
        throw core::ValueError(error.what());
    }
}

std::vector<std::filesystem::path> exportFrames(
    core::System& system, const std::vector<std::int64_t>& frames, const core::AtomGroup& atoms,
    const std::filesystem::path& outDir, const std::string& prefix,
    const std::function<void(double, const std::string&)>& progress, const std::string& stamp) {
    if (frames.empty()) {
        throw core::SelectionError("no frames selected for PDB export");
    }
    std::filesystem::create_directories(outDir);
    const std::int64_t count = system.trajectory->frameCount();
    const std::size_t width = std::max<std::size_t>(5, std::to_string(count).size());
    const core::TrajectoryInfo info = system.trajectory->info();
    std::vector<std::filesystem::path> paths;
    for (std::size_t k = 0; k < frames.size(); ++k) {
        const std::int64_t f = frames[k];
        const core::Frame& frame = system.trajectory->seek(f);
        const std::string name =
            prefix + "_" + zeroPadded(f, width) + (stamp.empty() ? "" : "_" + stamp) + ".pdb";
        const std::filesystem::path path = outDir / pathFromUtf8(name);
        try {
            writeText(path, frameText(atoms, frame, info));
        } catch (const PartialText& error) {
            if (error.fileRemoved) {
                std::error_code ignored;
                std::filesystem::remove(path, ignored);
            } else {
                writeText(path, error.text);
            }
            throw core::ValueError(error.what());
        }
        paths.push_back(path);
        if (progress) {
            progress(static_cast<double>(k + 1) / static_cast<double>(frames.size()),
                     "PDB frame " + std::to_string(f));
        }
    }
    return paths;
}

}  // namespace mdbridge::io
