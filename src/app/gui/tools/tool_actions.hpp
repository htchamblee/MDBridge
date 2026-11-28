#pragma once
// The work of the PDB clean-up and structure builder windows (milestone 6 W3-B), without Qt: each
// action calls the chem library as the command line does (src/app/commands/chem_commands.cpp:
// chem::clean, chem::annotate, chem::check, chem::buildHelix, chem::mutate, chem::fitToMap), logs
// the step in the project's step log with MDATK's window titles (mdatk/gui/pdbdialog.py,
// builderdialog.py: no command line in the entry) and writes the same files, so a window writes
// what `mdbridge pdb clean` / `build helix` / `build mutate` / `build fit` write. The windows turn
// their widgets into the requests below and show the texts these return.
//
// Errors propagate as the library throws them (io::PdbError, io::BuildError, core::ValueError,
// std::filesystem::filesystem_error, core::Cancelled); errorText() is what the windows show after
// "ERROR: ".
//
// The other pairings and triplexes (W2-A): buildHelix passes the pairing and triplex to
// chem::buildHelix (its fitted builder, chem/helix_fit.cpp); helixInBackground() sends those builds
// to the window's worker thread, as MDATK does.

#include <cstdint>
#include <exception>
#include <filesystem>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "mdbridge/chem/annotate.hpp"
#include "mdbridge/core/progress.hpp"

namespace mdbridge::app::tools {

// str() of an exception as MDATK's windows show it (app::pythonErrorText, as the command line
// prints errors): Python's OSError text ("[Errno 2] No such file or directory: 'x.pdb'") for file
// errors, the message otherwise.
std::string errorText(const std::exception& error);

// os.path.isfile(path) of a path as typed (relative to the current folder).
bool isFile(const std::string& path);
// os.path.abspath(a) == os.path.abspath(b)
bool samePath(const std::string& a, const std::string& b);
// The output the windows suggest for `input`: "<stem>_<suffix><extension, else .pdb>" (os.path.splitext;
// on Windows '\\' separates too); "" for "".
std::string suggestedOutput(const std::string& input, std::string_view suffix);

// What one write action did.
struct ToolResult {
    std::string output;                // the file written
    std::vector<std::string> notes;
    std::optional<std::string> fasta;  // the FASTA file written, when one was
    std::optional<double> cc;          // a density fit's correlation
};

// ---- PDB clean-up (PDBToolsWindow) ----------------------------------------------------------

// *Check input*: read the file, assign the chirality (detection on, `lRanges` forced L), and the
// report MDATK's run_check shows (the file's TER segments and CONECT records, per chain its
// residues, numbers, names, chirality calls and sequence, other residues, the notes).
std::string pdbCheckReport(const std::string& input, const std::string& lRanges);

// The window's options() and annotation_options(): pdbtools.clean's keyword arguments as typed.
struct PdbCleanRequest {
    std::string input;
    std::string output;
    std::string fasta;              // "" = no FASTA file
    std::string scheme;             // "" = keep the names
    std::string records;            // "" = follow the scheme
    bool stripConect = true;
    bool detectL = true;
    std::string lRanges;
    std::string chainReference;     // "" = none (chain mode other than "reference")
    std::optional<std::vector<std::string>> chainIds;
    bool assignChains = false;
    bool copyNumbers = false;
    bool splitSegments = false;
    std::string splitAt;
    std::vector<std::string> combine;
    std::optional<std::string> combineInto;
    bool sortChains = false;
    bool renumber = false;
    std::int64_t start = 1;
    bool perChain = false;
    bool keepGaps = false;
    bool seqres = false;
    bool terminalNames = true;
    bool strip5p = false;
    std::optional<chem::AnnotateOptions> annotation;  // nullopt: no REMARK 950 records
};

// *Write cleaned PDB*: step "pdb clean", "PDB clean-up of <name>": clean, then (when asked) the
// annotation of the cleaned text (its notes added unless already there), the output and the FASTA
// file written; the entry's files and summary as MDATK's window writes them.
ToolResult cleanPdb(const PdbCleanRequest& request, const std::filesystem::path& project);

// [c.strip() for c in text.split(",") if c.strip()], nullopt when that is empty.
std::optional<std::vector<std::string>> splitIds(const std::string& text);

// ---- structure builder (BuilderWindow) ------------------------------------------------------

// The *Ideal helix* page's helix_options().
struct HelixRequest {
    std::string sequence;           // as typed
    std::string partner = "auto";   // "auto", "none" or the partner sequence
    std::string chirality = "D";
    std::vector<std::string> chains{"A", "B"};       // strand 1, 2 (3 for a triplex)
    std::vector<std::int64_t> start{1, 1};
    std::string scheme;             // "" = the default names
    bool terminals = true;
    bool hydrogens = false;
    bool phosphate5 = false;
    std::string pairing = "wc";
    std::string triplex;            // "" = none
    std::string third = "auto";
    std::string output;
};

// The helix is built on the worker thread: the fitted pairings and triplexes take seconds
// (MDATK's `other`).
bool helixInBackground(const HelixRequest& request);
// *Check*: check_helix's messages, "LEVEL: text" per line.
std::string helixCheckReport(const HelixRequest& request);
// *Write PDB*: step "build helix", "Ideal <kind> <D|L>-RNA helix <sequence>", citations arnott-arna
// (+ pdb-4jrd for the fitted pairings, + pdb-triplexes for the Hoogsteen duplex and the triplex).
// `reporter`: the job's (progress, cancel): the fitted builds poll it, a build cancelled at its end
// writes nothing (core::Cancelled); the Watson-Crick builder does not use it.
ToolResult buildHelix(const HelixRequest& request, const std::filesystem::path& project,
                      const core::Reporter& reporter = core::Reporter());

// *Check* of the *Mutate* page: per mutation "<where> -> <base>: <residue name or NOT FOUND>".
std::string mutateCheckReport(const std::string& input, const std::string& mutations);
// *Write PDB* of the *Mutate* page: step "build mutate", "Mutation of <name>: <mutations>".
ToolResult mutate(const std::string& input, const std::string& mutations, bool compensatory, bool wobble,
                  const std::string& output, const std::filesystem::path& project);

// *Write PDB* of the *Fit to map* page: step "build fit", "Rigid-body fit of <name> into <map>";
// progress and Cancel through `reporter` (a cancelled fit throws core::Cancelled, and the step is
// logged as cancelled: MDATK 0.0.13).
ToolResult fitToMap(const std::string& input, const std::string& map, double resolution, const std::string& search,
                    double step, const std::string& output, const std::filesystem::path& project,
                    const core::Reporter& reporter);

}  // namespace mdbridge::app::tools
