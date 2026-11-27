#pragma once
// The PDB clean-up window (milestone 6 W3-B, owner decision D1): MDATK v0.0.4.1's PDBToolsWindow
// (mdatk/gui/pdbdialog.py) with 0.1.1's host (gui/host.py) and fixes, on the chem library through
// tools/tool_actions.hpp, so it writes what `mdbridge pdb clean` writes.
//
//   Files                 input PDB (the output suggested as <input>_clean.pdb while it is empty or
//                         still the suggestion), output PDB
//   Names and records     residue name scheme, 5'/3' terminal names, L detection, forced L ranges,
//                         ATOM / HETATM policy, CONECT, 5'-phosphates (ticked when Amber is chosen)
//   Chains                keep / reference PDB (+ its residue numbers) / ids per TER segment / A, B,
//                         C ...; split segments, split before, combine into, alphabetical order
//   Numbering and sequence  renumber (first number, per chain, keep gaps), SEQRES, FASTA file
//   Annotation            REMARK 950 records (engine, pairs, stacking, dot-bracket)
//   Check input           the report of the file (nothing is written)
//   Write cleaned PDB     the output (and the FASTA file), the step "pdb clean" in the project's log
//
// Differences from v0.0.4.1 (docs/DIFFERENCES.md TW-xx): the host interface (TW-01), *Check input*
// reports a file that is not UTF-8 text (TW-02, MDATK 0.1.0), the annotation engines as MDBridge
// has them (TW-05, CH-14).

#include <QString>

#include "gui/tools/tool_window.hpp"

class QCheckBox;
class QComboBox;
class QLineEdit;
class QSpinBox;

namespace mdbridge::app::tools {

struct PdbCleanRequest;

class PdbToolsWindow : public ToolWindow {
    Q_OBJECT

public:
    explicit PdbToolsWindow(Host& host, QWidget* parent = nullptr);

    // The buttons' actions (MDATK's run_check / run_clean).
    void runCheck();
    void runClean();

    // The window's options() and annotation_options() as a clean-up request for `input` / `output`.
    PdbCleanRequest request(const QString& input, const QString& output) const;

private:
    void sync();
    void onInput();
    // the input file, or "" after saying "Choose an existing input PDB file first."
    QString inputFile();

    QString m_autoOutput;
    PathEntry* m_input = nullptr;
    PathEntry* m_output = nullptr;
    QComboBox* m_scheme = nullptr;
    QCheckBox* m_terminals = nullptr;
    QCheckBox* m_detect = nullptr;
    QLineEdit* m_lRanges = nullptr;
    QComboBox* m_records = nullptr;
    QCheckBox* m_keepConect = nullptr;
    QCheckBox* m_strip5p = nullptr;
    QComboBox* m_chainMode = nullptr;
    PathEntry* m_reference = nullptr;
    QCheckBox* m_copyNumbers = nullptr;
    QLineEdit* m_chainIds = nullptr;
    QCheckBox* m_splitSegments = nullptr;
    QLineEdit* m_splitAt = nullptr;
    QLineEdit* m_combine = nullptr;
    QLineEdit* m_combineInto = nullptr;
    QCheckBox* m_sortChains = nullptr;
    QCheckBox* m_renumber = nullptr;
    QSpinBox* m_start = nullptr;
    QCheckBox* m_perChain = nullptr;
    QCheckBox* m_keepGaps = nullptr;
    QCheckBox* m_seqres = nullptr;
    PathEntry* m_fasta = nullptr;
    QCheckBox* m_annotate = nullptr;
    QComboBox* m_engine = nullptr;
    QCheckBox* m_annPairs = nullptr;
    QCheckBox* m_annStacking = nullptr;
    QCheckBox* m_annDbn = nullptr;
};

}  // namespace mdbridge::app::tools
