#pragma once
// The structure builder window (milestone 6 W3-B, owner decision D1): MDATK v0.0.4.1's
// BuilderWindow (mdatk/gui/builderdialog.py) with 0.1.1's host (gui/host.py) and fixes, on the chem
// library through tools/tool_actions.hpp, so it writes what `mdbridge build helix / mutate / fit`
// write.
//
//   Ideal helix   sequence, partner (complement / none / explicit), chirality, pairing, triplex
//                 and third strand, chains and first numbers per strand, residue names, terminal
//                 names, 5'-phosphate, hydrogens, output
//   Mutate        input PDB, mutations, compensatory partners (G·U kept), output (<input>_mut.pdb
//                 suggested)
//   Fit to map    model PDB, density map, resolution, global search or local refinement, angular
//                 step, output (<input>_fit.pdb suggested)
//   Check         the current page's inputs (nothing is written; the fit page has none)
//   Write PDB     the current page's build; the fitted pairings, triplexes and the density fit run
//                 on a worker thread (the window stays responsive) with the progress bar, the fit
//                 with *Cancel fit*
//
// The other pairings and triplexes are offered as in MDATK and built by the library's fitted helix
// builder (chem/helix_fit.cpp, milestone 6 W2-A) on the worker thread (tool_actions.hpp).
//
// Differences from v0.0.4.1 (docs/DIFFERENCES.md TW-xx): the host interface (TW-01), *Check* reports
// a file that is not UTF-8 text (TW-02, MDATK 0.1.0), a cancelled fit is logged as cancelled and
// reported as "Cancelled: the density fit stopped, nothing was written." (TW-03, MDATK 0.0.13),
// closing during a build asks and waits (TW-04, JB-01).

#include <memory>

#include <QString>

#include "gui/tools/tool_window.hpp"

class QCheckBox;
class QComboBox;
class QDoubleSpinBox;
class QLineEdit;
class QProgressBar;
class QPushButton;
class QSpinBox;
class QTabWidget;

namespace mdbridge::app::jobs {
class JobRunner;
}

namespace mdbridge::app::tools {

struct HelixRequest;
struct ToolResult;

class BuilderWindow : public ToolWindow {
    Q_OBJECT

public:
    explicit BuilderWindow(Host& host, QWidget* parent = nullptr);
    ~BuilderWindow() override;

    // "helix", "mutate", "fit" (MDATK's PAGES) and showing one.
    QString page() const;
    void setPage(const QString& name);

    // The buttons' actions (MDATK's check_current / run_current / _cancel).
    void checkCurrent();
    void runCurrent();
    void cancel();

    // The *Ideal helix* page's helix_options().
    HelixRequest helixRequest() const;

    // The window's worker (one build at a time).
    jobs::JobRunner& jobs() const { return *m_jobs; }
    bool busy() const;

protected:
    bool mayClose() override;

private:
    QWidget* helixPage();
    QWidget* mutatePage();
    QWidget* fitPage();
    void sync();
    void suggest(const QString& page);
    QString inputFile(PathEntry* entry, const QString& what);
    QString outputFile(PathEntry* entry, const QString& source = QString());
    void startJob(const QString& title, bool cancellable, bool pulse);
    void finishJob();
    void showWritten(const ToolResult& result);

    void runHelix();
    void runMutate();
    void runFit();

    jobs::JobRunner* m_jobs = nullptr;
    QString m_autoMutate;
    QString m_autoFit;
    QString m_cancelledText;  // the report after Cancel of the running job

    QTabWidget* m_tabs = nullptr;
    QProgressBar* m_progress = nullptr;
    QPushButton* m_cancel = nullptr;
    QPushButton* m_check = nullptr;
    QPushButton* m_write = nullptr;

    // helix
    QLineEdit* m_seq = nullptr;
    QComboBox* m_partnerMode = nullptr;
    QLineEdit* m_partner = nullptr;
    QComboBox* m_chirality = nullptr;
    QComboBox* m_pairing = nullptr;
    QComboBox* m_triplex = nullptr;
    QLineEdit* m_third = nullptr;
    QLineEdit* m_chain1 = nullptr;
    QLineEdit* m_chain2 = nullptr;
    QLineEdit* m_chain3 = nullptr;
    QSpinBox* m_start1 = nullptr;
    QSpinBox* m_start2 = nullptr;
    QSpinBox* m_start3 = nullptr;
    QComboBox* m_scheme = nullptr;
    QCheckBox* m_terminals = nullptr;
    QCheckBox* m_phosphate5 = nullptr;
    QCheckBox* m_hydrogens = nullptr;
    PathEntry* m_helixOut = nullptr;
    // mutate
    PathEntry* m_mutInput = nullptr;
    QLineEdit* m_mutations = nullptr;
    QCheckBox* m_compensatory = nullptr;
    QCheckBox* m_wobble = nullptr;
    PathEntry* m_mutOut = nullptr;
    // fit
    PathEntry* m_fitInput = nullptr;
    PathEntry* m_fitMap = nullptr;
    QDoubleSpinBox* m_resolution = nullptr;
    QComboBox* m_search = nullptr;
    QDoubleSpinBox* m_step = nullptr;
    PathEntry* m_fitOut = nullptr;
};

}  // namespace mdbridge::app::tools
