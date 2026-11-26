#pragma once
// The Input & cleaning page (MDATK v0.0.4.1 gui/app.py _tab_input): project notes, input files,
// frames / time / external data, trajectory cleaning with cpptraj.

#include <string>

#include "gui/pages/form_page.hpp"

class QPushButton;

namespace mdbridge::app::pages {

class InputPage : public FormPage {
    Q_OBJECT

public:
    explicit InputPage(forms::FormContext& context, QWidget* parent = nullptr);

    QPushButton* detectButton() const { return m_detect; }
    QPushButton* runCleaningButton() const { return m_runCleaning; }
    QPushButton* inspectButton() const { return m_inspect; }

    // *Detect*: the first conda environment that has bin/cpptraj goes into cleaning.conda_env (MDATK
    // detect_conda); without one the log says where it looked. Returns the environment chosen ("").
    QString detectConda();

signals:
    // A line for the window's log ("WARNING: ..." for problems), as MDATK's append_log.
    void logMessage(const QString& text);
    // The window's actions (milestone 5 W2-B): *Run cleaning now* and *Inspect system*.
    void runCleaningRequested();
    void inspectSystemRequested();
    // *Detect* found nothing: MDATK's info dialog ("No conda environment with cpptraj found", detail).
    void information(const QString& text, const QString& detail);

private:
    // field() with `label` (check box text, accessible name) and `hint` as tooltip when the spec has none.
    forms::FieldWidget* fieldAs(const std::string& path, const QString& label, const QString& hint = {});
    QPushButton* m_detect = nullptr;
    QPushButton* m_runCleaning = nullptr;
    QPushButton* m_inspect = nullptr;
};

}  // namespace mdbridge::app::pages
