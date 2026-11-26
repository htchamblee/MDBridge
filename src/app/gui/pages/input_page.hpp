#pragma once
// The Input & cleaning page (MDATK v0.0.4.1 gui/app.py _tab_input): project notes, input files,
// frames / time / external data, trajectory cleaning with cpptraj.

#include "gui/pages/settings_page.hpp"

class QPushButton;

namespace mdbridge::app::pages {

class InputPage : public SettingsPage {
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
    // The window's actions (milestone 5 W2-B): *Run cleaning now* and *Inspect system*.
    void runCleaningRequested();
    void inspectSystemRequested();
    // *Detect* found nothing: MDATK's info dialog ("No conda environment with cpptraj found", detail).
    void information(const QString& text, const QString& detail);

private:
    QPushButton* m_detect = nullptr;
    QPushButton* m_runCleaning = nullptr;
    QPushButton* m_inspect = nullptr;
};

}  // namespace mdbridge::app::pages
