#pragma once
// The "Outputs" settings page (milestone 5 W1-C; MDATK v0.0.4.1 gui/app.py:1336-1419
// _tab_outputs): Locations (output and cache folders, the parameter log), Video (MP4 on / off, file
// name, resolution, DPI, FPS / frame stride, codec / CRF, FFmpeg), Time charts (on / off, series,
// heatmaps, formats, size, DPI, style / line colour / width, title template, axis labels,
// smoothing / target line, prefix) and Frame PDBs (on / off, frames, times, atom subset, prefix,
// Export PDBs now, Re-render video + charts from cache), each a section of MDATK's frames.
//
// Differences from MDATK (docs/DIFFERENCES.md PG-xx): the chart formats are check boxes for every
// format the charts write (RN-11: PNG, SVG, PDF, JPG, TIF, EPS, PS), and entries the list has
// beyond them are kept with a box of their own (AP-04); the chart styles offered are the ones
// MDBridge draws ("mdatk", "default"; owner decision D5 of milestone 4) while any typed name is
// kept (AP-19). The two job buttons are the window's actions (milestone 5 W2-B): they stay
// disabled until the window connects them (setJobActionsAvailable).

#include "gui/pages/form_page.hpp"

class QPushButton;

namespace mdbridge::app::forms {
class ChoiceChecksField;
}

namespace mdbridge::app::pages {

class OutputsPage : public FormPage {
    Q_OBJECT

public:
    explicit OutputsPage(forms::FormContext& context, QWidget* parent = nullptr);

    forms::ChoiceChecksField* formats() const { return m_formats; }
    QPushButton* exportPdbsButton() const { return m_exportPdbs; }
    QPushButton* rerenderButton() const { return m_rerender; }

    // The window enables the job buttons when it handles exportPdbsRequested / rerenderRequested
    // (and disables them while a job runs).
    void setJobActionsAvailable(bool available);

signals:
    void exportPdbsRequested();
    void rerenderRequested();

private:
    void buildLocations();
    void buildVideo();
    void buildCharts();
    void buildPdbs();

    forms::ChoiceChecksField* m_formats = nullptr;
    QPushButton* m_exportPdbs = nullptr;
    QPushButton* m_rerender = nullptr;
};

}  // namespace mdbridge::app::pages
