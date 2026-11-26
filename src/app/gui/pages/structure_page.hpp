#pragma once
// The "Structure view" settings page (milestone 5 W1-C; MDATK v0.0.4.1 gui/app.py:758-1107
// _tab_visual, _representation_frame, _density_frame): Orientation (the three rotation sliders and
// the angle label, zoom, centre on, follow, reset buttons), Representations (3D) (one row per
// representation: style, groups, colouring, size, opacity, sugars; the other residues, backbone
// atom, hydrogens, renderer), Density map (file, show / colour / opacity / in front, level and
// resolution, fit frame / start / step, automatic fit, the fit buttons, the map-mode toggle, the
// placement status), Residue number tags and Appearance, each a section of MDATK's frames.
//
// Bound to the document: the page's fields show it and store edits as undo steps; the canvas
// (the window) follows the document. Global search and Refine are shown disabled until milestone 6
// brings the density fit (owner decision D2: "comes with milestone 6"); Centre on model needs the
// loaded reference frame and the map, which the window has: it is enabled once the window says so
// (setCentreOnModelAvailable) and asks it with centreOnModelRequested(). Choosing another map file
// forgets the placement of the previous map in the same undo step (MDATK _on_density_map_changed).

#include <QString>

#include "gui/pages/form_page.hpp"

class QLabel;
class QPushButton;

namespace mdbridge::app::forms {
class ListEditor;
class RotationField;
}  // namespace mdbridge::app::forms

namespace mdbridge::app::pages {

class StructurePage : public FormPage {
    Q_OBJECT

public:
    explicit StructurePage(forms::FormContext& context, QWidget* parent = nullptr);

    forms::RotationField* rotation() const { return m_rotation; }
    forms::ListEditor* representations() const { return m_reps; }
    QLabel* densityStatus() const { return m_densityStatus; }
    QPushButton* globalSearchButton() const { return m_globalSearch; }
    QPushButton* refineButton() const { return m_refine; }
    QPushButton* centreOnModelButton() const { return m_centre; }
    QPushButton* clearButton() const { return m_clear; }
    QPushButton* mapModeButton() const { return m_mapMode; }
    bool mapMode() const;

    // MDATK's density status line for the document's settings: "No map selected.", the fit info
    // (or "Placed (not fitted).") when the map has a placement, else "Not placed yet: run a global
    // search, or centre it and move it by hand."
    QString densityStatusText() const;

    // The window enables Centre on model once it has the reference frame shown on the canvas.
    void setCentreOnModelAvailable(bool available);

signals:
    // "Apply": redraw the preview now.
    void previewRequested();
    // "Centre on model" was pressed (the window places the map: core density centre transform).
    void centreOnModelRequested();
    // "Move the map with the mouse" toggled: the canvas's map mode.
    void mapModeChanged(bool on);

private:
    void buildOrientation();
    void buildRepresentations();
    void buildDensity();
    void buildTags();
    void buildAppearance();
    void documentChanged(const core::DocumentChange& change);

    forms::RotationField* m_rotation = nullptr;
    forms::ListEditor* m_reps = nullptr;
    QLabel* m_densityStatus = nullptr;
    QPushButton* m_globalSearch = nullptr;
    QPushButton* m_refine = nullptr;
    QPushButton* m_centre = nullptr;
    QPushButton* m_clear = nullptr;
    QPushButton* m_mapMode = nullptr;
    bool m_centreAvailable = false;
    bool m_fixingMap = false;
};

}  // namespace mdbridge::app::pages
