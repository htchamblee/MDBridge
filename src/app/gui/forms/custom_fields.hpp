#pragma once
// Field widgets the generic factory does not build (milestone 5 W1-C): the three orientation
// sliders of the structure view with MDATK's angle label (visual.rotation) and a row of check
// boxes for a list of choices (output.charts.formats, AP-04). Bound to the document like the other
// field widgets (field_widgets.hpp): they show its value, again after every change, and store what
// the user does as one undo step.

#include <functional>
#include <string>
#include <vector>

#include <QString>
#include <QStringList>

#include "gui/forms/field_widgets.hpp"

class QCheckBox;
class QGridLayout;
class QLabel;
class QSlider;

namespace mdbridge::app::forms {

// visual.rotation: one slider per axis (MDATK's Gtk.Scale, -180 ... 180 in steps of 0.5 degrees,
// the value drawn at its right with one decimal) and the angle label below them
// (render::formatAngles: "X/Y/Z = +15.5°, -30.0°, +120.2°  (... rad)", selectable).
//
// The sliders work in thousandths of a degree, so a value the user sets keeps MDATK's storage
// rounding (3 decimals). A stored angle outside -180 ... 180 is kept (AP-01): its slider shows the
// nearest end, the label and the value text show the stored angle, and moving another axis's slider
// stores the new angle with the others unchanged (MDATK stored all three sliders, clamped). A
// slider drag is one undo step ("Rotate structure", as the canvas's drag); each keyboard step is
// one.
//
// Names: the container is "visual.rotation", the sliders "visual.rotation.0" ... ".2" (with the
// property "mdbridgeScale" = 1000, the slider units per degree: the GUI harness sets them), their
// value texts "visual.rotation.N:value", the label "visual.rotation:angles".
class RotationField : public FieldWidget {
    Q_OBJECT

public:
    static constexpr int kScale = 1000;  // slider units per degree

    RotationField(FormContext& context, std::string path, FieldSpec spec, QWidget* parent = nullptr);

    QSlider* slider(int axis) const { return m_sliders[static_cast<std::size_t>(axis)]; }
    QLabel* valueLabel(int axis) const { return m_values[static_cast<std::size_t>(axis)]; }
    QLabel* angleLabel() const { return m_angles; }
    // The stored angles (three; zeros for what the value lacks).
    std::vector<double> angles() const;

    void refresh() override;
    QString displayedText() const override;

private:
    void sliderChanged(int axis, int value);

    std::vector<QSlider*> m_sliders;
    std::vector<QLabel*> m_values;
    QLabel* m_angles = nullptr;
    bool m_dragging = false;
};

// A list of text values chosen with check boxes (output.charts.formats): one box per offered
// choice (spec.choices, in order), then one per entry of the stored list that no box offers (kept,
// AP-04: "jpeg", "bogus"), labelled with the choice's label or the value in capitals, with "(not
// supported)" after values `supported` rejects. Ticking appends the value to the list as stored,
// unticking removes it (and every entry `normalize` maps to it); the list keeps its order and
// every other entry, and an empty list is stored as such (charts then write PNG, MDATK's
// fallback). The boxes of entries that are no longer stored stay until another project is opened
// (or the document replaced), so an entry unticked by mistake can be ticked again. `problem`
// gives a warning for the list ("" when there is none), shown when the document has no problem
// of the field.
//
// Names: the container is the field path, the boxes "<path>:<value>" (property "mdbridgeChoice"
// = the value; the GUI harness ticks them).
class ChoiceChecksField : public FieldWidget {
    Q_OBJECT

public:
    using Normalize = std::function<QString(const QString&)>;
    using Supported = std::function<bool(const QString&)>;
    using ListProblem = std::function<QString(const QStringList&)>;

    ChoiceChecksField(FormContext& context, std::string path, FieldSpec spec, Normalize normalize = {},
                      Supported supported = {}, ListProblem problem = {}, QWidget* parent = nullptr);

    // The boxes shown, in order.
    const std::vector<QCheckBox*>& boxes() const { return m_boxes; }
    QCheckBox* box(const QString& value) const;
    // The stored list.
    QStringList values() const;

    void refresh() override;
    QString displayedText() const override;

protected:
    void onDocumentChanged(const core::DocumentChange& change) override;

private:
    void toggled(const QString& value, bool on);
    QString normalized(const QString& value) const;
    QString boxText(const QString& value) const;
    void rebuild();

    Normalize m_normalize;
    Supported m_supported;
    ListProblem m_problem;
    QGridLayout* m_grid = nullptr;
    QStringList m_shown;  // the values of the boxes, in order
    std::vector<QCheckBox*> m_boxes;
};

}  // namespace mdbridge::app::forms
