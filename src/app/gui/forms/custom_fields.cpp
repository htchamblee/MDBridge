#include "gui/forms/custom_fields.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <utility>

#include <QCheckBox>
#include <QGridLayout>
#include <QLabel>
#include <QSlider>

#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/render/orientation.hpp"

namespace mdbridge::app::forms {

namespace {

constexpr int kAxes = 3;
constexpr double kLimit = 180.0;  // MDATK's slider range: -180 ... 180 degrees

int sliderUnits(double degrees) {
    if (!std::isfinite(degrees)) {
        return 0;
    }
    const double clamped = std::clamp(degrees, -kLimit, kLimit);
    return static_cast<int>(std::lround(clamped * RotationField::kScale));
}

}  // namespace

// ---- RotationField ---------------------------------------------------------------------------

RotationField::RotationField(FormContext& context, std::string path, FieldSpec spec, QWidget* parent)
    : FieldWidget(context, std::move(path), std::move(spec), parent) {
    auto* grid = new QGridLayout(this);
    grid->setContentsMargins(0, 0, 0, 0);
    grid->setHorizontalSpacing(12);
    auto* sliders = new QWidget(this);  // the editor named after the path
    auto* inner = new QGridLayout(sliders);
    inner->setContentsMargins(0, 0, 0, 0);
    inner->setHorizontalSpacing(12);
    static constexpr std::array<const char*, kAxes> kNames{"X", "Y", "Z"};
    for (int axis = 0; axis < kAxes; ++axis) {
        const QString name = tr("Rotation %1 (°)").arg(QLatin1String(kNames[static_cast<std::size_t>(axis)]));
        auto* axisLabel = new QLabel(name, sliders);
        auto* slider = new QSlider(Qt::Horizontal, sliders);
        slider->setRange(-static_cast<int>(kLimit) * kScale, static_cast<int>(kLimit) * kScale);
        slider->setSingleStep(kScale / 2);  // MDATK's step: 0.5 degrees
        slider->setPageStep(5 * kScale);
        slider->setMinimumWidth(220);
        slider->setObjectName(toQString(this->path()) + QStringLiteral(".%1").arg(axis));
        slider->setProperty("mdbridgeScale", kScale);
        slider->setAccessibleName(name);
        axisLabel->setBuddy(slider);
        guardWheel(slider);
        auto* value = new QLabel(sliders);
        value->setObjectName(slider->objectName() + QStringLiteral(":value"));
        value->setMinimumWidth(value->fontMetrics().horizontalAdvance(QStringLiteral("-180.0")));
        value->setAlignment(Qt::AlignRight | Qt::AlignVCenter);
        inner->addWidget(axisLabel, axis, 0);
        inner->addWidget(slider, axis, 1);
        inner->addWidget(value, axis, 2);
        m_sliders.push_back(slider);
        m_values.push_back(value);
        connect(slider, &QSlider::valueChanged, this, [this, axis](int v) {
            if (!refreshing()) {
                sliderChanged(axis, v);
            }
        });
        connect(slider, &QSlider::sliderPressed, this, [this] {
            if (!m_dragging) {
                m_dragging = true;
                document().beginGroup("Rotate structure");
            }
        });
        connect(slider, &QSlider::sliderReleased, this, [this] {
            if (m_dragging) {
                m_dragging = false;
                document().endGroup();
            }
        });
    }
    inner->setColumnStretch(1, 1);
    m_angles = new QLabel(this);
    m_angles->setObjectName(toQString(this->path()) + QStringLiteral(":angles"));
    m_angles->setTextInteractionFlags(Qt::TextSelectableByMouse | Qt::TextSelectableByKeyboard);
    grid->addWidget(sliders, 0, 0);
    grid->addWidget(m_angles, 1, 0);
    setEditor(sliders);
    setFocusProxy(m_sliders.front());
    refresh();
}

std::vector<double> RotationField::angles() const {
    std::vector<double> out(kAxes, 0.0);
    if (!fieldExists()) {
        return out;
    }
    const core::Value value = document().get(path());
    if (const core::Array* items = value.get<core::Array>()) {
        for (std::size_t k = 0; k < items->size() && k < out.size(); ++k) {
            if (const auto* d = (*items)[k].get<double>()) {
                out[k] = *d;
            } else if (const auto* i = (*items)[k].get<std::int64_t>()) {
                out[k] = static_cast<double>(*i);
            }
        }
    }
    return out;
}

void RotationField::sliderChanged(int axis, int value) {
    std::vector<double> a = angles();
    a[static_cast<std::size_t>(axis)] = static_cast<double>(value) / kScale;
    core::Array items;
    for (double v : a) {
        items.emplace_back(v);
    }
    try {
        document().set(path(), core::Value(std::move(items)), "Rotate structure");
        emit committed(true);
    } catch (const core::FieldError&) {
        emit committed(false);
    }
}

void RotationField::refresh() {
    const RefreshScope scope(*this);
    const std::vector<double> a = angles();
    for (std::size_t k = 0; k < m_sliders.size(); ++k) {
        if (m_sliders[k]->value() != sliderUnits(a[k])) {
            m_sliders[k]->setValue(sliderUnits(a[k]));
        }
        m_values[k]->setText(QString::number(a[k], 'f', 1));
        m_sliders[k]->setEnabled(fieldExists());
    }
    m_angles->setText(toQString(render::formatAngles({a[0], a[1], a[2]})));
    QString problem;
    for (double v : a) {
        if (std::isfinite(v) && std::abs(v) > kLimit) {
            problem = tr("%1 is outside the sliders' range -180 to 180 (kept as it is).")
                          .arg(toQString(core::pythonFloatRepr(v)));
            break;
        }
    }
    setLocalProblem(problem);
}

QString RotationField::displayedText() const {
    QStringList parts;
    for (const QSlider* s : m_sliders) {
        parts.append(toQString(core::pythonFloatRepr(static_cast<double>(s->value()) / kScale)));
    }
    return parts.join(QStringLiteral(", "));
}

// ---- ChoiceChecksField -----------------------------------------------------------------------

ChoiceChecksField::ChoiceChecksField(FormContext& context, std::string path, FieldSpec spec, Normalize normalize,
                                     Supported supported, ListProblem problem, QWidget* parent)
    : FieldWidget(context, std::move(path), std::move(spec), parent), m_normalize(std::move(normalize)),
      m_supported(std::move(supported)), m_problem(std::move(problem)) {
    auto* layout = new QGridLayout(this);
    layout->setContentsMargins(0, 0, 0, 0);
    auto* row = new QWidget(this);
    m_grid = new QGridLayout(row);
    m_grid->setContentsMargins(0, 0, 0, 0);
    m_grid->setHorizontalSpacing(12);
    layout->addWidget(row, 0, 0);
    setEditor(row);
    for (const Choice& c : this->spec().choices) {
        m_shown.append(c.value);
    }
    rebuild();
    refresh();
}

QString ChoiceChecksField::normalized(const QString& value) const { return m_normalize ? m_normalize(value) : value; }

QString ChoiceChecksField::boxText(const QString& value) const {
    for (const Choice& c : spec().choices) {
        if (c.value == value) {
            return c.label.isEmpty() ? value.toUpper() : c.label;
        }
    }
    QString text = value.toUpper();
    if (m_supported && !m_supported(value)) {
        text += tr(" (not supported)");
    }
    return text;
}

QCheckBox* ChoiceChecksField::box(const QString& value) const {
    for (std::size_t k = 0; k < m_boxes.size(); ++k) {
        if (m_shown[static_cast<qsizetype>(k)] == value) {
            return m_boxes[k];
        }
    }
    return nullptr;
}

QStringList ChoiceChecksField::values() const {
    QStringList out;
    if (!fieldExists()) {
        return out;
    }
    const core::Value value = document().get(path());
    if (const core::Array* items = value.get<core::Array>()) {
        for (const core::Value& item : *items) {
            out.append(toQString(item.stringOr()));
        }
    }
    return out;
}

void ChoiceChecksField::rebuild() {
    for (QCheckBox* b : m_boxes) {
        delete b;
    }
    m_boxes.clear();
    constexpr int kPerRow = 4;
    for (qsizetype k = 0; k < m_shown.size(); ++k) {
        const QString value = m_shown[k];
        auto* b = new QCheckBox(boxText(value), editor());
        b->setObjectName(toQString(path()) + QLatin1Char(':') + value);
        b->setProperty("mdbridgeChoice", value);
        b->setAccessibleName(label() + QStringLiteral(": ") + boxText(value));
        m_grid->addWidget(b, static_cast<int>(k / kPerRow), static_cast<int>(k % kPerRow));
        connect(b, &QCheckBox::clicked, this, [this, value](bool on) {
            if (!refreshing()) {
                toggled(value, on);
            }
        });
        m_boxes.push_back(b);
    }
    if (!m_boxes.empty()) {
        editor()->setFocusProxy(m_boxes.front());
    }
}

void ChoiceChecksField::refresh() {
    const RefreshScope scope(*this);
    const QStringList stored = values();
    QStringList shown = m_shown;
    for (const QString& entry : stored) {
        const bool matched = std::any_of(shown.begin(), shown.end(),
                                         [&](const QString& s) { return normalized(s) == normalized(entry); });
        if (!matched) {
            shown.append(entry);
        }
    }
    if (shown != m_shown) {
        m_shown = shown;
        rebuild();
    }
    for (std::size_t k = 0; k < m_boxes.size(); ++k) {
        const QString value = m_shown[static_cast<qsizetype>(k)];
        const bool on = std::any_of(stored.begin(), stored.end(),
                                    [&](const QString& e) { return normalized(e) == normalized(value); });
        m_boxes[k]->setChecked(on);
        m_boxes[k]->setEnabled(fieldExists());
    }
    setLocalProblem(m_problem ? m_problem(stored) : QString());
}

void ChoiceChecksField::toggled(const QString& value, bool on) {
    QStringList list = values();
    const auto matches = [&](const QString& e) { return normalized(e) == normalized(value); };
    if (on) {
        if (std::none_of(list.begin(), list.end(), matches)) {
            list.append(value);
        }
    } else {
        list.erase(std::remove_if(list.begin(), list.end(), matches), list.end());
    }
    core::Array items;
    for (const QString& v : list) {
        items.emplace_back(toStdString(v));
    }
    try {
        document().set(path(), core::Value(std::move(items)));
        emit committed(true);
    } catch (const core::FieldError&) {
        emit committed(false);
        refresh();
    }
}

void ChoiceChecksField::onDocumentChanged(const core::DocumentChange& change) {
    if (change.kind == core::DocumentChange::Kind::Replaced) {  // another project: its own boxes
        QStringList offered;
        for (const Choice& c : spec().choices) {
            offered.append(c.value);
        }
        if (offered != m_shown) {
            m_shown = offered;
            rebuild();
        }
    }
    FieldWidget::onDocumentChanged(change);
}

QString ChoiceChecksField::displayedText() const {
    QStringList on;
    for (std::size_t k = 0; k < m_boxes.size(); ++k) {
        if (m_boxes[k]->isChecked()) {
            on.append(m_shown[static_cast<qsizetype>(k)]);
        }
    }
    return on.join(QStringLiteral(", "));
}

}  // namespace mdbridge::app::forms
