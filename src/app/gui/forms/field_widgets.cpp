#include "gui/forms/field_widgets.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <exception>
#include <optional>
#include <string>
#include <utility>

#include <QAction>
#include <QCheckBox>
#include <QApplication>
#include <QColorDialog>
#include <QDir>
#include <QEvent>
#include <QFileDialog>
#include <QFileInfo>
#include <QHBoxLayout>
#include <QKeyEvent>
#include <QLineEdit>
#include <QPaintEvent>
#include <QPainter>
#include <QPen>
#include <QPalette>
#include <QPixmap>
#include <QSizePolicy>
#include <QStyle>
#include <QStyleOption>
#include <QStylePainter>
#include <QToolButton>

#include "gui/display_path.hpp"
#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/render/color.hpp"
#include "mdbridge/render/text_qt.hpp"

namespace mdbridge::app::forms {

namespace {

const char* const kGuardProperty = "mdbridgeWheelGuard";
const char* const kFlagMessage = "mdbridgeFlagMessage";
const char* const kFlagLevel = "mdbridgeFlagLevel";
const char* const kFlagAction = "mdbridgeFlagAction";

// ---- wheel ----

class WheelGuard : public QObject {
public:
    using QObject::QObject;

protected:
    bool eventFilter(QObject* watched, QEvent* event) override {
        if (event->type() == QEvent::Wheel) {
            auto* widget = qobject_cast<QWidget*>(watched);
            if (widget && !widget->hasFocus()) {
                event->ignore();  // the widget keeps its value; the event goes on to the page
                return true;
            }
        }
        return false;
    }
};

WheelGuard* wheelGuard() {
    static WheelGuard* guard = new WheelGuard(qApp);
    return guard;
}

// ---- flags (MDATK 0.1.1 widgets.FLAG_COLORS) ----

QColor flagColor(const QString& level) {
    return level == QLatin1String("warning") ? QColor(0xe5, 0xa5, 0x0a) : QColor(0xe0, 0x1b, 0x24);
}

QColor mix(const QColor& a, const QColor& b, double t) {
    return QColor::fromRgbF(static_cast<float>(a.redF() * (1 - t) + b.redF() * t),
                            static_cast<float>(a.greenF() * (1 - t) + b.greenF() * t),
                            static_cast<float>(a.blueF() * (1 - t) + b.blueF() * t));
}

QLineEdit* lineEditOf(QWidget* widget) {
    if (auto* edit = qobject_cast<QLineEdit*>(widget)) {
        return edit;
    }
    if (auto* combo = qobject_cast<QComboBox*>(widget)) {
        return combo->lineEdit();
    }
    if (qobject_cast<QAbstractSpinBox*>(widget)) {
        return widget->findChild<QLineEdit*>();
    }
    return nullptr;
}

QString withHint(const QString& message, const QString& hint) {
    if (message.isEmpty()) {
        return hint;
    }
    return hint.isEmpty() ? message : message + QStringLiteral("\n\n") + hint;
}

QString numberHintText(double v) { return QString::number(v, 'g', 12); }

// ---- dialogs ----

QStringList nativeFileDialog(QWidget* parent, const FileDialogRequest& r) {
    const QString filters = r.filters.join(QStringLiteral(";;"));
    switch (r.mode) {
    case PathMode::OpenFile: {
        const QString file = QFileDialog::getOpenFileName(parent, r.title, r.start, filters);
        return file.isEmpty() ? QStringList() : QStringList{file};
    }
    case PathMode::OpenFiles:
        return QFileDialog::getOpenFileNames(parent, r.title, r.start, filters);
    case PathMode::SaveFile: {
        const QString file = QFileDialog::getSaveFileName(
            parent, r.title, r.start, filters, nullptr,
            r.confirmOverwrite ? QFileDialog::Options() : QFileDialog::Options(QFileDialog::DontConfirmOverwrite));
        return file.isEmpty() ? QStringList() : QStringList{file};
    }
    case PathMode::Folder: {
        const QString dir = QFileDialog::getExistingDirectory(parent, r.title, r.start);
        return dir.isEmpty() ? QStringList() : QStringList{dir};
    }
    }
    return {};
}

QColor nativeColorDialog(QWidget* parent, const QColor& initial, const QString& title) {
    return QColorDialog::getColor(initial, parent, title);
}

FileDialogFunction& fileDialogFunction() {
    static FileDialogFunction f = nativeFileDialog;
    return f;
}

ColorDialogFunction& colorDialogFunction() {
    static ColorDialogFunction f = nativeColorDialog;
    return f;
}

QColor colourOfText(const QString& text) {
    try {
        const render::Rgba c = render::toRgba(toStdString(text));
        return QColor::fromRgbF(static_cast<float>(c.r), static_cast<float>(c.g), static_cast<float>(c.b),
                                static_cast<float>(c.a));
    } catch (const std::exception&) {
        return {};
    }
}

QIcon swatchIconOf(const QColor& colour) {
    QPixmap pixmap(26, 16);
    pixmap.fill(Qt::transparent);
    QPainter p(&pixmap);
    const QRect inner(2, 2, 22, 12);
    p.fillRect(pixmap.rect(), QColor(255, 255, 255));  // a light frame keeps the swatch distinct
    if (colour.isValid()) {
        for (int y = 0; y < 2; ++y) {  // a checkerboard shows transparency
            for (int x = 0; x < 4; ++x) {
                p.fillRect(QRect(inner.x() + x * 6, inner.y() + y * 6, 6, 6),
                           (x + y) % 2 ? QColor(204, 204, 204) : QColor(255, 255, 255));
            }
        }
        p.fillRect(inner, colour);
    } else {  // no colour: crossed out
        p.setPen(QPen(QColor(0xe0, 0x1b, 0x24), 2));
        p.drawLine(inner.topLeft(), inner.bottomRight());
        p.drawLine(inner.bottomLeft(), inner.topRight());
    }
    p.setPen(QColor(128, 128, 128));
    p.drawRect(inner.adjusted(0, 0, -1, -1));
    return QIcon(pixmap);
}

// The folder (or file) a dialog starts at: the field's path made absolute against the project's
// folder.
QString startOf(const core::ProjectDocument& document, const QString& text) {
    const QDir folder(toQString(document.folder()));
    QString t = text.trimmed();
    if (t.startsWith(QLatin1Char('~'))) {
        t = QDir::homePath() + t.mid(1);
    }
    return t.isEmpty() ? folder.absolutePath() : QDir::cleanPath(folder.absoluteFilePath(t));
}

// Text typed into a line edit and not stored yet. QLineEdit's own modified flag is not used: a
// spin box clears it when Return is pressed, before editingFinished.
const char* const kDirty = "mdbridgeDirty";

bool isDirty(const QLineEdit* edit) { return edit->property(kDirty).toBool(); }

void setDirty(QLineEdit* edit, bool dirty) { edit->setProperty(kDirty, dirty); }

void trackEdits(QLineEdit* edit) {
    QObject::connect(edit, &QLineEdit::textEdited, edit, [edit] { setDirty(edit, true); });
}

// Esc in a line edit shows the stored value again (the edit is dropped).
class EscapeReverts : public QObject {
public:
    EscapeReverts(QObject* parent, std::function<void()> revert) : QObject(parent), m_revert(std::move(revert)) {}

protected:
    bool eventFilter(QObject* watched, QEvent* event) override {
        if (event->type() == QEvent::KeyPress && static_cast<QKeyEvent*>(event)->key() == Qt::Key_Escape) {
            QLineEdit* edit = lineEditOf(qobject_cast<QWidget*>(watched));
            if (edit && isDirty(edit)) {
                m_revert();
                return true;
            }
        }
        return false;
    }

private:
    std::function<void()> m_revert;
};

}  // namespace

QColor matplotlibColor(const QString& text) { return colourOfText(text); }

QIcon colorSwatch(const QColor& colour) { return swatchIconOf(colour); }

void guardWheel(QWidget* widget) {
    if (!widget || hasWheelGuard(widget)) {
        return;
    }
    widget->installEventFilter(wheelGuard());
    widget->setProperty(kGuardProperty, true);
    if (widget->focusPolicy() == Qt::WheelFocus) {
        widget->setFocusPolicy(Qt::StrongFocus);
    }
}

bool hasWheelGuard(const QWidget* widget) { return widget && widget->property(kGuardProperty).toBool(); }

void setFlag(QWidget* widget, const QString& message, const QString& level, const QString& hint) {
    if (!widget) {
        return;
    }
    const QString lv = message.isEmpty() ? QString() : (level == QLatin1String("warning") ? level : QStringLiteral("error"));
    const bool same = widget->property(kFlagMessage).toString() == message && widget->property(kFlagLevel).toString() == lv;
    const QString shown = gui::displayPathsIn(message);  // paths in it shown natively (WN-12)
    widget->setToolTip(withHint(shown, hint));
    widget->setAccessibleDescription(withHint(shown, hint));
    if (same) {
        return;
    }
    widget->setProperty(kFlagMessage, message);
    widget->setProperty(kFlagLevel, lv);
    widget->setPalette(QPalette());  // the inherited palette (a theme change keeps working)
    QLineEdit* edit = lineEditOf(widget);
    if (edit) {
        if (auto* old = qobject_cast<QAction*>(edit->property(kFlagAction).value<QObject*>())) {
            edit->removeAction(old);
            delete old;
            edit->setProperty(kFlagAction, QVariant());
        }
    }
    if (message.isEmpty()) {
        return;
    }
    QPalette palette = widget->palette();
    for (const QPalette::ColorRole role : {QPalette::Base, QPalette::Button}) {
        palette.setColor(role, mix(palette.color(role), flagColor(lv), 0.25));
    }
    widget->setPalette(palette);
    if (edit) {
        const QIcon icon = widget->style()->standardIcon(lv == QLatin1String("warning") ? QStyle::SP_MessageBoxWarning
                                                                                        : QStyle::SP_MessageBoxCritical);
        QAction* action = edit->addAction(icon, QLineEdit::TrailingPosition);
        action->setToolTip(shown);
        action->setText(shown);
        edit->setProperty(kFlagAction, QVariant::fromValue<QObject*>(action));
    }
}

std::pair<QString, QString> flagOf(const QWidget* widget) {
    if (!widget) {
        return {};
    }
    return {widget->property(kFlagMessage).toString(), widget->property(kFlagLevel).toString()};
}

FileDialogFunction setFileDialogFunction(FileDialogFunction function) {
    return std::exchange(fileDialogFunction(), function ? std::move(function) : FileDialogFunction(nativeFileDialog));
}

ColorDialogFunction setColorDialogFunction(ColorDialogFunction function) {
    return std::exchange(colorDialogFunction(), function ? std::move(function) : ColorDialogFunction(nativeColorDialog));
}

QStringList runFileDialog(QWidget* parent, const FileDialogRequest& request) { return fileDialogFunction()(parent, request); }

QString pickedPathText(const core::ProjectDocument& document, const QString& picked) {
    const QString native = QDir::toNativeSeparators(picked);
    if (!document.hasPath() || picked.isEmpty()) {
        return native;
    }
    const QString absolute = QDir::cleanPath(QFileInfo(picked).absoluteFilePath());
    const std::optional<std::string> relative = core::relativePath(toStdString(absolute), document.folder());
    if (!relative || *relative == ".." || relative->rfind("../", 0) == 0) {
        return native;
    }
    const QString text = toQString(*relative);
    return text.startsWith(QLatin1Char('~')) ? QStringLiteral("./") + text : text;
}

QColor runColorDialog(QWidget* parent, const QColor& initial, const QString& title) {
    return colorDialogFunction()(parent, initial, title);
}

// ---- FieldWidget -----------------------------------------------------------------------------

FieldWidget::FieldWidget(FormContext& context, std::string path, FieldSpec spec, QWidget* parent)
    : QWidget(parent), m_context(context), m_path(std::move(path)), m_spec(std::move(spec)) {
    setObjectName(toQString(m_path) + QStringLiteral(":field"));
    m_context.addField(this);
    connect(&m_context, &FormContext::documentChanged, this,
            [this](const core::DocumentChange& change) { onDocumentChanged(change); });
    connect(&m_context, &FormContext::problemsChanged, this, &FieldWidget::updateFlag);
}

FieldWidget::~FieldWidget() {
    m_context.removeField(this);
    // QWidget's destructor takes the focus from a focused editor, whose focus-out (editingFinished)
    // would call into this half-destroyed widget: its editors' signals no longer reach it
    for (QObject* child : findChildren<QObject*>()) {
        QObject::disconnect(child, nullptr, this, nullptr);
    }
}

void FieldWidget::onDocumentChanged(const core::DocumentChange& change) {
    if (change.affects(m_path)) {
        refresh();
        updateFlag();
    }
}

QString FieldWidget::label() const {
    if (!m_spec.label.isEmpty()) {
        return m_spec.label;
    }
    const QString where = toQString(core::describeField(m_path, &document().project()));
    const qsizetype cut = where.lastIndexOf(QStringLiteral(" › "));
    return cut < 0 ? where : where.mid(cut + 3);
}

void FieldWidget::setEditor(QWidget* editor, QWidget* flagTarget) {
    m_editor = editor;
    m_flagTarget = flagTarget ? flagTarget : editor;
    editor->setObjectName(toQString(m_path));
    editor->setAccessibleName(label());
    if (QLineEdit* inner = lineEditOf(editor); inner && inner != editor) {
        inner->setAccessibleName(label());  // a spin box's or editable combo's text: what has the focus (D10)
    }
    editor->setToolTip(m_spec.hint);
    editor->setAccessibleDescription(m_spec.hint);
    setFocusProxy(editor);
    if (QLineEdit* edit = lineEditOf(editor)) {  // check as you type (W1-D)
        connect(edit, &QLineEdit::textEdited, &m_context, &FormContext::scheduleCheck);
    }
}

std::pair<QString, QString> FieldWidget::flag() const { return flagOf(m_flagTarget); }

bool FieldWidget::fieldExists() const { return core::fieldInfo(document().project(), m_path).has_value(); }

QString FieldWidget::documentText() const {
    if (!fieldExists()) {
        return {};
    }
    try {
        return toQString(document().text(m_path));
    } catch (const core::FieldError&) {
        return {};
    }
}

bool FieldWidget::commitText(const QString& text, const std::string& path) {
    const std::string field = path.empty() ? m_path : path;
    try {
        document().setText(field, toStdString(text));
    } catch (const core::FieldError& error) {
        const core::Problem problem = error.problem(&document().project());
        m_context.setFormError(field, &problem);
        updateFlag();
        emit committed(false);
        return false;
    } catch (const std::exception& error) {
        const core::Problem problem = core::locatedProblem(field, error.what(), &document().project());
        m_context.setFormError(field, &problem);
        updateFlag();
        emit committed(false);
        return false;
    }
    m_context.setFormError(field, nullptr);
    updateFlag();
    emit committed(true);
    return true;
}

void FieldWidget::setLocalProblem(const QString& message, const QString& level) {
    if (m_localMessage == message && m_localLevel == level) {
        return;
    }
    m_localMessage = message;
    m_localLevel = level;
    if (message.isEmpty()) {
        m_context.setLocalProblem(m_path, nullptr);
    } else {
        const core::Problem problem = core::locatedProblem(m_path, toStdString(message), &document().project(),
                                                           level == QLatin1String("warning") ? "warning" : "error");
        m_context.setLocalProblem(m_path, &problem);
    }
    updateFlag();
}

void FieldWidget::updateFlag() {
    if (!m_flagTarget) {
        return;
    }
    const std::vector<core::Problem> problems = m_context.problemsAt(m_path);
    QString message;
    QString level;
    if (!problems.empty()) {
        level = toQString(problems.front().level);
        QStringList lines;
        for (const core::Problem& p : problems) {
            const QString line = toQString(p.message);
            if (toQString(p.level) == level && !lines.contains(line)) {
                lines.append(line);
            }
        }
        message = lines.join(QLatin1Char('\n'));
    } else {
        message = m_localMessage;
        level = m_localLevel;
    }
    QString hint = m_spec.hint;
    if (m_spec.minimum || m_spec.maximum) {
        const QString range = QStringLiteral("Usual range: %1 to %2 (other values are kept)")
                                  .arg(m_spec.minimum ? numberHintText(*m_spec.minimum) : QStringLiteral("-"),
                                       m_spec.maximum ? numberHintText(*m_spec.maximum) : QStringLiteral("-"));
        hint = hint.isEmpty() ? range : hint + QLatin1Char('\n') + range;
    }
    setFlag(m_flagTarget, message, level, hint);
}

// ---- NumberField -----------------------------------------------------------------------------

NumberBox::NumberBox(QWidget* parent) : QAbstractSpinBox(parent) {
    setButtonSymbols(QAbstractSpinBox::UpDownArrows);
    setAccelerated(true);
    guardWheel(this);
}

void NumberBox::stepBy(int steps) {
    if (onStep) {
        onStep(steps);
    }
}

QAbstractSpinBox::StepEnabled NumberBox::stepEnabled() const {
    return onStepEnabled ? onStepEnabled() : StepEnabled(StepUpEnabled | StepDownEnabled);
}

NumberField::NumberField(FormContext& context, std::string path, FieldSpec spec, QWidget* parent)
    : FieldWidget(context, std::move(path), std::move(spec), parent) {
    auto* layout = new QHBoxLayout(this);
    layout->setContentsMargins(0, 0, 0, 0);
    m_box = new NumberBox(this);
    m_box->edit()->setPlaceholderText(this->spec().placeholder);
    trackEdits(m_box->edit());
    layout->addWidget(m_box);
    setEditor(m_box);
    m_box->onStep = [this](int steps) { step(steps); };
    m_box->onStepEnabled = [this]() -> QAbstractSpinBox::StepEnabled {
        const std::optional<double> v = currentNumber();
        QAbstractSpinBox::StepEnabled out = QAbstractSpinBox::StepNone;
        if (!fieldExists() || m_box->isReadOnly()) {
            return out;
        }
        if (!v || !this->spec().maximum || *v < *this->spec().maximum) {
            out |= QAbstractSpinBox::StepUpEnabled;
        }
        if (!v || !this->spec().minimum || *v > *this->spec().minimum) {
            out |= QAbstractSpinBox::StepDownEnabled;
        }
        return out;
    };
    connect(m_box, &QAbstractSpinBox::editingFinished, this, &NumberField::commit);
    connect(m_box->edit(), &QLineEdit::textEdited, this, &NumberField::checkRange);
    m_box->installEventFilter(new EscapeReverts(this, [this] {
        this->context().setFormError(this->path(), nullptr);
        refresh();
        updateFlag();
        this->context().scheduleCheck();  // the typed text is gone: check the value shown again
    }));
    refresh();
}

bool NumberField::isInteger() const {
    const auto info = core::fieldInfo(document().project(), path());
    return info && (info->kind == core::FieldKind::Integer || info->kind == core::FieldKind::OptionalInteger);
}

std::optional<double> NumberField::currentNumber() const {
    const std::string text = toStdString(m_box->edit()->text());
    if (const auto i = core::tryPythonInt(text)) {
        return static_cast<double>(*i);
    }
    const auto d = core::tryPythonFloat(text);
    return d && std::isfinite(*d) ? d : std::nullopt;
}

void NumberField::step(int steps) {
    if (!fieldExists()) {
        return;
    }
    const FieldSpec& s = spec();
    const std::optional<double> current = currentNumber();
    const double start = current ? *current : (s.minimum ? *s.minimum : 0.0);
    double next = start + steps * s.step;
    // the arrows stop at the range hint (a value already outside it is never moved into it by force)
    if (steps > 0 && s.maximum && start < *s.maximum) {
        next = std::min(next, *s.maximum);
    }
    if (steps < 0 && s.minimum && start > *s.minimum) {
        next = std::max(next, *s.minimum);
    }
    QString text;
    if (isInteger()) {
        const double rounded = std::round(next);
        if (!(std::abs(rounded) < 9.0e18)) {
            return;
        }
        text = QString::number(static_cast<std::int64_t>(rounded));
    } else {
        if (s.decimals >= 0) {
            const double scale = std::pow(10.0, s.decimals);
            next = std::round(next * scale) / scale;
        }
        text = toQString(core::pythonFloatRepr(next));
    }
    if (commitText(text)) {
        setDirty(m_box->edit(), false);
    }
}

void NumberField::refresh() {
    const RefreshScope scope(*this);
    const QString text = documentText();
    if (m_box->edit()->text() != text) {
        m_box->edit()->setText(text);
    }
    setDirty(m_box->edit(), false);
    m_box->setReadOnly(!fieldExists());
    checkRange(text);
}

void NumberField::checkRange(const QString& text) {
    QString problem;
    const std::optional<double> v = currentNumber();
    const FieldSpec& s = spec();
    if (v && ((s.minimum && *v < *s.minimum) || (s.maximum && *v > *s.maximum))) {
        problem = QStringLiteral("%1 is outside the usual range %2 to %3 (kept as it is).")
                      .arg(text.trimmed(), s.minimum ? numberHintText(*s.minimum) : QStringLiteral("-"),
                           s.maximum ? numberHintText(*s.maximum) : QStringLiteral("-"));
    }
    setLocalProblem(problem);
}

void NumberField::commit() {
    if (!isDirty(m_box->edit())) {
        return;
    }
    if (commitText(m_box->edit()->text())) {
        setDirty(m_box->edit(), false);
    }
}

QString NumberField::displayedText() const { return m_box->edit()->text(); }

std::optional<QString> NumberField::pendingText() const {
    return isDirty(m_box->edit()) ? std::optional<QString>(m_box->edit()->text()) : std::nullopt;
}

// ---- LineField -------------------------------------------------------------------------------


LineField::LineField(FormContext& context, std::string path, FieldSpec spec, QWidget* parent)
    : FieldWidget(context, std::move(path), std::move(spec), parent) {
    m_row = new QHBoxLayout(this);
    m_row->setContentsMargins(0, 0, 0, 0);
    m_row->setSpacing(4);
    m_edit = new QLineEdit(this);
    m_edit->setPlaceholderText(this->spec().placeholder);
    trackEdits(m_edit);
    m_row->addWidget(m_edit, 1);
    setEditor(m_edit);
    connect(m_edit, &QLineEdit::editingFinished, this, &LineField::commit);
    connect(m_edit, &QLineEdit::textEdited, this, [this](const QString& text) { textShown(text); });
    m_edit->installEventFilter(new EscapeReverts(this, [this] {
        this->context().setFormError(this->path(), nullptr);
        refresh();
        updateFlag();
        this->context().scheduleCheck();  // the typed text is gone: check the value shown again
    }));
    refresh();
}

void LineField::refresh() {
    const RefreshScope scope(*this);
    const QString text = documentText();
    if (m_edit->text() != text) {
        m_edit->setText(text);
        m_edit->setCursorPosition(0);
    }
    setDirty(m_edit, false);
    m_edit->setReadOnly(!fieldExists());
    textShown(text);
}

void LineField::textShown(const QString&) {}

void LineField::commit() {
    if (!isDirty(m_edit)) {
        return;
    }
    if (commitText(m_edit->text())) {
        setDirty(m_edit, false);
    }
}

QString LineField::displayedText() const { return m_edit->text(); }

std::optional<QString> LineField::pendingText() const {
    return isDirty(m_edit) ? std::optional<QString>(m_edit->text()) : std::nullopt;
}

// ---- PathField -------------------------------------------------------------------------------

PathField::PathField(FormContext& context, std::string path, FieldSpec spec, QWidget* parent)
    : LineField(context, std::move(path), std::move(spec), parent) {
    m_browse = new QToolButton(this);
    m_browse->setObjectName(toQString(this->path()) + QStringLiteral(":browse"));
    m_browse->setIcon(style()->standardIcon(this->spec().pathMode == PathMode::Folder ? QStyle::SP_DirOpenIcon
                                                                                      : QStyle::SP_DialogOpenButton));
    const QString title = this->spec().dialogTitle.isEmpty() ? label() : this->spec().dialogTitle;
    m_browse->setText(QStringLiteral("Browse…"));
    m_browse->setToolTip(QStringLiteral("Browse: ") + title);
    m_browse->setAccessibleName(QStringLiteral("Browse: ") + title);
    m_browse->setFocusPolicy(Qt::StrongFocus);
    m_row->addWidget(m_browse);
    connect(m_browse, &QToolButton::clicked, this, &PathField::browse);
}

void PathField::browse() {
    if (!fieldExists()) {
        return;
    }
    commit();
    FileDialogRequest request;
    request.mode = spec().pathMode == PathMode::OpenFiles ? PathMode::OpenFile : spec().pathMode;
    request.title = spec().dialogTitle.isEmpty() ? label() : spec().dialogTitle;
    request.start = startOf(document(), m_edit->text());
    request.filters = nameFilters(spec());
    const QStringList picked = runFileDialog(window(), request);
    if (!picked.isEmpty() && commitText(pickedPathText(document(), picked.front()))) {
        setDirty(m_edit, false);
    }
}

// ---- ColorField ------------------------------------------------------------------------------

ColorField::ColorField(FormContext& context, std::string path, FieldSpec spec, QWidget* parent)
    : LineField(context, std::move(path), std::move(spec), parent) {
    m_pick = new QToolButton(this);
    m_pick->setObjectName(toQString(this->path()) + QStringLiteral(":pick"));
    m_pick->setIconSize(QSize(26, 16));
    m_pick->setToolTip(QStringLiteral("Pick colour…"));
    m_pick->setAccessibleName(QStringLiteral("Pick colour: ") + label());
    m_pick->setFocusPolicy(Qt::StrongFocus);
    m_row->addWidget(m_pick);
    connect(m_pick, &QToolButton::clicked, this, &ColorField::pick);
    refresh();
}

void ColorField::textShown(const QString& text) {
    m_swatch = matplotlibColor(text);
    if (m_pick) {
        m_pick->setIcon(colorSwatch(m_swatch));
    }
    QString problem;
    if (!(spec().allowEmpty && text.trimmed().isEmpty())) {
        if (const auto message = render::colorProblem(toStdString(text))) {
            problem = toQString(*message);
        }
    }
    setLocalProblem(problem, QStringLiteral("error"));
}

void ColorField::pick() {
    if (!fieldExists()) {
        return;
    }
    commit();
    const QColor initial = m_swatch.isValid() ? m_swatch : QColor(Qt::black);
    const QColor chosen = runColorDialog(window(), initial, label());
    if (!chosen.isValid()) {
        return;  // cancelled
    }
    // the picker writes #rrggbb (AP-03); a colour typed as text stays as typed until then
    const render::Rgba rgba{chosen.redF(), chosen.greenF(), chosen.blueF(), 1.0};
    if (commitText(toQString(render::toHex(rgba)))) {
        setDirty(m_edit, false);
    }
}

// ---- MultiLineField --------------------------------------------------------------------------

MultiLineField::MultiLineField(FormContext& context, std::string path, FieldSpec spec, QWidget* parent)
    : FieldWidget(context, std::move(path), std::move(spec), parent) {
    auto* layout = new QHBoxLayout(this);
    layout->setContentsMargins(0, 0, 0, 0);
    m_edit = new QPlainTextEdit(this);
    m_edit->setPlaceholderText(this->spec().placeholder);
    m_edit->setTabChangesFocus(true);  // keyboard users leave the field with Tab
    m_edit->setMinimumHeight(54);
    if (this->spec().monospace) {
        m_edit->setFont(render::interfaceFixedFont(QApplication::font().pointSizeF()));
    } else {
        m_edit->setLineWrapMode(QPlainTextEdit::WidgetWidth);
    }
    layout->addWidget(m_edit);
    setEditor(m_edit);
    connect(m_edit, &QPlainTextEdit::textChanged, this, [this] {
        if (!refreshing()) {
            m_dirty = true;
            this->context().scheduleCheck();
        }
    });
    m_edit->installEventFilter(this);
    refresh();
}

bool MultiLineField::eventFilter(QObject* watched, QEvent* event) {
    if (watched == m_edit && event->type() == QEvent::FocusOut) {
        commit();
    }
    return FieldWidget::eventFilter(watched, event);
}

void MultiLineField::refresh() {
    const RefreshScope scope(*this);
    const QString text = documentText();
    if (m_edit->toPlainText() != text) {
        m_edit->setPlainText(text);
    }
    m_dirty = false;
    m_edit->setReadOnly(!fieldExists());
}

void MultiLineField::commit() {
    if (!m_dirty) {
        return;
    }
    if (commitText(m_edit->toPlainText())) {
        m_dirty = false;
    }
}

QString MultiLineField::displayedText() const { return m_edit->toPlainText(); }

std::optional<QString> MultiLineField::pendingText() const {
    return m_dirty ? std::optional<QString>(m_edit->toPlainText()) : std::nullopt;
}

// ---- ComboField ------------------------------------------------------------------------------

namespace {

QString choiceLabel(const Choice& c) { return c.label.isEmpty() ? c.value : c.label; }

QString unknownLabel(const QString& value) {
    return (value.isEmpty() ? QStringLiteral("\"\"") : value) + QStringLiteral(" (unknown)");
}

}  // namespace

ComboField::ComboField(FormContext& context, std::string path, FieldSpec spec, QWidget* parent)
    : FieldWidget(context, std::move(path), std::move(spec), parent) {
    auto* layout = new QHBoxLayout(this);
    layout->setContentsMargins(0, 0, 0, 0);
    m_combo = new QComboBox(this);
    m_combo->setSizeAdjustPolicy(QComboBox::AdjustToMinimumContentsLengthWithIcon);
    m_combo->setMinimumContentsLength(12);  // long labels are elided, not widening the side bar
    guardWheel(m_combo);
    layout->addWidget(m_combo);
    setEditor(m_combo);
    connect(m_combo, &QComboBox::activated, this, [this](int index) {
        if (index >= 0) {
            commitText(m_combo->itemData(index).toString());
        }
    });
    setChoices(this->spec().choices);
}

void ComboField::setChoices(const std::vector<Choice>& choices) {
    m_choices = choices;
    const RefreshScope scope(*this);
    m_combo->clear();
    for (const Choice& c : m_choices) {
        m_combo->addItem(choiceLabel(c), c.value);
    }
    refresh();
}

void ComboField::refresh() {
    const RefreshScope scope(*this);
    const QString value = documentText();
    // the unknown item (always the one after the choices) shows the value the document holds
    while (m_combo->count() > static_cast<int>(m_choices.size())) {
        m_combo->removeItem(m_combo->count() - 1);
    }
    int index = m_combo->findData(value);
    if (index < 0 && fieldExists()) {
        m_combo->addItem(unknownLabel(value), value);
        index = m_combo->count() - 1;
        setLocalProblem(QStringLiteral("%1 is not one of the choices (kept as it is).")
                            .arg(toQString(core::pythonRepr(core::Value(toStdString(value))))));
    } else {
        setLocalProblem({});
    }
    m_combo->setCurrentIndex(index);
    m_combo->setEnabled(fieldExists());
}

QString ComboField::displayedText() const { return m_combo->currentData().toString(); }

// ---- EditableComboField ----------------------------------------------------------------------

EditableComboField::EditableComboField(FormContext& context, std::string path, FieldSpec spec, QWidget* parent)
    : FieldWidget(context, std::move(path), std::move(spec), parent) {
    auto* layout = new QHBoxLayout(this);
    layout->setContentsMargins(0, 0, 0, 0);
    m_combo = new QComboBox(this);
    m_combo->setEditable(true);
    m_combo->setInsertPolicy(QComboBox::NoInsert);
    m_combo->setDuplicatesEnabled(true);  // Return never swaps the text for an item of other case
    m_combo->setCompleter(nullptr);       // nor completes it: the text stays as typed
    m_combo->setSizeAdjustPolicy(QComboBox::AdjustToMinimumContentsLengthWithIcon);
    m_combo->setMinimumContentsLength(8);
    m_combo->lineEdit()->setPlaceholderText(this->spec().placeholder);
    trackEdits(m_combo->lineEdit());
    guardWheel(m_combo);
    layout->addWidget(m_combo);
    setEditor(m_combo);
    connect(m_combo->lineEdit(), &QLineEdit::editingFinished, this, &EditableComboField::commit);
    connect(m_combo, &QComboBox::activated, this, [this](int index) {
        if (index >= 0 && commitText(m_combo->itemText(index))) {
            setDirty(m_combo->lineEdit(), false);
        }
    });
    m_combo->lineEdit()->installEventFilter(new EscapeReverts(this, [this] {
        this->context().setFormError(this->path(), nullptr);
        refresh();
        updateFlag();
        this->context().scheduleCheck();  // the typed text is gone: check the value shown again
    }));
    QStringList items;
    for (const Choice& c : this->spec().choices) {
        items.append(c.value);
    }
    setChoices(items);
}

void EditableComboField::setChoices(const QStringList& items) {
    const RefreshScope scope(*this);
    const QString text = m_combo->lineEdit()->text();
    const bool modified = isDirty(m_combo->lineEdit());
    m_combo->clear();
    m_combo->addItems(items);
    m_combo->setCurrentIndex(-1);
    m_combo->lineEdit()->setText(text);
    setDirty(m_combo->lineEdit(), modified);
    if (!modified) {
        refresh();
    }
}

void EditableComboField::refresh() {
    const RefreshScope scope(*this);
    const QString text = documentText();
    if (m_combo->lineEdit()->text() != text) {
        m_combo->lineEdit()->setText(text);
        m_combo->lineEdit()->setCursorPosition(0);
    }
    setDirty(m_combo->lineEdit(), false);
    m_combo->setEnabled(fieldExists());
}

void EditableComboField::commit() {
    if (!isDirty(m_combo->lineEdit())) {
        return;
    }
    if (commitText(m_combo->lineEdit()->text())) {
        setDirty(m_combo->lineEdit(), false);
    }
}

QString EditableComboField::displayedText() const { return m_combo->lineEdit()->text(); }

std::optional<QString> EditableComboField::pendingText() const {
    return isDirty(m_combo->lineEdit()) ? std::optional<QString>(m_combo->lineEdit()->text()) : std::nullopt;
}

// ---- CheckField ------------------------------------------------------------------------------

namespace {

// Partly checked only while an optional boolean is unset; a click sets true or false. A long text
// wraps onto further lines when the side bar is narrower than it (MDATK's GTK check buttons wrap),
// so a check box never pushes a page wider than the side bar.
class TwoStateClicks : public QCheckBox {
public:
    TwoStateClicks(const QString& text, QWidget* parent) : QCheckBox(text, parent) {
        QSizePolicy policy(QSizePolicy::Preferred, QSizePolicy::Minimum);  // may shrink below one line, then wraps
        policy.setHeightForWidth(true);
        setSizePolicy(policy);
    }

    QSize minimumSizeHint() const override {
        // the indicator and the longest word (or the whole text when it is shorter)
        int word = 0;
        for (const QString& part : text().split(QLatin1Char(' '), Qt::SkipEmptyParts)) {
            word = std::max(word, fontMetrics().horizontalAdvance(part));
        }
        const QSize full = QCheckBox::sizeHint();
        return {std::min(full.width(), textOffset() + word + 2), full.height()};
    }
    bool hasHeightForWidth() const override { return true; }
    int heightForWidth(int w) const override {
        const QRect bounds = fontMetrics().boundingRect(QRect(0, 0, std::max(w - textOffset(), 1), 100000),
                                                        Qt::TextWordWrap, text());
        return std::max(QCheckBox::sizeHint().height(), bounds.height() + 4);
    }

protected:
    void nextCheckState() override { setCheckState(checkState() == Qt::Checked ? Qt::Unchecked : Qt::Checked); }
    bool hitButton(const QPoint& pos) const override { return rect().contains(pos); }
    void paintEvent(QPaintEvent* event) override {
        if (QCheckBox::sizeHint().width() <= width()) {
            QCheckBox::paintEvent(event);  // fits on one line: the style draws it as usual
            return;
        }
        QStylePainter painter(this);
        QStyleOptionButton option;
        initStyleOption(&option);
        const QString caption = option.text;
        option.text.clear();
        const QRect indicator = style()->subElementRect(QStyle::SE_CheckBoxIndicator, &option, this);
        option.rect = QRect(rect().left(), rect().top(), indicator.width(), rect().height());
        painter.drawControl(QStyle::CE_CheckBox, option);
        const QRect area(textOffset(), 0, width() - textOffset(), height());
        painter.drawItemText(area, Qt::AlignLeft | Qt::AlignVCenter | Qt::TextWordWrap, palette(), isEnabled(), caption,
                             QPalette::WindowText);
        if (hasFocus()) {
            QStyleOptionFocusRect focus;
            focus.initFrom(this);
            focus.rect = area;
            painter.drawPrimitive(QStyle::PE_FrameFocusRect, focus);
        }
    }

private:
    int textOffset() const {
        return style()->pixelMetric(QStyle::PM_IndicatorWidth, nullptr, this) +
               style()->pixelMetric(QStyle::PM_CheckBoxLabelSpacing, nullptr, this);
    }
};

}  // namespace

CheckField::CheckField(FormContext& context, std::string path, FieldSpec spec, QWidget* parent)
    : FieldWidget(context, std::move(path), std::move(spec), parent) {
    auto* layout = new QHBoxLayout(this);
    layout->setContentsMargins(0, 0, 0, 0);
    m_check = new TwoStateClicks(label(), this);
    layout->addWidget(m_check);
    setEditor(m_check);
    connect(m_check, &QCheckBox::clicked, this, [this](bool checked) {
        if (!refreshing()) {
            commitText(checked ? QStringLiteral("true") : QStringLiteral("false"));
        }
    });
    refresh();
}

void CheckField::refresh() {
    const RefreshScope scope(*this);
    const QString text = documentText();
    m_check->setTristate(text.isEmpty());
    m_check->setCheckState(text == QLatin1String("true") ? Qt::Checked
                           : text.isEmpty()              ? Qt::PartiallyChecked
                                                         : Qt::Unchecked);
    m_check->setEnabled(fieldExists());
}

QString CheckField::displayedText() const {
    switch (m_check->checkState()) {
    case Qt::Checked:
        return QStringLiteral("true");
    case Qt::Unchecked:
        return QStringLiteral("false");
    case Qt::PartiallyChecked:
        break;
    }
    return {};
}

}  // namespace mdbridge::app::forms
