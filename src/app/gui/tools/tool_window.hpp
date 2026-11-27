#pragma once
// What the PDB clean-up and structure builder windows share (milestone 6 W3-B): MDATK's layout of
// both (gui/pdbdialog.py, gui/builderdialog.py): the form on the left (a scroll area 500 px wide at
// least, the splitter at 520 px), the report on the right (bold "Report", a read-only monospaced
// text), the buttons at the bottom right; and the form's controls, each named "<tool>.<key>" after
// MDATK's attribute (pdb_tools.input, builder.seq ...; the GUI cases and tests find them so), with
// an accessible name and MDATK's tooltip. Value widgets ignore the mouse wheel unless focused
// (forms::guardWheel); long drop-down labels are elided rather than widening the form.
//
// The windows are top-level windows above their host's window (MDATK's transient_for), deleted
// when closed (the host's tool registry forgets them, as MDATK's destroy handler did).

#include <cstdint>
#include <initializer_list>
#include <utility>
#include <vector>

#include <QPlainTextEdit>
#include <QString>
#include <QStringList>
#include <QWidget>

#include "gui/forms/field_spec.hpp"

class QCheckBox;
class QComboBox;
class QDoubleSpinBox;
class QHBoxLayout;
class QLineEdit;
class QPushButton;
class QScrollArea;
class QSpinBox;
class QSplitter;
class QToolButton;
class QVBoxLayout;

namespace mdbridge::app::tools {

class Host;

// A path typed or picked (MDATK's W.file_chooser_button): the line edit "<tool>.<key>" and its
// Browse button "<tool>.<key>:browse" (forms::runFileDialog, which tests replace).
class PathEntry : public QWidget {
    Q_OBJECT

public:
    PathEntry(const QString& name, const QString& dialogTitle, forms::PathMode mode, QStringList filters,
              QWidget* parent = nullptr);
    QLineEdit* edit() const { return m_edit; }
    QToolButton* browseButton() const { return m_browse; }
    // The text, stripped (MDATK's entry.get_text().strip()).
    QString path() const;
    void setPath(const QString& path);
    void browse();

signals:
    void pathChanged(const QString& text);

private:
    QLineEdit* m_edit = nullptr;
    QToolButton* m_browse = nullptr;
    QString m_title;
    forms::PathMode m_mode;
    QStringList m_filters;
};

class ToolWindow : public QWidget {
    Q_OBJECT

public:
    // `id`: "pdb_tools" or "builder", the prefix of the controls' names.
    ToolWindow(Host& host, QString id, const QString& title, QWidget* parent);
    ~ToolWindow() override;

    Host& host() const { return m_host; }
    QString toolId() const { return m_id; }
    // The control "<tool>.<key>" (null when there is none).
    QWidget* control(const QString& key) const;
    QPlainTextEdit* report() const { return m_report; }
    QString reportText() const { return m_report->toPlainText(); }
    // The report shows `text` (MDATK's _show).
    void showReport(const QString& text);

protected:
    // ---- the layout ----
    // The splitter; the form goes in first (addForm), then the report pane (addReport).
    QSplitter* splitter() const { return m_splitter; }
    // A scroll area for one form column (MDATK's W.scrolled), its body's layout.
    QScrollArea* makeFormScroll(const QString& key, const QString& accessibleName, QVBoxLayout** column, QWidget* parent);
    // The report pane (and returns its layout, for the builder's progress bar).
    QVBoxLayout* addReport();
    // The button bar's buttons, left to right (MDATK packs them from the right).
    QPushButton* addButton(const QString& key, const QString& text, const QString& tooltip = QString());
    // Done adding: sizes (MDATK's default size, the splitter at 520 px).
    void finishLayout(int width, int height);

    // ---- the controls (named "<tool>.<key>") ----
    QLineEdit* makeLine(const QString& key, const QString& accessibleName, const QString& placeholder = QString(),
                        const QString& tooltip = QString(), int widthChars = 0);
    QComboBox* makeCombo(const QString& key, const QString& accessibleName,
                         std::initializer_list<std::pair<const char*, const char*>> items, const char* active);
    QCheckBox* makeCheck(const QString& key, const QString& text, bool checked = false, const QString& tooltip = QString());
    QSpinBox* makeSpin(const QString& key, const QString& accessibleName, int lo, int hi, int value);
    QDoubleSpinBox* makeDoubleSpin(const QString& key, const QString& accessibleName, double lo, double hi, double step,
                                   int decimals, double value);
    PathEntry* makePath(const QString& key, const QString& accessibleName, const QString& dialogTitle,
                        forms::PathMode mode, const std::vector<forms::FileFilter>& filters);

    // The combo's current id (itemData).
    static QString currentId(const QComboBox* combo);

    // ---- closing ----
    // false keeps the window open (the builder: a job that should not be abandoned).
    virtual bool mayClose() { return true; }
    void closeEvent(QCloseEvent* event) override;

    QString nameOf(const QString& key) const { return m_id + QLatin1Char('.') + key; }

private:
    Host& m_host;
    QString m_id;
    QSplitter* m_splitter = nullptr;
    QPlainTextEdit* m_report = nullptr;
    QHBoxLayout* m_bar = nullptr;
};

}  // namespace mdbridge::app::tools
