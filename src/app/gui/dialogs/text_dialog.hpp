#pragma once
// A dialog showing a selectable text (MDATK's _show_text and About dialog; texts selectable as in
// MDATK 0.1.0's dialogtext.py): an optional heading (rich text), the text in a read-only view, extra
// buttons on the left and Close. The window opens it with open() (window-modal, not blocking), so
// tests can read and close it.

#include <functional>

#include <QDialog>
#include <QString>

class QPlainTextEdit;
class QPushButton;
class QHBoxLayout;

namespace mdbridge::app {

class TextDialog : public QDialog {
    Q_OBJECT

public:
    TextDialog(QWidget* parent, const QString& title, const QString& heading, const QString& text,
               bool monospace = true);

    QString text() const;
    QPlainTextEdit* view() const { return m_view; }
    // A button before Close, objectName "dialog.<name>".
    QPushButton* addButton(const QString& name, const QString& label, std::function<void()> onClick);

private:
    QPlainTextEdit* m_view = nullptr;
    QHBoxLayout* m_buttons = nullptr;
};

}  // namespace mdbridge::app
