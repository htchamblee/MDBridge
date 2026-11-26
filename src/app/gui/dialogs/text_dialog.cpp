#include "gui/dialogs/text_dialog.hpp"

#include <utility>

#include <QFontDatabase>
#include <QHBoxLayout>
#include <QLabel>
#include <QPlainTextEdit>
#include <QPushButton>
#include <QVBoxLayout>

namespace mdbridge::app {

TextDialog::TextDialog(QWidget* parent, const QString& title, const QString& heading, const QString& text,
                       bool monospace)
    : QDialog(parent) {
    setWindowTitle(title);
    setAttribute(Qt::WA_DeleteOnClose);
    auto* layout = new QVBoxLayout(this);
    if (!heading.isEmpty()) {
        auto* label = new QLabel(heading, this);
        label->setObjectName(QStringLiteral("dialog.heading"));
        label->setWordWrap(true);
        label->setTextFormat(Qt::RichText);
        label->setOpenExternalLinks(true);
        label->setTextInteractionFlags(Qt::TextBrowserInteraction);
        layout->addWidget(label);
    }
    m_view = new QPlainTextEdit(this);
    m_view->setObjectName(QStringLiteral("dialog.text"));
    m_view->setReadOnly(true);
    m_view->setPlainText(text);
    m_view->setAccessibleName(title);
    m_view->setLineWrapMode(QPlainTextEdit::WidgetWidth);
    if (monospace) {
        m_view->setFont(QFontDatabase::systemFont(QFontDatabase::FixedFont));
    }
    m_view->setMinimumSize(720, 420);  // MDATK: 720 x 420
    layout->addWidget(m_view, 1);
    m_buttons = new QHBoxLayout;
    m_buttons->addStretch(1);
    auto* close = new QPushButton(tr("Close"), this);
    close->setObjectName(QStringLiteral("dialog.close"));
    close->setDefault(true);
    connect(close, &QPushButton::clicked, this, &QDialog::accept);
    m_buttons->addWidget(close);
    layout->addLayout(m_buttons);
}

QString TextDialog::text() const { return m_view->toPlainText(); }

QPushButton* TextDialog::addButton(const QString& name, const QString& label, std::function<void()> onClick) {
    auto* button = new QPushButton(label, this);
    button->setObjectName(QStringLiteral("dialog.") + name);
    connect(button, &QPushButton::clicked, this, [fn = std::move(onClick)] { fn(); });
    m_buttons->insertWidget(m_buttons->count() - 2, button);  // before the stretch and Close
    return button;
}

}  // namespace mdbridge::app
