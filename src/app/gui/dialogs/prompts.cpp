#include "gui/dialogs/prompts.hpp"

#include <iostream>
#include <memory>
#include <string>

#include <QAbstractButton>
#include <QLabel>
#include <QMessageBox>
#include <QPushButton>

namespace mdbridge::app {

namespace {

// Details longer than this go into the message box's scrollable "Show Details" area.
constexpr qsizetype kInlineDetail = 600;

void setDetail(QMessageBox& box, const QString& detail) {
    if (detail.size() > kInlineDetail) {
        box.setInformativeText(detail.left(kInlineDetail) + QStringLiteral(" …"));
        box.setDetailedText(detail);
    } else {
        box.setInformativeText(detail);
    }
    // the texts can be selected and copied (MDATK 0.1.0, dialogtext.py)
    for (QLabel* label : box.findChildren<QLabel*>()) {
        label->setTextInteractionFlags(Qt::TextSelectableByMouse | Qt::TextSelectableByKeyboard);
    }
}

}  // namespace

QString Prompter::ask(QWidget* parent, const Question& question) {
    QMessageBox box(parent);
    box.setIcon(question.kind == QLatin1String("unsaved") || question.kind == QLatin1String("recover")
                    ? QMessageBox::Question
                    : QMessageBox::Warning);
    box.setWindowTitle(question.title);
    box.setText(question.text);
    setDetail(box, question.detail);
    std::vector<std::pair<QAbstractButton*, QString>> buttons;
    for (const auto& [id, label] : question.choices) {
        QMessageBox::ButtonRole role = QMessageBox::AcceptRole;
        if (id == QLatin1String("cancel")) role = QMessageBox::RejectRole;
        else if (id == QLatin1String("discard")) role = QMessageBox::DestructiveRole;
        QPushButton* button = box.addButton(label, role);
        button->setObjectName(QStringLiteral("answer.") + id);
        if (id == QLatin1String("cancel")) box.setEscapeButton(button);
        buttons.emplace_back(button, id);
    }
    for (const auto& [button, id] : buttons) {
        if (question.defaultId.isEmpty() ? button == buttons.front().first : id == question.defaultId) {
            box.setDefaultButton(static_cast<QPushButton*>(button));
            break;
        }
    }
    box.exec();
    for (const auto& [button, id] : buttons) {
        if (box.clickedButton() == button) return id;
    }
    for (const auto& [id, label] : question.choices) {  // closed by the window manager
        if (id == QLatin1String("cancel")) return id;
    }
    return question.choices.empty() ? QString() : question.choices.back().first;
}

void Prompter::message(QWidget* parent, MessageKind kind, const QString& text, const QString& detail) {
    messageBox(parent, kind, text, detail)->exec();
}

std::unique_ptr<QMessageBox> messageBox(QWidget* parent, MessageKind kind, const QString& text, const QString& detail) {
    auto box = std::make_unique<QMessageBox>(parent);
    box->setIcon(kind == MessageKind::Error ? QMessageBox::Critical : QMessageBox::Information);
    box->setWindowTitle(kind == MessageKind::Error ? QObject::tr("Error") : QObject::tr("Information"));
    box->setText(text);
    setDetail(*box, detail);
    box->setStandardButtons(QMessageBox::Ok);
    return box;
}

void Prompter::printError(const QString& line) { std::cerr << line.toStdString() << std::endl; }

Prompter& defaultPrompter() {
    static Prompter prompter;
    return prompter;
}

}  // namespace mdbridge::app
