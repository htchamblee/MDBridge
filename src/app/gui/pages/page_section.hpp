#pragma once
// A collapsible section of a settings page (MDATK's framed sections, app.py): a group box whose
// title toggles its content, objectName "section:<page>.<section>" (the window's openSection()
// and the GUI tests find it by that name). Pages put their rows into form() or content().

#include <QGroupBox>
#include <QString>

class QFormLayout;
class QWidget;

namespace mdbridge::app::pages {

class PageSection : public QGroupBox {
    Q_OBJECT

public:
    // `page` and `section` make the objectName ("section:analyses.contacts").
    PageSection(const QString& page, const QString& section, const QString& title, QWidget* parent = nullptr);

    QString sectionId() const { return m_section; }
    // The rows: a form layout (label, field) inside content().
    QFormLayout* form() const { return m_form; }
    QWidget* content() const { return m_content; }

    bool expanded() const;
    void setExpanded(bool expanded);

private:
    QString m_section;
    QWidget* m_content = nullptr;
    QFormLayout* m_form = nullptr;
};

// The section of `page` (a page widget or any widget inside one) named `section` (the part after
// "section:<page>."), or null.
PageSection* findSection(QWidget* page, const QString& pageId, const QString& section);

// Scrolls the scroll area that holds `widget` so that `widget` shows (with a margin), expanding
// the sections around it. Returns false when `widget` is in no scroll area (nothing to scroll).
bool scrollIntoView(QWidget* widget);

}  // namespace mdbridge::app::pages
