#pragma once
// A settings page built from the project model until its own page is delivered (milestone 5 W1-A;
// the pages of W1-B / W1-C / W2-A replace it in page_registry.cpp): every field the page shows
// (core::problemPage) with the field widget forms::createFieldWidget gives it, one collapsible
// section per part of the project ("inputs", "analyses.contacts", "output.video" ...), lists of
// records as list editors. It edits the document like the real page will, so the window is usable
// and the GUI cases that only need fields run, but it has none of MDATK's page layout, hints or
// buttons.

#include <string>
#include <vector>

#include <QString>
#include <QWidget>

#include "gui/forms/form_context.hpp"

class QVBoxLayout;

namespace mdbridge::app::pages {

class PlaceholderPage : public QWidget {
    Q_OBJECT

public:
    // `id`: one of core::kProblemPages.
    PlaceholderPage(forms::FormContext& context, const QString& id, QWidget* parent = nullptr);

    QString pageId() const { return m_id; }
    QString title() const;
    // Shows and focuses the field of `path` (or the list editor holding it); false when no field of
    // this page edits it.
    bool revealField(const QString& path);

private:
    void rebuild();
    std::vector<std::string> pagePaths() const;

    forms::FormContext& m_context;
    QString m_id;
    std::vector<std::string> m_paths;
    QWidget* m_body = nullptr;
    QVBoxLayout* m_layout = nullptr;
};

}  // namespace mdbridge::app::pages
