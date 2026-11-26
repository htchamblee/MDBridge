#include "gui/pages/page_section.hpp"

#include <algorithm>

#include <QFormLayout>
#include <QScrollBar>
#include <QScrollArea>
#include <QVBoxLayout>

namespace mdbridge::app::pages {

PageSection::PageSection(const QString& page, const QString& section, const QString& title, QWidget* parent)
    : QGroupBox(title, parent), m_section(section) {
    setObjectName(QStringLiteral("section:%1.%2").arg(page, section));
    setCheckable(true);  // the title's check box folds the section (MDATK's frames do not fold)
    setChecked(true);
    setAccessibleName(title);
    setToolTip(tr("Click the title to fold or unfold this section"));
    auto* outer = new QVBoxLayout(this);
    outer->setContentsMargins(8, 4, 8, 8);
    m_content = new QWidget(this);
    m_form = new QFormLayout(m_content);
    m_form->setContentsMargins(0, 0, 0, 0);
    m_form->setFieldGrowthPolicy(QFormLayout::AllNonFixedFieldsGrow);
    m_form->setRowWrapPolicy(QFormLayout::WrapLongRows);
    outer->addWidget(m_content);
    // a checkable group box disables its children when unchecked; here it folds them instead
    connect(this, &QGroupBox::toggled, this, [this](bool on) {
        m_content->setEnabled(true);
        m_content->setVisible(on);
    });
}

bool PageSection::expanded() const { return isChecked(); }

void PageSection::setExpanded(bool expanded) { setChecked(expanded); }

PageSection* findSection(QWidget* page, const QString& pageId, const QString& section) {
    if (!page) return nullptr;
    return page->findChild<PageSection*>(QStringLiteral("section:%1.%2").arg(pageId, section));
}

bool scrollIntoView(QWidget* widget) {
    if (!widget) return false;
    QScrollArea* area = nullptr;
    for (QWidget* w = widget->parentWidget(); w; w = w->parentWidget()) {
        if (auto* section = qobject_cast<PageSection*>(w); section && !section->expanded()) {
            section->setExpanded(true);
        }
        if (auto* a = qobject_cast<QScrollArea*>(w); a && !area) {
            area = a;
        }
    }
    if (!area) return false;
    if (QWidget* inner = area->widget()) {
        inner->adjustSize();  // the layout of sections just unfolded
    }
    if (qobject_cast<PageSection*>(widget)) {
        // a section: its top at the top of the view (MDATK's open_section: 6 px above it)
        const QPoint at = widget->mapTo(area->widget(), QPoint(0, 0));
        area->verticalScrollBar()->setValue(std::max(0, at.y() - 6));
    } else {
        area->ensureWidgetVisible(widget, 20, 40);
    }
    return true;
}

}  // namespace mdbridge::app::pages
