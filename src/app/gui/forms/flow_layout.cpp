#include "gui/forms/flow_layout.hpp"

#include <algorithm>
#include <vector>

#include <QLayoutItem>
#include <QMargins>
#include <QWidget>

namespace mdbridge::app::forms {

FlowLayout::FlowLayout(QWidget* parent, int horizontalSpacing, int verticalSpacing)
    : QLayout(parent), m_hSpacing(horizontalSpacing), m_vSpacing(verticalSpacing) {}

FlowLayout::~FlowLayout() {
    while (QLayoutItem* item = takeAt(0)) {
        delete item;
    }
}

void FlowLayout::addItem(QLayoutItem* item) { m_items.append(item); }

void FlowLayout::insertWidget(int index, QWidget* widget) {
    addChildWidget(widget);
    auto* item = new QWidgetItem(widget);
    if (index < 0 || index > m_items.size()) {
        m_items.append(item);
    } else {
        m_items.insert(index, item);
    }
    invalidate();
}

int FlowLayout::count() const { return static_cast<int>(m_items.size()); }

QLayoutItem* FlowLayout::itemAt(int index) const {
    return index >= 0 && index < m_items.size() ? m_items.at(index) : nullptr;
}

QLayoutItem* FlowLayout::takeAt(int index) {
    return index >= 0 && index < m_items.size() ? m_items.takeAt(index) : nullptr;
}

Qt::Orientations FlowLayout::expandingDirections() const { return {}; }

bool FlowLayout::hasHeightForWidth() const { return true; }

int FlowLayout::heightForWidth(int width) const { return arrange(QRect(0, 0, width, 0), false); }

QSize FlowLayout::minimumSize() const {
    QSize size;
    for (const QLayoutItem* item : m_items) {
        if (!item->isEmpty()) {
            size = size.expandedTo(item->minimumSize());
        }
    }
    const QMargins m = contentsMargins();
    return size + QSize(m.left() + m.right(), m.top() + m.bottom());
}

QSize FlowLayout::sizeHint() const {
    int lineWidth = 0;
    int lineHeight = 0;
    for (const QLayoutItem* item : m_items) {
        if (item->isEmpty()) {
            continue;
        }
        const QSize hint = item->sizeHint();
        lineWidth += (lineWidth > 0 ? m_hSpacing : 0) + hint.width();
        lineHeight = std::max(lineHeight, hint.height());
    }
    const QMargins m = contentsMargins();
    return QSize(lineWidth + m.left() + m.right(), lineHeight + m.top() + m.bottom());
}

void FlowLayout::setGeometry(const QRect& rect) {
    QLayout::setGeometry(rect);
    arrange(rect, true);
}

int FlowLayout::arrange(const QRect& rect, bool apply) const {
    const QRect area = rect.marginsRemoved(contentsMargins());
    const int available = std::max(area.width(), 1);
    struct Placed {
        QLayoutItem* item;
        int width;
        int height;
    };
    std::vector<Placed> line;
    int y = area.y();
    int lineWidth = 0;
    int lineHeight = 0;
    const auto finishLine = [&] {
        int x = area.x();
        for (const Placed& p : line) {
            if (apply) {
                p.item->setGeometry(QRect(x, y + (lineHeight - p.height) / 2, p.width, p.height));
            }
            x += p.width + m_hSpacing;
        }
        if (!line.empty()) {
            y += lineHeight + m_vSpacing;
        }
        line.clear();
        lineWidth = 0;
        lineHeight = 0;
    };
    for (QLayoutItem* item : m_items) {
        if (item->isEmpty()) {
            continue;
        }
        const QSize hint = item->sizeHint();
        const int width = std::max(std::min(hint.width(), available), item->minimumSize().width());
        const int height = std::max(hint.height(), item->minimumSize().height());
        if (!line.empty() && lineWidth + m_hSpacing + width > available) {
            finishLine();
        }
        lineWidth += (line.empty() ? 0 : m_hSpacing) + width;
        lineHeight = std::max(lineHeight, height);
        line.push_back({item, width, height});
    }
    finishLine();
    const int used = y > area.y() ? y - m_vSpacing - area.y() : 0;
    return used + contentsMargins().top() + contentsMargins().bottom();
}

}  // namespace mdbridge::app::forms
