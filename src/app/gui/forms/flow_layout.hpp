#pragma once
// A row of widgets that wraps onto further lines when the side bar is too narrow for all of them
// (GTK's FlowBox for MDATK's form rows): items keep their size hint (shrunk to the line width
// when one alone is wider, never below its minimum) and are placed left to right, vertically
// centred on each line. The layout reports its height for a width, so a scroll area sized to
// its content grows by the extra lines.

#include <QLayout>
#include <QList>
#include <QRect>
#include <QSize>

class QLayoutItem;
class QWidget;

namespace mdbridge::app::forms {

class FlowLayout : public QLayout {
public:
    explicit FlowLayout(QWidget* parent = nullptr, int horizontalSpacing = 6, int verticalSpacing = 4);
    ~FlowLayout() override;
    FlowLayout(const FlowLayout&) = delete;
    FlowLayout& operator=(const FlowLayout&) = delete;

    void addItem(QLayoutItem* item) override;
    // Inserts `widget` before item `index` (at the end for an index out of range).
    void insertWidget(int index, QWidget* widget);
    int count() const override;
    QLayoutItem* itemAt(int index) const override;
    QLayoutItem* takeAt(int index) override;
    Qt::Orientations expandingDirections() const override;
    bool hasHeightForWidth() const override;
    int heightForWidth(int width) const override;
    QSize minimumSize() const override;
    QSize sizeHint() const override;
    void setGeometry(const QRect& rect) override;

private:
    int arrange(const QRect& rect, bool apply) const;
    QList<QLayoutItem*> m_items;
    int m_hSpacing;
    int m_vSpacing;
};

}  // namespace mdbridge::app::forms
