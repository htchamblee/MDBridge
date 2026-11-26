#pragma once
// The still preview of *Render still* (MDATK's show_image: a dialog titled "<program> preview –
// <file name>" with the image and Close, not modal). MDATK scales the image to a fixed 1 200 px
// width; here it is scaled to fit the screen the window is on (AP-18): at most 90 % of the
// screen's available width and height, never enlarged, aspect kept. The dialog deletes itself
// when closed.

#include <QDialog>
#include <QImage>
#include <QString>

class QLabel;

namespace mdbridge::app {

class ImageDialog : public QDialog {
    Q_OBJECT

public:
    // Shows nothing until show(); `path` is read once (a file that cannot be read gives a dialog
    // saying so).
    ImageDialog(QWidget* parent, const QString& path);

    QString path() const { return m_path; }
    // The image as read (null when it could not be read) and the size it is shown at.
    const QImage& image() const { return m_image; }
    QSize shownSize() const;

    // The largest size within `room` that keeps the aspect of `image` without enlarging it.
    static QSize fitted(QSize image, QSize room);

private:
    QString m_path;
    QImage m_image;
    QLabel* m_label = nullptr;
};

}  // namespace mdbridge::app
