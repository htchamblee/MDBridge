#include "gui/dialogs/image_dialog.hpp"

#include <algorithm>

#include <QDialogButtonBox>
#include <QFileInfo>
#include <QGuiApplication>
#include <QImageReader>
#include <QLabel>
#include <QPixmap>
#include <QPushButton>
#include <QScreen>
#include <QVBoxLayout>

namespace mdbridge::app {

namespace {

constexpr double kScreenPart = 0.9;   // of the screen's available width and height
constexpr int kFrameRoom = 120;       // for the dialog's frame, margins and button row

}  // namespace

QSize ImageDialog::fitted(QSize image, QSize room) {
    if (image.isEmpty()) return {};
    if (room.width() < 1 || room.height() < 1) return {1, 1};
    if (image.width() <= room.width() && image.height() <= room.height()) return image;
    return image.scaled(room, Qt::KeepAspectRatio).expandedTo(QSize(1, 1));
}

ImageDialog::ImageDialog(QWidget* parent, const QString& path) : QDialog(parent), m_path(path) {
    setObjectName(QStringLiteral("stillDialog"));
    setAttribute(Qt::WA_DeleteOnClose);
    setModal(false);
    setWindowTitle(tr("MDBridge preview – %1").arg(QFileInfo(path).fileName()));
    auto* layout = new QVBoxLayout(this);
    m_label = new QLabel(this);
    m_label->setObjectName(QStringLiteral("still.image"));
    m_label->setAlignment(Qt::AlignCenter);
    QImageReader reader(path);
    m_image = reader.read();
    if (m_image.isNull()) {
        m_label->setText(tr("Cannot show %1: %2").arg(path, reader.errorString()));
        m_label->setTextInteractionFlags(Qt::TextSelectableByMouse);
    } else {
        const QScreen* screen = parent && parent->screen() ? parent->screen() : QGuiApplication::primaryScreen();
        QSize room(1200, 900);
        if (screen) {
            const QRect area = screen->availableGeometry();
            room = QSize(static_cast<int>(area.width() * kScreenPart), static_cast<int>(area.height() * kScreenPart) - kFrameRoom);
        }
        const QSize shown = fitted(m_image.size(), room);
        QPixmap pixmap = QPixmap::fromImage(m_image);
        if (shown != m_image.size()) pixmap = pixmap.scaled(shown, Qt::KeepAspectRatio, Qt::SmoothTransformation);
        m_label->setPixmap(pixmap);
        m_label->setFixedSize(pixmap.size());
    }
    layout->addWidget(m_label, 1, Qt::AlignCenter);
    auto* buttons = new QDialogButtonBox(QDialogButtonBox::Close, this);
    buttons->button(QDialogButtonBox::Close)->setObjectName(QStringLiteral("still.close"));
    connect(buttons, &QDialogButtonBox::rejected, this, &QDialog::close);
    layout->addWidget(buttons);
}

QSize ImageDialog::shownSize() const {
    const QPixmap pixmap = m_label->pixmap();
    return pixmap.isNull() ? QSize() : pixmap.size();
}

}  // namespace mdbridge::app
