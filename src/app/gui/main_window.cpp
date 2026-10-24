#include "main_window.hpp"

#include <iostream>

#include <QApplication>
#include <QDragEnterEvent>
#include <QDropEvent>
#include <QFileDialog>
#include <QFileInfo>
#include <QLabel>
#include <QMenuBar>
#include <QMessageBox>
#include <QMimeData>
#include <QSplitter>
#include <QStatusBar>

#include "../project_summary.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/render/vulkan_probe.hpp"
#include "mdbridge/version.hpp"
#include "project_panel.hpp"

namespace mdbridge::app {

namespace {

QString toQString(const std::string& text) { return QString::fromUtf8(text.data(), static_cast<qsizetype>(text.size())); }

std::filesystem::path toPath(const QString& file) { return core::pathFromUtf8(file.toStdString()); }

}  // namespace

MainWindow::MainWindow(render::Renderer renderer, QWidget* parent) : QMainWindow(parent) {
    setWindowTitle(QStringLiteral("MDBridge"));
    setAcceptDrops(true);

    m_splitter = new QSplitter(this);
    m_panel = new ProjectPanel(m_splitter);
    m_splitter->addWidget(m_panel);
    setCentralWidget(m_splitter);

    m_viewLabel = new QLabel(this);
    statusBar()->addPermanentWidget(m_viewLabel);
    createMenus();
    setView(render::createView(renderer));
    m_splitter->setStretchFactor(0, 0);
    m_splitter->setStretchFactor(1, 1);
    m_splitter->setSizes({440, 840});
}

MainWindow::~MainWindow() {
    delete m_container;  // deletes the view's window before the view itself
}

void MainWindow::createMenus() {
    QMenu* file = menuBar()->addMenu(tr("&File"));
    file->addAction(tr("&Open project..."), QKeySequence::Open, this, &MainWindow::chooseProject);
    file->addSeparator();
    file->addAction(tr("&Quit"), QKeySequence::Quit, qApp, &QApplication::quit);

    QMenu* view = menuBar()->addMenu(tr("&View"));
    view->addAction(tr("Use the &software 3D view"), this, [this] { useSoftwareView(tr("chosen in the View menu")); });
    view->addAction(tr("Use the &Vulkan 3D view"), this, [this] { setView(render::createView(render::Renderer::Vulkan)); });

    QMenu* help = menuBar()->addMenu(tr("&Help"));
    help->addAction(tr("&Graphics..."), this, &MainWindow::showGraphics);
    help->addAction(tr("&About MDBridge"), this, &MainWindow::showAbout);
}

void MainWindow::setView(render::ViewChoice choice) {
    // the old window goes with its container; then the old view
    delete m_container;
    m_container = nullptr;
    m_view = std::move(choice.view);
    m_fallbackReason = choice.fallbackReason;

    m_container = QWidget::createWindowContainer(m_view->window(), m_splitter);
    m_container->setMinimumSize(320, 240);
    m_container->setFocusPolicy(Qt::StrongFocus);
    m_splitter->insertWidget(1, m_container);
    connect(m_view.get(), &render::View::frameRendered, this, &MainWindow::viewFrameRendered);
    // queued: the failing view must not be deleted while it is still emitting
    connect(m_view.get(), &render::View::failed, this, &MainWindow::useSoftwareView, Qt::QueuedConnection);

    m_viewLabel->setText(m_view->description());
    if (!m_fallbackReason.isEmpty() && m_view->renderer() == render::Renderer::Software) {
        std::cerr << "mdbridge: 3D view: " << m_view->description().toStdString() << '\n';
    }
}

void MainWindow::useSoftwareView(const QString& reason) {
    setView({render::createSoftwareView(reason), reason});
}

bool MainWindow::openProject(const std::filesystem::path& path, bool interactive) {
    try {
        const OpenedProject project = app::openProject(path);
        m_panel->showProject(project);
        setWindowTitle(tr("%1 - MDBridge").arg(toQString(core::pathToUtf8(path.filename()))));
        statusBar()->showMessage(tr("Opened %1").arg(toQString(core::pathToUtf8(path))), 5000);
        return true;
    } catch (const core::ProjectError& error) {
        if (interactive) {
            QMessageBox::warning(this, tr("Cannot open the project"), toQString(error.what()));
        } else {
            std::cerr << "mdbridge: " << error.what() << '\n';
        }
        return false;
    }
}

void MainWindow::chooseProject() {
    const QString file = QFileDialog::getOpenFileName(
        this, tr("Open project"), QString(),
        tr("MDATK projects (*.mdatk.toml);;Parameter logs of runs (*params_log.toml);;TOML files (*.toml);;All files (*)"));
    if (!file.isEmpty()) {
        openProject(toPath(file));
    }
}

void MainWindow::dragEnterEvent(QDragEnterEvent* event) {
    const QList<QUrl> urls = event->mimeData()->urls();
    if (urls.size() == 1 && urls.front().isLocalFile()) {
        event->acceptProposedAction();
    }
}

void MainWindow::dropEvent(QDropEvent* event) {
    const QList<QUrl> urls = event->mimeData()->urls();
    if (urls.size() == 1 && urls.front().isLocalFile()) {
        openProject(toPath(urls.front().toLocalFile()));
        event->acceptProposedAction();
    }
}

void MainWindow::showGraphics() {
    const render::VulkanReport report = render::probeVulkan();
    QString text = tr("<p><b>3D view:</b> %1</p>").arg(m_view->description().toHtmlEscaped());
    text += tr("<p><b>Vulkan loader:</b> %1</p>")
                .arg(report.loaderFound ? tr("found, Vulkan %1").arg(toQString(render::vulkanVersionText(report.loaderVersion)))
                                        : tr("not found"));
    if (!report.devices.empty()) {
        text += QStringLiteral("<ul>");
        for (const render::VulkanDevice& d : report.devices) {
            text += tr("<li>%1 (%2, Vulkan %3)%4</li>")
                        .arg(toQString(d.name).toHtmlEscaped(), toQString(d.type),
                             toQString(render::vulkanVersionText(d.apiVersion)),
                             d.meetsMinimum() ? QString() : tr(": needs Vulkan %1.%2").arg(render::kVulkanMajor).arg(render::kVulkanMinor));
        }
        text += QStringLiteral("</ul>");
    }
    if (!report.problem.empty()) {
        text += tr("<p>%1</p>").arg(toQString(report.problem).toHtmlEscaped());
    }
    QMessageBox::information(this, tr("Graphics"), text);
}

void MainWindow::showAbout() {
    QMessageBox::about(
        this, tr("About MDBridge"),
        tr("<h3>MDBridge %1</h3>"
           "<p>Analysis and visualisation of molecular dynamics trajectories of RNA systems; the successor of "
           "MDATK. Opens MDATK %2 project files.</p>"
           "<p>Copyright © 2026 the MDBridge author. All rights reserved; the licence follows a legal review.</p>"
           "<p>Built with Qt %3 (LGPL v3), toml++ (MIT), CLI11 (BSD-3-Clause) and the Vulkan headers "
           "(Apache 2.0 / MIT); see the licenses folder of this program.</p>"
           "<p><a href=\"%4\">%4</a></p>")
            .arg(QString::fromLatin1(kVersion), QString::fromLatin1(kMdatkReference), QString::fromLatin1(qVersion()),
                 QString::fromLatin1(kHomepage)));
}

}  // namespace mdbridge::app
