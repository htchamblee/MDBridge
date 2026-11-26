#include "main_window.hpp"

#include <exception>
#include <iostream>
#include <optional>
#include <string>

#include <QApplication>
#include <QDragEnterEvent>
#include <QDropEvent>
#include <QFileDialog>
#include <QFileInfo>
#include <QLabel>
#include <QMenuBar>
#include <QMessageBox>
#include <QMimeData>
#include <QProgressBar>
#include <QSplitter>
#include <QStatusBar>

#include "../project_summary.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/project_io.hpp"
#include "mdbridge/render/example.hpp"
#include "preview_widget.hpp"
#include "project_view.hpp"
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
    m_right = new QSplitter(Qt::Vertical, m_splitter);
    m_splitter->addWidget(m_right);
    m_preview = new PreviewWidget(m_right);
    m_right->addWidget(m_preview);
    // object names for the GUI tests (tests/gui); field widgets will use their project key
    m_panel->setObjectName(QStringLiteral("projectPanel"));
    m_preview->setObjectName(QStringLiteral("preview"));
    setCentralWidget(m_splitter);

    m_progress = new QProgressBar(this);
    m_progress->setRange(0, 100);
    m_progress->setMaximumWidth(160);
    m_progress->hide();
    statusBar()->addPermanentWidget(m_progress);
    m_viewLabel = new QLabel(this);
    statusBar()->addPermanentWidget(m_viewLabel);
    createMenus();
    setView(render::createView(renderer));
    m_splitter->setStretchFactor(0, 0);
    m_splitter->setStretchFactor(1, 1);
    m_splitter->setSizes({440, 840});
    m_right->setStretchFactor(0, 3);
    m_right->setStretchFactor(1, 2);
    m_right->setSizes({460, 300});
}

MainWindow::~MainWindow() {
    cancelLoad();
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

    m_container = QWidget::createWindowContainer(m_view->window(), m_right);
    m_container->setObjectName(QStringLiteral("view3d"));
    m_container->setMinimumSize(320, 240);
    m_container->setFocusPolicy(Qt::StrongFocus);
    m_container->setToolTip(tr("Drag to rotate, scroll to zoom, double-click to reset the view"));
    m_right->insertWidget(0, m_container);
    if (m_molecule) {
        m_view->setMoleculeFrame(m_molecule);
    }
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
        m_projectPath = path;
        m_project.reset();
        try {
            m_project = core::loadProject(path);
        } catch (const std::exception& e) {  // shown, but cannot be saved
            std::cerr << "mdbridge: " << e.what() << '\n';
        }
        setWindowTitle(tr("%1 - MDBridge").arg(toQString(core::pathToUtf8(path.filename()))));
        statusBar()->showMessage(tr("Opened %1").arg(toQString(core::pathToUtf8(path))), 5000);
        startLoad(path);
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

bool MainWindow::saveProject(const std::filesystem::path& path, bool interactive) {
    const std::filesystem::path target = path.empty() ? m_projectPath : path;
    QString problem;
    if (!m_project) {
        problem = tr("No project with readable settings is open.");
    } else if (target.empty()) {
        problem = tr("The project has no file yet.");
    } else {
        try {
            core::saveProject(target, *m_project);
        } catch (const std::exception& e) {  // core::ProjectError, std::filesystem::filesystem_error
            problem = toQString(e.what());
        }
    }
    if (!problem.isEmpty()) {
        if (interactive) {
            QMessageBox::warning(this, tr("Cannot save the project"), problem);
        } else {
            std::cerr << "mdbridge: " << problem.toStdString() << '\n';
        }
        return false;
    }
    m_projectPath = target;
    setWindowTitle(tr("%1 - MDBridge").arg(toQString(core::pathToUtf8(target.filename()))));
    statusBar()->showMessage(tr("Saved %1").arg(toQString(core::pathToUtf8(target))), 5000);
    return true;
}

namespace {

// What the worker thread loads for an opened project.
struct Loaded {
    bool cancelled = false;
    std::shared_ptr<const render::MoleculeFrame> molecule;
    QString structureProblem;  // why the 3D view has no structure
    QString summary;
    QStringList notes;
    std::shared_ptr<const PreviewContent> preview;
    QString previewProblem;
};

QString describe(const std::exception& e) { return QString::fromUtf8(e.what()); }

Loaded loadStructure(const std::filesystem::path& path, const core::Reporter& reporter) {
    Loaded out;
    core::Project project;
    try {
        project = core::loadProject(path);
    } catch (const std::exception& e) {
        out.structureProblem = describe(e);
        out.previewProblem = describe(e);
        return out;
    }
    std::optional<render::StructurePreview> preview;
    try {
        const ReferenceFrame ref = loadReferenceFrame(project, reporter);
        for (const std::string& note : ref.notes) out.notes.push_back(toQString(note));
        try {
            preview = structurePreview(ref.dataset);
        } catch (const std::exception& e) {
            out.notes.push_back(MainWindow::tr("preview: %1").arg(describe(e)));
        }
        try {
            const StructureView view = structureView(project, ref.dataset);
            auto frame = std::make_shared<render::MoleculeFrame>();
            frame->items = view.items;
            frame->rotation = view.viewport.rotation;
            frame->xlim = view.viewport.xlim;
            frame->ylim = view.viewport.ylim;
            frame->zrange = view.viewport.zrange;
            frame->shading = view.shading;
            frame->background = view.background;
            out.molecule = std::move(frame);
            out.summary = MainWindow::tr("Reference frame %1: %2 residues, %3 atoms%4")
                              .arg(ref.frame)
                              .arg(view.residues)
                              .arg(view.atoms)
                              .arg(ref.cleaned ? MainWindow::tr(" (cleaned files)") : QString());
        } catch (const std::exception& e) {
            out.structureProblem = describe(e);
        }
    } catch (const core::Cancelled&) {
        out.cancelled = true;
        return out;
    } catch (const std::exception& e) {
        out.structureProblem = describe(e);
    }
    try {
        auto content = std::make_shared<PreviewContent>();
        content->example = render::exampleDataset(project, preview ? &*preview : nullptr);
        content->watermark = previewWatermark(preview.has_value());
        content->project = std::move(project);
        out.preview = std::move(content);
    } catch (const std::exception& e) {
        out.previewProblem = MainWindow::tr("Cannot draw this layout: %1").arg(describe(e));
    }
    return out;
}

}  // namespace

void MainWindow::cancelLoad() {
    if (m_cancel) {
        m_cancel->store(true);
        m_cancel.reset();
    }
    ++m_loadGeneration;
    m_loading = false;
    m_progress->hide();
}

void MainWindow::startLoad(const std::filesystem::path& path) {
    cancelLoad();
    const std::uint64_t generation = m_loadGeneration;
    m_cancel = std::make_shared<std::atomic<bool>>(false);
    m_loading = true;
    m_loadProblem.clear();
    m_structureSummary.clear();
    m_molecule.reset();
    m_view->setMoleculeFrame(nullptr);
    m_preview->setMessage(tr("Loading the project's structure ..."));
    m_progress->setValue(0);
    m_progress->show();
    statusBar()->showMessage(tr("Loading the structure ..."));

    const auto post = m_background.poster();
    const core::Reporter reporter(
        {},
        [this, post, generation](double fraction, const std::string& message) {
            post([this, generation, fraction, text = toQString(message)] {
                if (generation == m_loadGeneration) {
                    m_progress->setValue(static_cast<int>(fraction * 100));
                    if (!text.isEmpty()) statusBar()->showMessage(tr("Loading the structure: %1").arg(text));
                }
            });
        },
        m_cancel);
    m_background.run([path, reporter] { return loadStructure(path, reporter); },
                     [this, generation](Loaded loaded) {
                         if (generation != m_loadGeneration || loaded.cancelled) {
                             return;
                         }
                         m_loading = false;
                         m_cancel.reset();
                         m_progress->hide();
                         m_molecule = loaded.molecule;
                         m_structureSummary = loaded.summary;
                         m_loadProblem = loaded.structureProblem;
                         if (m_molecule) {
                             m_view->setMoleculeFrame(m_molecule);
                         }
                         if (loaded.preview) {
                             m_preview->setContent(loaded.preview);
                         } else {
                             m_preview->setMessage(loaded.previewProblem);
                         }
                         QString message = m_molecule ? m_structureSummary
                                                      : tr("3D view: the structure could not be loaded: %1").arg(m_loadProblem);
                         if (!loaded.notes.isEmpty()) {
                             message += QStringLiteral(" (") + loaded.notes.join(QStringLiteral("; ")) + QLatin1Char(')');
                         }
                         statusBar()->showMessage(message);
                         if (!m_molecule) {
                             std::cerr << "mdbridge: " << message.toStdString() << '\n';
                         }
                         emit projectLoaded();
                     });
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
