#include "main_window.hpp"

#include <algorithm>
#include <array>
#include <cerrno>
#include <chrono>
#include <exception>
#include <fstream>
#include <iostream>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <iterator>
#include <system_error>
#include <tuple>
#include <utility>
#include <vector>

#include <QAbstractButton>
#include <QAction>
#include <QApplication>
#include <QButtonGroup>
#include <QClipboard>
#include <QCloseEvent>
#include <QDir>
#include <QDragEnterEvent>
#include <QDropEvent>
#include <QFileInfo>
#include <QFileSystemWatcher>
#include <QFrame>
#include <QGuiApplication>
#include <QHBoxLayout>
#include <QLabel>
#include <QMenu>
#include <QMimeData>
#include <QMouseEvent>
#include <QPalette>
#include <QProgressBar>
#include <QPushButton>
#include <QScrollArea>
#include <QScrollBar>
#include <QStackedWidget>
#include <QStatusBar>
#include <QThreadPool>
#include <QToolBar>
#include <QToolButton>
#include <QVBoxLayout>
#include <QWidgetAction>

#include "citations.hpp"
#include "gui/actions/window_actions.hpp"
#include "demo.hpp"
#include "gui/dialogs/prompts.hpp"
#include "gui/dialogs/text_dialog.hpp"
#include "gui/forms/field_spec.hpp"
#include "gui/forms/field_widgets.hpp"
#include "gui/jobs/job_runner.hpp"
#include "gui/jobs/log_view.hpp"
#include "gui/pages/form_page.hpp"
#include "gui/pages/layout_page.hpp"
#include "gui/pages/structure_page.hpp"
#include "gui/problems/problems_panel.hpp"
#include "mdbridge/core/citations.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/problems.hpp"
#include "mdbridge/core/project_file.hpp"
#include "mdbridge/core/project_io.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/core/settings.hpp"
#include "mdbridge/render/panels/structure.hpp"
#include "mdbridge/render/panels/vulkan_stage.hpp"
#include "mdbridge/render/shading.hpp"
#include "mdbridge/render/vulkan_probe.hpp"
#include "mdbridge/version.hpp"
#include "log_capture.hpp"
#include "preview_widget.hpp"
#include "project_view.hpp"

namespace mdbridge::app {

namespace fs = std::filesystem;

namespace {

constexpr int kSideWidth = 520;               // MDATK's SIDE_WIDTH
constexpr int kRecoveryMs = 2 * 60 * 1000;    // D6: a recovery copy every 2 minutes
constexpr int kPreviewDelayMs = 120;          // MDATK's canvas debounce
constexpr int kReloadDelayMs = 600;           // input files typed: load the reference frame again
constexpr int kDiskCheckDelayMs = 300;        // after a file-watcher event
// a recovery copy written this recently may belong to another MDBridge window that is still open
constexpr auto kLiveRecovery = std::chrono::seconds(150);
constexpr const char* kProjectSuffix = ".mdatk.toml";

QString qs(const std::string& text) { return QString::fromUtf8(text.data(), static_cast<qsizetype>(text.size())); }
std::string utf8(const QString& text) { return forms::toStdString(text); }
QString display(const fs::path& path) { return qs(core::pathToUtf8(path)); }
fs::path toPath(const QString& file) { return core::pathFromUtf8(utf8(file)); }

// MDATK's name for a new file: the plan's project name unless it is the default.
std::string suggestedName(const core::Project& project) {
    return project.workflow.project_name != "project" ? project.workflow.project_name : std::string("MDATK");
}

// Python's str(datetime): the isoformat with a space.
std::string pythonTime(const core::Value& value) {
    if (const auto* t = value.get<core::DateTime>()) {
        std::string text = t->isoformat();
        if (text.size() > 10 && text[10] == 'T') text[10] = ' ';
        return text;
    }
    if (const auto* s = value.get<std::string>()) return *s;
    return "None";
}

std::string entryText(const core::Value& entry, std::string_view key, const char* missing) {
    const core::Value* v = entry.find(key);
    if (!v) return missing;
    if (const auto* s = v->get<std::string>()) return *s;
    return core::toJson(*v, -1);
}

// What the worker thread loads for a project.
struct Loaded {
    bool cancelled = false;
    std::shared_ptr<const render::MoleculeFrame> molecule;
    std::shared_ptr<const render::StructurePreview> preview;
    QString structureProblem;
    QString summary;
    QStringList notes;
};

QString describe(const std::exception& e) { return QString::fromUtf8(e.what()); }

Loaded loadStructure(const core::Project& project, const core::Reporter& reporter) {
    Loaded out;
    try {
        const ReferenceFrame ref = loadReferenceFrame(project, reporter);
        for (const std::string& note : ref.notes) out.notes.push_back(qs(note));
        try {
            out.preview = std::make_shared<render::StructurePreview>(structurePreview(ref.dataset));
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
    } catch (const std::exception& e) {
        out.structureProblem = describe(e);
    }
    return out;
}

}  // namespace

struct MainWindow::CitationsOf {
    core::Project project;
    ProjectCitations citations;
};

// ---- construction ------------------------------------------------------------------------------

MainWindow::MainWindow(render::Renderer renderer, QWidget* parent)
    : QMainWindow(parent), m_renderer(renderer), m_prompter(&defaultPrompter()) {
    setObjectName(QStringLiteral("mainWindow"));
    setAcceptDrops(true);
    resize(1320, 920);  // MDATK's default size
    m_context = std::make_unique<forms::FormContext>(m_document);

    m_log = new jobs::LogView;
    m_log->setObjectName(QStringLiteral("log"));
    m_log->textView()->setMinimumHeight(120);  // MDATK: min content height 120
    m_log->textView()->setMaximumHeight(200);  // the canvas keeps the room
    m_log->setExpanded(false);  // MDATK's log expander starts folded (Ctrl+L)
    m_jobs = new jobs::JobRunner(m_log->sink(), this);

    // WARNING and above of the libraries' messages go to the log too (MDATK's GuiLogHandler)
    {
        auto sink = m_log->sink(jobs::LogView::Kind::Plain);
        auto previous = std::make_shared<core::LogSink>();
        m_previousSink = core::setLogSink({});
        *previous = m_previousSink;
        core::setLogSink([sink, previous](core::LogLevel level, std::string_view message) {
            if (std::vector<std::string>* captured = ThreadLogCapture::current()) {  // the preview's drawing
                if (level >= core::LogLevel::Warning) {
                    captured->push_back(std::string(core::logLevelName(level)) + ": " + std::string(message));
                }
                return;
            }
            if (level >= core::LogLevel::Warning) {
                sink(qs(std::string(core::logLevelName(level))) + QStringLiteral(": ") + qs(std::string(message)));
            }
            if (*previous) {
                (*previous)(level, message);
            } else if (level >= core::LogLevel::Warning) {
                std::cerr << core::logLevelName(level) << ": " << message << '\n';
            }
        });
    }

    createActions();
    addToolBar(Qt::TopToolBarArea, static_cast<QToolBar*>(createHeader()));

    auto* central = new QWidget(this);
    central->setObjectName(QStringLiteral("central"));
    auto* outer = new QVBoxLayout(central);
    outer->setContentsMargins(0, 0, 0, 0);
    outer->setSpacing(0);
    outer->addWidget(createExternalBar());
    outer->addWidget(createPageBar());
    auto* body = new QHBoxLayout;
    body->setSpacing(0);
    body->addWidget(createSideBar());
    auto* line = new QFrame(central);
    line->setFrameShape(QFrame::VLine);
    line->setFrameShadow(QFrame::Sunken);
    body->addWidget(line);
    m_preview = new PreviewWidget(central);
    m_preview->setObjectName(QStringLiteral("preview"));
    body->addWidget(m_preview, 1);
    outer->addLayout(body, 1);
    m_log->setParent(central);
    outer->addWidget(m_log);
    outer->addWidget(createActionBar());
    setCentralWidget(central);
    addAction(m_log->foldAction());  // Ctrl+L (as RNAdo)

    connect(m_context.get(), &forms::FormContext::documentChanged, this, &MainWindow::onDocumentChanged);
    // a drawing's warnings are logged when they change after an edit; those of a project as opened
    // are not (MDATK draws it before its log handler is there)
    connect(m_preview, &PreviewWidget::drawn, this, [this] {
        const std::vector<std::string>& now = m_preview->warnings();
        if (m_warningsKnown && now != m_lastWarnings) {
            for (const std::string& warning : now) {
                if (std::find(m_lastWarnings.begin(), m_lastWarnings.end(), warning) == m_lastWarnings.end()) {
                    appendLog(qs(warning));
                }
            }
        }
        if (m_preview->current() || !m_preview->problem().isEmpty()) {
            m_lastWarnings = now;
            m_warningsKnown = true;
        }
    });
    connectCanvas();
    connect(m_jobs, &jobs::JobRunner::errorReported, this, [this](const QString& text, const QString& detail) {
        prompter().message(this, MessageKind::Error, text, detail);
    });
    connect(m_jobs, &jobs::JobRunner::statusChanged, this, [this](double fraction, const QString& text) {
        m_progress->setValue(static_cast<int>(fraction * 1000.0));
        m_progress->setFormat(text);
    });
    m_windowActions = new WindowActions(*this);  // the action bar's and the pages' jobs (W2-B)

    m_watcher = new QFileSystemWatcher(this);
    connect(m_watcher, &QFileSystemWatcher::fileChanged, this,
            [this] { QTimer::singleShot(kDiskCheckDelayMs, this, &MainWindow::checkDisk); });

    m_recoveryTimer.setInterval(kRecoveryMs);
    connect(&m_recoveryTimer, &QTimer::timeout, this, &MainWindow::writeRecoveryNow);
    m_recoveryTimer.start();
    m_previewTimer.setSingleShot(true);
    m_previewTimer.setInterval(kPreviewDelayMs);
    connect(&m_previewTimer, &QTimer::timeout, this, &MainWindow::updatePreview);
    m_reloadTimer.setSingleShot(true);
    m_reloadTimer.setInterval(kReloadDelayMs);
    connect(&m_reloadTimer, &QTimer::timeout, this, &MainWindow::startLoad);

    showPage(QStringLiteral("input"));
    statusBar()->showMessage(tr("Open a project (Ctrl+O), or start from these default settings."));
    afterReplace({});
}

MainWindow::~MainWindow() {
    cancelLoad();
    if (m_gpu) {  // the GPU check and the canvas's drawing end before the device goes
        m_preview->setFastStageRenderer(nullptr);
        QThreadPool::globalInstance()->waitForDone();
        m_gpu.reset();
    }
    m_previewTimer.stop();
    m_reloadTimer.stop();
    // the field widgets go before their FormContext (a member, destroyed before QWidget's children)
    delete takeCentralWidget();
    core::setLogSink(m_previousSink);
}

QAction* MainWindow::makeAction(const QString& name, const QString& text, const QKeySequence& key, const QString& tip) {
    auto* a = new QAction(text, this);
    a->setObjectName(QStringLiteral("action.") + name);
    if (!key.isEmpty()) a->setShortcut(key);
    if (!tip.isEmpty()) a->setToolTip(tip);
    a->setStatusTip(tip);
    addAction(a);  // the shortcut works wherever the action is shown (or not)
    m_actions[a->objectName()] = a;
    return a;
}

QAction* MainWindow::action(const QString& name) const {
    const auto it = m_actions.find(name.startsWith(QLatin1String("action.")) ? name : QStringLiteral("action.") + name);
    return it == m_actions.end() ? nullptr : it->second;
}

void MainWindow::createActions() {
    const QString later6 = tr("Comes with milestone 6 of MDBridge");
    const QString later7 = tr("Comes with milestone 7 of MDBridge");
    connect(makeAction(QStringLiteral("open"), tr("Open project"), QKeySequence::Open, tr("Open a project file (Ctrl+O)")),
            &QAction::triggered, this, &MainWindow::chooseAndOpen);
    connect(makeAction(QStringLiteral("browse"), tr("Browse…"), {}, tr("Choose a project file")), &QAction::triggered,
            this, &MainWindow::chooseAndOpen);
    connect(makeAction(QStringLiteral("demo_project"), tr("Demo project…"), {},
                       tr("Write the synthetic aptamer/target demo system into a folder and open it")),
            &QAction::triggered, this, &MainWindow::createDemo);
    connect(makeAction(QStringLiteral("save"), tr("Save"), QKeySequence(Qt::CTRL | Qt::Key_S), tr("Save the project (Ctrl+S)")),
            &QAction::triggered, this, [this] { save(); });
    connect(makeAction(QStringLiteral("save_as"), tr("Save as…"), QKeySequence(Qt::CTRL | Qt::SHIFT | Qt::Key_S),
                       tr("Save the project to a new file (Ctrl+Shift+S)")),
            &QAction::triggered, this, [this] { saveAs(); });
    connect(makeAction(QStringLiteral("new"), tr("New project (defaults)"), QKeySequence(Qt::CTRL | Qt::Key_N),
                       tr("Start from the default settings (Ctrl+N)")),
            &QAction::triggered, this, &MainWindow::newProject);
    connect(makeAction(QStringLiteral("close_project"), tr("Close project"), QKeySequence(Qt::CTRL | Qt::Key_W),
                       tr("Close the project file and return to the default settings, tied to no file (Ctrl+W)")),
            &QAction::triggered, this, &MainWindow::closeProject);
    connect(makeAction(QStringLiteral("create_demo"), tr("Create demo system…"), {},
                       tr("Write the synthetic aptamer/target demo system into a folder and open it")),
            &QAction::triggered, this, &MainWindow::createDemo);
    connect(makeAction(QStringLiteral("project_log"), tr("Project log…"), {}, tr("The steps run in this project")),
            &QAction::triggered, this, &MainWindow::showProjectLog);
    connect(makeAction(QStringLiteral("export_citations"), tr("Export citations (BibTeX)…"), {},
                       tr("Write the references of what the project uses as a BibTeX file")),
            &QAction::triggered, this, &MainWindow::exportCitations);
    connect(makeAction(QStringLiteral("citations"), tr("Citations…"), {},
                       tr("The references of the resources this project uses")),
            &QAction::triggered, this, &MainWindow::showCitations);
    connect(makeAction(QStringLiteral("graphics"), tr("Graphics…"), {}, tr("The 3D renderer and the Vulkan devices found")),
            &QAction::triggered, this, &MainWindow::showGraphics);
    connect(makeAction(QStringLiteral("about"), tr("About and citations"), {}, tr("About MDBridge, and what to cite")),
            &QAction::triggered, this, &MainWindow::showAbout);
    connect(makeAction(QStringLiteral("quit"), tr("Quit"), QKeySequence(Qt::CTRL | Qt::Key_Q), tr("Close the window (Ctrl+Q)")),
            &QAction::triggered, this, &QWidget::close);

    QAction* undo = makeAction(QStringLiteral("undo"), tr("Undo"), QKeySequence(Qt::CTRL | Qt::Key_Z), tr("Undo (Ctrl+Z)"));
    connect(undo, &QAction::triggered, this, [this] {
        m_context->commitPendingEdits();
        if (m_document.canUndo()) m_document.undo();
    });
    QAction* redo = makeAction(QStringLiteral("redo"), tr("Redo"), QKeySequence(Qt::CTRL | Qt::SHIFT | Qt::Key_Z),
                               tr("Redo (Ctrl+Shift+Z)"));
    redo->setShortcuts({QKeySequence(Qt::CTRL | Qt::SHIFT | Qt::Key_Z), QKeySequence(Qt::CTRL | Qt::Key_Y)});
    connect(redo, &QAction::triggered, this, [this] {
        m_context->commitPendingEdits();
        if (m_document.canRedo()) m_document.redo();
    });
    connect(makeAction(QStringLiteral("refresh"), tr("Refresh"), QKeySequence(Qt::Key_F5),
                       tr("Redraw the layout preview with the current settings (F5)")),
            &QAction::triggered, this, &MainWindow::refreshPreview);

    // the workflow programs (D1, D2): shells until their milestones
    for (const auto& [name, text, key, tip] :
         {std::tuple<QString, QString, QKeySequence, QString>{QStringLiteral("planner"), tr("Project Planner…"),
                                                             QKeySequence(Qt::CTRL | Qt::Key_P), later7},
          {QStringLiteral("wizard"), tr("Script wizard…"), QKeySequence(Qt::CTRL | Qt::SHIFT | Qt::Key_W), later7},
          {QStringLiteral("pdb_tools"), tr("PDB clean-up…"), QKeySequence(), later6},
          {QStringLiteral("builder"), tr("Structure builder…"), QKeySequence(), later6}}) {
        makeAction(name, text, key, tip)->setEnabled(false);
    }
    // the action bar's jobs (connected by WindowActions, W2-B)
    makeAction(QStringLiteral("render"), tr("Render"), QKeySequence(Qt::CTRL | Qt::Key_R),
               tr("Run cleaning, analyses and all enabled outputs (Ctrl+R)"));
    QAction* cancel = makeAction(QStringLiteral("cancel"), tr("Cancel"), {}, tr("Stop the running job"));
    connect(cancel, &QAction::triggered, this, [this] { m_jobs->cancel(); });
    cancel->setEnabled(false);
    m_jobs->addStoppable(cancel);
    makeAction(QStringLiteral("open_output"), tr("Open output folder"), {}, tr("Show the output folder in the file manager"));
    makeAction(QStringLiteral("open_video"), tr("Open video"), {}, tr("Play the most recent video"));
    m_jobs->addBusy(action(QStringLiteral("create_demo")));
    m_jobs->addBusy(action(QStringLiteral("demo_project")));

    for (std::size_t k = 0; k < std::size(core::kProblemPages); ++k) {
        const QString id = qs(std::string(core::kProblemPages[k].id));
        const QString title = qs(std::string(core::kProblemPages[k].title));
        QAction* page = makeAction(QStringLiteral("page.") + id, title,
                                   QKeySequence(Qt::CTRL | static_cast<Qt::Key>(Qt::Key_1 + static_cast<int>(k))),
                                   tr("%1 (Ctrl+%2)").arg(title).arg(k + 1));
        connect(page, &QAction::triggered, this, [this, id] { showPage(id); });
    }
}

QWidget* MainWindow::createHeader() {
    auto* bar = new QToolBar(tr("Project"), this);
    bar->setObjectName(QStringLiteral("header"));
    bar->setMovable(false);
    bar->setFloatable(false);
    bar->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);

    auto* open = new QToolButton(bar);
    open->setObjectName(QStringLiteral("header.open"));
    open->setDefaultAction(action(QStringLiteral("open")));
    m_openMenu = new QMenu(open);
    m_openMenu->setObjectName(QStringLiteral("menu.open"));
    m_openMenu->setToolTipsVisible(true);
    connect(m_openMenu, &QMenu::aboutToShow, this, &MainWindow::fillOpenMenu);
    open->setMenu(m_openMenu);
    open->setPopupMode(QToolButton::MenuButtonPopup);
    open->setToolTip(tr("Open a project file (Ctrl+O); the arrow: recent projects, the demo project, or browse"));
    bar->addWidget(open);
    bar->addAction(action(QStringLiteral("save")));
    bar->addAction(action(QStringLiteral("save_as")));

    auto* workflow = new QToolButton(bar);
    workflow->setObjectName(QStringLiteral("header.workflow"));
    workflow->setText(tr("Workflow"));
    workflow->setToolTip(tr("Project Planner, script wizard, PDB clean-up and structure builder (milestones 6 and 7)"));
    workflow->setPopupMode(QToolButton::InstantPopup);
    auto* workflowMenu = new QMenu(workflow);
    workflowMenu->setObjectName(QStringLiteral("menu.workflow"));
    workflowMenu->setToolTipsVisible(true);
    for (const char* name : {"planner", "wizard", "pdb_tools", "builder"}) {
        workflowMenu->addAction(action(QString::fromLatin1(name)));
    }
    workflow->setMenu(workflowMenu);
    bar->addWidget(workflow);
    bar->addSeparator();
    bar->addAction(action(QStringLiteral("undo")));
    bar->addAction(action(QStringLiteral("redo")));

    auto* spacer = new QWidget(bar);
    spacer->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Preferred);
    bar->addWidget(spacer);

    // the problems list (W1-D; MDATK 0.1.1's problems button and popover)
    auto* problemsButton = new QToolButton(bar);
    problemsButton->setObjectName(QStringLiteral("header.problems"));
    problemsButton->setText(tr("No problems"));
    problemsButton->setPopupMode(QToolButton::InstantPopup);
    problemsButton->setToolTip(tr("The problems of the settings; click one to go to its field"));
    auto* problemsMenu = new QMenu(problemsButton);
    problemsMenu->setObjectName(QStringLiteral("menu.problems"));
    auto* panel = new problems::ProblemsPanel(*m_context, problemsMenu);
    m_problems = panel;
    panel->setMinimumSize(460, 220);
    panel->setRevealFunction([this, problemsMenu](const QString& field) {
        problemsMenu->close();
        return revealField(field);
    });
    auto* holder = new QWidgetAction(problemsMenu);
    holder->setDefaultWidget(panel);
    problemsMenu->addAction(holder);
    problemsButton->setMenu(problemsMenu);
    connect(panel, &problems::ProblemsPanel::problemsShown, problemsButton, [problemsButton](int errors, int warnings) {
        const problems::ProblemCount count{errors, warnings};
        problemsButton->setText(count.total() == 0 ? tr("No problems") : problems::countText(count));
    });
    bar->addWidget(problemsButton);

    auto* menuButton = new QToolButton(bar);
    menuButton->setObjectName(QStringLiteral("header.menu"));
    menuButton->setText(tr("Menu"));
    menuButton->setIcon(QIcon::fromTheme(QStringLiteral("open-menu-symbolic")));
    menuButton->setAccessibleName(tr("Main menu"));
    menuButton->setPopupMode(QToolButton::InstantPopup);
    auto* menu = new QMenu(menuButton);
    menu->setObjectName(QStringLiteral("menu.main"));
    menu->setToolTipsVisible(true);
    for (const char* name : {"new", "close_project", "", "create_demo", "project_log", "", "export_citations",
                             "citations", "graphics", "about", "", "quit"}) {
        if (*name == '\0') menu->addSeparator();
        else menu->addAction(action(QString::fromLatin1(name)));
    }
    menuButton->setMenu(menu);
    bar->addWidget(menuButton);
    return bar;
}

QWidget* MainWindow::createPageBar() {
    auto* bar = new QWidget(this);
    bar->setObjectName(QStringLiteral("pagebar"));
    auto* layout = new QHBoxLayout(bar);
    layout->setContentsMargins(4, 4, 4, 4);
    layout->setSpacing(2);
    m_pageButtons = new QButtonGroup(bar);
    m_pageButtons->setExclusive(true);
    for (std::size_t k = 0; k < std::size(core::kProblemPages); ++k) {
        const QString id = qs(std::string(core::kProblemPages[k].id));
        const QString title = qs(std::string(core::kProblemPages[k].title));
        auto* button = new QToolButton(bar);
        button->setObjectName(QStringLiteral("page:") + id);
        button->setText(QString(title).replace(QLatin1Char('&'), QStringLiteral("&&")));  // no mnemonic
        button->setCheckable(true);
        button->setAutoRaise(true);
        button->setToolTip(tr("%1 (Ctrl+%2)").arg(title).arg(k + 1));
        button->setAccessibleName(title);
        m_pageButtons->addButton(button, static_cast<int>(k));
        connect(button, &QToolButton::clicked, this, [this, id] { showPage(id); });
        layout->addWidget(button);
        auto* badge = new problems::ProblemBadge(*m_context, id, bar);
        badge->setCursor(Qt::PointingHandCursor);
        badge->installEventFilter(this);
        layout->addWidget(badge);
        layout->addSpacing(4);
    }
    layout->addStretch(1);
    auto* hint = new QLabel(tr("Drag panels to swap them · drag the structure to rotate it · double-click a panel to "
                               "edit it"),
                            bar);
    hint->setObjectName(QStringLiteral("pagebar.hint"));
    hint->setEnabled(false);
    layout->addWidget(hint);
    auto* refresh = new QToolButton(bar);
    refresh->setDefaultAction(action(QStringLiteral("refresh")));
    layout->addWidget(refresh);
    return bar;
}

QWidget* MainWindow::createSideBar() {
    auto* side = new QWidget(this);
    side->setObjectName(QStringLiteral("sidebar"));
    side->setFixedWidth(kSideWidth);
    auto* layout = new QVBoxLayout(side);
    layout->setContentsMargins(0, 0, 0, 0);
    layout->setSpacing(0);
    m_sideTitle = new QLabel(side);
    m_sideTitle->setObjectName(QStringLiteral("side.title"));
    m_sideTitle->setContentsMargins(6, 6, 6, 6);
    layout->addWidget(m_sideTitle);
    auto* line = new QFrame(side);
    line->setFrameShape(QFrame::HLine);
    line->setFrameShadow(QFrame::Sunken);
    layout->addWidget(line);
    m_stack = new QStackedWidget(side);
    m_stack->setObjectName(QStringLiteral("side.pages"));
    for (QWidget* widget : pages::createPages(*m_context)) {
        m_pages.push_back(Page{pages::pageIdOf(widget), pages::pageTitleOf(widget), widget, pages::isPlaceholder(widget)});
    }
    for (const Page& page : m_pages) {
        page.widget->setParent(m_stack);  // a page scrolls itself (pages/form_page.hpp)
        m_stack->addWidget(page.widget);
    }
    layout->addWidget(m_stack, 1);
    return side;
}

QWidget* MainWindow::createActionBar() {
    auto* bar = new QWidget(this);
    bar->setObjectName(QStringLiteral("actionbar"));
    auto* layout = new QHBoxLayout(bar);
    layout->setContentsMargins(8, 8, 8, 8);
    layout->setSpacing(8);
    const auto button = [&](const char* name) {
        auto* b = new QToolButton(bar);
        b->setObjectName(QStringLiteral("actionbar.") + QString::fromLatin1(name));
        b->setDefaultAction(action(QString::fromLatin1(name)));
        b->setToolButtonStyle(Qt::ToolButtonTextOnly);
        return b;
    };
    QToolButton* render = button("render");
    QFont bold = render->font();
    bold.setBold(true);
    render->setFont(bold);  // MDATK's suggested action
    layout->addWidget(render);
    layout->addWidget(button("cancel"));
    m_progress = new QProgressBar(bar);
    m_progress->setObjectName(QStringLiteral("progress"));
    m_progress->setRange(0, 1000);
    m_progress->setValue(0);
    m_progress->setTextVisible(true);
    m_progress->setFormat(m_jobs->statusText());
    layout->addWidget(m_progress, 1);
    layout->addWidget(button("open_output"));
    layout->addWidget(button("open_video"));
    return bar;
}

QWidget* MainWindow::createExternalBar() {
    m_externalBar = new QFrame(this);
    m_externalBar->setObjectName(QStringLiteral("external"));
    m_externalBar->setFrameShape(QFrame::StyledPanel);
    m_externalBar->setAutoFillBackground(true);
    QPalette palette = m_externalBar->palette();
    palette.setColor(QPalette::Window, palette.color(QPalette::ToolTipBase));
    palette.setColor(QPalette::WindowText, palette.color(QPalette::ToolTipText));
    m_externalBar->setPalette(palette);
    auto* layout = new QHBoxLayout(m_externalBar);
    layout->setContentsMargins(8, 4, 8, 4);
    m_externalText = new QLabel(m_externalBar);
    m_externalText->setObjectName(QStringLiteral("external.text"));
    m_externalText->setWordWrap(true);
    layout->addWidget(m_externalText, 1);
    auto* reload = new QPushButton(tr("Reload"), m_externalBar);
    reload->setObjectName(QStringLiteral("external.reload"));
    reload->setToolTip(tr("Read the file again (your unsaved changes are lost)"));
    connect(reload, &QPushButton::clicked, this, &MainWindow::reloadFromDisk);
    layout->addWidget(reload);
    auto* keep = new QPushButton(tr("Keep mine"), m_externalBar);
    keep->setObjectName(QStringLiteral("external.keep"));
    keep->setToolTip(tr("Go on with the settings shown; saving asks before it replaces the file"));
    connect(keep, &QPushButton::clicked, m_externalBar, &QWidget::hide);
    layout->addWidget(keep);
    m_externalBar->hide();
    return m_externalBar;
}

void MainWindow::fillOpenMenu() {
    m_openMenu->clear();
    const std::vector<fs::path> recent = core::recentProjects();
    if (!recent.empty()) {
        QAction* head = m_openMenu->addAction(tr("Recent projects"));
        head->setEnabled(false);
        int n = 0;
        for (const fs::path& path : recent) {
            QAction* item = m_openMenu->addAction(display(path.filename()));
            item->setObjectName(QStringLiteral("recent:%1").arg(n++));
            item->setToolTip(display(path));
            item->setData(display(path));
            connect(item, &QAction::triggered, this, [this, path] { openPath(path); });
        }
        m_openMenu->addSeparator();
    }
    m_openMenu->addAction(action(QStringLiteral("demo_project")));
    m_openMenu->addAction(action(QStringLiteral("browse")));
}

void MainWindow::setPrompter(Prompter* prompter) { m_prompter = prompter ? prompter : &defaultPrompter(); }

// ---- pages ---------------------------------------------------------------------------------------

QAbstractButton* MainWindow::pageButton(const QString& id) const {
    for (QAbstractButton* b : m_pageButtons->buttons()) {
        if (b->objectName() == QStringLiteral("page:") + id) return b;
    }
    return nullptr;
}

QString MainWindow::currentPage() const {
    const int index = m_stack->currentIndex();
    return index >= 0 && static_cast<std::size_t>(index) < m_pages.size() ? m_pages[static_cast<std::size_t>(index)].id
                                                                           : QString();
}

bool MainWindow::showPage(const QString& id) {
    for (std::size_t k = 0; k < m_pages.size(); ++k) {
        if (m_pages[k].id != id) continue;
        const bool changed = m_stack->currentIndex() != static_cast<int>(k);
        m_stack->setCurrentIndex(static_cast<int>(k));
        m_sideTitle->setText(QStringLiteral("<b>%1</b>").arg(m_pages[k].title.toHtmlEscaped()));
        if (QAbstractButton* b = pageButton(id); b && !b->isChecked()) b->setChecked(true);
        if (changed) emit pageChanged(id);
        return true;
    }
    return false;
}

bool MainWindow::openSection(const QString& page, const QString& section) {
    if (!showPage(page)) return false;
    QWidget* shown = m_stack->currentWidget();
    auto* area = shown->findChild<QScrollArea*>();
    if (section.isEmpty()) {
        if (area) area->verticalScrollBar()->setValue(0);
        return true;
    }
    auto* target = shown->findChild<QWidget*>(QStringLiteral("section:%1.%2").arg(page, section));
    if (!target) return false;
    if (auto* s = qobject_cast<pages::Section*>(target)) s->setExpanded(true);
    // MDATK's open_section: the section's top 6 px below the top, once the page's layout is done
    const auto scroll = [area, target] {
        if (!area || !area->widget()) return;
        const QPoint at = target->mapTo(area->widget(), QPoint(0, 0));
        area->verticalScrollBar()->setValue(std::max(0, at.y() - 6));
    };
    scroll();
    QTimer::singleShot(0, target, scroll);
    return true;
}

bool MainWindow::revealField(const QString& path) {
    const QString page = qs(std::string(core::problemPage(utf8(path))));
    if (page.isEmpty() || !showPage(page)) return false;
    for (const Page& entry : m_pages) {
        if (entry.id == page) return pages::revealFieldOn(entry.widget, path);
    }
    return false;
}

bool MainWindow::eventFilter(QObject* watched, QEvent* event) {
    if (event->type() == QEvent::MouseButtonRelease) {
        if (auto* badge = qobject_cast<problems::ProblemBadge*>(watched)) {
            // a badge click goes to the page's next problem (MDATK 0.1.1)
            showPage(badge->page());
            if (const core::Problem* problem = m_problems->nextOnPage(badge->page())) {
                revealField(qs(problem->field));
            }
            return true;
        }
    }
    return QMainWindow::eventFilter(watched, event);
}

// ---- log and messages ------------------------------------------------------------------------------

void MainWindow::appendLog(const QString& text) { m_log->feed(text); }

void MainWindow::showError(const QString& text, const QString& detail) {
    appendLog(QStringLiteral("ERROR: ") + text + (detail.isEmpty() ? QString() : QStringLiteral(": ") + detail));
    prompter().message(this, MessageKind::Error, text, detail);
}

void MainWindow::showInfo(const QString& text, const QString& detail) {
    appendLog(text + (detail.isEmpty() ? QString() : QStringLiteral(": ") + detail));
    prompter().message(this, MessageKind::Info, text, detail);
}

QString MainWindow::ask(const QString& kind, const QString& path, const QString& title, const QString& text,
                        const QString& detail, const std::vector<std::pair<QString, QString>>& choices,
                        const QString& defaultId) {
    Question q;
    q.kind = kind;
    q.path = path;
    q.title = title;
    q.text = text;
    q.detail = detail;
    q.choices = choices;
    q.defaultId = defaultId;
    return prompter().ask(this, q);
}

// ---- the document ----------------------------------------------------------------------------------

fs::path MainWindow::projectFile() const {
    return m_document.isParamsLog() ? fs::path() : m_document.path();
}

void MainWindow::updateTitle() {
    QString name;
    if (m_document.isParamsLog()) {
        name = tr("%1 (settings of a run, unsaved)").arg(display(m_document.path().filename()));
    } else if (m_document.hasPath()) {
        name = display(m_document.path());
    } else {
        name = tr("unsaved project");
    }
    setWindowTitle(tr("%1[*] - MDBridge").arg(name));
    setWindowModified(m_document.isModified());
}

void MainWindow::updateUndo() {
    QAction* undo = action(QStringLiteral("undo"));
    QAction* redo = action(QStringLiteral("redo"));
    undo->setEnabled(m_document.canUndo());
    redo->setEnabled(m_document.canRedo());
    undo->setToolTip(m_document.canUndo() ? tr("Undo %1 (Ctrl+Z)").arg(qs(m_document.undoLabel())) : tr("Nothing to undo"));
    redo->setToolTip(m_document.canRedo() ? tr("Redo %1 (Ctrl+Shift+Z)").arg(qs(m_document.redoLabel()))
                                          : tr("Nothing to redo"));
}

void MainWindow::onDocumentChanged(const core::DocumentChange& change) {
    updateTitle();
    updateUndo();
    using Kind = core::DocumentChange::Kind;
    if (change.kind == Kind::Replaced) {
        afterReplace(m_document.path());
        return;
    }
    if (change.kind == Kind::Saved) {
        m_externalNoticed = false;
        m_externalBar->hide();
        watchFile();
        if (!change.fields.empty()) schedulePreview();  // Save as rebased paths
        return;
    }
    if (m_canvasEdit) {
        // the canvas shows it already (AP-20: no rebuild of the figure for a drag)
    } else if (!change.fields.empty() && std::all_of(change.fields.begin(), change.fields.end(), [](const std::string& f) {
                   return f.rfind("visual.rotation", 0) == 0;
               })) {
        // the sliders (or undo): the structure panel turns, the figure is not rebuilt (AP-20)
        const std::vector<double>& r = m_document.project().visual.rotation;
        if (r.size() == 3) m_preview->setRotation({r[0], r[1], r[2]});
        else schedulePreview();
    } else {
        schedulePreview();
    }
    for (const char* part : {"inputs", "cleaning", "groups"}) {
        if (change.affects(part)) {
            m_reloadTimer.start();
            break;
        }
    }
}

void MainWindow::afterReplace(const fs::path& /*opened*/) {
    render::resetWarnOnce();  // a project shown anew warns again about what its drawing finds
    m_warningsKnown = false;
    m_lastWarnings.clear();
    m_externalNoticed = false;
    m_externalBar->hide();
    watchFile();
    updateTitle();
    updateUndo();
    startLoad();
    updatePreview();
}

void MainWindow::watchFile() {
    if (!m_watcher->files().isEmpty()) m_watcher->removePaths(m_watcher->files());
    const fs::path file = projectFile();
    std::error_code ec;
    if (!file.empty() && fs::exists(file, ec)) m_watcher->addPath(display(file));
}

void MainWindow::checkDisk() {
    const fs::path file = projectFile();
    if (file.empty()) return;
    std::error_code ec;
    if (fs::exists(file, ec) && !m_watcher->files().contains(display(file))) {
        m_watcher->addPath(display(file));  // an atomic save by another program replaced the file watched
    }
    if (m_externalNoticed || !m_document.changedOnDisk()) return;
    m_externalNoticed = true;
    const QString what = fs::exists(file, ec) ? tr("%1 was changed by another program.").arg(display(file))
                                              : tr("%1 was deleted or moved by another program.").arg(display(file));
    m_externalText->setText(what + QLatin1Char(' ') + tr("Reload it, or keep the settings shown (Save then asks before "
                                                         "it replaces the file)."));
    m_externalBar->show();
    appendLog(what);
}

bool MainWindow::externalChangeShown() const { return m_externalBar->isVisible(); }

QString MainWindow::formErrorText() {
    std::vector<const core::Problem*> errors;
    for (const auto& [path, problem] : m_context->formErrors()) errors.push_back(&problem);
    const auto pageIndex = [](const core::Problem* p) {
        const std::string_view page = p->page();
        for (std::size_t k = 0; k < std::size(core::kProblemPages); ++k) {
            if (core::kProblemPages[k].id == page) return k;
        }
        return std::size(core::kProblemPages);
    };
    std::stable_sort(errors.begin(), errors.end(),
                     [&](const core::Problem* a, const core::Problem* b) { return pageIndex(a) < pageIndex(b); });
    QStringList lines;
    for (const core::Problem* p : errors) lines.push_back(qs(p->text()));
    return lines.join(QLatin1Char('\n'));
}

bool MainWindow::confirmUnsaved(const QString& doing) {
    m_context->commitPendingEdits();
    if (!m_document.isModified()) return true;
    const fs::path file = projectFile();
    const QString text = file.empty()
                             ? tr("Save the changes to this project before you %1?").arg(doing)
                             : tr("Save the changes to %1 before you %2?").arg(display(file.filename()), doing);
    const QString answer =
        ask(QStringLiteral("unsaved"), display(file), tr("Unsaved changes"), text,
            tr("Your changes are lost if you do not save them."),
            {{QStringLiteral("save"), tr("Save")}, {QStringLiteral("discard"), tr("Discard")}, {QStringLiteral("cancel"), tr("Cancel")}});
    if (answer == QLatin1String("save")) return save();
    if (answer == QLatin1String("discard")) {
        dropRecovery();
        return true;
    }
    return false;
}

// ---- start -------------------------------------------------------------------------------------------

bool MainWindow::startWith(const fs::path& path) {
    std::vector<std::string> warnings;
    try {
        m_document.open(path, &warnings);
    } catch (const std::exception& e) {
        const QString reason = describe(e);
        prompter().printError(QStringLiteral("warning: could not load %1: %2").arg(display(path), reason));
        // AP-12 (MDATK 0.1.0): started from a menu no terminal shows that warning
        appendLog(QStringLiteral("ERROR: Could not open %1: %2\nThese are the default settings instead (unsaved).")
                      .arg(display(path), reason));
        m_log->setExpanded(true);
        statusBar()->showMessage(tr("Could not open %1; these are the default settings (unsaved)").arg(display(path)));
        return false;
    }
    // on the terminal, as MDATK's start (its log handler is not there yet when the file is read)
    for (const std::string& warning : warnings) std::cerr << "WARNING: " << warning << '\n';
    statusBar()->showMessage(tr("Opened %1").arg(display(path)), 5000);
    return true;
}

void MainWindow::offerRecoveries() {
    const auto now = fs::file_time_type::clock::now();
    for (const core::RecoveryFile& r : core::findRecoveries()) {
        if (r.file == m_recoveryFile) continue;
        std::error_code ec;
        const auto written = fs::last_write_time(r.file, ec);
        if (!ec && now - written < kLiveRecovery) continue;  // probably another MDBridge window's, still open
        const QString what = r.original.empty() ? tr("a new project") : display(r.original);
        const QString when = qs(pythonTime(core::Value(r.saved)));
        const QString answer = ask(QStringLiteral("recover"), display(r.original), tr("Recover unsaved changes"),
                                   tr("MDBridge kept a copy of the unsaved changes of %1.").arg(what),
                                   tr("The copy was written at %1 by a window that was not closed normally.").arg(when),
                                   {{QStringLiteral("recover"), tr("Recover")},
                                    {QStringLiteral("discard"), tr("Discard the copy")},
                                    {QStringLiteral("later"), tr("Not now")}});
        if (answer == QLatin1String("discard")) {
            core::removeRecovery(r.file);
            continue;
        }
        if (answer != QLatin1String("recover")) continue;
        if (!confirmUnsaved(tr("recover the copy"))) return;
        try {
            std::vector<std::string> warnings;
            core::openRecovery(m_document, r, &warnings);
            for (const std::string& warning : warnings) appendLog(QStringLiteral("WARNING: ") + qs(warning));
        } catch (const std::exception& e) {
            showError(tr("Could not recover %1").arg(display(r.file)), describe(e));
            continue;
        }
        core::removeRecovery(r.file);
        m_recoveryFile.clear();
        writeRecoveryNow();
        appendLog(tr("Recovered the unsaved changes of %1 (copy of %2).").arg(what, when));
        return;
    }
}

void MainWindow::writeRecoveryNow() {
    if (!m_document.isModified()) return;
    try {
        const fs::path file = core::writeRecovery(m_document);
        if (!m_recoveryFile.empty() && m_recoveryFile != file) core::removeRecovery(m_recoveryFile);
        m_recoveryFile = file;
    } catch (const std::exception& e) {
        core::logWarning(std::string("could not write a recovery copy: ") + e.what());
    }
}

void MainWindow::dropRecovery() {
    if (!m_recoveryFile.empty()) {
        core::removeRecovery(m_recoveryFile);
        m_recoveryFile.clear();
    }
}

// ---- new, open, close -------------------------------------------------------------------------------

void MainWindow::newProject() {
    if (!confirmUnsaved(tr("start a new project"))) return;
    dropRecovery();
    m_document.newProject();
    appendLog(QStringLiteral("New project with default settings (unsaved)."));
}

void MainWindow::closeProject() {
    const fs::path file = projectFile();
    if (!confirmUnsaved(tr("close it"))) return;
    dropRecovery();
    m_document.newProject();
    appendLog(file.empty() ? QStringLiteral("Closed the project.")
                           : QStringLiteral("Closed %1: default settings (unsaved).").arg(display(file.filename())));
}

void MainWindow::chooseAndOpen() {
    if (!confirmUnsaved(tr("open another project"))) return;
    forms::FileDialogRequest request;
    request.mode = forms::PathMode::OpenFile;
    request.title = tr("Open project");
    request.start = qs(m_document.folder());
    request.filters = {forms::nameFilter(forms::projectFiles()), forms::nameFilter(forms::runLogFiles()),
                       tr("All files (*)")};
    const QStringList picked = forms::runFileDialog(this, request);
    if (!picked.isEmpty()) openPath(toPath(picked.front()), false);
}

bool MainWindow::openPath(const fs::path& path, bool askUnsaved) {
    if (askUnsaved && !confirmUnsaved(tr("open another project"))) return false;
    m_context->commitPendingEdits();
    const fs::path previousRecovery = m_recoveryFile;
    std::vector<std::string> warnings;
    try {
        m_document.open(path, &warnings);
    } catch (const std::exception& e) {
        showError(tr("Could not open %1").arg(display(path)), describe(e));
        return false;
    }
    if (!previousRecovery.empty()) {
        core::removeRecovery(previousRecovery);
        if (m_recoveryFile == previousRecovery) m_recoveryFile.clear();
    }
    for (const std::string& warning : warnings) appendLog(QStringLiteral("WARNING: ") + qs(warning));
    if (m_document.isParamsLog()) {
        appendLog(QStringLiteral("Loaded the settings of the run logged in %1 (save them as a project to keep them)")
                      .arg(display(path)));
    } else {
        core::addRecent(path);
        appendLog(QStringLiteral("Opened project %1").arg(display(path)));
    }
    statusBar()->showMessage(tr("Opened %1").arg(display(path)), 5000);
    return true;
}

void MainWindow::reloadFromDisk() {
    const fs::path file = projectFile();
    if (file.empty()) return;
    m_context->commitPendingEdits();
    if (m_document.isModified() &&
        ask(QStringLiteral("reload"), display(file), tr("Reload the project"),
            tr("Reload %1 and discard your unsaved changes?").arg(display(file.filename())),
            tr("The file on disk was changed by another program."),
            {{QStringLiteral("reload"), tr("Reload")}, {QStringLiteral("cancel"), tr("Cancel")}},
            QStringLiteral("cancel")) != QLatin1String("reload")) {
        return;
    }
    dropRecovery();
    try {
        std::vector<std::string> warnings;
        m_document.open(file, &warnings);
        for (const std::string& warning : warnings) appendLog(QStringLiteral("WARNING: ") + qs(warning));
    } catch (const std::exception& e) {
        showError(tr("Could not open %1").arg(display(file)), describe(e));
        return;
    }
    appendLog(tr("Reloaded project %1").arg(display(file)));
}

void MainWindow::closeEvent(QCloseEvent* event) {
    if (!confirmUnsaved(tr("close the window"))) {
        event->ignore();
        return;
    }
    if (!m_jobs->confirmClose(this)) {  // AP-14: a job running is cancelled and waited for
        event->ignore();
        return;
    }
    dropRecovery();
    cancelLoad();
    event->accept();
}

void MainWindow::changeEvent(QEvent* event) {
    QMainWindow::changeEvent(event);
    if (event->type() == QEvent::ActivationChange && isActiveWindow()) {
        QTimer::singleShot(0, this, &MainWindow::checkDisk);
    }
}

// ---- save --------------------------------------------------------------------------------------------

fs::path MainWindow::askSavePath() {
    forms::FileDialogRequest request;
    request.mode = forms::PathMode::SaveFile;
    request.title = tr("Save project");
    request.start = qs(core::pathToUtf8(core::pathFromUtf8(m_document.folder()) /
                                        core::pathFromUtf8(suggestedName(m_document.project()) + kProjectSuffix)));
    request.filters = {forms::nameFilter(forms::projectFiles())};
    request.confirmOverwrite = false;  // asked below, about the final name (AP-11)
    const QStringList picked = forms::runFileDialog(this, request);
    if (picked.isEmpty() || picked.front().isEmpty()) return {};
    QString file = picked.front();
    if (!file.endsWith(QLatin1String(".toml"))) file += QLatin1String(kProjectSuffix);  // before the question
    const fs::path path = toPath(file);
    std::error_code ec;
    if (fs::exists(path, ec) &&
        ask(QStringLiteral("overwrite"), display(path), tr("Replace the file"),
            tr("Replace %1?").arg(display(path.filename())), tr("%1 exists already. Saving replaces it.").arg(display(path)),
            {{QStringLiteral("replace"), tr("Replace")}, {QStringLiteral("cancel"), tr("Cancel")}},
            QStringLiteral("cancel")) != QLatin1String("replace")) {
        return {};
    }
    return path;
}

bool MainWindow::save() {
    m_context->commitPendingEdits();
    if (const QString errors = formErrorText(); !errors.isEmpty()) {
        showError(tr("Some settings are invalid"), errors);
        return false;
    }
    if (projectFile().empty()) {
        const fs::path path = askSavePath();
        return !path.empty() && saveTo(path, true);
    }
    const fs::path file = projectFile();
    std::error_code ec;
    if (m_document.changedOnDisk() && fs::exists(file, ec) &&
        ask(QStringLiteral("external"), display(file), tr("The file was changed"),
            tr("Another program changed %1 since it was opened here. Replace its version with yours?")
                .arg(display(file.filename())),
            display(file), {{QStringLiteral("replace"), tr("Replace it")}, {QStringLiteral("cancel"), tr("Cancel")}},
            QStringLiteral("cancel")) != QLatin1String("replace")) {
        return false;
    }
    return saveTo(file, false);
}

bool MainWindow::saveAs() {
    m_context->commitPendingEdits();
    if (const QString errors = formErrorText(); !errors.isEmpty()) {
        showError(tr("Some settings are invalid"), errors);
        return false;
    }
    const fs::path path = askSavePath();
    return !path.empty() && saveTo(path, true);
}

bool MainWindow::saveTo(const fs::path& target, bool asNew) {
    try {
        if (asNew) {
            std::error_code ec;
            if (target.has_parent_path()) fs::create_directories(target.parent_path(), ec);
            m_document.saveAs(target);
        } else {
            m_document.save();
        }
    } catch (const std::exception& e) {  // core::ProjectError, std::filesystem::filesystem_error
        showError(tr("Could not save the project"), describe(e));
        return false;
    }
    dropRecovery();
    core::addRecent(target);
    appendLog(QStringLiteral("Saved project %1").arg(display(target)));
    statusBar()->showMessage(tr("Saved %1").arg(display(target)), 5000);
    updateTitle();
    return true;
}

// ---- dialogs -------------------------------------------------------------------------------------

TextDialog* MainWindow::openTextDialog(const QString& title, const QString& heading, const QString& text, bool monospace) {
    if (m_dialog) m_dialog->close();
    auto* dialog = new TextDialog(this, title, heading, text, monospace);
    m_dialog = dialog;
    return dialog;
}

void MainWindow::showProjectLog() {
    m_context->commitPendingEdits();
    const fs::path file = projectFile();
    if (file.empty()) {
        showInfo(QStringLiteral("The project has no file yet"),
                 QStringLiteral("Save the project; the steps you run are then logged in it."));
        return;
    }
    std::string warning;
    const core::Array entries = core::readLog(file, &warning);
    QStringList lines;
    for (const core::Value& e : entries) {
        const core::Value* time = e.find("time");
        lines.push_back(qs((time ? pythonTime(*time) : std::string("None")) + "  " + entryText(e, "step", "?") + "  [" +
                           entryText(e, "status", "?") + "]  " + entryText(e, "title", "")));
        const std::string summary = entryText(e, "summary", "");
        if (!summary.empty()) lines.push_back(QStringLiteral("    ") + qs(summary));
    }
    if (!warning.empty()) lines.push_back(QStringLiteral("WARNING: ") + qs(warning));
    TextDialog* dialog = openTextDialog(tr("Project log"), QString(),
                                        lines.isEmpty() ? QStringLiteral("Nothing has been run in this project yet.")
                                                        : lines.join(QLatin1Char('\n')),
                                        true);
    dialog->open();
}

namespace {

QString aboutHeading() {
    return MainWindow::tr("<h3>MDBridge %1</h3>"
                          "<p>Offline analysis and multi-panel video rendering of molecular dynamics trajectories of RNA "
                          "systems; the successor of MDATK. Opens MDATK %2 project files.</p>"
                          "<p>Copyright © 2026 the MDBridge author. All rights reserved; the licence follows a legal "
                          "review.</p>"
                          "<p>Built with Qt %3 (LGPL v3), toml++ (MIT), CLI11 (BSD-3-Clause) and the Vulkan headers "
                          "(Apache 2.0 / MIT); see the licenses folder of this program. "
                          "<a href=\"%4\">%4</a></p>")
        .arg(QString::fromLatin1(kVersion), QString::fromLatin1(kMdatkReference), QString::fromLatin1(qVersion()),
             QString::fromLatin1(kHomepage));
}

}  // namespace

void MainWindow::showCitations() {
    m_context->commitPendingEdits();
    auto of = std::make_shared<CitationsOf>();
    try {
        of->project = m_document.resolvedProject();
    } catch (const std::exception&) {
        of->project = m_document.project();
    }
    of->citations = projectCitations(of->project);
    const QString text = qs(citationText(of->citations));
    TextDialog* dialog = openTextDialog(tr("Citations"), tr("<p>What to cite for the resources this project uses "
                                                           "(and, after a run, the ones its results record).</p>"),
                                        text, false);
    dialog->addButton(QStringLiteral("copy"), tr("Copy"), [text] { QGuiApplication::clipboard()->setText(text); });
    dialog->addButton(QStringLiteral("copy_json"), tr("Copy as JSON"), [of] {
        QGuiApplication::clipboard()->setText(qs(core::toJson(citationJson(of->citations), 2)));
    });
    dialog->addButton(QStringLiteral("bibtex"), tr("Export BibTeX…"), [this] { exportCitations(); });
    dialog->open();
}

void MainWindow::showAbout() {
    m_context->commitPendingEdits();
    core::Project project;
    try {
        project = m_document.resolvedProject();
    } catch (const std::exception&) {
        project = m_document.project();
    }
    const QString text = qs(citationText(projectCitations(project)));
    TextDialog* dialog = openTextDialog(tr("About MDBridge"), aboutHeading(), text, false);
    dialog->addButton(QStringLiteral("copy"), tr("Copy citations"), [text] { QGuiApplication::clipboard()->setText(text); });
    dialog->addButton(QStringLiteral("bibtex"), tr("Export BibTeX…"), [this] { exportCitations(); });
    dialog->addButton(QStringLiteral("graphics"), tr("Graphics…"), [this] { showGraphics(); });
    dialog->open();
}

void MainWindow::exportCitations() {
    m_context->commitPendingEdits();
    core::Project project;
    try {
        project = m_document.resolvedProject();
    } catch (const std::exception&) {
        project = m_document.project();
    }
    const ProjectCitations citations = projectCitations(project);
    if (citations.keys.empty()) {
        showInfo(QStringLiteral("Nothing to cite yet"), QStringLiteral("The project uses no external resource with a reference."));
        return;
    }
    std::string base = m_document.hasPath() ? m_document.folder() : project.output.out_dir;
    std::error_code ec;
    if (base.empty() || !fs::is_directory(core::pathFromUtf8(base), ec)) base = m_document.folder();
    const std::string name = project.workflow.project_name.empty() ? "MDATK" : project.workflow.project_name;
    forms::FileDialogRequest request;
    request.mode = forms::PathMode::SaveFile;
    request.title = tr("Export citations (BibTeX)");
    request.start = qs(core::pathToUtf8(core::pathFromUtf8(base) / core::pathFromUtf8(name + "_citations.bib")));
    request.filters = {tr("BibTeX files (*.bib)")};
    const QStringList picked = forms::runFileDialog(this, request);
    if (picked.isEmpty() || picked.front().isEmpty()) return;
    const fs::path path = toPath(picked.front());
    {
        std::ofstream out(path, std::ios::binary);
        out << core::citations::bibtex(citations.keys);
        if (!out) {
            showError(QStringLiteral("Could not write the BibTeX file"),
                      tr("%1: %2").arg(display(path), qs(std::generic_category().message(errno))));
            return;
        }
    }
    appendLog(QStringLiteral("Wrote %1 references to %2").arg(citations.keys.size()).arg(display(path)));
}

void MainWindow::showGraphics() {
    const render::VulkanReport report = render::probeVulkan();
    QString text = tr("3D renderer chosen: %1 (the structure panel draws with it while you rotate it)\n")
                       .arg(QString::fromLatin1(render::rendererName(m_renderer)));
    text += tr("Vulkan loader: %1\n")
                .arg(report.loaderFound ? tr("found, Vulkan %1").arg(qs(render::vulkanVersionText(report.loaderVersion)))
                                        : tr("not found"));
    for (const render::VulkanDevice& d : report.devices) {
        text += tr("  %1 (%2, Vulkan %3)%4\n")
                    .arg(qs(d.name), qs(d.type), qs(render::vulkanVersionText(d.apiVersion)),
                         d.meetsMinimum() ? QString() : tr(": needs Vulkan %1.%2").arg(render::kVulkanMajor).arg(render::kVulkanMinor));
    }
    if (!report.problem.empty()) text += qs(report.problem) + QLatin1Char('\n');
    openTextDialog(tr("Graphics"), QString(), text, true)->open();
}

void MainWindow::createDemo() {
    forms::FileDialogRequest request;
    request.mode = forms::PathMode::Folder;
    request.title = tr("Folder for the synthetic demo system");
    request.start = qs(m_document.folder());
    const QStringList picked = forms::runFileDialog(this, request);
    if (picked.isEmpty() || picked.front().isEmpty()) return;
    const std::string folder = utf8(picked.front());
    m_jobs->start(
        tr("Creating demo"),
        [folder](const core::Reporter& reporter) {
            reporter.log("Writing synthetic RNA aptamer/target system to " + folder + " ...");
            return buildDemo(folder, DemoOptions{}, reporter);
        },
        [this](DemoFiles files) {
            appendLog(QStringLiteral("Demo system ready: %1").arg(qs(files.project)));
            openPath(core::pathFromUtf8(files.project), true);
        });
}

// ---- drag and drop -------------------------------------------------------------------------------

void MainWindow::dragEnterEvent(QDragEnterEvent* event) {
    const QList<QUrl> urls = event->mimeData()->urls();
    if (urls.size() == 1 && urls.front().isLocalFile()) event->acceptProposedAction();
}

void MainWindow::dropEvent(QDropEvent* event) {
    const QList<QUrl> urls = event->mimeData()->urls();
    if (urls.size() == 1 && urls.front().isLocalFile()) {
        event->acceptProposedAction();
        const fs::path path = toPath(urls.front().toLocalFile());
        // after the drop has finished: the question must not run inside the drag's event
        QTimer::singleShot(0, this, [this, path] { openPath(path); });
    }
}

// ---- the reference frame and the preview ---------------------------------------------------------------

void MainWindow::refreshPreview() {
    m_context->commitPendingEdits();
    updatePreview();
}

void MainWindow::schedulePreview() { m_previewTimer.start(); }

bool MainWindow::previewSettled() const { return !m_previewTimer.isActive() && m_preview->settled(); }

void MainWindow::updatePreview() {
    m_previewTimer.stop();
    core::Project project;
    try {
        project = m_document.resolvedProject();
    } catch (const std::exception& e) {
        m_preview->setMessage(tr("Cannot draw this layout: %1").arg(describe(e)));
        return;
    }
    auto content = std::make_shared<PreviewContent>();
    try {
        content->example = render::exampleDataset(project, m_structurePreview.get());
    } catch (const std::exception& e) {
        m_preview->setMessage(tr("Cannot draw this layout: %1").arg(describe(e)));
        return;
    }
    content->watermark = previewWatermark(m_structurePreview != nullptr);
    content->project = std::move(project);
    m_preview->setContent(std::move(content));
}

void MainWindow::cancelLoad() {
    if (m_cancel) {
        m_cancel->store(true);
        m_cancel.reset();
    }
    ++m_loadGeneration;
    m_loading = false;
}

void MainWindow::startLoad() {
    cancelLoad();
    m_reloadTimer.stop();
    const std::uint64_t generation = m_loadGeneration;
    m_loadProblem.clear();
    m_structureSummary.clear();
    m_molecule.reset();
    m_structurePreview.reset();
    core::Project project;
    try {
        project = m_document.resolvedProject();
    } catch (const std::exception& e) {
        m_loadProblem = describe(e);
        emit projectLoaded();
        return;
    }
    m_cancel = std::make_shared<std::atomic<bool>>(false);
    m_loading = true;
    statusBar()->showMessage(tr("Loading the reference frame ..."));
    const auto post = m_background.poster();
    const core::Reporter reporter(
        {},
        [this, post, generation](double fraction, const std::string& message) {
            post([this, generation, fraction, text = qs(message)] {
                if (generation == m_loadGeneration && !text.isEmpty()) {
                    statusBar()->showMessage(tr("Loading the reference frame: %1 (%2 %)")
                                                 .arg(text)
                                                 .arg(static_cast<int>(fraction * 100)));
                }
            });
        },
        m_cancel);
    m_background.run([project = std::move(project), reporter] { return loadStructure(project, reporter); },
                     [this, generation](Loaded loaded) {
                         if (generation != m_loadGeneration || loaded.cancelled) return;
                         m_loading = false;
                         m_cancel.reset();
                         m_molecule = loaded.molecule;
                         m_structurePreview = loaded.preview;
                         m_structureSummary = loaded.summary;
                         m_loadProblem = loaded.structureProblem;
                         updatePreview();
                         QString message = m_molecule ? m_structureSummary
                                                      : tr("No reference frame: %1").arg(m_loadProblem);
                         if (!loaded.notes.isEmpty()) {
                             message += QStringLiteral(" (") + loaded.notes.join(QStringLiteral("; ")) + QLatin1Char(')');
                         }
                         statusBar()->showMessage(message);
                         emit projectLoaded();
                     });
}

// ---- the canvas (W2-A) ------------------------------------------------------------------------------

namespace {

// Sets m_canvasEdit while the canvas's own change is stored.
struct CanvasEdit {
    explicit CanvasEdit(bool& target) : flag(target), was(target) { target = true; }
    ~CanvasEdit() { flag = was; }
    CanvasEdit(const CanvasEdit&) = delete;
    CanvasEdit& operator=(const CanvasEdit&) = delete;
    bool& flag;
    bool was;
};

constexpr const char* kPlacedByHand = "Placed by hand (not fitted) - 'Refine' fits from here.";

}  // namespace

void MainWindow::connectCanvas() {
    for (const Page& page : m_pages) {
        if (auto* layout = qobject_cast<pages::LayoutPage*>(page.widget)) m_layoutPage = layout;
        if (auto* structure = qobject_cast<pages::StructurePage*>(page.widget)) m_structurePage = structure;
    }
    // D4: the GPU renderer draws the structure panel's 3D layer while dragging when there is a GPU
    if (m_renderer != render::Renderer::Software) {
        auto gpu = std::make_shared<render::panels::GpuStageRenderers>(render::panels::GpuStageOptions{});
        m_gpu = gpu;
        m_background.run(
            [gpu]() -> std::shared_ptr<const render::panels::StageRenderer> {
                try {
                    if (gpu->usable() && !gpu->cpuDevice()) return gpu->factory()->choose("vulkan").renderer;
                } catch (const std::exception&) {
                }
                return nullptr;
            },
            [this](std::shared_ptr<const render::panels::StageRenderer> renderer) {
                if (renderer) m_preview->setFastStageRenderer(std::move(renderer));
            });
    }
    // a drag is one undo step (D5)
    connect(m_preview, &PreviewWidget::dragStarted, this, [this](const QString& mode) {
        m_document.beginGroup(mode == QStringLiteral("rotate") ? "Rotate structure"
                              : mode == QStringLiteral("map") ? "Move density map"
                                                               : "Rotate density map");
    });
    connect(m_preview, &PreviewWidget::dragFinished, this, [this] {
        const CanvasEdit edit(m_canvasEdit);
        if (m_document.inGroup()) m_document.endGroup();
    });
    connect(m_preview, &PreviewWidget::dragCancelled, this, [this] {
        if (m_document.inGroup()) m_document.cancelGroup();  // Esc: the project as before the drag
    });
    connect(m_preview, &PreviewWidget::rotated, this, [this](const std::array<double, 3>& angles) {
        const CanvasEdit edit(m_canvasEdit);
        core::Array stored;  // the sliders' 3 decimals (MDATK)
        for (const double a : angles) stored.emplace_back(core::pythonRound(a, 3));
        m_document.set("visual.rotation", core::Value(std::move(stored)), "Rotate structure");
    });
    connect(m_preview, &PreviewWidget::mapMoved, this, [this](const std::vector<double>& transform) {
        // MDATK's _on_map_moved: the placement belongs to the frame shown on the canvas
        const CanvasEdit edit(m_canvasEdit);
        m_document.edit("Move density map", [&](core::Project& p) {
            render::panels::storeDensityPlacement(p.density, transform, m_preview->structureModel(), kPlacedByHand);
        });
    });
    connect(m_preview, &PreviewWidget::zoomRequested, this, [this](double factor) {
        // the stored zoom (AP-01: not the one a spin box shows), within MDATK's wheel limits
        const double zoom = std::min(50.0, std::max(0.05, m_document.project().visual.zoom * factor));
        m_preview->setZoom(zoom);
        // a burst of notches is one undo step (LE-10): a pause ends it (any other change too)
        if (!m_wheelClock.isValid() || m_wheelClock.elapsed() > kWheelBurstMs) m_document.endMerge();
        m_wheelClock.start();
        m_document.setMerged("wheel", "visual.zoom", core::Value(zoom), "Zoom the structure");
    });
    connect(m_preview, &PreviewWidget::swapRequested, this, [this](int from, int to) {
        if (m_layoutPage) m_layoutPage->swapPanels(from, to);
    });
    connect(m_preview, &PreviewWidget::slotActivated, this,
            [this](int slot) { openSection(QStringLiteral("layout"), QStringLiteral("slot%1").arg(slot + 1)); });
    if (m_layoutPage) {
        connect(m_layoutPage, &pages::LayoutPage::logLine, this, &MainWindow::appendLog);
        connect(m_layoutPage, &pages::LayoutPage::previewNotUpdated, this, [this](const QString& problem) {
            appendLog(tr("Preview not updated: %1").arg(problem));
            // the next drawing's warnings are news again (MDATK logs those of each redraw)
            m_lastWarnings.clear();
            m_warningsKnown = true;
        });
        connect(m_layoutPage, &pages::LayoutPage::previewRequested, this, &MainWindow::refreshPreview);
    }
    if (m_structurePage) {
        connect(m_structurePage, &pages::StructurePage::mapModeChanged, m_preview, &PreviewWidget::setMapMode);
    }
}

}  // namespace mdbridge::app
