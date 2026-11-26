// mdbridge_gui_canvas_tests: the figure canvas as layout editor and 3D view (src/app/gui/preview_widget.cpp,
// milestone 5 W2-A) and the Layout & panels page (gui/pages/layout_page.cpp): the counterparts of
// MDATK's tests/test_gui.py test_canvas_shows_layout_with_example_data,
// test_drag_rotates_structure_and_swaps_panels and test_slot_defaults_for_residue_panels_and_wheel_zoom,
// one undo step per drag, Esc, AP-20 (rotation without rebuilding the figure) and the automatic
// titles after a group rename (MDATK 0.1.1). The saved files of the same mouse paths are compared
// with MDATK's by the canvas cases of mdbridge_gui_tests (tests/data/gui/canvas).

#include <algorithm>
#include <array>
#include <filesystem>
#include <string>
#include <vector>

#include <QApplication>
#include <QComboBox>
#include <QCoreApplication>
#include <QLineEdit>
#include <QLocale>
#include <QMouseEvent>
#include <QSignalSpy>
#include <QTemporaryDir>
#include <QTest>
#include <QWheelEvent>

#include "gui/main_window.hpp"
#include "gui/pages/layout_page.hpp"
#include "gui/preview_widget.hpp"
#include "mdbridge/analysis/registry.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/project.hpp"
#include "mdbridge/core/project_document.hpp"

namespace fs = std::filesystem;
using namespace mdbridge;

namespace {

const fs::path kData = core::pathFromUtf8(MDBRIDGE_TEST_DATA);

void send(QWidget* w, QEvent::Type type, const QPointF& at, Qt::MouseButton button, Qt::MouseButtons buttons,
          Qt::KeyboardModifiers mods = Qt::NoModifier) {
    QMouseEvent e(type, at, w->mapToGlobal(at), button, buttons, mods);
    QCoreApplication::sendEvent(w, &e);
}

void drag(QWidget* w, const QPointF& from, const std::vector<QPointF>& path, Qt::KeyboardModifiers mods = Qt::NoModifier) {
    send(w, QEvent::MouseButtonPress, from, Qt::LeftButton, Qt::LeftButton, mods);
    for (const QPointF& p : path) send(w, QEvent::MouseMove, p, Qt::NoButton, Qt::LeftButton, mods);
    send(w, QEvent::MouseButtonRelease, path.empty() ? from : path.back(), Qt::LeftButton, Qt::NoButton, mods);
}

void wheel(QWidget* w, const QPointF& at, int notches) {
    QWheelEvent e(at, w->mapToGlobal(at), QPoint(), QPoint(0, 120 * notches), Qt::NoButton, Qt::NoModifier,
                  Qt::NoScrollPhase, false);
    QCoreApplication::sendEvent(w, &e);
}

QPointF centre(const QRectF& r) { return r.center(); }

app::pages::LayoutPage* layoutPage(app::MainWindow& window) {
    for (const app::MainWindow::Page& p : window.pages()) {
        if (auto* l = qobject_cast<app::pages::LayoutPage*>(p.widget)) return l;
    }
    return nullptr;
}

}  // namespace

class CanvasTests : public QObject {
    Q_OBJECT

private slots:
    void init();
    void cleanup();
    void canvasShowsLayout();
    void dragRotatesAndSwapsPanels();
    void dragIsOneUndoStepAndEscCancels();
    void rotationRedrawsOnlyTheStructure();
    void slotDefaultsAndWheelZoom();
    void titlesFollowAGroupRename();

private:
    bool settled() {
        return QTest::qWaitFor([this] { return !m_window->loading() && m_window->previewSettled() && m_window->preview()->current(); },
                               60000);
    }
    QTemporaryDir* m_dir = nullptr;
    app::MainWindow* m_window = nullptr;
};

// The demo project of the canvas cases (a structure panel, two time series, a heat map; 24 residues).
void CanvasTests::init() {
    m_dir = new QTemporaryDir;
    QVERIFY(m_dir->isValid());
    const fs::path dir = core::pathFromUtf8(m_dir->path().toStdString());
    for (const auto& [name, source] : std::vector<std::pair<const char*, const char*>>{
             {"demo_dry.parm7", "topology/demo_dry.parm7"},
             {"demo_dry.nc", "io/trajectories/demo_dry.nc"},
             {"p.mdatk.toml", "gui/canvas/drag_rotate_swap/in/p.mdatk.toml"}}) {
        fs::copy_file(kData / source, dir / name);
    }
    m_window = new app::MainWindow(render::Renderer::Software);
    m_window->preview()->setFixedSize(960, 540);
    m_window->resize(1600, 1000);
    m_window->show();
    QVERIFY(QTest::qWaitForWindowExposed(m_window));
    QVERIFY(m_window->startWith(dir / "p.mdatk.toml"));
    QVERIFY(settled());
}

void CanvasTests::cleanup() {
    m_window->cancelLoad();
    delete m_window;
    m_window = nullptr;
    delete m_dir;
    m_dir = nullptr;
}

void CanvasTests::canvasShowsLayout() {
    app::PreviewWidget* canvas = m_window->preview();
    QVERIFY(!canvas->image().isNull());
    QCOMPARE(canvas->slotRects().size(), m_window->document().project().layout.panels.size());
    QCOMPARE(canvas->structureSlot(), 0);
    QVERIFY(canvas->structureCentre().has_value());
    QCOMPARE(canvas->structureModel().size(), std::size_t{24});
    QCOMPARE(m_window->currentPage(), QStringLiteral("input"));
    QVERIFY(m_window->openSection(QStringLiteral("layout"), QStringLiteral("slot2")));
    QCOMPARE(m_window->currentPage(), QStringLiteral("layout"));
    // a double-click on a panel opens its slot's editor
    m_window->showPage(QStringLiteral("input"));
    const QPointF at = centre(canvas->slotRects()[2]);
    send(canvas, QEvent::MouseButtonPress, at, Qt::LeftButton, Qt::LeftButton);
    send(canvas, QEvent::MouseButtonRelease, at, Qt::LeftButton, Qt::NoButton);
    send(canvas, QEvent::MouseButtonDblClick, at, Qt::LeftButton, Qt::LeftButton);
    send(canvas, QEvent::MouseButtonRelease, at, Qt::LeftButton, Qt::NoButton);
    QCOMPARE(m_window->currentPage(), QStringLiteral("layout"));
}

void CanvasTests::dragRotatesAndSwapsPanels() {
    app::PreviewWidget* canvas = m_window->preview();
    const core::ProjectDocument& doc = m_window->document();
    const double y0 = doc.project().visual.rotation[1];
    const QPointF c = centre(canvas->slotRects()[0]);
    drag(canvas, c, {c + QPointF(40, 0)});
    QVERIFY(settled());
    QCOMPARE(doc.project().visual.rotation[1], y0 + 20.0);  // 0.5 degrees per pixel
    // a press moved less than 4 pixels is no drag
    drag(canvas, c, {c + QPointF(1, 1), c + QPointF(2, 0)});
    QVERIFY(settled());
    QCOMPARE(doc.project().visual.rotation[1], y0 + 20.0);
    // panel 2 onto panel 3: swapped; the drag shows the target and the panel's type meanwhile
    const std::vector<core::PanelSpec> before = doc.project().layout.panels;
    const std::vector<QRectF> rects = canvas->slotRects();
    send(canvas, QEvent::MouseButtonPress, centre(rects[1]), Qt::LeftButton, Qt::LeftButton);
    send(canvas, QEvent::MouseMove, centre(rects[2]), Qt::NoButton, Qt::LeftButton);
    QCOMPARE(canvas->dragMode(), QStringLiteral("move"));
    QCOMPARE(canvas->dragTarget(), 2);
    send(canvas, QEvent::MouseButtonRelease, centre(rects[2]), Qt::LeftButton, Qt::NoButton);
    QVERIFY(settled());
    const std::vector<core::PanelSpec>& after = doc.project().layout.panels;
    QCOMPARE(after[1].series, before[2].series);
    QCOMPARE(after[2].series, before[1].series);
    // the structure panel moves by its title strip (or with Shift)
    drag(canvas, QPointF(c.x(), rects[0].top() + 0.05 * rects[0].height()), {centre(rects[1])});
    QVERIFY(settled());
    QCOMPARE(doc.project().layout.panels[1].type, std::string("structure"));
    QCOMPARE(canvas->structureSlot(), 1);
}

void CanvasTests::dragIsOneUndoStepAndEscCancels() {
    app::PreviewWidget* canvas = m_window->preview();
    core::ProjectDocument& doc = m_window->document();
    const std::vector<double> start = doc.project().visual.rotation;
    const QPointF c = centre(canvas->slotRects()[0]);
    std::vector<QPointF> path;
    for (int k = 1; k <= 10; ++k) path.push_back(c + QPointF(6.0 * k, 2.0 * k));
    drag(canvas, c, path);
    QVERIFY(settled());
    QVERIFY(doc.project().visual.rotation != start);
    QCOMPARE(QString::fromStdString(doc.undoLabel()), QStringLiteral("Rotate structure"));
    doc.undo();  // ten motions, one step
    QVERIFY(settled());
    QCOMPARE(doc.project().visual.rotation, start);
    QVERIFY(!doc.canUndo());
    // Esc during a drag: the project as before it, no step
    send(canvas, QEvent::MouseButtonPress, c, Qt::LeftButton, Qt::LeftButton);
    send(canvas, QEvent::MouseMove, c + QPointF(30, 0), Qt::NoButton, Qt::LeftButton);
    QVERIFY(doc.project().visual.rotation != start);
    QTest::keyClick(canvas, Qt::Key_Escape);
    send(canvas, QEvent::MouseButtonRelease, c + QPointF(30, 0), Qt::LeftButton, Qt::NoButton);
    QVERIFY(settled());
    QCOMPARE(doc.project().visual.rotation, start);
    QVERIFY(!doc.canUndo());
}

// AP-20: turning the model (sliders, undo) redraws the structure panel of the figure already built
// on the worker; the figure is not built again.
void CanvasTests::rotationRedrawsOnlyTheStructure() {
    app::PreviewWidget* canvas = m_window->preview();
    QSignalSpy built(canvas, &app::PreviewWidget::drawn);
    const QImage before = canvas->image();
    m_window->document().set("visual.rotation", core::Value(core::Array{30.0, 10.0, 0.0}));
    QVERIFY(settled());
    QCOMPARE(built.count(), 0);
    QVERIFY(canvas->image() != before);
    QCOMPARE(canvas->angles(), (std::array<double, 3>{30.0, 10.0, 0.0}));
    // any other setting rebuilds it
    m_window->document().set("visual.title", core::Value(std::string("Demo")));
    QVERIFY(settled());
    QCOMPARE(built.count(), 1);
}

void CanvasTests::slotDefaultsAndWheelZoom() {
    app::pages::LayoutPage* page = layoutPage(*m_window);
    QVERIFY(page);
    core::ProjectDocument& doc = m_window->document();
    auto* type = m_window->findChild<QComboBox*>(QStringLiteral("layout.panels.1.type"));
    QVERIFY(type);
    const int secondary = type->findData(QStringLiteral("secondary"));
    QVERIFY(secondary >= 0);
    type->setCurrentIndex(secondary);
    emit type->activated(secondary);
    QVERIFY(page->slotState(1).options.contains(QStringLiteral("\"source\": \"dynamic\"")));
    QVERIFY(doc.project().layout.panels[1].series.empty());  // residue panels take no data keys
    QCOMPARE(doc.project().layout.panels[1].type, std::string("secondary"));
    QVERIFY(settled());
    // the wheel over the structure zooms by 1.15 per notch; over another panel nothing
    app::PreviewWidget* canvas = m_window->preview();
    const double z0 = doc.project().visual.zoom;
    wheel(canvas, centre(canvas->slotRects()[0]), 1);
    QVERIFY(settled());
    QCOMPARE(doc.project().visual.zoom, z0 * 1.15);
    wheel(canvas, centre(canvas->slotRects()[2]), 1);
    QVERIFY(settled());
    QCOMPARE(doc.project().visual.zoom, z0 * 1.15);
    // MDATK's wheel limits
    for (int k = 0; k < 40; ++k) wheel(canvas, centre(canvas->slotRects()[0]), 1);
    QVERIFY(settled());
    QCOMPARE(doc.project().visual.zoom, 50.0);
}

// MDATK 0.1.1 _adopt_project: a renamed group renames the data keys, and the titles that were their
// data's labels follow; default titles and typed ones stay.
void CanvasTests::titlesFollowAGroupRename() {
    core::Project p;
    p.layout.panels = {
        {"timeseries", {"contacts@Aptamer+Target/atom_pairs"}, "", {}, {}},
        {"timeseries", {"rmsd/Complex"}, "RMSD", {}, {}},
        {"timeseries", {"hbonds@Aptamer+Target/count"}, "My own", {}, {}},
    };
    p.layout.panels[0].title = analysis::keyLabel(p.layout.panels[0].series[0], p);
    const core::Project before = p;
    const std::vector<std::string> changed = core::renameGroup(p, "Aptamer", "Apt");
    const std::vector<std::string> titles = app::pages::retitleAutomaticPanels(p, before, changed);
    QCOMPARE(titles, std::vector<std::string>{"layout.panels.0.title"});
    QCOMPARE(p.layout.panels[0].series[0], std::string("contacts@Apt+Target/atom_pairs"));
    QCOMPARE(p.layout.panels[0].title, analysis::keyLabel("contacts@Apt+Target/atom_pairs", p));
    QVERIFY(p.layout.panels[0].title != before.layout.panels[0].title);
    QCOMPARE(p.layout.panels[1].title, std::string("RMSD"));
    QCOMPARE(p.layout.panels[2].title, std::string("My own"));
}

int main(int argc, char** argv) {
    qputenv("QT_QPA_PLATFORM", "offscreen");
    QLocale::setDefault(QLocale::c());
    QApplication app(argc, argv);
    CanvasTests tests;
    return QTest::qExec(&tests, argc, argv);
}

#include "test_canvas.moc"
