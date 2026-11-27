// mdbridge_gui_canvas_tests: the figure canvas as layout editor and 3D view (src/app/gui/preview_widget.cpp,
// milestone 5 W2-A) and the Layout & panels page (gui/pages/layout_page.cpp): the counterparts of
// MDATK's tests/test_gui.py test_canvas_shows_layout_with_example_data,
// test_drag_rotates_structure_and_swaps_panels and test_slot_defaults_for_residue_panels_and_wheel_zoom,
// one undo step per drag, Esc, AP-20 (rotation without rebuilding the figure) and the automatic
// titles after a group rename (MDATK 0.1.1), and AP-08's Remove buttons (MDATK 0.1.2). The saved
// files of the same mouse paths are compared with MDATK's by the canvas cases of mdbridge_gui_tests
// (tests/data/gui/canvas).

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <string>
#include <vector>

#include <QApplication>
#include <QComboBox>
#include <QCoreApplication>
#include <QLabel>
#include <QLineEdit>
#include <QLocale>
#include <QMouseEvent>
#include <QPushButton>
#include <QSignalSpy>
#include <QTemporaryDir>
#include <QTest>
#include <QWheelEvent>
#include <QKeyEvent>

#include "gui/main_window.hpp"
#include "gui/pages/layout_page.hpp"
#include "gui/preview_widget.hpp"
#include "mdbridge/analysis/registry.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/project.hpp"
#include "mdbridge/core/project_document.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/render/orientation.hpp"

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
    void wheelBurstIsOneUndoStep();
    void keyboardRotatesAndZooms();
    void titlesFollowAGroupRename();
    void removePanelsAndPanelsWithoutSlot();

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

// LE-10: notches at most MainWindow::kWheelBurstMs apart are one undo step; after a pause, a new one
void CanvasTests::wheelBurstIsOneUndoStep() {
    app::PreviewWidget* canvas = m_window->preview();
    core::ProjectDocument& doc = m_window->document();
    const double z0 = doc.project().visual.zoom;
    const QPointF at = centre(canvas->slotRects()[0]);
    for (int k = 0; k < 5; ++k) wheel(canvas, at, 1);
    QVERIFY(settled());
    QTest::qWait(app::MainWindow::kWheelBurstMs + 200);
    for (int k = 0; k < 3; ++k) wheel(canvas, at, -1);
    QVERIFY(settled());
    QCOMPARE(QString::fromStdString(doc.undoLabel()), QStringLiteral("Zoom the structure"));
    doc.undo();
    QVERIFY(settled());
    QVERIFY(std::abs(doc.project().visual.zoom - z0 * std::pow(1.15, 5)) < 1e-12);
    doc.undo();
    QVERIFY(settled());
    QCOMPARE(doc.project().visual.zoom, z0);
    QVERIFY(!doc.canUndo());
}

// D10 (WN-13; MDATK's canvas has no keyboard use): the canvas takes the focus with Tab, shows it,
// has an accessible name; the arrow keys turn the structure as a 10-pixel drag (Shift: 2 pixels), +
// and - zoom it as a wheel notch; a burst of keys is one undo step; Ctrl+arrow does nothing here.
void CanvasTests::keyboardRotatesAndZooms() {
    app::PreviewWidget* canvas = m_window->preview();
    core::ProjectDocument& doc = m_window->document();
    QVERIFY(canvas->focusPolicy() & Qt::TabFocus);
    QVERIFY(!canvas->accessibleName().isEmpty());
    QVERIFY(canvas->accessibleDescription().contains(QStringLiteral("arrow keys")));
    canvas->setFocus(Qt::TabFocusReason);
    QTRY_VERIFY(canvas->hasFocus());
    const std::vector<double> r0 = doc.project().visual.rotation;
    const double z0 = doc.project().visual.zoom;
    QTest::keyClick(canvas, Qt::Key_Right);
    QTest::keyClick(canvas, Qt::Key_Right);
    QTest::keyClick(canvas, Qt::Key_Up, Qt::ShiftModifier);
    QVERIFY(settled());
    // the angles of the same mouse path: right 10 + 10 pixels, then up 2
    render::Vec3d expected = render::drag({r0[0], r0[1], r0[2]}, 10.0, 0.0, 0.5);
    expected = render::drag(expected, 10.0, 0.0, 0.5);
    expected = render::drag(expected, 0.0, -2.0, 0.5);
    const std::vector<double>& r = doc.project().visual.rotation;
    QCOMPARE(r.size(), std::size_t{3});
    for (std::size_t k = 0; k < 3; ++k) QVERIFY(std::abs(r[k] - core::pythonRound(expected[k], 3)) < 1e-9);
    QCOMPARE(canvas->angles()[0], expected[0]);
    QCOMPARE(QString::fromStdString(doc.undoLabel()), QStringLiteral("Rotate structure"));
    QTest::keyClick(canvas, Qt::Key_Right, Qt::ControlModifier);  // not the canvas's
    QCOMPARE(doc.project().visual.rotation, r);
    QTest::qWait(app::MainWindow::kWheelBurstMs + 200);
    QTest::keyClick(canvas, Qt::Key_Plus);
    QTest::keyClick(canvas, Qt::Key_Plus);
    QTest::keyClick(canvas, Qt::Key_Minus);
    QVERIFY(settled());
    QVERIFY(std::abs(doc.project().visual.zoom - z0 * 1.15) < 1e-12);
    QCOMPARE(QString::fromStdString(doc.undoLabel()), QStringLiteral("Zoom the structure"));
    doc.undo();  // the zoom keys
    QVERIFY(settled());
    QCOMPARE(doc.project().visual.zoom, z0);
    doc.undo();  // the three arrow keys
    QVERIFY(settled());
    QCOMPARE(doc.project().visual.rotation, r0);
    QVERIFY(!doc.canUndo());
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

// AP-08 (MDATK 0.1.2): "Remove this panel" deletes a slot's panel (the panels after it move up, one
// undo step; an empty slot has nothing to remove); the panels a smaller preset has no slot for are
// listed with a note and a Remove each; a removal that would move unparsed options is refused.
void CanvasTests::removePanelsAndPanelsWithoutSlot() {
    app::pages::LayoutPage* page = layoutPage(*m_window);
    QVERIFY(page);
    core::ProjectDocument& doc = m_window->document();
    const std::vector<core::PanelSpec> panels = doc.project().layout.panels;
    QCOMPARE(panels.size(), std::size_t{4});
    const auto button = [this](const QString& name) { return m_window->findChild<QPushButton*>(name); };
    QVERIFY(page->slotlessPanels().empty());
    QVERIFY(!m_window->findChild<QLabel*>(QStringLiteral("layout.slotless.note")));

    // a smaller preset: three panels without a slot, listed in order
    doc.set("layout.preset", core::Value(std::string("single")), "Set Grid preset");
    QCOMPARE(page->slotCount(), 1);
    const std::vector<QString> listed = page->slotlessPanels();
    QCOMPARE(listed.size(), std::size_t{3});
    QCOMPARE(listed[0], QStringLiteral("Panel 2: Time series \u2013 RMSD"));  // its title
    QVERIFY(listed[1].startsWith(QStringLiteral("Panel 3: Time series \u2013 ")));  // its data's label
    QVERIFY(listed[2].startsWith(QStringLiteral("Panel 4: Residue heatmap over time \u2013 ")));
    QLabel* note = m_window->findChild<QLabel*>(QStringLiteral("layout.slotless.note"));
    QVERIFY(note);
    QVERIFY(note->text().startsWith(QStringLiteral("3 more panels of this project have no slot in this grid")));
    QPushButton* third = button(QStringLiteral("layout.panels.2.remove"));
    QVERIFY(third);
    third->click();
    QCOMPARE(doc.project().layout.panels.size(), std::size_t{3});
    QCOMPARE(doc.project().layout.panels[2].type, panels[3].type);
    QCOMPARE(QString::fromStdString(doc.undoLabel()), QStringLiteral("Remove panel 3"));
    QCOMPARE(page->slotlessPanels().size(), std::size_t{2});
    doc.undo();
    QCOMPARE(doc.project().layout.panels, panels);
    doc.undo();  // the preset
    QCOMPARE(doc.project().layout.preset, std::string("large_plus_3"));
    QVERIFY(page->slotlessPanels().empty());
    QVERIFY(QTest::qWaitFor([this] { return !m_window->findChild<QLabel*>(QStringLiteral("layout.slotless.note")); }));

    // a slot's panel: the panels after it move up, the last slot is empty and has nothing to remove
    QPushButton* second = button(QStringLiteral("layout.panels.1.remove"));
    QVERIFY(second && second->isEnabled());
    second->click();
    QCOMPARE(doc.project().layout.panels.size(), std::size_t{3});
    QCOMPARE(doc.project().layout.panels[1], panels[2]);
    QCOMPARE(doc.project().layout.panels[2], panels[3]);
    QCOMPARE(page->slotState(1).title, QString());  // panel 3's (no title of its own)
    QCOMPARE(page->slotState(2).type, QString::fromStdString(panels[3].type));
    QVERIFY(!button(QStringLiteral("layout.panels.3.remove"))->isEnabled());
    QVERIFY(button(QStringLiteral("layout.panels.2.remove"))->isEnabled());
    doc.undo();
    QCOMPARE(doc.project().layout.panels, panels);
    QVERIFY(button(QStringLiteral("layout.panels.3.remove"))->isEnabled());

    // options that do not parse in slot 3: removing a panel before it is refused, removing it is not
    auto* options = m_window->findChild<QLineEdit*>(QStringLiteral("layout.panels.2.options"));
    QVERIFY(options);
    options->setText(QStringLiteral("{bad"));
    emit options->textEdited(options->text());
    emit options->editingFinished();
    QVERIFY(!page->optionsProblem().isEmpty());
    QSignalSpy logged(page, &app::pages::LayoutPage::logLine);
    QVERIFY(!page->removePanel(0));
    QCOMPARE(logged.count(), 1);
    QVERIFY(logged.at(0).at(0).toString().startsWith(
        QStringLiteral("Panel not removed; fix this first: Layout & panels")));
    QCOMPARE(doc.project().layout.panels, panels);
    QVERIFY(page->removePanel(3));  // after it: nothing moves
    QVERIFY(page->removePanel(2));  // the panel itself: its options go with it
    QVERIFY(page->optionsProblem().isEmpty());
    QCOMPARE(doc.project().layout.panels.size(), std::size_t{2});
    QVERIFY(settled());
}

int main(int argc, char** argv) {
    qputenv("QT_QPA_PLATFORM", "offscreen");
    QLocale::setDefault(QLocale::c());
    QApplication app(argc, argv);
    CanvasTests tests;
    return QTest::qExec(&tests, argc, argv);
}

#include "test_canvas.moc"
