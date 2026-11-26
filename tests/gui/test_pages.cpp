// mdbridge_gui_pages_tests: the Structure view and Outputs settings pages (milestone 5 W1-C,
// src/app/gui/pages): MDATK's test_representation_editor_and_3d_preview (the editor part; the 3D
// preview is the canvas's, W2-A) and test_value_widgets_ignore_the_mouse_wheel, the custom editors
// (rotation sliders and angle label, chart format check boxes, AP-04), the density buttons (D2),
// the map placement rule, sections and revealField, and a screenshot of each page. Qt Test on the
// offscreen platform.

#include <cmath>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include <QAbstractSlider>
#include <QAbstractSpinBox>
#include <QApplication>
#include <QCheckBox>
#include <QComboBox>
#include <QLabel>
#include <QLocale>
#include <QPushButton>
#include <QScrollArea>
#include <QScrollBar>
#include <QSignalSpy>
#include <QSlider>
#include <QTableView>
#include <QTest>
#include <QToolButton>
#include <QWheelEvent>

#include "gui/forms/custom_fields.hpp"
#include "gui/forms/field_widgets.hpp"
#include "gui/forms/form_context.hpp"
#include "gui/forms/list_editor.hpp"
#include "gui/pages/form_page.hpp"
#include "gui/pages/outputs_page.hpp"
#include "gui/pages/page_registry.hpp"
#include "gui/pages/structure_page.hpp"
#include "gui_harness.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/project_document.hpp"

namespace fs = std::filesystem;
using namespace mdbridge;
using app::forms::ChoiceChecksField;
using app::forms::RotationField;
using app::pages::OutputsPage;
using app::pages::StructurePage;

namespace {

const fs::path kData = core::pathFromUtf8(MDBRIDGE_TEST_DATA);

QString text(const core::ProjectDocument& doc, const char* path) { return QString::fromStdString(doc.text(path)); }

template <class Page>
struct Shown {
    core::ProjectDocument doc;
    app::forms::FormContext context{doc};
    std::unique_ptr<Page> page;
    Shown() {
        page = std::make_unique<Page>(context);
        page->resize(560, 900);
        page->show();
        QVERIFY(QTest::qWaitForWindowExposed(page.get()));
    }
    ~Shown() { page.reset(); }
};

void wheelOver(QWidget* target) {
    QWheelEvent event(QPointF(5, 5), target->mapToGlobal(QPointF(5, 5)), QPoint(), QPoint(0, 120), Qt::NoButton,
                      Qt::NoModifier, Qt::NoScrollPhase, false);
    QApplication::sendEvent(target, &event);
}

}  // namespace

class PagesTest : public QObject {
    Q_OBJECT

private slots:
    void registryInMdatksOrder();
    void sectionsAreNamedAndFold();
    void representationEditorRowsRoundTrip();
    void valueWidgetsIgnoreTheMouseWheel();
    void rotationSlidersAndAngleLabel();
    void resetButtons();
    void densityButtonsAndStatus();
    void anotherMapForgetsThePlacementInOneStep();
    void chartFormatsKeepEveryEntry();
    void chartStylesAreMdbridgesAndTypedNamesStay();
    void revealFieldOpensTheSection();
    void screenshots();
};

void PagesTest::registryInMdatksOrder() {
    core::ProjectDocument doc;
    app::forms::FormContext context(doc);
    std::vector<QWidget*> pages = app::pages::createPages(context);
    QStringList ids;
    for (QWidget* p : pages) {
        ids.append(app::pages::pageIdOf(p));
        QVERIFY(!app::pages::pageTitleOf(p).isEmpty());
    }
    // MDATK's order (core::kProblemPages), the pages there are so far
    QStringList order;
    for (const core::ProblemPage& p : core::kProblemPages) {
        const QString id = QString::fromUtf8(p.id.data(), static_cast<qsizetype>(p.id.size()));
        if (ids.contains(id)) order.append(id);
    }
    QCOMPARE(ids, order);
    QVERIFY(ids.contains(QStringLiteral("structure")));
    QVERIFY(ids.contains(QStringLiteral("outputs")));
    for (QWidget* p : pages) {
        if (app::pages::pageIdOf(p) == QStringLiteral("structure")) {
            QCOMPARE(app::pages::pageTitleOf(p), QStringLiteral("Structure view"));
        }
        delete p;
    }
}

void PagesTest::sectionsAreNamedAndFold() {
    Shown<StructurePage> s;
    Shown<OutputsPage> o;
    QStringList names;
    for (const auto* sec : s.page->sections()) names.append(sec->objectName());
    for (const auto* sec : o.page->sections()) names.append(sec->objectName());
    QCOMPARE(names, QStringList({"section:structure.orientation", "section:structure.representations",
                                 "section:structure.density", "section:structure.tags", "section:structure.appearance",
                                 "section:outputs.locations", "section:outputs.video", "section:outputs.charts",
                                 "section:outputs.pdb"}));
    QCOMPARE(s.page->section(QStringLiteral("density"))->title(), QStringLiteral("Density map"));
    auto* tags = s.page->section(QStringLiteral("tags"));
    QVERIFY(tags->form()->isVisible());
    QTest::mouseClick(tags->header(), Qt::LeftButton);
    QVERIFY(!tags->isExpanded());
    QVERIFY(!tags->form()->isVisible());
    QTest::mouseClick(tags->header(), Qt::LeftButton);
    QVERIFY(tags->form()->isVisible());
}

// MDATK test_representation_editor_and_3d_preview, the editor: rows added and filled, read back,
// one removed, the rest reloaded unchanged.
void PagesTest::representationEditorRowsRoundTrip() {
    Shown<StructurePage> s;
    QWidget* root = s.page.get();
    QCOMPARE(s.page->representations()->addButton()->text(), QStringLiteral("Add representation"));
    guitest::tableAction(root, QStringLiteral("visual.representations"), "add", -1, 0);
    guitest::setKey(root, QStringLiteral("visual.representations.0.groups"), "Aptamer");
    guitest::setKey(root, QStringLiteral("visual.representations.0.sugars"), true);
    guitest::tableAction(root, QStringLiteral("visual.representations"), "add", -1, 0);
    guitest::setKey(root, QStringLiteral("visual.representations.1.style"), "licorice");
    guitest::setKey(root, QStringLiteral("visual.representations.1.groups"), "Target");
    guitest::setKey(root, QStringLiteral("visual.representations.1.color"), "element");
    guitest::setKey(root, QStringLiteral("visual.representations.1.size"), 0.8);
    const std::vector<core::RepSpec>& reps = s.doc.project().visual.representations;
    QCOMPARE(reps.size(), std::size_t{2});
    QCOMPARE(QString::fromStdString(reps[0].style), QStringLiteral("cartoon"));
    QCOMPARE(reps[0].groups, std::vector<std::string>{"Aptamer"});
    QCOMPARE(QString::fromStdString(reps[0].color), QStringLiteral("auto"));
    QVERIFY(reps[0].sugars);
    QCOMPARE(QString::fromStdString(reps[1].style), QStringLiteral("licorice"));
    QCOMPARE(reps[1].groups, std::vector<std::string>{"Target"});
    QCOMPARE(QString::fromStdString(reps[1].color), QStringLiteral("element"));
    QCOMPARE(reps[1].size, 0.8);
    guitest::tableAction(root, QStringLiteral("visual.representations"), "remove", 1, 0);
    QCOMPARE(s.doc.project().visual.representations.size(), std::size_t{1});
    QCOMPARE(QString::fromStdString(s.doc.project().visual.representations[0].style), QStringLiteral("cartoon"));
    // load the project again: the rows show it, and reading them back changes nothing
    const core::Project p = s.doc.project();
    s.doc.replace(p);
    QCOMPARE(s.page->representations()->rowCount(), 1);
    s.context.commitPendingEdits();
    QVERIFY(s.doc.project() == p);
    QVERIFY(!s.doc.canUndo());
}

// MDATK test_value_widgets_ignore_the_mouse_wheel on these pages: every spin box, drop-down and
// slider (not scroll bars) ignores the wheel without the focus; a wheel over one changes nothing.
void PagesTest::valueWidgetsIgnoreTheMouseWheel() {
    core::ProjectDocument doc;
    app::forms::FormContext context(doc);
    std::vector<std::unique_ptr<QWidget>> pages;
    for (QWidget* p : app::pages::createPages(context)) pages.emplace_back(p);
    QStringList missing;
    int seen = 0;
    for (const auto& page : pages) {
        page->resize(560, 900);
        page->show();
        for (QWidget* w : page->findChildren<QWidget*>()) {
            const bool value = qobject_cast<QAbstractSpinBox*>(w) || qobject_cast<QComboBox*>(w) ||
                               (qobject_cast<QAbstractSlider*>(w) && !qobject_cast<QScrollBar*>(w));
            if (value) {
                ++seen;
                if (!app::forms::hasWheelGuard(w)) missing.append(app::pages::pageIdOf(page.get()) + ": " + w->objectName());
            }
        }
    }
    QVERIFY2(seen >= 30, qPrintable(QString::number(seen)));
    QVERIFY2(missing.isEmpty(), qPrintable(missing.join(QStringLiteral(", "))));
    auto* structure = qobject_cast<StructurePage*>(pages[0].get());
    QVERIFY(structure);
    QVERIFY(QTest::qWaitForWindowExposed(structure));
    auto* zoom = structure->findChild<QAbstractSpinBox*>(QStringLiteral("visual.zoom"));
    auto* slider = structure->findChild<QSlider*>(QStringLiteral("visual.rotation.1"));
    auto* theme = structure->findChild<QComboBox*>(QStringLiteral("visual.theme"));
    QVERIFY(zoom && slider && theme);
    for (QWidget* w : {static_cast<QWidget*>(zoom), static_cast<QWidget*>(slider), static_cast<QWidget*>(theme)}) {
        wheelOver(w);
    }
    QCOMPARE(text(doc, "visual.zoom"), QStringLiteral("1.0"));
    QCOMPARE(text(doc, "visual.rotation"), QStringLiteral("0.0, 0.0, 0.0"));
    QCOMPARE(text(doc, "visual.theme"), QStringLiteral("light"));
    QVERIFY(!doc.canUndo());
}

void PagesTest::rotationSlidersAndAngleLabel() {
    Shown<StructurePage> s;
    RotationField* r = s.page->rotation();
    QVERIFY(r);
    QCOMPARE(r->angleLabel()->text(), QStringLiteral("X/Y/Z = +0.0°, +0.0°, +0.0°  (+0.0000, +0.0000, +0.0000 rad)"));
    // a stored angle beyond the sliders is kept (AP-01): the slider shows the end, the texts the value
    s.doc.set("visual.rotation", core::Value(core::Array{12.34567, -200.0, 0.0}));
    QCOMPARE(r->slider(0)->value(), 12346);
    QCOMPARE(r->slider(1)->value(), -180000);
    QCOMPARE(r->valueLabel(1)->text(), QStringLiteral("-200.0"));
    QVERIFY(r->angleLabel()->text().startsWith(QStringLiteral("X/Y/Z = +12.3°, -200.0°, +0.0°")));
    QCOMPARE(r->flag().second, QStringLiteral("warning"));
    // moving one slider stores that angle (3 decimals) and keeps the others as stored
    guitest::setField(r->slider(2), 90.5);
    QCOMPARE(text(s.doc, "visual.rotation"), QStringLiteral("12.34567, -200.0, 90.5"));
    // the harness's list form: the three sliders
    guitest::setKey(s.page.get(), QStringLiteral("visual.rotation"), core::Value(core::Array{15.5, -30.0, 120.25}));
    QCOMPARE(text(s.doc, "visual.rotation"), QStringLiteral("15.5, -30.0, 120.25"));
    // a drag is one undo step
    s.doc.clearUndo();
    r->slider(0)->setSliderDown(true);
    for (int v : {16000, 18000, 21500}) r->slider(0)->setValue(v);
    r->slider(0)->setSliderDown(false);
    QCOMPARE(text(s.doc, "visual.rotation"), QStringLiteral("21.5, -30.0, 120.25"));
    QVERIFY(s.doc.canUndo());
    s.doc.undo();
    QVERIFY(!s.doc.canUndo());
    QCOMPARE(text(s.doc, "visual.rotation"), QStringLiteral("15.5, -30.0, 120.25"));
    QCOMPARE(r->slider(0)->value(), 15500);
    // keyboard: one step of 0.5 degrees
    r->slider(1)->setFocus();
    QTest::keyClick(r->slider(1), Qt::Key_Right);
    QCOMPARE(text(s.doc, "visual.rotation"), QStringLiteral("15.5, -29.5, 120.25"));
}

void PagesTest::resetButtons() {
    Shown<StructurePage> s;
    s.doc.set("visual.rotation", core::Value(core::Array{10.0, 20.0, 30.0}));
    s.doc.set("visual.zoom", 2.5);
    s.doc.set("visual.center_on", "B:12");
    QTest::mouseClick(guitest::findField(s.page.get(), QStringLiteral("structure.reset_orientation")), Qt::LeftButton);
    QCOMPARE(text(s.doc, "visual.rotation"), QStringLiteral("0.0, 0.0, 0.0"));
    QTest::mouseClick(guitest::findField(s.page.get(), QStringLiteral("structure.reset_view")), Qt::LeftButton);
    QCOMPARE(text(s.doc, "visual.zoom"), QStringLiteral("1.0"));
    QCOMPARE(text(s.doc, "visual.center_on"), QString());
    s.doc.undo();  // zoom and centre were one step
    QCOMPARE(text(s.doc, "visual.zoom"), QStringLiteral("2.5"));
    QCOMPARE(text(s.doc, "visual.center_on"), QStringLiteral("B:12"));
}

void PagesTest::densityButtonsAndStatus() {
    Shown<StructurePage> s;
    // owner decision D2: the fit is milestone 6's
    for (QPushButton* b : {s.page->globalSearchButton(), s.page->refineButton()}) {
        QVERIFY(!b->isEnabled());
        QVERIFY2(b->toolTip().contains(QStringLiteral("comes with milestone 6")), qPrintable(b->toolTip()));
    }
    QCOMPARE(s.page->globalSearchButton()->objectName(), QStringLiteral("density.global_search"));
    QVERIFY(!s.page->centreOnModelButton()->isEnabled());
    QSignalSpy centre(s.page.get(), &StructurePage::centreOnModelRequested);
    s.page->setCentreOnModelAvailable(true);
    QTest::mouseClick(s.page->centreOnModelButton(), Qt::LeftButton);
    QCOMPARE(centre.count(), 1);
    QSignalSpy mode(s.page.get(), &StructurePage::mapModeChanged);
    guitest::setKey(s.page.get(), QStringLiteral("density.move_with_mouse"), true);
    QCOMPARE(mode.count(), 1);
    QVERIFY(s.page->mapMode());

    QLabel* status = s.page->densityStatus();
    QCOMPARE(status->text(), QStringLiteral("No map selected."));
    s.doc.edit("map", [](core::Project& p) { p.density.map_path = "map.mrc"; });
    QCOMPARE(status->text(), QStringLiteral("Not placed yet: run a global search, or centre it and move it by hand."));
    s.doc.edit("place", [](core::Project& p) {
        p.density.transform = {1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 2, 3};
        p.density.anchor = {{0, 0, 0}};
    });
    QCOMPARE(status->text(), QStringLiteral("Placed (not fitted)."));
    s.doc.edit("info", [](core::Project& p) { p.density.fit_info = "Placed by hand (not fitted) - 'Refine' fits from here."; });
    QCOMPARE(status->text(), QStringLiteral("Placed by hand (not fitted) - 'Refine' fits from here."));
    QTest::mouseClick(s.page->clearButton(), Qt::LeftButton);
    const core::DensityParams& d = s.doc.project().density;
    QVERIFY(d.transform.empty() && d.anchor.empty() && d.fit_info.empty() && d.fit_cc == 0.0);
    QCOMPARE(QString::fromStdString(d.map_path), QStringLiteral("map.mrc"));
    QCOMPARE(status->text(), QStringLiteral("Not placed yet: run a global search, or centre it and move it by hand."));
}

// MDATK _on_density_map_changed: a placement belongs to one map. MDBridge forgets it in the same
// undo step as the new file name.
void PagesTest::anotherMapForgetsThePlacementInOneStep() {
    Shown<StructurePage> s;
    s.doc.edit("placed", [](core::Project& p) {
        p.density.map_path = "a.mrc";
        p.density.transform = {1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0};
        p.density.fit_cc = 0.7;
        p.density.fit_info = "fitted";
    });
    s.doc.clearUndo();
    guitest::typeKey(s.page.get(), QStringLiteral("density.map_path"), QStringLiteral("b.mrc"));
    QTRY_VERIFY(s.doc.project().density.transform.empty());
    const core::DensityParams& d = s.doc.project().density;
    QCOMPARE(QString::fromStdString(d.map_path), QStringLiteral("b.mrc"));
    QVERIFY(d.fit_info.empty() && d.fit_cc == 0.0);
    s.doc.undo();
    QVERIFY(!s.doc.canUndo());  // one step
    QCOMPARE(QString::fromStdString(s.doc.project().density.map_path), QStringLiteral("a.mrc"));
    QCOMPARE(s.doc.project().density.transform.size(), std::size_t{12});
    QCOMPARE(QString::fromStdString(s.doc.project().density.fit_info), QStringLiteral("fitted"));
    // other density settings leave the placement alone
    guitest::setKey(s.page.get(), QStringLiteral("density.level"), 2.5);
    QTest::qWait(20);
    QCOMPARE(s.doc.project().density.transform.size(), std::size_t{12});
}

// AP-04: the formats as set, every format the charts write offered, others kept with a box.
void PagesTest::chartFormatsKeepEveryEntry() {
    Shown<OutputsPage> o;
    ChoiceChecksField* f = o.page->formats();
    QStringList offered;
    for (QCheckBox* b : f->boxes()) offered.append(b->property("mdbridgeChoice").toString());
    QCOMPARE(offered, QStringList({"png", "svg", "pdf", "jpg", "tif", "eps", "ps"}));
    QVERIFY(f->box(QStringLiteral("png"))->isChecked());
    QCOMPARE(f->box(QStringLiteral("png"))->text(), QStringLiteral("PNG"));
    // MDATK's lossy case: ["jpg", "svg"] was saved as ["svg"]
    o.doc.set("output.charts.formats", core::Value(core::Array{"jpg", "svg"}));
    QVERIFY(f->box(QStringLiteral("jpg"))->isChecked() && f->box(QStringLiteral("svg"))->isChecked());
    QVERIFY(!f->box(QStringLiteral("png"))->isChecked());
    QCOMPARE(f->boxes().size(), std::size_t{7});
    // entries no box offers get one of their own; unsupported ones are named and flagged
    core::Project p = o.doc.project();
    p.output.charts.formats = {"png", "bogus", "jpeg", ".PDF"};
    o.doc.replace(p);
    QCOMPARE(f->boxes().size(), std::size_t{9});
    QCOMPARE(f->box(QStringLiteral("bogus"))->text(), QStringLiteral("BOGUS (not supported)"));
    QCOMPARE(f->box(QStringLiteral("jpeg"))->text(), QStringLiteral("JPEG"));
    QVERIFY(f->box(QStringLiteral("pdf"))->isChecked());  // ".PDF" is PDF
    QCOMPARE(f->flag().second, QStringLiteral("warning"));
    QVERIFY(f->flag().first.contains(QStringLiteral("chart format 'bogus' is not supported")));
    // unticking removes the entry (and its spellings), ticking appends; the order stays
    QTest::keyClick(f->box(QStringLiteral("pdf")), Qt::Key_Space);
    QCOMPARE(text(o.doc, "output.charts.formats"), QStringLiteral("png, bogus, jpeg"));
    QTest::keyClick(f->box(QStringLiteral("svg")), Qt::Key_Space);
    QCOMPARE(text(o.doc, "output.charts.formats"), QStringLiteral("png, bogus, jpeg, svg"));
    QTest::keyClick(f->box(QStringLiteral("bogus")), Qt::Key_Space);
    QCOMPARE(text(o.doc, "output.charts.formats"), QStringLiteral("png, jpeg, svg"));
    QVERIFY(f->box(QStringLiteral("bogus")));  // stays, to be ticked again
    QCOMPARE(f->flag().first, QString());
    QTest::keyClick(f->box(QStringLiteral("bogus")), Qt::Key_Space);
    QCOMPARE(text(o.doc, "output.charts.formats"), QStringLiteral("png, jpeg, svg, bogus"));
    // none ticked is stored as such (charts write PNG then, MDATK's fallback)
    guitest::setKey(o.page.get(), QStringLiteral("output.charts.formats"), core::Value(core::Array{}));
    QCOMPARE(text(o.doc, "output.charts.formats"), QString());
    for (int k = 0; k < 4; ++k) o.doc.undo();  // a box at a time: one step each
    QCOMPARE(text(o.doc, "output.charts.formats"), QStringLiteral("png, jpeg, svg, bogus"));
    // another project: its own boxes
    o.doc.replace(core::Project{});
    QCOMPARE(f->boxes().size(), std::size_t{7});
}

void PagesTest::chartStylesAreMdbridgesAndTypedNamesStay() {
    Shown<OutputsPage> o;
    auto* style = o.page->findChild<QComboBox*>(QStringLiteral("output.charts.style"));
    QVERIFY(style && style->isEditable());
    QStringList items;
    for (int k = 0; k < style->count(); ++k) items.append(style->itemText(k));
    QCOMPARE(items, QStringList({"mdatk", "default"}));
    guitest::typeKey(o.page.get(), QStringLiteral("output.charts.style"), QStringLiteral("ggplot"));
    QCOMPARE(text(o.doc, "output.charts.style"), QStringLiteral("ggplot"));
    QVERIFY(!o.page->exportPdbsButton()->isEnabled());
    o.page->setJobActionsAvailable(true);
    QSignalSpy exported(o.page.get(), &OutputsPage::exportPdbsRequested);
    QTest::mouseClick(o.page->exportPdbsButton(), Qt::LeftButton);
    QCOMPARE(exported.count(), 1);
}

void PagesTest::revealFieldOpensTheSection() {
    Shown<StructurePage> s;
    auto* appearance = s.page->section(QStringLiteral("appearance"));
    appearance->setExpanded(false);
    QVERIFY(s.page->revealField(QStringLiteral("visual.glow_color")));
    QVERIFY(appearance->isExpanded());
    QTRY_VERIFY(QApplication::focusWidget() &&
                QApplication::focusWidget()->objectName() == QStringLiteral("visual.glow_color"));
    QVERIFY(!s.page->revealField(QStringLiteral("analyses.contacts.cutoff")));
    QVERIFY(!s.page->revealField(QStringLiteral("output.video.fps")));
    s.doc.insertItem("visual.representations", 0, s.page->representations()->newItem());
    QVERIFY(s.page->revealField(QStringLiteral("visual.representations.0.size")));
    QCOMPARE(s.page->representations()->view()->currentIndex(),
             s.page->representations()->cellIndex("visual.representations.0.size"));
    QVERIFY(s.page->revealField(QStringLiteral("density.transform")));  // a placement problem: the map's fields
    Shown<OutputsPage> o;
    QVERIFY(app::pages::revealFieldOn(o.page.get(), QStringLiteral("output.video.fps")));
    QTRY_VERIFY(QApplication::focusWidget() && QApplication::focusWidget()->objectName() == QStringLiteral("output.video.fps"));
}

// A screenshot of each page, whole, with the demo-like full project (CI artifacts, never compared).
void PagesTest::screenshots() {
    core::ProjectDocument doc;
    doc.open(kData / "projects" / "full.mdatk.toml");
    app::forms::FormContext context(doc);
    const fs::path dir = guitest::screenshotDir(core::pathFromUtf8(MDBRIDGE_GUI_SCREENSHOT_DIR));
    for (QWidget* raw : app::pages::createPages(context)) {
        std::unique_ptr<QWidget> page(raw);
        page->resize(540, 900);
        page->show();
        QVERIFY(QTest::qWaitForWindowExposed(page.get()));
        auto* form = qobject_cast<app::pages::FormPage*>(page.get());
        QVERIFY(form);
        const fs::path shot = guitest::saveScreenshot(form->scrollArea()->widget(), dir,
                                                      "page-" + app::pages::pageIdOf(page.get()).toStdString());
        qInfo("screenshot %s", core::pathToUtf8(shot).c_str());
        QVERIFY(fs::file_size(shot) > 1000);
    }
}

int main(int argc, char** argv) {
    const QByteArray platform = qgetenv("MDBRIDGE_GUI_TEST_PLATFORM");
    qputenv("QT_QPA_PLATFORM", platform.isEmpty() ? QByteArray("offscreen") : platform);
    QLocale::setDefault(QLocale::c());
    QApplication app(argc, argv);
    PagesTest tests;
    return QTest::qExec(&tests, argc, argv);
}

#include "test_pages.moc"
