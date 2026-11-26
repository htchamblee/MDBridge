// The settings pages Input & cleaning, Residue groups and Analyses (milestone 5 W1-B (src/app/gui/pages): Input & cleaning, Residue groups,
// Analyses. Qt Test on the offscreen platform: MDATK's order, sections and fields; the analysis
// group drop-downs follow the groups (MDATK test_panel_group_dropdowns_and_auto_title, the group
// part); renaming a group renames its uses (MDATK 0.1.1 tests/test_groups_gui.py); removing a used
// group says where; Pick colour…; Groups from chains (D15, MDATK 0.1.1 chaingroups: the expected
// groups are what MDATK 0.1.1 computes for the same files, checked once with its own code); reveal;
// a screenshot of each page (CI artifacts, never compared).

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include <QApplication>
#include <QCheckBox>
#include <QComboBox>
#include <QRegularExpression>
#include <QSet>
#include <QToolButton>
#include <QLabel>
#include <QLineEdit>
#include <QLocale>
#include <QPushButton>
#include <QSignalSpy>
#include <QTableView>
#include <QTemporaryDir>
#include <QTest>

#include "chain_groups.hpp"
#include "gui/forms/field_widgets.hpp"
#include "gui/forms/form_context.hpp"
#include "gui/forms/list_editor.hpp"
#include "gui/pages/analyses_page.hpp"
#include "gui/pages/groups_page.hpp"
#include "gui/pages/input_page.hpp"
#include "gui/pages/form_page.hpp"
#include "gui/pages/page_registry.hpp"
#include "gui_harness.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/project_document.hpp"
#include "page_host.hpp"
#include <QTabWidget>

namespace fs = std::filesystem;
using namespace mdbridge;
using namespace mdbridge::app;

namespace {

fs::path screenshots() { return guitest::screenshotDir(core::pathFromUtf8(MDBRIDGE_GUI_SCREENSHOT_DIR)); }

// Two RNA strands of four residues (chains A and B, both numbered 1-4 unless `renumber` is false:
// then 1-4 and 5-8), consecutive residues 1.6 Å apart (O3'-P), a water and a Mg2+ (chains W, M;
// residue 1); no CONECT records. Written as MDATK 0.1.1's tests/test_groups.py write_pdb does.
std::string chainsPdb(bool chainIds, bool renumber = true) {
    struct Atom {
        char chain;
        int number;
        std::string resname, name, element;
        double x, y, z;
    };
    std::vector<Atom> atoms;
    const char* chains[] = {"A", "B"};
    for (int c = 0; c < 2; ++c) {
        for (int k = 0; k < 4; ++k) {
            const int number = renumber ? k + 1 : 4 * c + k + 1;
            const std::string resname(1, "GCGC"[k]);
            const double x = 20.0 * c;
            atoms.push_back({chains[c][0], number, resname, "P", "P", x, 0.0, 4.6 * k});
            atoms.push_back({chains[c][0], number, resname, "C4'", "C", x, 0.0, 4.6 * k + 1.5});
            atoms.push_back({chains[c][0], number, resname, "O3'", "O", x, 0.0, 4.6 * k + 3.0});
        }
    }
    atoms.push_back({'W', 1, "HOH", "O", "O", 40, 0, 0});
    atoms.push_back({'W', 1, "HOH", "H1", "H", 40.96, 0, 0});
    atoms.push_back({'W', 1, "HOH", "H2", "H", 39.76, 0.93, 0});
    atoms.push_back({'M', 1, "MG", "MG", "MG", 50, 0, 0});
    std::string out;
    for (std::size_t k = 0; k < atoms.size(); ++k) {
        const Atom& a = atoms[k];
        const char* record = a.resname == "HOH" || a.resname == "MG" ? "HETATM" : "ATOM  ";
        char name[8];
        std::snprintf(name, sizeof name, a.name.size() < 4 ? " %-3s" : "%s", a.name.c_str());
        char line[128];
        std::snprintf(line, sizeof line, "%s%5d %s %3s %c%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s\n", record,
                      static_cast<int>(k + 1), name, a.resname.c_str(), chainIds ? a.chain : ' ', a.number, a.x, a.y,
                      a.z, a.element.c_str());
        out += line;
    }
    return out + "END\n";
}

fs::path writeFile(const fs::path& path, const std::string& text) {
    std::ofstream(path, std::ios::binary) << text;
    return path;
}

template <class Page>
Page* findPage(guitest::PageHost& host) {
    for (QWidget* w : host.pages()) {
        if (auto* p = qobject_cast<Page*>(w)) return p;
    }
    return nullptr;
}

QStringList comboItems(QComboBox* box) {
    QStringList out;
    for (int k = 0; k < box->count(); ++k) out.append(box->itemText(k));
    return out;
}

std::string text(core::ProjectDocument& doc, const char* path) { return doc.text(path); }

// Shows the page of `key` and reveals it there, as the window's side bar will.
void reveal(guitest::PageHost& host, const QString& key) {
    for (QWidget* page : host.pages()) {
        if (pages::revealFieldOn(page, key)) {
            host.tabs()->setCurrentWidget(page);
            pages::revealFieldOn(page, key);
            return;
        }
    }
}

}  // namespace

class GroupPagesTests : public QObject {
    Q_OBJECT

private slots:
    void pagesInMdatksOrder();
    void sectionsAndFields();
    void groupDropdownsFollowTheGroups();
    void renamingAGroupRenamesItEverywhere();
    void badNamesAreRefused();
    void removingAUsedGroupSaysWhere();
    void pickColourWritesHex();
    void chainGroupsAsMdatkLatest();
    void groupsFromChainsButton();
    void choosingATopologyReplacesUntouchedDefaults();
    void loadingAProjectDropsPendingReads();
    void revealFindsFieldsCellsAndRows();
    void screenshotsOfEachPage();
};

void GroupPagesTests::pagesInMdatksOrder() {
    core::ProjectDocument doc;
    forms::FormContext context(doc);
    std::vector<QWidget*> pages = pages::createPages(context);
    QCOMPARE(pages.size(), std::size_t{6});  // layout: a placeholder until W2-A
    QVERIFY(pages::isPlaceholder(pages[4]));
    const char* ids[] = {"input", "groups", "analyses"};
    const char* titles[] = {"Input & cleaning", "Residue groups", "Analyses"};
    for (std::size_t k = 0; k < 3; ++k) {
        QString id;
        QString title;
        QVERIFY(QMetaObject::invokeMethod(pages[k], "pageId", Q_RETURN_ARG(QString, id)));
        QVERIFY(QMetaObject::invokeMethod(pages[k], "title", Q_RETURN_ARG(QString, title)));
        QCOMPARE(id, QString::fromLatin1(ids[k]));
        QCOMPARE(title, QString::fromUtf8(titles[k]));
        QCOMPARE(QString::fromUtf8(core::kProblemPages[k].title.data(), static_cast<qsizetype>(core::kProblemPages[k].title.size())), title);
    }
    for (QWidget* p : pages) delete p;
}

void GroupPagesTests::sectionsAndFields() {
    guitest::PageHost host;
    const auto sectionTitles = [](pages::FormPage* page) {
        QStringList out;
        for (pages::Section* s : page->sections()) {
            out.append(s->objectName().startsWith(QStringLiteral("section:%1.").arg(page->pageId()))
                           ? s->title()
                           : QStringLiteral("(misnamed) ") + s->objectName());
        }
        return out;
    };
    auto* input = findPage<pages::InputPage>(host);
    auto* groups = findPage<pages::GroupsPage>(host);
    auto* analyses = findPage<pages::AnalysesPage>(host);
    QVERIFY(input && groups && analyses);
    QCOMPARE(sectionTitles(input), (QStringList{"Project notes", "Input files", "Frames, time and external data",
                                                "Trajectory cleaning (cpptraj)"}));
    QCOMPARE(sectionTitles(groups), QStringList{"Residue groups"});
    QCOMPARE(sectionTitles(analyses),
             (QStringList{"General", "RMSD", "Pairwise contacts", "Hydrogen bonds", "O–Mg–O contacts and Mg occupancy",
                          "π–π stacking", "Base pairs (secondary structure)", "Intermolecular energy"}));
    // every field MDATK's pages bind, in its order, once
    const std::vector<std::pair<pages::FormPage*, std::vector<std::string>>> fields{
        {input,
         {"notes", "inputs.topology", "inputs.trajectories", "inputs.start", "inputs.stop", "inputs.step", "inputs.dt_ps",
          "inputs.time_unit", "inputs.external", "cleaning.enabled", "cleaning.cpptraj_path", "cleaning.conda_env",
          "cleaning.conda_exe", "cleaning.solvent_mask", "cleaning.ion_mask", "cleaning.extra_strip_mask",
          "cleaning.autoimage", "cleaning.autoimage_anchor", "cleaning.fit", "cleaning.fit_mask",
          "cleaning.fit_reference", "cleaning.reference_frame", "cleaning.reference_file",
          "cleaning.remove_box_after_fit", "cleaning.start", "cleaning.stop", "cleaning.offset",
          "cleaning.extra_commands"}},
        {groups, {"groups"}},
        {analyses,
         {"analyses.use_pbc", "analyses.proximity_search_radius", "analyses.rmsd.enabled", "analyses.rmsd.targets",
          "analyses.rmsd.reference_frame", "analyses.rmsd.superposition", "analyses.rmsd.mass_weighted",
          "analyses.contacts.enabled", "analyses.contacts.group_a", "analyses.contacts.group_b",
          "analyses.contacts.cutoff", "analyses.contacts.heavy_only", "analyses.hbonds.enabled",
          "analyses.hbonds.group_a", "analyses.hbonds.group_b", "analyses.hbonds.distance_cutoff",
          "analyses.hbonds.use_angle", "analyses.hbonds.angle_cutoff", "analyses.hbonds.donor_elements",
          "analyses.hbonds.acceptor_elements", "analyses.mg.enabled", "analyses.mg.mg_selection",
          "analyses.mg.ligand_selection", "analyses.mg.ligand_elements", "analyses.mg.group_a", "analyses.mg.group_b",
          "analyses.mg.inner_cutoff", "analyses.mg.outer_cutoff", "analyses.mg.bridge_shell", "analyses.mg.geometry",
          "analyses.mg.angle_tolerance", "analyses.mg.occupancy_selection", "analyses.mg.occupancy_elements",
          "analyses.stacking.enabled", "analyses.stacking.selection", "analyses.stacking.scope",
          "analyses.stacking.group_a", "analyses.stacking.group_b", "analyses.stacking.distance_cutoff",
          "analyses.stacking.angle_cutoff", "analyses.stacking.offset_cutoff", "analyses.stacking.planarity_cutoff",
          "analyses.basepairs.engine", "analyses.basepairs.stride", "analyses.basepairs.cutoff",
          "analyses.basepairs.wobble", "analyses.basepairs.noncanonical", "analyses.energy.enabled",
          "analyses.energy.group_a", "analyses.energy.group_b", "analyses.energy.cutoff", "analyses.energy.dielectric",
          "analyses.energy.distance_dependent", "analyses.energy.use_parm7", "analyses.energy.exclude_bonded",
          "analyses.energy.target_value"}},
    };
    for (const auto& [page, paths] : fields) {
        std::vector<std::string> shown;
        for (forms::FieldWidget* f : page->findChildren<forms::FieldWidget*>()) {
            if (!qobject_cast<forms::ListEditor*>(f->parentWidget()) && f->path().find(':') == std::string::npos) {
                shown.push_back(f->path());
            }
        }
        QCOMPARE(shown, paths);
    }
    // MDATK's check box texts and tooltips
    auto* clean = qobject_cast<QCheckBox*>(guitest::findField(input, QStringLiteral("cleaning.enabled")));
    QVERIFY(clean);
    QCOMPARE(clean->text(), QStringLiteral("Strip and clean the trajectory with cpptraj before analysis"));
    QCOMPARE(guitest::findField(input, QStringLiteral("inputs.dt_ps"))->toolTip(),
             QStringLiteral("Leave blank to use the trajectory's time stamps"));
    QCOMPARE(guitest::findField(analyses, QStringLiteral("analyses.stacking.distance_cutoff"))->toolTip().left(30),
             QStringLiteral("Ring-centroid distance and ang"));
    QVERIFY(guitest::findField(input, QStringLiteral("cleaning.detect")));
    QVERIFY(guitest::findField(input, QStringLiteral("cleaning.run")));
    QVERIFY(guitest::findField(input, QStringLiteral("cleaning.inspect")));
    // a section folds and unfolds
    pages::Section* rmsd = analyses->section(QStringLiteral("rmsd"));
    rmsd->header()->click();
    QVERIFY(!rmsd->isExpanded() && rmsd->form()->isHidden());
    rmsd->header()->click();
    QVERIFY(rmsd->isExpanded());
}

// MDATK's _groups_changed: every group drop-down lists the current group names (the text is kept).
void GroupPagesTests::groupDropdownsFollowTheGroups() {
    guitest::PageHost host;
    host.show();
    QVERIFY(QTest::qWaitForWindowExposed(&host));
    auto* analyses = findPage<pages::AnalysesPage>(host);
    core::ProjectDocument& doc = host.document();
    QCOMPARE(analyses->groupCombos().size(), std::size_t{10});
    for (forms::EditableComboField* c : analyses->groupCombos()) {
        QCOMPARE(comboItems(c->comboBox()), (QStringList{"Aptamer", "Target"}));
    }
    reveal(host, QStringLiteral("groups"));
    guitest::tableAction(&host, QStringLiteral("groups"), "add", -1, 0);
    guitest::setKey(&host, QStringLiteral("groups.2.name"), "Loop");
    for (forms::EditableComboField* c : analyses->groupCombos()) {
        QCOMPARE(comboItems(c->comboBox()), (QStringList{"Aptamer", "Target", "Loop"}));
    }
    QCOMPARE(analyses->groupCombos()[0]->comboBox()->currentText(), QStringLiteral("Aptamer"));
    reveal(host, QStringLiteral("analyses.contacts.group_b"));
    guitest::setKey(&host, QStringLiteral("analyses.contacts.group_b"), "Loop");
    QCOMPARE(QString::fromStdString(text(doc, "analyses.contacts.group_b")), QStringLiteral("Loop"));
    guitest::tableAction(&host, QStringLiteral("groups"), "remove", 0, 0);
    for (forms::EditableComboField* c : analyses->groupCombos()) {
        QCOMPARE(comboItems(c->comboBox()), (QStringList{"Target", "Loop"}));
    }
    // a group that is gone stays named where it was used (validation flags it)
    QCOMPARE(analyses->groupCombos()[0]->comboBox()->currentText(), QStringLiteral("Aptamer"));
    doc.undo();
    QCOMPARE(comboItems(analyses->groupCombos()[3]->comboBox()), (QStringList{"Aptamer", "Target", "Loop"}));
    doc.newProject();
    QCOMPARE(comboItems(analyses->groupCombos()[3]->comboBox()), (QStringList{"Aptamer", "Target"}));
}

// MDATK 0.1.1 test_renaming_a_group_renames_it_everywhere (the settings MDBridge's pages have so far).
void GroupPagesTests::renamingAGroupRenamesItEverywhere() {
    guitest::PageHost host;
    host.show();
    QVERIFY(QTest::qWaitForWindowExposed(&host));
    core::ProjectDocument& doc = host.document();
    auto* groups = findPage<pages::GroupsPage>(host);
    auto* analyses = findPage<pages::AnalysesPage>(host);
    doc.setText("visual.structure_groups", "Aptamer, Target");
    doc.setText("output.pdb.selection", "Aptamer, Target");
    reveal(host, QStringLiteral("groups.0.name"));
    guitest::setKey(&host, QStringLiteral("groups.0.name"), "Binder");  // the cell, Enter
    QCOMPARE(QString::fromStdString(text(doc, "groups.0.name")), QStringLiteral("Binder"));
    QVERIFY2(QRegularExpression(QStringLiteral("^Renamed Aptamer → Binder in \\d+ places$"))
                 .match(groups->statusText())
                 .hasMatch(),
             qPrintable(groups->statusText()));
    QVERIFY(!groups->statusIsError());
    QVERIFY(groups->statusDetail().contains(QStringLiteral("analyses.contacts.group_a")));
    for (const char* p : {"analyses.contacts.group_a", "analyses.hbonds.group_a", "analyses.mg.group_a",
                          "analyses.stacking.group_a", "analyses.energy.group_a"}) {
        QCOMPARE(QString::fromStdString(text(doc, p)), QStringLiteral("Binder"));
    }
    QCOMPARE(QString::fromStdString(text(doc, "visual.structure_groups")), QStringLiteral("Binder, Target"));
    QCOMPARE(QString::fromStdString(text(doc, "output.pdb.selection")), QStringLiteral("Binder, Target"));
    QVERIFY(QString::fromStdString(text(doc, "analyses.rmsd.targets.0.selection"))
                .startsWith(QStringLiteral("(group:Binder or group:Target)")));
    QCOMPARE(analyses->groupCombos()[0]->comboBox()->currentText(), QStringLiteral("Binder"));
    QCOMPARE(comboItems(analyses->groupCombos()[0]->comboBox()), (QStringList{"Binder", "Target"}));
    // one undo step: undone, every use comes back
    doc.undo();
    QCOMPARE(QString::fromStdString(text(doc, "groups.0.name")), QStringLiteral("Aptamer"));
    QCOMPARE(QString::fromStdString(text(doc, "analyses.contacts.group_a")), QStringLiteral("Aptamer"));
    // a new row's first name has no uses to follow
    guitest::tableAction(&host, QStringLiteral("groups"), "add", -1, 0);
    guitest::setKey(&host, QStringLiteral("groups.2.name"), "Loop");
    QCOMPARE(groups->statusText(), QStringLiteral("Renamed NewGroup → Loop in 0 places"));
}

// MDATK 0.1.1 test_group_names_are_checked_and_leaving_the_cell_commits (the refusals).
void GroupPagesTests::badNamesAreRefused() {
    guitest::PageHost host;
    host.show();
    QVERIFY(QTest::qWaitForWindowExposed(&host));
    core::ProjectDocument& doc = host.document();
    auto* groups = findPage<pages::GroupsPage>(host);
    reveal(host, QStringLiteral("groups.1.name"));
    guitest::setKey(&host, QStringLiteral("groups.1.name"), "two words");
    QCOMPARE(QString::fromStdString(text(doc, "groups.1.name")), QStringLiteral("Target"));
    QVERIFY(groups->statusText().startsWith(QStringLiteral("Group Target not renamed: ")));
    QVERIFY(groups->statusText().contains(QStringLiteral("no spaces")));
    QVERIFY(groups->statusIsError());
    guitest::setKey(&host, QStringLiteral("groups.1.name"), "Aptamer");
    QCOMPARE(QString::fromStdString(text(doc, "groups.1.name")), QStringLiteral("Target"));
    QVERIFY(groups->statusText().contains(QStringLiteral("already")));
    QCOMPARE(QString::fromStdString(text(doc, "analyses.contacts.group_b")), QStringLiteral("Target"));
    QVERIFY(!doc.isModified());
    guitest::setKey(&host, QStringLiteral("groups.1.name"), "Ligand");
    QCOMPARE(QString::fromStdString(text(doc, "analyses.contacts.group_b")), QStringLiteral("Ligand"));
    QCOMPARE(QString::fromStdString(text(doc, "visual.proximity_set")), QStringLiteral("Ligand"));
    QVERIFY(groups->statusText().startsWith(QStringLiteral("Renamed Target → Ligand")));
}

// MDATK 0.1.1 test_deleting_a_used_group_says_where_it_is_used.
void GroupPagesTests::removingAUsedGroupSaysWhere() {
    guitest::PageHost host;
    host.show();
    QVERIFY(QTest::qWaitForWindowExposed(&host));
    auto* groups = findPage<pages::GroupsPage>(host);
    QSignalSpy logged(groups, &pages::GroupsPage::logMessage);
    reveal(host, QStringLiteral("groups"));
    guitest::tableAction(&host, QStringLiteral("groups"), "remove", 1, 0);
    QVERIFY(groups->statusText().startsWith(QStringLiteral("Removed group still in use")));
    for (const char* where : {"analyses.contacts.group_b", "visual.proximity_set", "analyses.rmsd.targets.0.selection"}) {
        QVERIFY2(groups->statusText().contains(QLatin1String(where)), where);
    }
    QVERIFY(groups->statusIsError());
    QVERIFY(!logged.isEmpty());
    QVERIFY(logged.last().at(0).toString().startsWith(QStringLiteral("WARNING: Removed group still in use")));
    // an unused group goes silently
    guitest::tableAction(&host, QStringLiteral("groups"), "add", -1, 0);
    const QString before = groups->statusText();
    guitest::tableAction(&host, QStringLiteral("groups"), "remove", 1, 0);
    QCOMPARE(groups->statusText(), before);
}

void GroupPagesTests::pickColourWritesHex() {
    guitest::PageHost host;
    host.show();
    QVERIFY(QTest::qWaitForWindowExposed(&host));
    core::ProjectDocument& doc = host.document();
    auto* groups = findPage<pages::GroupsPage>(host);
    QColor offered;
    const auto previous = forms::setColorDialogFunction([&](QWidget*, const QColor& start, const QString& title) {
        offered = start;
        return title == QStringLiteral("Group colour") ? QColor(0xaa, 0x22, 0xcc) : QColor();
    });
    groups->groupEditor()->view()->clearSelection();
    groups->groupEditor()->view()->setCurrentIndex({});
    QTest::mouseClick(groups->pickColourButton(), Qt::LeftButton);
    QCOMPARE(groups->statusText(), QStringLiteral("Select a group first."));
    groups->groupEditor()->selectRow(1);
    QTest::mouseClick(groups->pickColourButton(), Qt::LeftButton);
    QCOMPARE(offered, QColor(0xeb, 0x68, 0x34));
    QCOMPARE(QString::fromStdString(text(doc, "groups.1.color")), QStringLiteral("#aa22cc"));
    // the swatch of the cell shows it
    const QModelIndex cell = groups->groupEditor()->cellIndex("groups.1.color");
    QVERIFY(!cell.data(Qt::DecorationRole).value<QIcon>().isNull());
    forms::setColorDialogFunction(previous);
}

// chain_groups.cpp against MDATK 0.1.1's chaingroups.topology_chain_groups / merge_chain_groups of
// the same files (run once with MDATK 0.1.1 at the commit regenerate_all.py names).
void GroupPagesTests::chainGroupsAsMdatkLatest() {
    QTemporaryDir tmp;
    const fs::path dir = core::pathFromUtf8(tmp.path().toStdString());
    using Row = std::tuple<std::string, std::string, std::string, std::string, std::string, std::size_t>;
    const auto rows = [](const std::vector<ChainGroup>& chains) {
        std::vector<Row> out;
        for (const ChainGroup& c : chains) {
            out.emplace_back(c.group.name, c.group.residues, c.group.selection, c.group.color, c.source, c.residues);
        }
        return out;
    };
    const core::Project defaults;
    const auto ids = topologyChainGroups(writeFile(dir / "c11.pdb", chainsPdb(true)), defaults.groups);
    QCOMPARE(rows(ids), (std::vector<Row>{{"Chain_A", "A:1-4", "", "#2f9e6e", "chainID", 4},
                                          {"Chain_B", "B:1-4", "", "#c2408f", "chainID", 4}}));
    const auto noIds = topologyChainGroups(writeFile(dir / "c01.pdb", chainsPdb(false)), defaults.groups);
    QCOMPARE(rows(noIds), (std::vector<Row>{{"Chain_1", "", "index 0:11", "#2f9e6e", "fragment", 4},
                                            {"Chain_2", "", "index 12:23", "#c2408f", "fragment", 4}}));
    const auto numbered = topologyChainGroups(writeFile(dir / "c00.pdb", chainsPdb(false, false)), defaults.groups);
    QCOMPARE(rows(numbered), (std::vector<Row>{{"Chain_1", "", "index 0:11", "#2f9e6e", "fragment", 4},
                                               {"Chain_2", "5-8", "", "#c2408f", "fragment", 4}}));

    core::Project p;
    ChainMerge merge = mergeChainGroups(p, ids, true);
    QCOMPARE(QString::fromStdString(merge.summary("x.pdb")),
             QStringLiteral("Default groups replaced by the chains of x.pdb (Aptamer → Chain_A, Target → Chain_B; 15 "
                            "reference(s) updated)."));
    QCOMPARE(p.groups.size(), std::size_t{2});
    QCOMPARE(p.groups[0].name, std::string("Chain_A"));
    QCOMPARE(p.groups[0].residues, std::string("A:1-4"));
    QCOMPARE(p.groups[0].color, std::string("#2a78d6"));  // the default's colour
    QCOMPARE(p.groups[1].color, std::string("#eb6834"));
    QCOMPARE(p.analyses.hbonds.group_a, std::string("Chain_A"));
    QCOMPARE(p.analyses.mg.ligand_selection, std::string("group:Chain_A or group:Chain_B"));
    core::Project q;
    merge = mergeChainGroups(q, ids, false);  // the button keeps the defaults
    QCOMPARE(QString::fromStdString(merge.summary("x.pdb")),
             QStringLiteral("Added 2 group(s) from the chains of x.pdb: Chain_A, Chain_B."));
    QCOMPARE(core::groupNames(q), (std::vector<std::string>{"Aptamer", "Target", "Chain_A", "Chain_B"}));
    QCOMPARE(q.analyses.hbonds.group_a, std::string("Aptamer"));
    // again: every chain has a group now
    merge = mergeChainGroups(q, topologyChainGroups(dir / "c11.pdb", q.groups), false);
    QCOMPARE(QString::fromStdString(merge.summary("x.pdb")), QStringLiteral("Every chain of x.pdb already has a group."));
    // user-made groups are never replaced; the chains are only named
    core::Project r;
    r.groups[0].residues = "A:1-2";
    merge = mergeChainGroups(r, topologyChainGroups(dir / "c11.pdb", r.groups), true);
    QVERIFY(!merge.changed());
    QCOMPARE(merge.pending, (std::vector<std::string>{"Chain_A", "Chain_B"}));
    QVERIFY(QString::fromStdString(merge.summary("x.pdb")).endsWith(QStringLiteral("'Groups from chains' adds them.")));
    // names and colours avoid existing groups
    const std::vector<core::ResidueGroup> existing{{"Chain_A", "B:1-2", "", "#2F9E6E"}, {"Strand", "B:1-4", "", "#c2408f"}};
    const auto avoid = topologyChainGroups(dir / "c11.pdb", existing);
    QCOMPARE(avoid.size(), std::size_t{2});
    QCOMPARE(avoid[0].group.name, std::string("Chain_A_2"));
    QCOMPARE(avoid[1].coveredBy, std::string("Strand"));
    QCOMPARE(avoid[0].group.color, std::string("#2a78d6"));
    // polymer mask: water and ions left out
    QVERIFY(untouchedDefaults(core::Project().groups));
}

// MDATK 0.1.1 test_groups_from_chains: the button adds the chains, once.
void GroupPagesTests::groupsFromChainsButton() {
    QTemporaryDir tmp;
    const fs::path dir = core::pathFromUtf8(tmp.path().toStdString());
    writeFile(dir / "two.pdb", chainsPdb(true));
    guitest::PageHost host;
    host.show();
    QVERIFY(QTest::qWaitForWindowExposed(&host));
    core::ProjectDocument& doc = host.document();
    auto* groups = findPage<pages::GroupsPage>(host);
    QTest::mouseClick(groups->chainsButton(), Qt::LeftButton);
    QCOMPARE(groups->statusText(), QStringLiteral("Choose a topology file on the Input page first."));
    QVERIFY(groups->statusIsError());
    doc.saveAs(dir / "p.mdatk.toml");
    doc.setText("groups.1.residues", "13-23");  // no longer the defaults
    doc.setText("inputs.topology", "two.pdb");  // relative to the project's folder
    QTRY_VERIFY_WITH_TIMEOUT(groups->statusText().contains(QStringLiteral("without a group")), 10000);
    QCOMPARE(core::groupNames(doc.project()), (std::vector<std::string>{"Aptamer", "Target"}));
    QTest::mouseClick(groups->chainsButton(), Qt::LeftButton);
    QCOMPARE(groups->statusText(), QStringLiteral("Reading the chains of two.pdb …"));
    QTRY_VERIFY_WITH_TIMEOUT(groups->statusText().startsWith(QStringLiteral("Added")), 10000);
    QCOMPARE(core::groupNames(doc.project()), (std::vector<std::string>{"Aptamer", "Target", "Chain_A", "Chain_B"}));
    QCOMPARE(QString::fromStdString(text(doc, "analyses.contacts.group_a")), QStringLiteral("Aptamer"));
    QCOMPARE(groups->groupEditor()->rowCount(), 4);
    QTest::mouseClick(groups->chainsButton(), Qt::LeftButton);
    QTRY_VERIFY_WITH_TIMEOUT(groups->statusText().contains(QStringLiteral("already has a group")), 10000);
    QCOMPARE(doc.project().groups.size(), std::size_t{4});
    doc.undo();  // the chains came as one step
    QCOMPARE(doc.project().groups.size(), std::size_t{2});
}

// MDATK 0.1.1 test_choosing_a_topology_replaces_untouched_defaults.
void GroupPagesTests::choosingATopologyReplacesUntouchedDefaults() {
    QTemporaryDir tmp;
    const fs::path dir = core::pathFromUtf8(tmp.path().toStdString());
    const fs::path pdb = writeFile(dir / "two.pdb", chainsPdb(true));
    guitest::PageHost host;
    host.show();
    QVERIFY(QTest::qWaitForWindowExposed(&host));
    core::ProjectDocument& doc = host.document();
    auto* groups = findPage<pages::GroupsPage>(host);
    auto* analyses = findPage<pages::AnalysesPage>(host);
    reveal(host, QStringLiteral("inputs.topology"));
    guitest::typeKey(&host, QStringLiteral("inputs.topology"), QString::fromStdString(core::pathToUtf8(pdb)));
    QTRY_VERIFY_WITH_TIMEOUT(groups->statusText().startsWith(QStringLiteral("Default groups replaced")), 10000);
    QCOMPARE(core::groupNames(doc.project()), (std::vector<std::string>{"Chain_A", "Chain_B"}));
    QCOMPARE(QString::fromStdString(text(doc, "analyses.hbonds.group_a")), QStringLiteral("Chain_A"));
    QCOMPARE(QString::fromStdString(text(doc, "analyses.hbonds.group_b")), QStringLiteral("Chain_B"));
    QCOMPARE(QString::fromStdString(text(doc, "analyses.mg.ligand_selection")),
             QStringLiteral("group:Chain_A or group:Chain_B"));
    QSet<QString> combos;
    for (forms::EditableComboField* c : analyses->groupCombos()) combos.insert(c->comboBox()->currentText());
    QCOMPARE(combos, (QSet<QString>{"Chain_A", "Chain_B"}));
    QVERIFY(groups->statusDetail().contains(QStringLiteral("analyses.contacts.group_a")));
    // the same topology again offers nothing; a new project offers it again
    doc.newProject();
    QCOMPARE(groups->statusText(), QString());
    doc.setText("inputs.topology", core::pathToUtf8(pdb));
    QTRY_VERIFY_WITH_TIMEOUT(groups->statusText().startsWith(QStringLiteral("Default groups replaced")), 10000);
}

// MDATK 0.1.1 test_loading_a_project_drops_pending_chain_reads: a project opened meanwhile keeps
// its groups.
void GroupPagesTests::loadingAProjectDropsPendingReads() {
    QTemporaryDir tmp;
    const fs::path dir = core::pathFromUtf8(tmp.path().toStdString());
    const fs::path pdb = writeFile(dir / "two.pdb", chainsPdb(true));
    guitest::PageHost host;
    core::ProjectDocument& doc = host.document();
    auto* groups = findPage<pages::GroupsPage>(host);
    {
        core::ProjectDocument other;
        other.setText("inputs.topology", core::pathToUtf8(pdb));
        other.saveAs(dir / "old.mdatk.toml");
    }
    doc.setText("inputs.topology", core::pathToUtf8(pdb));  // the 500 ms timer runs ...
    doc.open(dir / "old.mdatk.toml");                      // ... when a project opens: loading never fills groups
    QTest::qWait(pages::GroupsPage::kTopologySettleMs + 300);
    QVERIFY(!groups->chainsPending());
    QCOMPARE(core::groupNames(doc.project()), (std::vector<std::string>{"Aptamer", "Target"}));
    doc.setText("groups.1.residues", "13-23");
    groups->groupsFromChains();  // read for this project ...
    doc.open(dir / "old.mdatk.toml");  // ... answered after another one opened: dropped
    QTRY_VERIFY_WITH_TIMEOUT(!groups->chainsPending(), 10000);
    QCOMPARE(core::groupNames(doc.project()), (std::vector<std::string>{"Aptamer", "Target"}));
    QCOMPARE(groups->statusText(), QString());
}

void GroupPagesTests::revealFindsFieldsCellsAndRows() {
    guitest::PageHost host;
    host.show();
    QVERIFY(QTest::qWaitForWindowExposed(&host));
    auto* input = findPage<pages::InputPage>(host);
    auto* groups = findPage<pages::GroupsPage>(host);
    auto* analyses = findPage<pages::AnalysesPage>(host);
    host.tabs()->setCurrentWidget(analyses);
    pages::Section* energy = analyses->section(QStringLiteral("energy"));
    energy->setExpanded(false);
    QVERIFY(analyses->revealField(QStringLiteral("analyses.energy.target_value")));
    QVERIFY(energy->isExpanded());
    QCOMPARE(QApplication::focusWidget(), guitest::findField(analyses, QStringLiteral("analyses.energy.target_value")));
    QVERIFY(!analyses->revealField(QStringLiteral("inputs.stop")));
    host.tabs()->setCurrentWidget(groups);
    QVERIFY(groups->revealField(QStringLiteral("groups.1.residues")));
    QCOMPARE(groups->groupEditor()->view()->currentIndex().row(), 1);
    QVERIFY(groups->revealField(QStringLiteral("groups.0")));
    QCOMPARE(groups->groupEditor()->view()->currentIndex().row(), 0);
    host.tabs()->setCurrentWidget(input);
    QVERIFY(input->revealField(QStringLiteral("cleaning.extra_commands")));
    QVERIFY(input->revealField(QStringLiteral("inputs.external.0.path")));
    QVERIFY(!input->revealField(QStringLiteral("visual.zoom")));
}

void GroupPagesTests::screenshotsOfEachPage() {
    const fs::path data = core::pathFromUtf8(MDBRIDGE_TEST_DATA);
    guitest::PageHost host;
    host.resize(520, 1400);
    host.show();
    QVERIFY(QTest::qWaitForWindowExposed(&host));
    host.document().open(data / "projects" / "full.mdatk.toml");
    for (QWidget* page : host.pages()) {
        QString id = pages::pageIdOf(page);
        if (id != QStringLiteral("input") && id != QStringLiteral("groups") && id != QStringLiteral("analyses")) {
            continue;  // the other pages' tests take theirs
        }
        host.tabs()->setCurrentWidget(page);
        QApplication::processEvents();
        const fs::path shot = guitest::saveScreenshot(&host, screenshots(), "page-" + id.toStdString());
        qInfo("screenshot %s", core::pathToUtf8(shot).c_str());
        QVERIFY(fs::file_size(shot) > 1000);
    }
}

int main(int argc, char** argv) {
    const QByteArray platform = qgetenv("MDBRIDGE_GUI_TEST_PLATFORM");
    qputenv("QT_QPA_PLATFORM", platform.isEmpty() ? QByteArray("offscreen") : platform);
    QLocale::setDefault(QLocale::c());
    QApplication app(argc, argv);
    guitest::useBundledFont();
    GroupPagesTests tests;
    return QTest::qExec(&tests, argc, argv);
}

#include "test_group_pages.moc"
