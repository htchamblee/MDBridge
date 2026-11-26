#include "gui/forms/field_spec.hpp"

#include <functional>
#include <map>
#include <string>
#include <utility>

#include <array>

namespace mdbridge::app::forms {

FileFilter topologyFiles() {
    return {QStringLiteral("Topology files"),
            {"*.parm7", "*.prmtop", "*.top", "*.pdb", "*.psf", "*.gro"}};
}

FileFilter trajectoryFiles() {
    return {QStringLiteral("Trajectory files"),
            {"*.nc", "*.mdcrd", "*.crd", "*.trj", "*.dcd", "*.trr", "*.gro", "*.xtc", "*.pdb", "*.xyz", "*.binpos",
             "*.ncdf"}};
}

FileFilter referenceFiles() {
    return {QStringLiteral("Structure files"), {"*.pdb", "*.rst7", "*.ncrst", "*.inpcrd", "*.crd", "*.gro", "*.cif"}};
}

FileFilter dataFiles() { return {QStringLiteral("Data files"), {"*.dat", "*.csv", "*.tsv", "*.txt"}}; }

FileFilter mapFiles() {
    return {QStringLiteral("MRC / CCP4 maps"), {"*.mrc", "*.map", "*.ccp4", "*.map.gz", "*.mrc.gz"}};
}

FileFilter projectFiles() { return {QStringLiteral("MDATK projects"), {"*.mdatk.toml"}}; }

FileFilter runLogFiles() { return {QStringLiteral("Parameter logs of runs"), {"*_params_log.toml"}}; }

QString nameFilter(const FileFilter& filter) {
    QStringList patterns;
    const auto add = [&](const QString& p) {
        if (!patterns.contains(p)) {
            patterns.append(p);
        }
    };
    for (const QString& p : filter.patterns) {
        add(p);
    }
    for (const QString& p : filter.patterns) {
        add(p.toLower());
        add(p.toUpper());
    }
    return filter.name + QStringLiteral(" (") + patterns.join(QLatin1Char(' ')) + QLatin1Char(')');
}

QStringList nameFilters(const FieldSpec& spec) {
    QStringList out;
    for (const FileFilter& f : spec.filters) {
        out.append(nameFilter(f));
    }
    if (!out.isEmpty()) {
        out.append(QStringLiteral("All files (*)"));
    }
    return out;
}

namespace {

// ---- builders of the table below ----

FieldSpec number(double lo, double hi, double step = 1.0, int decimals = -1) {
    FieldSpec s;
    s.editor = Editor::Number;
    s.minimum = lo;
    s.maximum = hi;
    s.step = step;
    s.decimals = decimals;
    return s;
}

FieldSpec line(QString placeholder = {}, QString hint = {}) {
    FieldSpec s;
    s.editor = Editor::Line;
    s.placeholder = std::move(placeholder);
    s.hint = std::move(hint);
    return s;
}

FieldSpec combo(std::vector<Choice> choices) {
    FieldSpec s;
    s.editor = Editor::Combo;
    s.choices = std::move(choices);
    return s;
}

FieldSpec editableCombo(const QStringList& items = {}) {
    FieldSpec s;
    s.editor = Editor::EditableCombo;
    for (const QString& item : items) {
        s.choices.push_back({item, {}});
    }
    return s;
}

FieldSpec path(QString title, std::vector<FileFilter> filters = {}, PathMode mode = PathMode::OpenFile,
               QString hint = {}) {
    FieldSpec s;
    s.editor = Editor::Path;
    s.dialogTitle = std::move(title);
    s.filters = std::move(filters);
    s.pathMode = mode;
    s.hint = std::move(hint);
    return s;
}

FieldSpec color(bool allowEmpty = false, QString placeholder = {}) {
    FieldSpec s;
    s.editor = Editor::Color;
    s.allowEmpty = allowEmpty;
    s.placeholder = std::move(placeholder);
    return s;
}

FieldSpec multiLine(bool monospace, QString hint = {}) {
    FieldSpec s;
    s.editor = Editor::MultiLine;
    s.monospace = monospace;
    s.hint = std::move(hint);
    return s;
}

FieldSpec table(std::vector<Column> columns, core::Value newItem = {}) {
    FieldSpec s;
    s.editor = Editor::Table;
    s.columns = std::move(columns);
    s.newItem = std::move(newItem);
    return s;
}

FieldSpec withHint(FieldSpec s, QString hint) {
    s.hint = std::move(hint);
    return s;
}

std::vector<Choice> plainChoices(const auto& values) {
    std::vector<Choice> out;
    for (std::string_view v : values) {
        out.push_back({QString::fromUtf8(v.data(), static_cast<qsizetype>(v.size())), {}});
    }
    return out;
}

const QStringList kXUnits{"ps", "ns", "us", "frame", "frame1"};

// mdatk/gui/app.py of MDATK v0.0.4.1: the widget of every field its pages bind (the AP items of
// docs/specs/m5-plan.md section 4 change what the widgets do with values, not these settings)
std::map<std::string, FieldSpec, std::less<>> makeSpecs() {
    std::map<std::string, FieldSpec, std::less<>> m;
    // ---- Input & cleaning (app.py:422-559)
    m["notes"] = multiLine(false, QStringLiteral("Free text saved in the project file (aims, decisions, to-dos)."));
    m["inputs.topology"] = path(QStringLiteral("Select topology"), {topologyFiles()}, PathMode::OpenFile,
                                QStringLiteral("Amber .parm7 (preferred), cleaned PDB or another topology"));
    {
        FieldSpec s = table({{"", QStringLiteral("Trajectory file(s), in order")}}, core::Value(std::string()));
        s.browse = true;
        s.pathMode = PathMode::OpenFiles;
        s.dialogTitle = QStringLiteral("Select trajectories");
        s.filters = {trajectoryFiles()};
        s.hint = QStringLiteral("Amber NetCDF .nc (preferred), DCD, XTC, ... (several files are concatenated)");
        m["inputs.trajectories"] = s;
    }
    m["inputs.start"] = number(0, 1e9);
    m["inputs.stop"] = line(QStringLiteral("last"));
    m["inputs.step"] = number(1, 1e6);
    m["inputs.dt_ps"] = line(QStringLiteral("from file"), QStringLiteral("Leave blank to use the trajectory's time stamps"));
    m["inputs.time_unit"] = combo({{"ps", "ps"}, {"ns", "ns"}, {"us", QStringLiteral("µs")}, {"frame", "frames"}});
    {
        core::Table item{{"x_unit", core::Value(std::string("ps"))}};
        FieldSpec s = table({{"path", QStringLiteral("Path")},
                             {"label", QStringLiteral("Label")},
                             {"x_column", QStringLiteral("x col")},
                             {"x_unit", QStringLiteral("x unit")},
                             {"columns", QStringLiteral("Columns")}},
                            core::Value(item));
        s.browse = true;
        s.browseKey = "path";
        s.pathMode = PathMode::OpenFiles;
        s.dialogTitle = QStringLiteral("Select data files");
        s.filters = {dataFiles()};
        s.hint = QStringLiteral("Extra traces (e.g. cpptraj .dat files) interpolated onto the analysed frames");
        m["inputs.external"] = s;
    }
    m["inputs.external.x_unit"] = editableCombo(kXUnits);
    m["inputs.external.path"] = path(QStringLiteral("Select a data file"), {dataFiles()});
    m["cleaning.cpptraj_path"] = path(QStringLiteral("Select cpptraj"), {}, PathMode::OpenFile,
                                      QStringLiteral("Name on PATH or absolute path"));
    m["cleaning.conda_env"] = line(QStringLiteral("e.g. AmberTools26 (blank: cpptraj on PATH)"));
    m["cleaning.conda_exe"] = path(QStringLiteral("Select conda / mamba / micromamba"), {}, PathMode::OpenFile,
                                   QStringLiteral("Blank: found automatically (conda, mamba or micromamba)"));
    m["cleaning.extra_strip_mask"] = line(QStringLiteral("optional"));
    m["cleaning.autoimage_anchor"] = line(QStringLiteral("first molecule"));
    m["cleaning.fit_reference"] = combo({{"first", "first frame"}, {"frame", "frame #"}, {"file", "file"}});
    m["cleaning.reference_frame"] = number(1, 1e9);
    m["cleaning.reference_file"] = path(QStringLiteral("Select reference structure"), {referenceFiles()});
    m["cleaning.start"] = number(1, 1e9);
    m["cleaning.stop"] = line(QStringLiteral("last"));
    m["cleaning.offset"] = number(1, 1e6);
    m["cleaning.extra_commands"] = multiLine(true, QStringLiteral("Advanced: inserted before trajout"));

    // ---- Residue groups (app.py:561-637)
    m["groups"] = table({{"name", QStringLiteral("Name")},
                         {"residues", QStringLiteral("Residues")},
                         {"selection", QStringLiteral("Selection (overrides)")},
                         {"color", QStringLiteral("Colour")}},
                        core::Value(core::Table{{"name", core::Value(std::string("NewGroup"))},
                                                {"residues", core::Value(std::string())},
                                                {"selection", core::Value(std::string())},
                                                {"color", core::Value(std::string("#6250d6"))}}));
    m["groups.color"] = color();

    // ---- Analyses (app.py:639-756)
    m["analyses.proximity_search_radius"] = number(1, 30, 0.5, 1);
    m["analyses.rmsd.targets"] = table({{"label", QStringLiteral("Label")},
                                        {"selection", QStringLiteral("Selection")},
                                        {"fit_selection", QStringLiteral("Fit selection (optional)")}},
                                       core::Value(core::Table{
                                           {"label", core::Value(std::string("Target"))},
                                           {"selection", core::Value(std::string("group:Aptamer and not name H*"))},
                                           {"fit_selection", core::Value(std::string())}}));
    m["analyses.rmsd.reference_frame"] = number(0, 1e9);
    for (const char* analysis : {"contacts", "hbonds", "mg", "stacking", "energy"}) {
        for (const char* side : {"group_a", "group_b"}) {
            m[std::string("analyses.") + analysis + "." + side] = editableCombo();  // the group names
        }
    }
    m["analyses.contacts.cutoff"] = number(0.5, 20, 0.1, 2);
    m["analyses.hbonds.distance_cutoff"] = number(2.0, 5.0, 0.05, 2);
    m["analyses.hbonds.angle_cutoff"] = number(90, 180, 1, 0);
    m["analyses.mg.inner_cutoff"] = number(1.5, 4.0, 0.05, 2);
    m["analyses.mg.outer_cutoff"] = number(2.0, 8.0, 0.1, 2);
    m["analyses.mg.bridge_shell"] = combo({{"inner", "inner sphere"}, {"outer", "inner + outer sphere"}});
    m["analyses.mg.angle_tolerance"] = number(5, 45, 1, 0);
    m["analyses.mg.occupancy_selection"] = line(QStringLiteral("default: ligand selection"));
    m["analyses.stacking.scope"] = combo({{"all", "all ring pairs"}, {"inter", "group A vs group B"}});
    m["analyses.stacking.distance_cutoff"] = number(2.5, 8, 0.1, 2);
    m["analyses.stacking.angle_cutoff"] = number(5, 90, 1, 0);
    m["analyses.stacking.offset_cutoff"] = number(0.5, 5, 0.1, 2);
    m["analyses.stacking.planarity_cutoff"] = number(0.01, 2, 0.01, 2);
    // AP-19: Barnaba is never bundled, so only the built-in engine is offered ("barnaba" in a file
    // is kept and shown as unknown, AP-05)
    m["analyses.basepairs.engine"] = combo({{"builtin", "built-in (canonical pairs from H-bond geometry)"}});
    m["analyses.basepairs.stride"] = number(1, 100000);
    m["analyses.basepairs.cutoff"] = number(2.0, 5.0, 0.05, 2);
    m["analyses.energy.cutoff"] = number(4, 50, 0.5, 1);
    m["analyses.energy.dielectric"] = number(1, 100, 0.5, 1);
    m["analyses.energy.exclude_bonded"] = number(0, 3);
    m["analyses.energy.target_value"] =
        line(QStringLiteral("none"), QStringLiteral("Drawn as a reference line on charts and panels"));

    // ---- Structure view (app.py:758-1107)
    m["visual.zoom"] = withHint(number(0.05, 50, 0.1, 2),
                                QStringLiteral("1 fits the whole model; 2 shows it twice as large"));
    m["visual.center_on"] = line(QStringLiteral("model centre; 12, B:12 or origin"));
    m["visual.label_every"] = withHint(number(0, 10000), QStringLiteral("0 turns interval labels off"));
    m["visual.label_residues"] = line(QStringLiteral("e.g. 1, 12, 24, 40-45"));
    m["visual.structure_groups"] = line(QStringLiteral("all groups"), QStringLiteral("Comma separated group names"));
    m["visual.color_mode"] = combo(plainChoices(std::array<std::string_view, 3>{"chain", "groups", "proximity"}));
    m["visual.base_color"] = color();
    m["visual.ion_color"] = color();
    m["visual.chain_palette"] = editableCombo({"tab10", "Set2", "Dark2", "viridis"});
    m["visual.proximity_set"] = line({}, QStringLiteral("Group name, or sel:<selection>"));
    m["visual.proximity_cutoff"] = number(1, 30, 0.1, 1);
    m["visual.proximity_color"] = color();
    m["visual.glow_mode"] = combo({{"none", "none"},
                                   {"ever", "whole movie if ever in proximity"},
                                   {"fade", "while in proximity, then fade"},
                                   {"current", "only while in proximity"}});
    m["visual.glow_color"] = color();
    m["visual.fade_frames"] = number(0, 10000);
    m["visual.glow_strength"] = number(0, 3, 0.1, 1);
    m["visual.line_width"] = number(0.2, 10, 0.1, 1);
    m["visual.point_size"] = number(1, 200, 1, 0);
    m["visual.theme"] = combo({{"light", "light"}, {"dark", "dark"}});
    m["visual.font_size"] = number(6, 30, 0.5, 1);
    m["visual.title"] = line(QStringLiteral("figure title (optional)"));
    m["visual.representations"] = table({{"style", QStringLiteral("Style")},
                                         {"groups", QStringLiteral("Groups")},
                                         {"color", QStringLiteral("Colouring")},
                                         {"size", QStringLiteral("Size")},
                                         {"opacity", QStringLiteral("Opacity")},
                                         {"sugars", QStringLiteral("Sugars")}});
    m["visual.representations.style"] = combo({{"cartoon", "New Cartoon (bases as slabs)"},
                                               {"ribbon", "Ribbon"},
                                               {"tube", "Tube"},
                                               {"licorice", "Licorice"},
                                               {"vdw", "VDW spheres"},
                                               {"surf", "Surface"},
                                               {"trace", "Trace (balls and sticks)"}});
    m["visual.representations.groups"] = line(QStringLiteral("all groups"));
    m["visual.representations.color"] = combo({{"auto", "auto"},
                                               {"residue", "residue colours"},
                                               {"element", "elements"},
                                               {"nucleotide", "nucleotides (A/G/C/U)"}});
    m["visual.representations.size"] = number(0.2, 5, 0.1, 1);
    m["visual.representations.opacity"] = number(0.05, 1, 0.05, 2);
    m["visual.rep_rest"] = combo({{"trace", "as a trace"}, {"hide", "hidden"}});
    m["visual.backbone_atom"] = combo({{"P", "P"}, {"C4'", "C4'"}});
    // AP-19 / RN-20: MDBridge's renderers ("vtk" means the GPU renderer, as "vulkan")
    m["visual.renderer"] = combo({{"auto", "auto"},
                                  {"software", "built-in (software)"},
                                  {"vulkan", "GPU (Vulkan)"},
                                  {"vtk", "GPU (vtk: as Vulkan)"}});
    m["density.map_path"] = path(QStringLiteral("Select a density map"), {mapFiles()}, PathMode::OpenFile,
                                 QStringLiteral("MRC2014 / CCP4 map (.mrc, .map)"));
    m["density.color"] = color();
    m["density.opacity"] = number(0, 1, 0.05, 2);
    m["density.level"] = number(-5, 20, 0.1, 1);
    m["density.resolution"] = number(0.5, 50, 0.1, 1);
    m["density.fit_frame"] = line(QStringLiteral("ref"));
    m["density.search"] = combo({{"global", "global search"}, {"manual", "my placement"}});
    m["density.angular_step"] = number(5, 45, 1, 0);

    // ---- Outputs (app.py:1336-1419)
    m["output.out_dir"] = path(QStringLiteral("Select output directory"), {}, PathMode::Folder);
    m["output.cache_dir"] = path(QStringLiteral("Select cache directory"), {}, PathMode::Folder,
                                 QStringLiteral("Blank: <output>/MDATK_cache (analysis .npz + JSON, cleaned trajectories)"));
    m["output.video.width"] = number(64, 7680, 2);
    m["output.video.height"] = number(64, 4320, 2);
    m["output.video.dpi"] = number(10, 600);
    m["output.video.fps"] = number(1, 240);
    m["output.video.frame_stride"] = number(1, 10000);
    m["output.video.codec"] = editableCombo({"libx264", "libx265", "mpeg4"});
    m["output.video.crf"] = number(0, 51);
    m["output.video.ffmpeg_path"] = path(QStringLiteral("Select ffmpeg"));
    m["output.charts.series"] =
        line(QStringLiteral("all series"), QStringLiteral("Comma separated data keys, blank = every time series"));
    m["output.charts.width"] = number(2, 30, 0.5, 1);
    m["output.charts.height"] = number(2, 30, 0.5, 1);
    m["output.charts.dpi"] = number(10, 1200, 10);
    m["output.charts.style"] = editableCombo({"mdatk"});
    m["output.charts.color"] = color(true, QStringLiteral("palette"));
    m["output.charts.line_width"] = number(0.2, 8, 0.1, 1);
    m["output.charts.title"] = line(QStringLiteral("{label}"));
    m["output.charts.xlabel"] = line(QStringLiteral("auto"));
    m["output.charts.ylabel"] = line(QStringLiteral("auto"));
    m["output.charts.smoothing"] = number(0, 100000);
    m["output.pdb.frames"] = line(QStringLiteral("e.g. 0, 10, 20-50:10, last"));
    m["output.pdb.times"] = line(QStringLiteral("e.g. 1.5ns, 2000ps"));
    m["output.pdb.selection"] = line(QStringLiteral("all atoms; or e.g. Aptamer, Target"));
    return m;
}

}  // namespace

FieldSpec defaultFieldSpec(std::string_view modelPath) {
    static const std::map<std::string, FieldSpec, std::less<>> specs = makeSpecs();
    const auto it = specs.find(modelPath);
    return it == specs.end() ? FieldSpec{} : it->second;
}

}  // namespace mdbridge::app::forms
