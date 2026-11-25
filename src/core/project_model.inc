// Descriptor tables of the parameter model (mdatk/config.py of MDATK v0.0.4.1).
// Included by src/core/project.cpp inside its anonymous namespace.
//
// fieldsOf(T*): the fields of a model struct in declaration order, by project-file key.
// kModelComments: projectfile.config_doc of every path of the model (list indices left out), in
// declaration order; tests/core/test_project_model.cpp checks them against
// tests/data/project-model/config_doc.json (tests/parity/make_model_cases.py).

template <class S, class M>
struct Field {
    std::string_view name;
    M S::*member;
    // written only when set (an optional with a value): the keys MDATK 0.1.x added, which a file
    // of MDATK 0.0.4.1 does not have and MDBridge does not write unless the file had them (D9)
    bool onlyWhenSet = false;
};

template <class S, class M>
constexpr Field<S, M> field(std::string_view name, M S::*member) {
    return {name, member};
}

template <class S, class M>
constexpr Field<S, std::optional<M>> fieldWhenSet(std::string_view name, std::optional<M> S::*member) {
    return {name, member, true};
}

constexpr auto fieldsOf(const ResidueGroup*) {
    using S = ResidueGroup;
    return std::make_tuple(field("name", &S::name), field("residues", &S::residues),
                           field("selection", &S::selection), field("color", &S::color));
}

constexpr auto fieldsOf(const ExternalSeries*) {
    using S = ExternalSeries;
    return std::make_tuple(field("path", &S::path), field("label", &S::label),
                           field("x_column", &S::x_column), field("x_unit", &S::x_unit),
                           field("columns", &S::columns), field("delimiter", &S::delimiter));
}

constexpr auto fieldsOf(const InputParams*) {
    using S = InputParams;
    return std::make_tuple(field("topology", &S::topology), field("trajectories", &S::trajectories),
                           field("start", &S::start), field("stop", &S::stop),
                           field("step", &S::step), field("dt_ps", &S::dt_ps),
                           field("time_unit", &S::time_unit), field("external", &S::external));
}

constexpr auto fieldsOf(const CleaningParams*) {
    using S = CleaningParams;
    return std::make_tuple(
        field("enabled", &S::enabled), field("cpptraj_path", &S::cpptraj_path),
        field("conda_env", &S::conda_env), field("conda_exe", &S::conda_exe),
        field("solvent_mask", &S::solvent_mask), field("ion_mask", &S::ion_mask),
        field("extra_strip_mask", &S::extra_strip_mask), field("autoimage", &S::autoimage),
        field("autoimage_anchor", &S::autoimage_anchor), field("fit", &S::fit),
        field("fit_mask", &S::fit_mask), field("fit_reference", &S::fit_reference),
        field("remove_box_after_fit", &S::remove_box_after_fit),
        field("reference_frame", &S::reference_frame), field("reference_file", &S::reference_file),
        field("start", &S::start), field("stop", &S::stop), field("offset", &S::offset),
        field("extra_commands", &S::extra_commands), field("output_dir", &S::output_dir));
}

constexpr auto fieldsOf(const RMSDTarget*) {
    using S = RMSDTarget;
    return std::make_tuple(field("label", &S::label), field("selection", &S::selection),
                           field("fit_selection", &S::fit_selection));
}

constexpr auto fieldsOf(const RMSDParams*) {
    using S = RMSDParams;
    return std::make_tuple(field("enabled", &S::enabled), field("targets", &S::targets),
                           field("reference_frame", &S::reference_frame),
                           field("superposition", &S::superposition),
                           field("mass_weighted", &S::mass_weighted));
}

constexpr auto fieldsOf(const ContactParams*) {
    using S = ContactParams;
    return std::make_tuple(field("enabled", &S::enabled), field("group_a", &S::group_a),
                           field("group_b", &S::group_b), field("cutoff", &S::cutoff),
                           field("heavy_only", &S::heavy_only));
}

constexpr auto fieldsOf(const HBondParams*) {
    using S = HBondParams;
    return std::make_tuple(field("enabled", &S::enabled), field("group_a", &S::group_a),
                           field("group_b", &S::group_b),
                           field("distance_cutoff", &S::distance_cutoff),
                           field("use_angle", &S::use_angle),
                           field("angle_cutoff", &S::angle_cutoff),
                           field("donor_elements", &S::donor_elements),
                           field("acceptor_elements", &S::acceptor_elements));
}

constexpr auto fieldsOf(const MgParams*) {
    using S = MgParams;
    return std::make_tuple(
        field("enabled", &S::enabled), field("mg_selection", &S::mg_selection),
        field("ligand_selection", &S::ligand_selection),
        field("ligand_elements", &S::ligand_elements), field("group_a", &S::group_a),
        field("group_b", &S::group_b), field("inner_cutoff", &S::inner_cutoff),
        field("outer_cutoff", &S::outer_cutoff), field("bridge_shell", &S::bridge_shell),
        field("geometry", &S::geometry), field("angle_tolerance", &S::angle_tolerance),
        field("occupancy_selection", &S::occupancy_selection),
        field("occupancy_elements", &S::occupancy_elements));
}

constexpr auto fieldsOf(const StackingParams*) {
    using S = StackingParams;
    return std::make_tuple(
        field("enabled", &S::enabled), field("selection", &S::selection), field("scope", &S::scope),
        field("group_a", &S::group_a), field("group_b", &S::group_b),
        field("distance_cutoff", &S::distance_cutoff), field("angle_cutoff", &S::angle_cutoff),
        field("offset_cutoff", &S::offset_cutoff), field("planarity_cutoff", &S::planarity_cutoff));
}

constexpr auto fieldsOf(const BasePairParams*) {
    using S = BasePairParams;
    return std::make_tuple(field("engine", &S::engine), field("stride", &S::stride),
                           field("cutoff", &S::cutoff), field("wobble", &S::wobble),
                           field("noncanonical", &S::noncanonical));
}

constexpr auto fieldsOf(const EnergyParams*) {
    using S = EnergyParams;
    return std::make_tuple(
        field("enabled", &S::enabled), field("group_a", &S::group_a), field("group_b", &S::group_b),
        field("cutoff", &S::cutoff), field("dielectric", &S::dielectric),
        field("distance_dependent", &S::distance_dependent), field("use_parm7", &S::use_parm7),
        field("target_value", &S::target_value), field("exclude_bonded", &S::exclude_bonded));
}

constexpr auto fieldsOf(const AnalysisParams*) {
    using S = AnalysisParams;
    return std::make_tuple(
        field("rmsd", &S::rmsd), field("contacts", &S::contacts), field("hbonds", &S::hbonds),
        field("mg", &S::mg), field("stacking", &S::stacking), field("energy", &S::energy),
        field("basepairs", &S::basepairs), field("use_pbc", &S::use_pbc),
        field("proximity_search_radius", &S::proximity_search_radius));
}

constexpr auto fieldsOf(const RepSpec*) {
    using S = RepSpec;
    return std::make_tuple(field("style", &S::style), field("groups", &S::groups),
                           field("color", &S::color), field("size", &S::size),
                           field("opacity", &S::opacity), field("sugars", &S::sugars));
}

constexpr auto fieldsOf(const VisualParams*) {
    using S = VisualParams;
    return std::make_tuple(
        field("structure_groups", &S::structure_groups), field("trace_atoms", &S::trace_atoms),
        field("representations", &S::representations), field("renderer", &S::renderer),
        field("rep_rest", &S::rep_rest), field("backbone_atom", &S::backbone_atom),
        field("rep_hydrogens", &S::rep_hydrogens), field("show_ions", &S::show_ions),
        field("align_frames", &S::align_frames), field("rotation", &S::rotation),
        field("zoom", &S::zoom), field("center_on", &S::center_on),
        field("center_follow", &S::center_follow), field("color_mode", &S::color_mode),
        field("base_color", &S::base_color), field("chain_palette", &S::chain_palette),
        field("proximity_set", &S::proximity_set), field("proximity_cutoff", &S::proximity_cutoff),
        field("proximity_color", &S::proximity_color), field("glow_mode", &S::glow_mode),
        field("glow_color", &S::glow_color), field("fade_frames", &S::fade_frames),
        field("glow_strength", &S::glow_strength), field("ion_color", &S::ion_color),
        field("line_width", &S::line_width), field("point_size", &S::point_size),
        field("depth_cue", &S::depth_cue), field("theme", &S::theme),
        field("font_size", &S::font_size), field("show_time_label", &S::show_time_label),
        field("label_every", &S::label_every), field("label_residues", &S::label_residues),
        field("title", &S::title));
}

constexpr auto fieldsOf(const DensityParams*) {
    using S = DensityParams;
    return std::make_tuple(
        field("map_path", &S::map_path), field("show", &S::show),
        field("resolution", &S::resolution), field("level", &S::level), field("color", &S::color),
        field("opacity", &S::opacity),
        field("front", &S::front), field("fit_frame", &S::fit_frame), field("search", &S::search),
        field("angular_step", &S::angular_step), field("auto_fit", &S::auto_fit),
        field("transform", &S::transform), field("anchor", &S::anchor), field("fit_cc", &S::fit_cc),
        field("fit_info", &S::fit_info));
}

constexpr auto fieldsOf(const PanelSpec*) {
    using S = PanelSpec;
    return std::make_tuple(field("type", &S::type), field("series", &S::series),
                           field("title", &S::title), field("options", &S::options));
}

constexpr auto fieldsOf(const LayoutParams*) {
    using S = LayoutParams;
    return std::make_tuple(field("preset", &S::preset), field("panels", &S::panels));
}

constexpr auto fieldsOf(const VideoParams*) {
    using S = VideoParams;
    return std::make_tuple(
        field("enabled", &S::enabled), field("filename", &S::filename), field("fps", &S::fps),
        field("width", &S::width), field("height", &S::height), field("dpi", &S::dpi),
        field("codec", &S::codec), field("crf", &S::crf), field("frame_stride", &S::frame_stride),
        field("ffmpeg_path", &S::ffmpeg_path));
}

constexpr auto fieldsOf(const ChartParams*) {
    using S = ChartParams;
    return std::make_tuple(
        field("enabled", &S::enabled), field("series", &S::series),
        field("include_heatmaps", &S::include_heatmaps), field("formats", &S::formats),
        field("dpi", &S::dpi), field("width", &S::width), field("height", &S::height),
        field("style", &S::style), field("color", &S::color), field("line_width", &S::line_width),
        field("title", &S::title), field("xlabel", &S::xlabel), field("ylabel", &S::ylabel),
        field("smoothing", &S::smoothing), field("show_target", &S::show_target),
        field("prefix", &S::prefix));
}

constexpr auto fieldsOf(const PDBParams*) {
    using S = PDBParams;
    return std::make_tuple(field("enabled", &S::enabled), field("frames", &S::frames),
                           field("times", &S::times), field("selection", &S::selection),
                           field("prefix", &S::prefix));
}

constexpr auto fieldsOf(const OutputParams*) {
    using S = OutputParams;
    return std::make_tuple(field("out_dir", &S::out_dir), field("cache_dir", &S::cache_dir),
                           field("video", &S::video), field("charts", &S::charts),
                           field("pdb", &S::pdb), field("write_params_log", &S::write_params_log),
                           field("timestamp_files", &S::timestamp_files));
}

constexpr auto fieldsOf(const PlanNode*) {
    using S = PlanNode;
    return std::make_tuple(field("id", &S::id), field("tool", &S::tool), field("step", &S::step),
                           field("title", &S::title), field("x", &S::x), field("y", &S::y),
                           field("settings", &S::settings), field("files", &S::files),
                           fieldWhenSet("adopted", &S::adopted), fieldWhenSet("folder_tag", &S::folder_tag));
}

constexpr auto fieldsOf(const PlanEdge*) {
    using S = PlanEdge;
    return std::make_tuple(field("source", &S::source), field("output", &S::output),
                           field("target", &S::target), field("input", &S::input));
}

constexpr auto fieldsOf(const WorkflowParams*) {
    using S = WorkflowParams;
    return std::make_tuple(field("project_name", &S::project_name), field("stamp", &S::stamp),
                           field("use_stamp", &S::use_stamp), field("out_dir", &S::out_dir),
                           field("step_folders", &S::step_folders),
                           fieldWhenSet("supersede_runs", &S::supersede_runs),
                           fieldWhenSet("parallel", &S::parallel), fieldWhenSet("cpu_cores", &S::cpu_cores),
                           fieldWhenSet("gpus", &S::gpus), fieldWhenSet("last_node", &S::last_node),
                           field("nodes", &S::nodes),
                           field("edges", &S::edges), field("scratch", &S::scratch));
}

constexpr auto fieldsOf(const Project*) {
    using S = Project;
    return std::make_tuple(field("notes", &S::notes), field("inputs", &S::inputs),
                           field("cleaning", &S::cleaning), field("groups", &S::groups),
                           field("analyses", &S::analyses), field("visual", &S::visual),
                           field("density", &S::density), field("layout", &S::layout),
                           field("output", &S::output), field("workflow", &S::workflow));
}

// Characters for which Python's str.isdigit() holds but \d (str.isdecimal) does not.
inline constexpr std::pair<char32_t, char32_t> kDigitNotDecimal[] = {
    {0xB2, 0xB3}, {0xB9, 0xB9}, {0x1369, 0x1371}, {0x19DA, 0x19DA},
    {0x2070, 0x2070}, {0x2074, 0x2079}, {0x2080, 0x2089}, {0x2460, 0x2468},
    {0x2474, 0x247C}, {0x2488, 0x2490}, {0x24EA, 0x24EA}, {0x24F5, 0x24FD},
    {0x24FF, 0x24FF}, {0x2776, 0x277E}, {0x2780, 0x2788}, {0x278A, 0x2792},
    {0x10A40, 0x10A43}, {0x10E60, 0x10E68}, {0x11052, 0x1105A}, {0x1F100, 0x1F10A},
};

inline constexpr std::pair<std::string_view, std::string_view> kModelComments[] = {
    {"notes", "Free-text notes about the project (what, why, decisions, to-dos)"},
    {"inputs", "Input files and frame range"},
    {"inputs.topology", "Topology: Amber .parm7 / .prmtop (preferred), PDB, PSF, GRO, ..."},
    {"inputs.trajectories", "Trajectory files in order; several files are concatenated"},
    {"inputs.start",
     "First analysed frame (0-based; of the cleaned trajectory when cleaning is on)"},
    {"inputs.stop", "Stop before this frame (0-based, exclusive; not set: up to the last frame)"},
    {"inputs.step", "Analyse every n-th frame"},
    {"inputs.dt_ps",
     "Time between stored frames in ps (not set: from the trajectory's time stamps)"},
    {"inputs.time_unit", "Unit shown on time axes: \"ps\", \"ns\", \"us\" or \"frame\""},
    {"inputs.external",
     "External data files shown as extra traces (interpolated onto the analysed frames)"},
    {"inputs.external.path", "CSV or whitespace-delimited text file (e.g. a cpptraj .dat file)"},
    {"inputs.external.label", "Name of the trace in panels and charts (\"\": the file name)"},
    {"inputs.external.x_column",
     "Column holding the time axis: name or 0-based index (\"\": the first column)"},
    {"inputs.external.x_unit",
     "Unit of that column: \"ps\", \"ns\", \"us\", \"frame\" (0-based) or \"frame1\" (1-based)"},
    {"inputs.external.columns", "Columns to show, by name or index (empty: all other columns)"},
    {"inputs.external.delimiter", "Column separator (\"\": detect comma, tab or whitespace)"},
    {"cleaning", "Trajectory cleaning with cpptraj"},
    {"cleaning.enabled", "Strip, image and fit the trajectory with cpptraj before the analyses"},
    {"cleaning.cpptraj_path", "cpptraj executable: a name on PATH or an absolute path"},
    {"cleaning.conda_env",
     "Conda environment providing cpptraj (a name such as AmberTools26 or a prefix path); every "
     "call runs through 'conda run -n <env>'"},
    {"cleaning.conda_exe", "conda / mamba / micromamba executable (\"\": found automatically)"},
    {"cleaning.solvent_mask", "cpptraj mask of the solvent residues to strip"},
    {"cleaning.ion_mask", "cpptraj mask of the ions to strip (Mg2+ is kept by default)"},
    {"cleaning.extra_strip_mask", "Further atoms to strip (cpptraj mask; optional)"},
    {"cleaning.autoimage", "Re-image molecules into the primary box (cpptraj autoimage)"},
    {"cleaning.autoimage_anchor",
     "Anchor mask for autoimage (\"\": cpptraj's default, the first molecule)"},
    {"cleaning.fit", "RMS-fit every frame onto the reference"},
    {"cleaning.fit_mask", "cpptraj mask of the atoms used for the RMS fit"},
    {"cleaning.fit_reference",
     "Fit reference: the \"first\" frame, a given \"frame\" number or a \"file\""},
    {"cleaning.remove_box_after_fit",
     "Remove the box after fitting (a rotated box cannot be stored in NetCDF)"},
    {"cleaning.reference_frame",
     "Reference frame for fit_reference = \"frame\" (1-based, cpptraj convention)"},
    {"cleaning.reference_file", "Reference structure for fit_reference = \"file\""},
    {"cleaning.start", "First frame read by cpptraj (1-based trajin argument)"},
    {"cleaning.stop", "Last frame read by cpptraj (1-based, inclusive; not set: the last frame)"},
    {"cleaning.offset", "Read every n-th frame (cpptraj trajin offset)"},
    {"cleaning.extra_commands", "Advanced: cpptraj commands inserted before trajout"},
    {"cleaning.output_dir", "Folder for the cleaned files (\"\": <cache_dir>/cleaned)"},
    {"groups", "Residue groups"},
    {"groups.name", "Group name, used in selections (group:<name>), panels and legends"},
    {"groups.residues",
     "Residue numbers or ranges, e.g. \"1-40, 45\"; may be chain-qualified (\"A:1-40\")"},
    {"groups.selection", "MDAnalysis selection that replaces 'residues' when given"},
    {"groups.color", "Colour of the group in the structure view, legends and charts (#rrggbb)"},
    {"analyses", "Analyses"},
    {"analyses.rmsd", "RMSD time series"},
    {"analyses.rmsd.enabled", "Compute RMSD time series"},
    {"analyses.rmsd.targets", "One RMSD series per target"},
    {"analyses.rmsd.targets.label", "Name of the series (data key rmsd/<label>)"},
    {"analyses.rmsd.targets.selection",
     "Atoms whose RMSD is computed (selection; group:<name> allowed)"},
    {"analyses.rmsd.targets.fit_selection",
     "Atoms superposed before computing the RMSD (\"\": the same atoms)"},
    {"analyses.rmsd.reference_frame",
     "Reference structure: index into the analysed frames (0: the first)"},
    {"analyses.rmsd.superposition",
     "Superpose each frame onto the reference before computing the RMSD"},
    {"analyses.rmsd.mass_weighted", "Mass-weighted superposition and RMSD"},
    {"analyses.contacts", "Heavy-atom contacts between two groups"},
    {"analyses.contacts.enabled", "Count heavy-atom contacts between two groups"},
    {"analyses.contacts.group_a", "First residue group"},
    {"analyses.contacts.group_b", "Second residue group"},
    {"analyses.contacts.cutoff", "Contact distance cutoff (Å)"},
    {"analyses.contacts.heavy_only", "Ignore hydrogens"},
    {"analyses.hbonds", "Hydrogen bonds between two groups"},
    {"analyses.hbonds.enabled", "Count hydrogen bonds between two groups"},
    {"analyses.hbonds.group_a", "First residue group"},
    {"analyses.hbonds.group_b", "Second residue group"},
    {"analyses.hbonds.distance_cutoff", "Maximum donor-acceptor distance (Å)"},
    {"analyses.hbonds.use_angle",
     "Also require a minimum D-H...A angle (needs hydrogens in the topology)"},
    {"analyses.hbonds.angle_cutoff", "Minimum D-H...A angle (degrees)"},
    {"analyses.hbonds.donor_elements", "Elements that can donate (comma separated)"},
    {"analyses.hbonds.acceptor_elements", "Elements that can accept (comma separated)"},
    {"analyses.mg", "Mg2+ coordination, bridges and occupancy"},
    {"analyses.mg.enabled",
     "Analyse Mg2+ coordination: inner/outer-sphere contacts, bridges and occupancy"},
    {"analyses.mg.mg_selection", "Selection of the Mg2+ ions"},
    {"analyses.mg.ligand_selection", "Atoms that can coordinate Mg2+"},
    {"analyses.mg.ligand_elements", "Elements of the coordinating atoms (comma separated)"},
    {"analyses.mg.group_a", "First group of a Mg2+ bridge"},
    {"analyses.mg.group_b", "Second group of a Mg2+ bridge"},
    {"analyses.mg.inner_cutoff", "Inner-sphere Mg-O distance cutoff (Å)"},
    {"analyses.mg.outer_cutoff", "Outer-sphere (water-mediated) Mg-O distance cutoff (Å)"},
    {"analyses.mg.bridge_shell", "Shell that counts as a bridge: \"inner\" or \"outer\""},
    {"analyses.mg.geometry", "Require octahedral O-Mg-O angles (about 90 or 180 degrees)"},
    {"analyses.mg.angle_tolerance", "Allowed deviation from the octahedral angles (degrees)"},
    {"analyses.mg.occupancy_selection",
     "Atoms whose Mg2+ occupancy is reported (\"\": the ligand selection)"},
    {"analyses.mg.occupancy_elements", "Elements counted for the occupancy (comma separated)"},
    {"analyses.stacking", "Base stacking"},
    {"analyses.stacking.enabled", "Detect base stacking from ring geometry"},
    {"analyses.stacking.selection", "Residues whose rings are considered"},
    {"analyses.stacking.scope", "\"all\" ring pairs, or \"inter\" (group A against group B only)"},
    {"analyses.stacking.group_a", "First group for scope = \"inter\""},
    {"analyses.stacking.group_b", "Second group for scope = \"inter\""},
    {"analyses.stacking.distance_cutoff", "Maximum ring-centroid distance (Å)"},
    {"analyses.stacking.angle_cutoff", "Maximum angle between the ring normals (degrees)"},
    {"analyses.stacking.offset_cutoff", "Maximum lateral displacement of the centroids (Å)"},
    {"analyses.stacking.planarity_cutoff",
     "Maximum RMS deviation of ring atoms from their plane (Å)"},
    {"analyses.energy", "Non-bonded interaction energy"},
    {"analyses.energy.enabled", "Non-bonded interaction energy between two groups"},
    {"analyses.energy.group_a", "First residue group"},
    {"analyses.energy.group_b", "Second residue group"},
    {"analyses.energy.cutoff", "Pair distance cutoff (Å)"},
    {"analyses.energy.dielectric", "Relative dielectric constant"},
    {"analyses.energy.distance_dependent",
     "Distance-dependent dielectric (epsilon = dielectric x r)"},
    {"analyses.energy.use_parm7", "Use the exact Lennard-Jones parameters of a parm7 topology"},
    {"analyses.energy.target_value",
     "Reference energy drawn on charts and panels (kcal/mol; not set: none)"},
    {"analyses.energy.exclude_bonded", "Exclude pairs within n bonds (0: exclude none)"},
    {"analyses.basepairs", "Base pairs (secondary structure)"},
    {"analyses.basepairs.engine",
     "Base pairs of the residues in the structure view, for the secondary-structure panel and "
     "PDB annotation; computed whenever a panel needs them. Engine: \"builtin\" (canonical "
     "pairs from H-bond geometry) or \"barnaba\" (optional)"},
    {"analyses.basepairs.stride",
     "Re-detect the pairs every n-th analysed frame (pairs hold in between)"},
    {"analyses.basepairs.cutoff",
     "Builtin: maximum donor-acceptor distance of a base-pair H-bond (Å)"},
    {"analyses.basepairs.wobble", "Builtin: also count G-U wobble pairs"},
    {"analyses.basepairs.noncanonical",
     "Barnaba: also keep non-canonical (non Watson-Crick) pairs"},
    {"analyses.use_pbc", "Minimum-image distances (cleaned, imaged systems rarely need it)"},
    {"analyses.proximity_search_radius",
     "Largest proximity cutoff usable later without re-analysing (Å)"},
    {"visual", "Structure view: colouring, glow, 3D representations, labels"},
    {"visual.structure_groups", "Residue groups shown in the structure view (empty: all)"},
    {"visual.trace_atoms",
     "Atom used for each residue's 2D trace point, first one found wins (P, C4', C1', CA ...)"},
    {"visual.representations", "3D representations (empty: the 2D backbone trace)"},
    {"visual.representations.style",
     "Representation: cartoon, ribbon, tube, licorice, vdw, surf or trace"},
    {"visual.representations.groups", "Residue groups drawn (empty: every shown group)"},
    {"visual.representations.color", "Colouring: auto, residue, element or nucleotide"},
    {"visual.representations.size", "Scales widths and radii"},
    {"visual.representations.opacity", "Opacity, 0 (invisible) to 1"},
    {"visual.representations.sugars", "Cartoon only: draw the ribose rings (shows the pucker)"},
    {"visual.renderer", "3D renderer: \"auto\" / \"software\" (built-in) or \"vtk\" (optional)"},
    {"visual.rep_rest", "Residues no representation draws: \"trace\" (3D trace) or \"hide\""},
    {"visual.backbone_atom",
     "Spline control atom of ribbons, cartoons and traces: \"P\" or \"C4'\""},
    {"visual.rep_hydrogens", "Include hydrogens in all-atom representations"},
    {"visual.show_ions", "Show Mg2+ ions in the structure view"},
    {"visual.align_frames", "Superpose each frame's trace onto the reference frame"},
    {"visual.rotation", "View rotation in degrees about X, Y and Z (R = Rz Ry Rx)"},
    {"visual.zoom",
     "Magnification of the structure view (1: the whole model fits; 2: twice as large)"},
    {"visual.center_on",
     "View centre: \"\" (centre of the model), \"origin\" (0, 0, 0 of the coordinates) or a "
     "residue, e.g. \"12\" or \"B:12\""},
    {"visual.center_follow",
     "Keep the centre residue in the middle in every frame (not only the first)"},
    {"visual.color_mode", "Residue colouring: \"chain\", \"groups\" or \"proximity\""},
    {"visual.base_color", "Colour of residues no group or proximity colour applies to"},
    {"visual.chain_palette", "Matplotlib palette for chain colouring"},
    {"visual.proximity_set",
     "Group (or sel:<selection>) whose proximity colours and glows residues"},
    {"visual.proximity_cutoff", "Distance that counts as proximity (Å)"},
    {"visual.proximity_color", "Colour of residues in proximity"},
    {"visual.glow_mode",
     "Glow of residues near the proximity set: \"none\", \"ever\" (whole movie), \"fade\" or "
     "\"current\""},
    {"visual.glow_color", "Glow colour"},
    {"visual.fade_frames", "Frames over which the glow fades (glow_mode = \"fade\")"},
    {"visual.glow_strength", "Glow intensity factor"},
    {"visual.ion_color", "Colour of Mg2+ ions"},
    {"visual.line_width", "Line width of the 2D trace (points)"},
    {"visual.point_size", "Marker size of the 2D trace points"},
    {"visual.depth_cue", "Fade distant parts toward the background colour"},
    {"visual.theme", "\"light\" or \"dark\""},
    {"visual.font_size", "Base font size (points)"},
    {"visual.show_time_label", "Show the time label in the structure view"},
    {"visual.label_every", "Residue-number tags on every n-th residue (0: off)"},
    {"visual.label_residues", "Explicit residue numbers / ranges to tag, e.g. \"1, 12, 24\""},
    {"visual.title", "Figure title (\"\": none)"},
    {"density", "Density map shadow and fit"},
    {"density.map_path", "Density map file (.mrc / .map / .ccp4)"},
    {"density.show", "Draw the map shadow"},
    {"density.resolution", "Resolution in Å; width of the simulated map used for fitting"},
    {"density.level", "Contour level in map standard deviations above the mean"},
    {"density.color", "Shadow colour"},
    {"density.opacity", "Shadow opacity, 0 (invisible) to 1"},
    {"density.front", "Draw the shadow in front of 3D representations (else behind them)"},
    {"density.fit_frame",
     "\"ref\" (the first analysed frame) or a trajectory frame number to fit on"},
    {"density.search",
     "Automatic fit: \"global\" search or \"manual\" (refine the stored placement)"},
    {"density.angular_step", "Degrees between sampled orientations of the global search"},
    {"density.auto_fit", "Fit before rendering when no placement is stored"},
    {"density.transform", "Placement: 3x3 rotation (row-major) + translation, model -> map"},
    {"density.anchor",
     "Backbone points of the reference frame the placement was made on (the map follows the "
     "molecule when the reference changes: other start frame, cleaned vs raw files)"},
    {"density.fit_cc", "Correlation of the stored placement (0: not fitted)"},
    {"density.fit_info", "Summary of the last fit"},
    {"layout", "Figure layout: grid and panels"},
    {"layout.preset", "Grid layout (see 'mdatk list')"},
    {"layout.panels", "Panels in slot order"},
    {"layout.panels.type", "Panel type: structure, timeseries, dual_axis, heatmap or contact_map"},
    {"layout.panels.series", "Data keys shown, e.g. rmsd/Complex or contacts/per_residue_a"},
    {"layout.panels.title", "Panel title (\"\": automatic)"},
    {"layout.panels.options", "Panel-specific options, e.g. axis limits or the contact-map mode"},
    {"output", "Outputs: folders, video, charts, PDB frames"},
    {"output.out_dir", "Output folder (a relative path is relative to the project file)"},
    {"output.cache_dir",
     "Cache of analyses and cleaned trajectories (\"\": <out_dir>/MDATK_cache)"},
    {"output.video", "MP4 video"},
    {"output.video.enabled", "Write the MP4 video"},
    {"output.video.filename", "Video file name (time-stamped when output.timestamp_files is on)"},
    {"output.video.fps", "Frames per second"},
    {"output.video.width", "Width in pixels"},
    {"output.video.height", "Height in pixels"},
    {"output.video.dpi", "Figure resolution; sets text and line sizes relative to the frame"},
    {"output.video.codec", "FFmpeg video codec"},
    {"output.video.crf",
     "Constant rate factor: lower means better quality and larger files (libx264: 18-28)"},
    {"output.video.frame_stride", "Use every n-th analysed frame"},
    {"output.video.ffmpeg_path", "FFmpeg executable"},
    {"output.charts", "Time charts"},
    {"output.charts.enabled", "Write high-resolution time charts"},
    {"output.charts.series", "Data keys to chart (empty: every scalar series)"},
    {"output.charts.include_heatmaps", "Also chart per-residue heatmaps"},
    {"output.charts.formats", "File formats, e.g. png, pdf, svg"},
    {"output.charts.dpi", "Raster resolution"},
    {"output.charts.width", "Width in inches"},
    {"output.charts.height", "Height in inches"},
    {"output.charts.style", "\"mdatk\" (built-in look) or any matplotlib style name"},
    {"output.charts.color", "Line colour (\"\": the first palette colour)"},
    {"output.charts.line_width", "Line width (points)"},
    {"output.charts.title",
     "Chart title (\"\": the series label; \"{label}\" placeholders allowed)"},
    {"output.charts.xlabel", "x-axis label (\"\": automatic)"},
    {"output.charts.ylabel", "y-axis label (\"\": automatic)"},
    {"output.charts.smoothing", "Running-average window in frames (0: off)"},
    {"output.charts.show_target", "Draw the energy target value when one is set"},
    {"output.charts.prefix", "File name prefix"},
    {"output.pdb", "PDB export of selected frames"},
    {"output.pdb.enabled", "Export selected frames as PDB files"},
    {"output.pdb.frames", "Frames to export, e.g. \"0, 10, 20-50:10, last\""},
    {"output.pdb.times", "Times to export, e.g. \"1.5ns, 2000ps\" (no unit: inputs.time_unit)"},
    {"output.pdb.selection", "Group names (comma separated) or a selection (\"\": everything)"},
    {"output.pdb.prefix", "File name prefix"},
    {"output.write_params_log",
     "Write MDATK_params_log.toml (settings, versions, cache keys, citations) next to the "
     "outputs"},
    {"output.timestamp_files",
     "Time-stamp output file names, e.g. MDATK_video_2026-10-15_143205.mp4 (cache and log keep "
     "fixed names)"},
    {"workflow", "Project Planner and script wizard"},
    {"workflow.project_name", "Project name used in every generated file name"},
    {"workflow.stamp", "Time stamp in generated file names (fixed once files are written)"},
    {"workflow.use_stamp", "Put the time stamp into file names"},
    {"workflow.out_dir",
     "Folder of the generated files (\"\": <output.out_dir>/<project_name>_workflow)"},
    {"workflow.step_folders", "One sub-folder per step, e.g. 03_amber-build/"},
    // the keys MDATK 0.1.x added (its config.py; config_doc.json, of v0.0.4.1, does not list them)
    {"workflow.supersede_runs",
     "Before a step's files are written again, move its completed run of the plan's stamp (outputs, the "
     "scripts and inputs that made them, logs) into <step folder>/superseded/<stamp>_<time>/"},
    {"workflow.parallel",
     "The most scripts MDATK runs at the same time (Run all / Run selected), within the CPU and GPU budgets "
     "below; steps wait for their inputs (0: no limit; a project file without it: 1, as before 0.0.12)"},
    {"workflow.cpu_cores",
     "CPU cores the scripts MDATK runs at the same time may use together (0: half of them)"},
    {"workflow.gpus", "GPUs they may use, each script its own (-1: the GPUs found; 0: none)"},
    {"workflow.last_node",
     "The highest step id number given out so far (n1, n2, ...): a new step gets the next one, never the id "
     "of a deleted step (0: none recorded yet; ids then follow the highest one in the plan)"},
    {"workflow.nodes", "Planner steps"},
    {"workflow.nodes.id", "Unique id of the step"},
    {"workflow.nodes.tool", "Step type: the script generator (see 'mdatk workflow tools')"},
    {"workflow.nodes.step", "Step number; orders the steps and their file names"},
    {"workflow.nodes.title", "Box title (\"\": the step type's label)"},
    {"workflow.nodes.x", "Position of the box on the planner canvas"},
    {"workflow.nodes.y", "Position of the box on the planner canvas"},
    {"workflow.nodes.settings", "The step's options (each one is explained by its comment)"},
    {"workflow.nodes.files",
     "Hand-picked input files for ports that no connection feeds ({port = path})"},
    {"workflow.nodes.adopted",
     "Time stamp of an earlier run whose output files the plan uses for this step, which is then neither "
     "written nor run again (\"\": none; set by 'Check previous runs' or 'mdatk workflow status --adopt')"},
    {"workflow.nodes.folder_tag",
     "Added to the step's folder name, e.g. \"model12\": 02_amber-build_model12 instead of 02_amber-build, "
     "so the folder says which branch or model it belongs to (ASCII letters and digits, '-' between them; "
     "file names do not change)"},
    {"workflow.edges", "Planner connections"},
    {"workflow.edges.source", "Id of the producing step"},
    {"workflow.edges.output", "Output port of the producing step"},
    {"workflow.edges.target", "Id of the receiving step"},
    {"workflow.edges.input", "Input port of the receiving step"},
    {"workflow.scratch", "Stand-alone script wizard settings, one entry per step type"},
    {"workflow.scratch.id", "Unique id of the step"},
    {"workflow.scratch.tool", "Step type: the script generator (see 'mdatk workflow tools')"},
    {"workflow.scratch.step", "Step number; orders the steps and their file names"},
    {"workflow.scratch.title", "Box title (\"\": the step type's label)"},
    {"workflow.scratch.x", "Position of the box on the planner canvas"},
    {"workflow.scratch.y", "Position of the box on the planner canvas"},
    {"workflow.scratch.settings", "The step's options (each one is explained by its comment)"},
    {"workflow.scratch.files",
     "Hand-picked input files for ports that no connection feeds ({port = path})"},
    {"workflow.scratch.adopted",
     "Time stamp of an earlier run whose output files the plan uses for this step, which is then neither "
     "written nor run again (\"\": none; set by 'Check previous runs' or 'mdatk workflow status --adopt')"},
    {"workflow.scratch.folder_tag",
     "Added to the step's folder name, e.g. \"model12\": 02_amber-build_model12 instead of 02_amber-build, "
     "so the folder says which branch or model it belongs to (ASCII letters and digits, '-' between them; "
     "file names do not change)"},
};
