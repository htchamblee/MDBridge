#pragma once
// The parameter model of a project: MDATK's mdatk/config.py (v0.0.4.1), field for field.
//
// Every user-adjustable setting (GUI, CLI, project file, per-run parameter log) is a field of the
// tree rooted at Project. A project file is toPlain(project) plus a few meta entries, every entry
// preceded by its field's comment (the comments above the fields below are those texts).
// Reading goes TOML -> plain data -> projectFromPlain (type-checked; unknown keys are warned about
// and ignored, missing keys keep their defaults) -> optional path resolution.
//
// Field names are the project-file keys, hence snake_case. Defaults that contain "MDATK"
// (MDATK_output, MDATK_video.mp4, ...) are kept so that MDBridge writes the same files as MDATK.

#include <array>
#include <cstdint>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "mdbridge/core/value.hpp"

namespace mdbridge::core {

inline constexpr std::array<std::string_view, 7> kRepStyles{"cartoon", "ribbon", "tube", "licorice",
                                                            "vdw",     "surf",   "trace"};
inline constexpr std::array<std::string_view, 2> kRepRest{"trace", "hide"};
inline constexpr std::array<std::string_view, 4> kRepColors{"auto", "residue", "element", "nucleotide"};
// Panel titles that count as automatic: a panel that changes type gets the new type's title.
inline constexpr std::array<std::string_view, 5> kDefaultPanelTitles{
    "Structure", "RMSD", "Inter-strand H-bonds", "Aptamer residue contacts", "Mg occupancy"};

// [groups[]]
struct ResidueGroup {
    // Group name, used in selections (group:<name>), panels and legends
    std::string name = "Group";
    // Residue numbers or ranges, e.g. "1-40, 45"; may be chain-qualified ("A:1-40")
    std::string residues;
    // MDAnalysis selection that replaces 'residues' when given
    std::string selection;
    // Colour of the group in the structure view, legends and charts (#rrggbb)
    std::string color = "#2a78d6";

    bool operator==(const ResidueGroup&) const = default;
};

// [inputs.external[]]
struct ExternalSeries {
    // CSV or whitespace-delimited text file (e.g. a cpptraj .dat file)
    std::string path;
    // Name of the trace in panels and charts ("": the file name)
    std::string label;
    // Column holding the time axis: name or 0-based index ("": the first column)
    std::string x_column;
    // Unit of that column: "ps", "ns", "us", "frame" (0-based) or "frame1" (1-based)
    std::string x_unit = "ps";
    // Columns to show, by name or index (empty: all other columns)
    std::vector<std::string> columns;
    // Column separator ("": detect comma, tab or whitespace)
    std::string delimiter;

    bool operator==(const ExternalSeries&) const = default;
};

// [inputs]
struct InputParams {
    // Topology: Amber .parm7 / .prmtop (preferred), PDB, PSF, GRO, ...
    std::string topology;
    // Trajectory files in order; several files are concatenated
    std::vector<std::string> trajectories;
    // First analysed frame (0-based; of the cleaned trajectory when cleaning is on)
    std::int64_t start = 0;
    // Stop before this frame (0-based, exclusive; not set: up to the last frame)
    std::optional<std::int64_t> stop;
    // Analyse every n-th frame
    std::int64_t step = 1;
    // Time between stored frames in ps (not set: from the trajectory's time stamps)
    std::optional<double> dt_ps;
    // Unit shown on time axes: "ps", "ns", "us" or "frame"
    std::string time_unit = "ns";
    // External data files shown as extra traces (interpolated onto the analysed frames)
    std::vector<ExternalSeries> external;

    bool operator==(const InputParams&) const = default;
};

// [cleaning]
struct CleaningParams {
    // Strip, image and fit the trajectory with cpptraj before the analyses
    bool enabled = false;
    // cpptraj executable: a name on PATH or an absolute path
    std::string cpptraj_path = "cpptraj";
    // Conda environment providing cpptraj (a name such as AmberTools26 or a prefix path); every call
    // runs through 'conda run -n <env>'
    std::string conda_env;
    // conda / mamba / micromamba executable ("": found automatically)
    std::string conda_exe;
    // cpptraj mask of the solvent residues to strip
    std::string solvent_mask = ":WAT,HOH,TIP3,SOL,T3P,T4P,SPC";
    // cpptraj mask of the ions to strip (Mg2+ is kept by default)
    std::string ion_mask = ":Na+,Cl-,K+,NA,CL,K,Na,Cl,SOD,CLA,POT";
    // Further atoms to strip (cpptraj mask; optional)
    std::string extra_strip_mask;
    // Re-image molecules into the primary box (cpptraj autoimage)
    bool autoimage = true;
    // Anchor mask for autoimage ("": cpptraj's default, the first molecule)
    std::string autoimage_anchor;
    // RMS-fit every frame onto the reference
    bool fit = true;
    // cpptraj mask of the atoms used for the RMS fit
    std::string fit_mask = "@P,C4'";
    // Fit reference: the "first" frame, a given "frame" number or a "file"
    std::string fit_reference = "first";
    // Remove the box after fitting (a rotated box cannot be stored in NetCDF)
    bool remove_box_after_fit = true;
    // Reference frame for fit_reference = "frame" (1-based, cpptraj convention)
    std::int64_t reference_frame = 1;
    // Reference structure for fit_reference = "file"
    std::string reference_file;
    // First frame read by cpptraj (1-based trajin argument)
    std::int64_t start = 1;
    // Last frame read by cpptraj (1-based, inclusive; not set: the last frame)
    std::optional<std::int64_t> stop;
    // Read every n-th frame (cpptraj trajin offset)
    std::int64_t offset = 1;
    // Advanced: cpptraj commands inserted before trajout
    std::string extra_commands;
    // Folder for the cleaned files ("": <cache_dir>/cleaned)
    std::string output_dir;

    bool operator==(const CleaningParams&) const = default;
};

// [analyses.rmsd.targets[]]
struct RMSDTarget {
    // Name of the series (data key rmsd/<label>)
    std::string label = "Complex";
    // Atoms whose RMSD is computed (selection; group:<name> allowed)
    std::string selection = "(group:Aptamer or group:Target) and not name H*";
    // Atoms superposed before computing the RMSD ("": the same atoms)
    std::string fit_selection;

    bool operator==(const RMSDTarget&) const = default;
};

// [analyses.rmsd]
struct RMSDParams {
    // Compute RMSD time series
    bool enabled = true;
    // One RMSD series per target
    std::vector<RMSDTarget> targets{{"Complex", "(group:Aptamer or group:Target) and not name H*", ""}, {"Backbone", "(group:Aptamer or group:Target) and name C4'", ""}};
    // Reference structure: index into the analysed frames (0: the first)
    std::int64_t reference_frame = 0;
    // Superpose each frame onto the reference before computing the RMSD
    bool superposition = true;
    // Mass-weighted superposition and RMSD
    bool mass_weighted = false;

    bool operator==(const RMSDParams&) const = default;
};

// [analyses.contacts]
struct ContactParams {
    // Count heavy-atom contacts between two groups
    bool enabled = true;
    // First residue group
    std::string group_a = "Aptamer";
    // Second residue group
    std::string group_b = "Target";
    // Contact distance cutoff (Å)
    double cutoff = 4.5;
    // Ignore hydrogens
    bool heavy_only = true;

    bool operator==(const ContactParams&) const = default;
};

// [analyses.hbonds]
struct HBondParams {
    // Count hydrogen bonds between two groups
    bool enabled = true;
    // First residue group
    std::string group_a = "Aptamer";
    // Second residue group
    std::string group_b = "Target";
    // Maximum donor-acceptor distance (Å)
    double distance_cutoff = 3.5;
    // Also require a minimum D-H...A angle (needs hydrogens in the topology)
    bool use_angle = true;
    // Minimum D-H...A angle (degrees)
    double angle_cutoff = 135.0;
    // Elements that can donate (comma separated)
    std::string donor_elements = "N,O";
    // Elements that can accept (comma separated)
    std::string acceptor_elements = "N,O";

    bool operator==(const HBondParams&) const = default;
};

// [analyses.mg]
struct MgParams {
    // Analyse Mg2+ coordination: inner/outer-sphere contacts, bridges and occupancy
    bool enabled = true;
    // Selection of the Mg2+ ions
    std::string mg_selection = "resname MG Mg MG2 or name MG Mg";
    // Atoms that can coordinate Mg2+
    std::string ligand_selection = "group:Aptamer or group:Target";
    // Elements of the coordinating atoms (comma separated)
    std::string ligand_elements = "O";
    // First group of a Mg2+ bridge
    std::string group_a = "Aptamer";
    // Second group of a Mg2+ bridge
    std::string group_b = "Target";
    // Inner-sphere Mg-O distance cutoff (Å)
    double inner_cutoff = 2.6;
    // Outer-sphere (water-mediated) Mg-O distance cutoff (Å)
    double outer_cutoff = 4.5;
    // Shell that counts as a bridge: "inner" or "outer"
    std::string bridge_shell = "inner";
    // Require octahedral O-Mg-O angles (about 90 or 180 degrees)
    bool geometry = false;
    // Allowed deviation from the octahedral angles (degrees)
    double angle_tolerance = 25.0;
    // Atoms whose Mg2+ occupancy is reported ("": the ligand selection)
    std::string occupancy_selection;
    // Elements counted for the occupancy (comma separated)
    std::string occupancy_elements = "N,O";

    bool operator==(const MgParams&) const = default;
};

// [analyses.stacking]
struct StackingParams {
    // Detect base stacking from ring geometry
    bool enabled = true;
    // Residues whose rings are considered
    std::string selection = "group:Aptamer or group:Target";
    // "all" ring pairs, or "inter" (group A against group B only)
    std::string scope = "all";
    // First group for scope = "inter"
    std::string group_a = "Aptamer";
    // Second group for scope = "inter"
    std::string group_b = "Target";
    // Maximum ring-centroid distance (Å)
    double distance_cutoff = 4.5;
    // Maximum angle between the ring normals (degrees)
    double angle_cutoff = 30.0;
    // Maximum lateral displacement of the centroids (Å)
    double offset_cutoff = 2.5;
    // Maximum RMS deviation of ring atoms from their plane (Å)
    double planarity_cutoff = 0.25;

    bool operator==(const StackingParams&) const = default;
};

// [analyses.basepairs]
struct BasePairParams {
    // Base pairs of the residues in the structure view, for the secondary-structure panel and PDB
    // annotation; computed whenever a panel needs them. Engine: "builtin" (canonical pairs from H-bond
    // geometry) or "barnaba" (optional)
    std::string engine = "builtin";
    // Re-detect the pairs every n-th analysed frame (pairs hold in between)
    std::int64_t stride = 10;
    // Builtin: maximum donor-acceptor distance of a base-pair H-bond (Å)
    double cutoff = 3.5;
    // Builtin: also count G-U wobble pairs
    bool wobble = true;
    // Barnaba: also keep non-canonical (non Watson-Crick) pairs
    bool noncanonical = false;

    bool operator==(const BasePairParams&) const = default;
};

// [analyses.energy]
struct EnergyParams {
    // Non-bonded interaction energy between two groups
    bool enabled = true;
    // First residue group
    std::string group_a = "Aptamer";
    // Second residue group
    std::string group_b = "Target";
    // Pair distance cutoff (Å)
    double cutoff = 12.0;
    // Relative dielectric constant
    double dielectric = 1.0;
    // Distance-dependent dielectric (epsilon = dielectric x r)
    bool distance_dependent = false;
    // Use the exact Lennard-Jones parameters of a parm7 topology
    bool use_parm7 = true;
    // Reference energy drawn on charts and panels (kcal/mol; not set: none)
    std::optional<double> target_value;
    // Exclude pairs within n bonds (0: exclude none)
    std::int64_t exclude_bonded = 3;

    bool operator==(const EnergyParams&) const = default;
};

// [analyses]
struct AnalysisParams {
    // RMSD time series
    RMSDParams rmsd;
    // Heavy-atom contacts between two groups
    ContactParams contacts;
    // Hydrogen bonds between two groups
    HBondParams hbonds;
    // Mg2+ coordination, bridges and occupancy
    MgParams mg;
    // Base stacking
    StackingParams stacking;
    // Non-bonded interaction energy
    EnergyParams energy;
    // Base pairs (secondary structure)
    BasePairParams basepairs;
    // Minimum-image distances (cleaned, imaged systems rarely need it)
    bool use_pbc = false;
    // Largest proximity cutoff usable later without re-analysing (Å)
    double proximity_search_radius = 8.0;

    bool operator==(const AnalysisParams&) const = default;
};

// [visual.representations[]]
struct RepSpec {
    // Representation: cartoon, ribbon, tube, licorice, vdw, surf or trace
    std::string style = "cartoon";
    // Residue groups drawn (empty: every shown group)
    std::vector<std::string> groups;
    // Colouring: auto, residue, element or nucleotide
    std::string color = "auto";
    // Scales widths and radii
    double size = 1.0;
    // Opacity, 0 (invisible) to 1
    double opacity = 1.0;
    // Cartoon only: draw the ribose rings (shows the pucker)
    bool sugars = false;

    bool operator==(const RepSpec&) const = default;
};

// [visual]
struct VisualParams {
    // Residue groups shown in the structure view (empty: all)
    std::vector<std::string> structure_groups;
    // Atom used for each residue's 2D trace point, first one found wins (P, C4', C1', CA ...)
    std::string trace_atoms = "C4',CA";
    // 3D representations (empty: the 2D backbone trace)
    std::vector<RepSpec> representations;
    // 3D renderer: "auto" / "software" (built-in) or "vtk" (optional)
    std::string renderer = "auto";
    // Residues no representation draws: "trace" (3D trace) or "hide"
    std::string rep_rest = "trace";
    // Spline control atom of ribbons, cartoons and traces: "P" or "C4'"
    std::string backbone_atom = "P";
    // Include hydrogens in all-atom representations
    bool rep_hydrogens = false;
    // Show Mg2+ ions in the structure view
    bool show_ions = true;
    // Superpose each frame's trace onto the reference frame
    bool align_frames = true;
    // View rotation in degrees about X, Y and Z (R = Rz Ry Rx)
    std::vector<double> rotation{0.0, 0.0, 0.0};
    // Magnification of the structure view (1: the whole model fits; 2: twice as large)
    double zoom = 1.0;
    // View centre: "" (centre of the model), "origin" (0, 0, 0 of the coordinates) or a residue, e.g.
    // "12" or "B:12"
    std::string center_on;
    // Keep the centre residue in the middle in every frame (not only the first)
    bool center_follow = false;
    // Residue colouring: "chain", "groups" or "proximity"
    std::string color_mode = "groups";
    // Colour of residues no group or proximity colour applies to
    std::string base_color = "#b5b3ad";
    // Matplotlib palette for chain colouring
    std::string chain_palette = "tab10";
    // Group (or sel:<selection>) whose proximity colours and glows residues
    std::string proximity_set = "Target";
    // Distance that counts as proximity (Å)
    double proximity_cutoff = 4.0;
    // Colour of residues in proximity
    std::string proximity_color = "#e34948";
    // Glow of residues near the proximity set: "none", "ever" (whole movie), "fade" or "current"
    std::string glow_mode = "fade";
    // Glow colour
    std::string glow_color = "#ffb000";
    // Frames over which the glow fades (glow_mode = "fade")
    std::int64_t fade_frames = 15;
    // Glow intensity factor
    double glow_strength = 1.0;
    // Colour of Mg2+ ions
    std::string ion_color = "#008300";
    // Line width of the 2D trace (points)
    double line_width = 2.2;
    // Marker size of the 2D trace points
    double point_size = 22.0;
    // Fade distant parts toward the background colour
    bool depth_cue = true;
    // "light" or "dark"
    std::string theme = "light";
    // Base font size (points)
    double font_size = 11.0;
    // Show the time label in the structure view
    bool show_time_label = true;
    // Residue-number tags on every n-th residue (0: off)
    std::int64_t label_every = 0;
    // Explicit residue numbers / ranges to tag, e.g. "1, 12, 24"
    std::string label_residues;
    // Figure title ("": none)
    std::string title;

    bool operator==(const VisualParams&) const = default;
};

// [density]
struct DensityParams {
    // Density map file (.mrc / .map / .ccp4)
    std::string map_path;
    // Draw the map shadow
    bool show = true;
    // Resolution in Å; width of the simulated map used for fitting
    double resolution = 5.0;
    // Contour level in map standard deviations above the mean
    double level = 2.0;
    // Shadow colour
    std::string color = "#2a78d6";
    // Shadow opacity, 0 (invisible) to 1
    double opacity = 0.45;
    // Draw the shadow in front of 3D representations (else behind them)
    bool front = false;
    // "ref" (the first analysed frame) or a trajectory frame number to fit on
    std::string fit_frame = "ref";
    // Automatic fit: "global" search or "manual" (refine the stored placement)
    std::string search = "global";
    // Degrees between sampled orientations of the global search
    double angular_step = 30.0;
    // Fit before rendering when no placement is stored
    bool auto_fit = true;
    // Placement: 3x3 rotation (row-major) + translation, model -> map
    std::vector<double> transform;
    // Backbone points of the reference frame the placement was made on (the map follows the molecule
    // when the reference changes: other start frame, cleaned vs raw files)
    std::vector<std::vector<double>> anchor;
    // Correlation of the stored placement (0: not fitted)
    double fit_cc = 0.0;
    // Summary of the last fit
    std::string fit_info;

    bool operator==(const DensityParams&) const = default;
};

// [layout.panels[]]
struct PanelSpec {
    // Panel type: structure, timeseries, dual_axis, heatmap or contact_map
    std::string type = "structure";
    // Data keys shown, e.g. rmsd/Complex or contacts/per_residue_a
    std::vector<std::string> series;
    // Panel title ("": automatic)
    std::string title;
    // Panel-specific options, e.g. axis limits or the contact-map mode
    Table options;

    bool operator==(const PanelSpec&) const = default;
};

// [layout]
struct LayoutParams {
    // Grid layout (see 'mdatk list')
    std::string preset = "large_plus_3";
    // Panels in slot order
    std::vector<PanelSpec> panels{{"structure", {}, "Structure", {}}, {"timeseries", {"rmsd/Complex"}, "RMSD", {}}, {"timeseries", {"hbonds/count"}, "Inter-strand H-bonds", {}}, {"contact_map", {"contacts/map"}, "", Table{{"mode", "instant"}, {"transpose", true}}}};

    bool operator==(const LayoutParams&) const = default;
};

// [output.video]
struct VideoParams {
    // Write the MP4 video
    bool enabled = true;
    // Video file name (time-stamped when output.timestamp_files is on)
    std::string filename = "MDATK_video.mp4";
    // Frames per second
    std::int64_t fps = 24;
    // Width in pixels
    std::int64_t width = 1920;
    // Height in pixels
    std::int64_t height = 1080;
    // Figure resolution; sets text and line sizes relative to the frame
    std::int64_t dpi = 120;
    // FFmpeg video codec
    std::string codec = "libx264";
    // Constant rate factor: lower means better quality and larger files (libx264: 18-28)
    std::int64_t crf = 20;
    // Use every n-th analysed frame
    std::int64_t frame_stride = 1;
    // FFmpeg executable
    std::string ffmpeg_path = "ffmpeg";

    bool operator==(const VideoParams&) const = default;
};

// [output.charts]
struct ChartParams {
    // Write high-resolution time charts
    bool enabled = true;
    // Data keys to chart (empty: every scalar series)
    std::vector<std::string> series;
    // Also chart per-residue heatmaps
    bool include_heatmaps = false;
    // File formats, e.g. png, pdf, svg
    std::vector<std::string> formats{"png"};
    // Raster resolution
    std::int64_t dpi = 300;
    // Width in inches
    double width = 8.0;
    // Height in inches
    double height = 4.5;
    // "mdatk" (built-in look) or any matplotlib style name
    std::string style = "mdatk";
    // Line colour ("": the first palette colour)
    std::string color;
    // Line width (points)
    double line_width = 1.5;
    // Chart title ("": the series label; "{label}" placeholders allowed)
    std::string title;
    // x-axis label ("": automatic)
    std::string xlabel;
    // y-axis label ("": automatic)
    std::string ylabel;
    // Running-average window in frames (0: off)
    std::int64_t smoothing = 0;
    // Draw the energy target value when one is set
    bool show_target = true;
    // File name prefix
    std::string prefix = "MDATK_chart";

    bool operator==(const ChartParams&) const = default;
};

// [output.pdb]
struct PDBParams {
    // Export selected frames as PDB files
    bool enabled = false;
    // Frames to export, e.g. "0, 10, 20-50:10, last"
    std::string frames;
    // Times to export, e.g. "1.5ns, 2000ps" (no unit: inputs.time_unit)
    std::string times;
    // Group names (comma separated) or a selection ("": everything)
    std::string selection;
    // File name prefix
    std::string prefix = "MDATK_frame";

    bool operator==(const PDBParams&) const = default;
};

// [output]
struct OutputParams {
    // Output folder (a relative path is relative to the project file)
    std::string out_dir = "MDATK_output";
    // Cache of analyses and cleaned trajectories ("": <out_dir>/MDATK_cache)
    std::string cache_dir;
    // MP4 video
    VideoParams video;
    // Time charts
    ChartParams charts;
    // PDB export of selected frames
    PDBParams pdb;
    // Write MDATK_params_log.toml (settings, versions, cache keys, citations) next to the outputs
    bool write_params_log = true;
    // Time-stamp output file names, e.g. MDATK_video_2026-10-15_143205.mp4 (cache and log keep fixed
    // names)
    bool timestamp_files = true;

    bool operator==(const OutputParams&) const = default;
};

// [workflow.nodes[], workflow.scratch[]]
struct PlanNode {
    // Unique id of the step
    std::string id;
    // Step type: the script generator (see 'mdatk workflow tools')
    std::string tool = "note";
    // Step number; orders the steps and their file names
    std::int64_t step = 1;
    // Box title ("": the step type's label)
    std::string title;
    // Position of the box on the planner canvas
    double x = 40.0;
    // Position of the box on the planner canvas
    double y = 40.0;
    // The step's options (each one is explained by its comment)
    Table settings;
    // Hand-picked input files for ports that no connection feeds ({port = path})
    Table files;

    bool operator==(const PlanNode&) const = default;
};

// [workflow.edges[]]
struct PlanEdge {
    // Id of the producing step
    std::string source;
    // Output port of the producing step
    std::string output;
    // Id of the receiving step
    std::string target;
    // Input port of the receiving step
    std::string input;

    bool operator==(const PlanEdge&) const = default;
};

// [workflow]
struct WorkflowParams {
    // Project name used in every generated file name
    std::string project_name = "project";
    // Time stamp in generated file names (fixed once files are written)
    std::string stamp;
    // Put the time stamp into file names
    bool use_stamp = true;
    // Folder of the generated files ("": <output.out_dir>/<project_name>_workflow)
    std::string out_dir;
    // One sub-folder per step, e.g. 03_amber-build/
    bool step_folders = true;
    // Planner steps
    std::vector<PlanNode> nodes;
    // Planner connections
    std::vector<PlanEdge> edges;
    // Stand-alone script wizard settings, one entry per step type
    std::vector<PlanNode> scratch;

    bool operator==(const WorkflowParams&) const = default;
};

// [<root>]
struct Project {
    // Free-text notes about the project (what, why, decisions, to-dos)
    std::string notes;
    // Input files and frame range
    InputParams inputs;
    // Trajectory cleaning with cpptraj
    CleaningParams cleaning;
    // Residue groups
    std::vector<ResidueGroup> groups{{"Aptamer", "1-12", "", "#2a78d6"}, {"Target", "13-24", "", "#eb6834"}};
    // Analyses
    AnalysisParams analyses;
    // Structure view: colouring, glow, 3D representations, labels
    VisualParams visual;
    // Density map shadow and fit
    DensityParams density;
    // Figure layout: grid and panels
    LayoutParams layout;
    // Outputs: folders, video, charts, PDB frames
    OutputParams output;
    // Project Planner and script wizard
    WorkflowParams workflow;

    bool operator==(const Project&) const = default;
};

// ---- plain data ------------------------------------------------------------------------------

// to_dict(): fields in declaration order, unset optionals as Null; Null entries of the free-form
// tables (options, settings, files) dropped at every depth.
Value toPlain(const Project& project);

// Project.from_dict(data, base_dir): the meta entries (format, format_version, mdatk_version,
// saved, log, log_file) are ignored, the rest is converted with MDATK's rules and messages, and,
// with a non-empty baseDir, relative paths are resolved against it (resolvePaths).
//
// Throws ProjectError with MDATK's message for the first entry (in the order of `data`) that has
// the wrong type, e.g. "project.analyses.contacts.cutoff: expected a number, got 'far'".
// Unknown keys are reported to `warnings` (or, when it is null, logged as warnings), e.g.
// "project.inputs: ignoring unknown key 'bogus'". Aliases: [input] for [inputs], and
// inputs.trajectory (a string or a list) for inputs.trajectories.
//
// Deviations, where MDATK has no usable behaviour: an integer field given an integral number
// beyond the 64-bit range is an "expected an integer" error (Python keeps a big integer);
// NaN / infinity in an integer field give ProjectError("cannot convert float NaN to integer") /
// ("cannot convert float infinity to integer") (Python raises ValueError / OverflowError).
Project projectFromPlain(const Value& data, std::string_view baseDir = {},
                         std::vector<std::string>* warnings = nullptr);

// The comment of a model entry in project files: `path` is the dot-joined path from the root
// with list indices left out ("analyses.contacts.cutoff", "groups.name", "layout.panels").
// std::nullopt for anything that is not a field of the model.
std::optional<std::string_view> modelComment(std::string_view path);

// ---- paths -----------------------------------------------------------------------------------

// Project.resolve_paths(base_dir): makes the file paths of the project absolute against baseDir,
// in place and in MDATK's order: inputs.topology, inputs.trajectories, inputs.external[].path,
// cleaning.reference_file, cleaning.output_dir, cleaning.cpptraj_path (only when it contains a
// "/"), output.out_dir, output.cache_dir, density.map_path, workflow.out_dir, then the files and
// path settings of every planner node and scratch node (resolveNodePaths). "~" is expanded;
// absolute paths are left as they are; relative ones become normpath(join(baseDir, path)).
// POSIX path rules (posixpath), as MDATK on Linux.
void resolvePaths(Project& project, std::string_view baseDir);

// resolve_node_paths: each files value is split at ";", the pieces stripped, empty ones dropped,
// the rest resolved and joined with ";"; settings named by the step type's pathSettings are
// resolved when they are non-empty strings. A files value that is not a string throws
// ProjectError("<where>: expected a file path, got <repr>") (MDATK crashes there).
void resolveNodePaths(PlanNode& node, std::string_view baseDir, std::string_view where);

// fix(p) of resolve_paths: "" stays "", "~" is expanded, an absolute path is returned unchanged, a
// relative one becomes normpath(join(baseDir, p)).
std::string resolvePath(std::string_view path, std::string_view baseDir);

// Folders, absolute and normalised against the current directory when relative:
// output.out_dir or "MDATK_output"; output.cache_dir or <out_dir>/MDATK_cache;
// workflow.out_dir or <out_dir>/<sanitizeProject(workflow.project_name)>_workflow.
std::string outDir(const Project& project);
std::string cacheDir(const Project& project);
std::string workflowDir(const Project& project);

// ---- groups and checks -----------------------------------------------------------------------

std::vector<std::string> groupNames(const Project& project);

// The first group with exactly this name. Throws std::out_of_range("no residue group named 'X'
// (defined: Aptamer, Target)") ("(defined: none)" without groups).
const ResidueGroup& group(const Project& project, std::string_view name);

// A density map is set, shown, and a structure panel exists.
bool densityActive(const Project& project);

// Every problem of the project in MDATK's order and words (empty when it can be run). Files are
// looked up relative to the current directory. The defaults give
// {"No topology file selected.", "No trajectory file selected."}.
std::vector<std::string> validate(const Project& project);

// The density-map part of validate().
std::vector<std::string> densityProblems(const Project& project);

// What is wrong with a stored rigid transform (12 numbers: 3x3 rotation row-major + translation),
// in the words of MDATK's RigidTransform.from_list, or std::nullopt for a valid one.
std::optional<std::string> rigidTransformProblem(std::span<const double> numbers);

// Hints about the figure layout (more than one structure panel).
std::vector<std::string> layoutWarnings(const Project& project);

}  // namespace mdbridge::core
