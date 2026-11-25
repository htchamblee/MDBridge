#include "mdbridge/core/project.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <stdexcept>
#include <system_error>
#include <tuple>
#include <type_traits>
#include <utility>

#include "mdbridge/core/log.hpp"
#include "mdbridge/core/naming.hpp"
#include "mdbridge/core/ntpath.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/project_file.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/core/rigid_transform.hpp"
#include "mdbridge/core/step_types.hpp"
#include "mdbridge/core/text.hpp"

namespace mdbridge::core {

namespace {

#include "project_model.inc"

template <class T>
concept Model = requires(const T* p) { fieldsOf(p); };

template <class S, class F>
void forEachField(S& s, F&& f) {
    using T = std::remove_const_t<S>;
    std::apply([&](const auto&... fd) { (f(fd.name, s.*(fd.member)), ...); },
               fieldsOf(static_cast<const T*>(nullptr)));
}

// Calls f(member) for the field named `key`; false when there is none.
template <class S, class F>
bool withField(S& s, std::string_view key, F&& f) {
    return std::apply(
        [&](const auto&... fd) {
            return ((fd.name == key ? (f(s.*(fd.member)), true) : false) || ...);
        },
        fieldsOf(static_cast<const S*>(nullptr)));
}

std::string repr(std::string_view text) { return pythonRepr(Value(std::string(text))); }

// _to_plain of a free-form table: entries holding None are left out, at every depth.
Value withoutNulls(const Value& value) {
    if (const Table* table = value.get<Table>()) {
        Table out;
        for (const auto& [key, item] : *table) {
            if (!item.isNull()) {
                out.emplace_back(key, withoutNulls(item));
            }
        }
        return out;
    }
    if (const Array* array = value.get<Array>()) {
        Array out;
        out.reserve(array->size());
        for (const Value& item : *array) {
            out.push_back(withoutNulls(item));
        }
        return out;
    }
    return value;
}

struct Writer {
    static Value plain(const std::string& v) { return v; }
    static Value plain(std::int64_t v) { return v; }
    static Value plain(double v) { return v; }
    static Value plain(bool v) { return v; }
    static Value plain(const Table& v) { return withoutNulls(Value(v)); }

    template <class T>
    static Value plain(const std::optional<T>& v) {
        return v ? plain(*v) : Value();
    }

    template <class T>
    static Value plain(const std::vector<T>& v) {
        Array out;
        out.reserve(v.size());
        for (const T& item : v) {
            out.push_back(plain(item));
        }
        return out;
    }

    template <Model T>
    static Value plain(const T& s) {
        Table out;
        forEachField(s, [&](std::string_view name, const auto& member) {
            out.emplace_back(std::string(name), plain(member));
        });
        return out;
    }
};

// _from_plain: the type checks and messages of config.py.
class Reader {
public:
    explicit Reader(std::vector<std::string>* warnings) : m_warnings(warnings) {}

    void read(const Value& v, const std::string& /*where*/, std::string& out) {
        out = v.isNull() ? std::string() : pythonStr(v);
    }

    void read(const Value& v, const std::string& where, double& out) {
        if (const auto* i = v.get<std::int64_t>()) {
            out = static_cast<double>(*i);
        } else if (const auto* d = v.get<double>()) {
            out = *d;
        } else {
            fail(where + ": expected a number, got " + pythonRepr(v));
        }
    }

    void read(const Value& v, const std::string& where, std::int64_t& out) {
        if (const auto* i = v.get<std::int64_t>()) {
            out = *i;
            return;
        }
        const double* d = v.get<double>();
        if (!d) {
            fail(where + ": expected an integer, got " + pythonRepr(v));
        }
        if (std::isnan(*d)) {
            fail("cannot convert float NaN to integer");
        }
        if (std::isinf(*d)) {
            fail("cannot convert float infinity to integer");
        }
        // [-2^63, 2^63): the integral doubles an int64 holds
        if (std::trunc(*d) != *d || *d < -9223372036854775808.0 || *d >= 9223372036854775808.0) {
            fail(where + ": expected an integer, got " + pythonRepr(v));
        }
        out = static_cast<std::int64_t>(*d);
    }

    void read(const Value& v, const std::string& where, bool& out) {
        const bool* b = v.get<bool>();
        if (!b) {
            fail(where + ": expected true/false, got " + pythonRepr(v));
        }
        out = *b;
    }

    void read(const Value& v, const std::string& where, Table& out) {
        const Table* table = v.get<Table>();
        if (!table) {
            fail(where + ": expected an object");
        }
        out = *table;
    }

    template <class T>
    void read(const Value& v, const std::string& where, std::optional<T>& out) {
        if (v.isNull()) {
            out.reset();
            return;
        }
        T value{};
        read(v, where, value);
        out = std::move(value);
    }

    template <class T>
    void read(const Value& v, const std::string& where, std::vector<T>& out) {
        if (v.isNull()) {
            out.clear();
            return;
        }
        const Array* array = v.get<Array>();
        if (!array) {
            fail(where + ": expected a list");
        }
        std::vector<T> items(array->size());
        for (std::size_t i = 0; i < array->size(); ++i) {
            read((*array)[i], where + "[" + std::to_string(i) + "]", items[i]);
        }
        out = std::move(items);
    }

    template <Model T>
    void read(const Value& v, const std::string& where, T& out) {
        const Table* table = v.get<Table>();
        if (!table) {
            fail(where + ": expected an object, got " + std::string(pythonTypeName(v)));
        }
        T result{};
        for (const auto& [given, value] : *table) {
            const std::string_view key = alias<T>(given);
            const std::string at = where + "." + std::string(key);
            bool oneTrajectory = false;  // a string for inputs.trajectories: a one-item list
            if constexpr (std::is_same_v<T, InputParams>) {
                oneTrajectory = key == "trajectories" && value.isString();
            }
            const bool known = withField(result, key, [&](auto& member) {
                read(oneTrajectory ? Value(Array{value}) : value, at, member);
            });
            if (!known) {
                warn(where + ": ignoring unknown key " + repr(key) + unknownKeyHint<T>(key));
            }
        }
        out = std::move(result);
    }

private:
    // The folders of Project that are properties, not fields (PRJ-07): name the field to set.
    template <class T>
    static std::string unknownKeyHint([[maybe_unused]] std::string_view key) {
        if constexpr (std::is_same_v<T, Project>) {
            if (key == "out_dir" || key == "cache_dir") {
                return " (did you mean output." + std::string(key) + "?)";
            }
            if (key == "workflow_dir") {
                return " (did you mean workflow.out_dir?)";
            }
        }
        return {};
    }

    template <class T>
    static std::string_view alias(std::string_view key) {
        if constexpr (std::is_same_v<T, InputParams>) {
            return key == "trajectory" ? "trajectories" : key;
        } else if constexpr (std::is_same_v<T, Project>) {
            return key == "input" ? "inputs" : key;
        } else {
            return key;
        }
    }

    [[noreturn]] static void fail(const std::string& message) { throw ProjectError(message); }

    void warn(const std::string& message) const {
        if (m_warnings) {
            m_warnings->push_back(message);
        } else {
            logWarning(message);
        }
    }

    std::vector<std::string>* m_warnings;
};

constexpr std::array<std::string_view, 6> kMetaKeys{"format", "format_version", "mdatk_version",
                                                    "saved",  "log",            "log_file"};

// os.path.exists, relative to the current directory
bool exists(const std::string& path) {
    if (path.empty() || path.find('\0') != std::string::npos) {
        return false;
    }
    std::error_code ec;
    return std::filesystem::exists(pathFromUtf8(path), ec);
}

template <std::size_t N>
bool oneOf(std::string_view value, const std::array<std::string_view, N>& options) {
    return std::find(options.begin(), options.end(), value) != options.end();
}

bool oneOf(std::string_view value, std::initializer_list<std::string_view> options) {
    return std::find(options.begin(), options.end(), value) != options.end();
}

bool contains(const std::vector<std::string>& names, const std::string& name) {
    return std::find(names.begin(), names.end(), name) != names.end();
}

std::string joined(const std::vector<std::string>& items, std::string_view separator) {
    std::string out;
    for (std::size_t i = 0; i < items.size(); ++i) {
        if (i) {
            out += separator;
        }
        out += items[i];
    }
    return out;
}

template <std::size_t N>
std::string joined(const std::array<std::string_view, N>& items, std::string_view separator) {
    return joined(std::vector<std::string>(items.begin(), items.end()), separator);
}

// re.fullmatch(r"[A-Za-z0-9_.\-]+", name)
bool isPlainName(std::string_view name) {
    return !name.empty() && std::all_of(name.begin(), name.end(), [](char c) {
        return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9') ||
               c == '_' || c == '.' || c == '-';
    });
}

// str.isdigit(): Unicode digits, decimal or not ("²", "①" ...)
bool isDigitString(std::string_view text) {
    const std::u32string codes = decodeUtf8(text);
    return !codes.empty() && std::all_of(codes.begin(), codes.end(), [](char32_t c) {
        const auto inRange = [c](const auto& range) {
            return range.first <= c && c <= range.second;
        };
        return isPythonDigit(c) ||
               std::any_of(std::begin(kDigitNotDecimal), std::end(kDigitNotDecimal), inRange);
    });
}

std::string asciiLower(std::string text) {
    for (char& c : text) {
        if (c >= 'A' && c <= 'Z') {
            c = static_cast<char>(c - 'A' + 'a');
        }
    }
    return text;
}

}  // namespace

// ---- plain data ------------------------------------------------------------------------------

Value toPlain(const Project& project) { return Writer::plain(project); }

Project projectFromPlain(const Value& data, std::string_view baseDir,
                         std::vector<std::string>* warnings) {
    Value input = data;
    if (Table* table = input.get<Table>()) {
        std::erase_if(*table, [](const auto& entry) { return oneOf(entry.first, kMetaKeys); });
    }
    Project project;
    Reader(warnings).read(input, "project", project);
    if (!baseDir.empty()) {
        resolvePaths(project, baseDir);
    }
    return project;
}

std::optional<std::string_view> modelComment(std::string_view path) {
    for (const auto& [key, text] : kModelComments) {
        if (key == path) {
            return text;
        }
    }
    return std::nullopt;
}

// ---- paths -----------------------------------------------------------------------------------

namespace {

// os.path.abspath against the current directory: POSIX rules (as MDATK); on Windows MDBridge's
// Windows rules (ntpath.hpp), in generic form.
std::string absoluteFolder(const std::string& path) {
#ifdef _WIN32
    return ntpath::abspath(path, currentDirectory());
#else
    return posixpath::abspath(path, currentDirectory());
#endif
}

}  // namespace

std::string resolvePath(std::string_view path, std::string_view baseDir) {
#ifdef _WIN32
    return ntpath::resolve(path, baseDir);
#else
    if (path.empty()) {
        return {};
    }
    std::string expanded = posixpath::expanduser(path);
    if (posixpath::isabs(expanded)) {
        return expanded;
    }
    return posixpath::normpath(posixpath::join(baseDir, expanded));
#endif
}

void resolveNodePaths(PlanNode& node, std::string_view baseDir, std::string_view where) {
    Table files;
    for (const auto& [port, value] : node.files) {
        const std::string* text = value.get<std::string>();
        if (!text) {
            throw ProjectError(std::string(where) + ".files." + port +
                               ": expected a file path, got " + pythonRepr(value));
        }
        std::vector<std::string> pieces;
        std::size_t start = 0;
        for (;;) {
            const std::size_t end = text->find(';', start);
            const std::string piece = pythonStrip(std::string_view(*text).substr(
                start, end == std::string::npos ? std::string::npos : end - start));
            if (!piece.empty()) {
                pieces.push_back(resolvePath(piece, baseDir));
            }
            if (end == std::string::npos) {
                break;
            }
            start = end + 1;
        }
        files.emplace_back(port, joined(pieces, ";"));
    }
    node.files = std::move(files);
    if (const StepType* type = findStepType(node.tool)) {
        for (const std::string& key : type->pathSettings) {
            for (auto& [name, value] : node.settings) {
                if (name == key) {
                    const std::string* text = value.get<std::string>();
                    if (text && !text->empty()) {
                        value = resolvePath(*text, baseDir);
                    }
                    break;
                }
            }
        }
    }
}

void resolvePaths(Project& p, std::string_view baseDir) {
    const auto fix = [baseDir](std::string& path) { path = resolvePath(path, baseDir); };
    fix(p.inputs.topology);
    for (std::string& path : p.inputs.trajectories) {
        fix(path);
    }
    for (ExternalSeries& series : p.inputs.external) {
        fix(series.path);
    }
    fix(p.cleaning.reference_file);
    fix(p.cleaning.output_dir);
#ifdef _WIN32
    const bool cpptrajIsAPath = p.cleaning.cpptraj_path.find_first_of("/\\") != std::string::npos;
#else
    const bool cpptrajIsAPath = p.cleaning.cpptraj_path.find('/') != std::string::npos;
#endif
    if (cpptrajIsAPath) {
        fix(p.cleaning.cpptraj_path);
    }
    fix(p.output.out_dir);
    fix(p.output.cache_dir);
    fix(p.density.map_path);
    fix(p.workflow.out_dir);
    for (std::size_t i = 0; i < p.workflow.nodes.size(); ++i) {
        resolveNodePaths(p.workflow.nodes[i], baseDir,
                         "project.workflow.nodes[" + std::to_string(i) + "]");
    }
    for (std::size_t i = 0; i < p.workflow.scratch.size(); ++i) {
        resolveNodePaths(p.workflow.scratch[i], baseDir,
                         "project.workflow.scratch[" + std::to_string(i) + "]");
    }
}

std::string outDir(const Project& project) {
    const std::string& dir = project.output.out_dir;
    return absoluteFolder(dir.empty() ? std::string("MDATK_output") : dir);
}

std::string cacheDir(const Project& project) {
    if (!project.output.cache_dir.empty()) {
        return absoluteFolder(project.output.cache_dir);
    }
    return posixpath::join(outDir(project), "MDATK_cache");
}

std::string workflowDir(const Project& project) {
    if (!project.workflow.out_dir.empty()) {
        return absoluteFolder(project.workflow.out_dir);
    }
    return posixpath::join(outDir(project),
                           sanitizeProject(project.workflow.project_name) + "_workflow");
}

// ---- groups and checks -----------------------------------------------------------------------

std::vector<std::string> groupNames(const Project& project) {
    std::vector<std::string> names;
    names.reserve(project.groups.size());
    for (const ResidueGroup& g : project.groups) {
        names.push_back(g.name);
    }
    return names;
}

const ResidueGroup& group(const Project& project, std::string_view name) {
    for (const ResidueGroup& g : project.groups) {
        if (g.name == name) {
            return g;
        }
    }
    const std::string defined = joined(groupNames(project), ", ");
    throw std::out_of_range("no residue group named " + repr(name) + " (defined: " +
                            (defined.empty() ? std::string("none") : defined) + ")");
}

bool densityActive(const Project& project) {
    const DensityParams& d = project.density;
    return !d.map_path.empty() && d.show &&
           std::any_of(project.layout.panels.begin(), project.layout.panels.end(),
                       [](const PanelSpec& panel) { return panel.type == "structure"; });
}

std::optional<std::string> rigidTransformProblem(std::span<const double> numbers) {
    try {
        RigidTransform::fromList(numbers);
    } catch (const ValueError& error) {
        return std::string(error.what());
    }
    return std::nullopt;
}

std::vector<std::string> densityProblems(const Project& project) {
    const DensityParams& d = project.density;
    std::vector<std::string> out;
    if (d.map_path.empty()) {
        return out;
    }
    if (densityActive(project) && !exists(d.map_path)) {
        out.push_back("Density map not found: " + d.map_path +
                      " (untick 'show' to render without it)");
    }
    if (d.resolution <= 0) {
        out.push_back("Density map resolution must be positive.");
    }
    if (!(0.0 <= d.opacity && d.opacity <= 1.0)) {
        out.push_back("Density map opacity must be between 0 and 1.");
    }
    if (!(5.0 <= d.angular_step && d.angular_step <= 45.0)) {
        out.push_back("Density search angular step must be between 5 and 45 degrees.");
    }
    if (!oneOf(d.search, {"global", "manual"})) {
        out.push_back("Unknown density fit start " + repr(d.search) + " (global or manual).");
    }
    if (!d.transform.empty()) {
        if (const auto problem = rigidTransformProblem(d.transform)) {
            out.push_back("The stored density map placement is invalid (" + *problem +
                          "); clear it or fit again.");
        }
    }
    const std::string frame = asciiLower(pythonStrip(d.fit_frame));
    if (!frame.empty() && frame != "ref" && !isDigitString(frame)) {
        out.push_back("Density fit frame must be 'ref' or a frame number, not " +
                      repr(d.fit_frame) + ".");
    }
    return out;
}

namespace {

// _out_key of the RMSD analysis: re.sub(r"[^A-Za-z0-9_.+\-]+", "_", label).strip("_") or "target<k>"
std::string rmsdDataKey(std::string_view label, std::size_t k) {
    std::string key;
    bool inRun = false;
    for (const char c : label) {
        const bool keep = (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9') ||
                          c == '_' || c == '.' || c == '+' || c == '-';
        if (keep) {
            key += c;
        } else if (!inRun) {
            key += '_';
        }
        inRun = !keep;
    }
    const std::size_t first = key.find_first_not_of('_');
    if (first == std::string::npos) {
        return "target" + std::to_string(k);
    }
    return key.substr(first, key.find_last_not_of('_') - first + 1);
}

}  // namespace

std::vector<std::string> analysisParameterProblems(const Project& project) {
    std::vector<std::string> out;
    const AnalysisParams& a = project.analyses;
    // a cutoff, radius or tolerance: a finite number > 0 (NaN fails every comparison)
    const auto positive = [&out](double value, std::string_view what) {
        if (!(value > 0)) {
            out.push_back(std::string(what) + " must be positive.");
        } else if (!std::isfinite(value)) {
            out.push_back(std::string(what) + " must be finite.");
        }
        return value > 0 && std::isfinite(value);
    };
    const auto angle = [&out](double value, std::string_view what) {
        if (!(0 <= value && value <= 180)) {
            out.push_back(std::string(what) + " must be between 0 and 180 degrees.");
        }
    };
    positive(a.contacts.cutoff, "Contact cutoff");
    positive(a.hbonds.distance_cutoff, "H-bond distance cutoff");
    angle(a.hbonds.angle_cutoff, "H-bond angle cutoff");
    const bool inner = positive(a.mg.inner_cutoff, "Mg inner-sphere cutoff");
    const bool outer = positive(a.mg.outer_cutoff, "Mg outer-sphere cutoff");
    if (inner && outer && a.mg.inner_cutoff > a.mg.outer_cutoff) {
        out.push_back("Mg inner-sphere cutoff must not exceed the outer-sphere cutoff.");
    }
    positive(a.mg.angle_tolerance, "Mg angle tolerance");
    positive(a.stacking.distance_cutoff, "Stacking distance cutoff");
    angle(a.stacking.angle_cutoff, "Stacking angle cutoff");
    positive(a.stacking.offset_cutoff, "Stacking offset cutoff");
    positive(a.stacking.planarity_cutoff, "Stacking planarity cutoff");
    positive(a.energy.cutoff, "Energy cutoff");
    positive(a.energy.dielectric, "Energy dielectric constant");
    positive(a.proximity_search_radius, "Proximity search radius");
    const std::vector<RMSDTarget>& targets = a.rmsd.targets;
    for (std::size_t j = 1; j < targets.size(); ++j) {
        const std::string key = rmsdDataKey(targets[j].label, j);
        for (std::size_t i = 0; i < j; ++i) {
            if (rmsdDataKey(targets[i].label, i) == key) {
                out.push_back("RMSD targets " + std::to_string(i + 1) + " and " + std::to_string(j + 1) +
                              " give the same data key " + repr(key) + "; rename one.");
                break;
            }
        }
    }
    return out;
}

std::vector<std::string> validate(const Project& project) {
    std::vector<std::string> problems;
    const InputParams& in = project.inputs;
    if (in.topology.empty()) {
        problems.push_back("No topology file selected.");
    } else if (!exists(in.topology)) {
        problems.push_back("Topology file not found: " + in.topology);
    }
    if (in.trajectories.empty()) {
        problems.push_back("No trajectory file selected.");
    }
    for (const std::string& t : in.trajectories) {
        if (!exists(t)) {
            problems.push_back("Trajectory file not found: " + t);
        }
    }
    if (in.step < 1) {
        problems.push_back("Frame step must be >= 1.");
    }
    // LD-04 (docs/DIFFERENCES.md): MDATK lets NaN and negative values through to the readers
    if (in.dt_ps && !(*in.dt_ps >= 0.0 && std::isfinite(*in.dt_ps))) {
        problems.push_back("inputs.dt_ps must be positive (0: the trajectory's own time step).");
    }
    const std::vector<std::string> names = groupNames(project);
    std::vector<std::string> sorted = names;
    std::sort(sorted.begin(), sorted.end());
    if (std::adjacent_find(sorted.begin(), sorted.end()) != sorted.end()) {
        problems.push_back("Residue group names must be unique.");
    }
    for (const ResidueGroup& g : project.groups) {
        if (pythonStrip(g.name).empty()) {
            problems.push_back("Residue groups need a name.");
        } else if (!isPlainName(g.name)) {
            problems.push_back("Group name " + repr(g.name) +
                               " may only contain letters, digits, '_', '.' and '-' (no spaces).");
        }
        if (pythonStrip(g.residues).empty() && pythonStrip(g.selection).empty()) {
            problems.push_back("Group " + repr(g.name) + " has neither residues nor a selection.");
        }
    }
    const AnalysisParams& a = project.analyses;
    const std::array<std::tuple<std::string_view, bool, const std::string*, const std::string*>, 5>
        pairs{{{"contacts", a.contacts.enabled, &a.contacts.group_a, &a.contacts.group_b},
               {"H-bonds", a.hbonds.enabled, &a.hbonds.group_a, &a.hbonds.group_b},
               {"Mg", a.mg.enabled, &a.mg.group_a, &a.mg.group_b},
               {"stacking", a.stacking.enabled, &a.stacking.group_a, &a.stacking.group_b},
               {"energy", a.energy.enabled, &a.energy.group_a, &a.energy.group_b}}};
    for (const auto& [label, enabled, groupA, groupB] : pairs) {
        if (!enabled) {
            continue;
        }
        for (const std::string* name : {groupA, groupB}) {
            if (!name->empty() && !contains(names, *name)) {
                problems.push_back(std::string(label) + ": unknown residue group " + repr(*name) +
                                   ".");
            }
        }
    }
    for (std::string& problem : analysisParameterProblems(project)) {
        problems.push_back(std::move(problem));
    }
    const VisualParams& v = project.visual;
    if (!oneOf(v.color_mode, {"chain", "groups", "proximity"})) {
        problems.push_back("Unknown colour mode " + repr(v.color_mode) + ".");
    }
    if (!oneOf(v.glow_mode, {"none", "ever", "fade", "current"})) {
        problems.push_back("Unknown glow mode " + repr(v.glow_mode) + ".");
    }
    for (std::size_t k = 0; k < v.representations.size(); ++k) {
        const RepSpec& rs = v.representations[k];
        const std::string head = "Representation " + std::to_string(k + 1) + ": ";
        if (!oneOf(rs.style, kRepStyles)) {
            problems.push_back(head + "unknown style " + repr(rs.style) + " (" +
                               joined(kRepStyles, ", ") + ").");
        }
        if (!oneOf(rs.color, kRepColors)) {
            problems.push_back(head + "unknown colouring " + repr(rs.color) + ".");
        }
        for (const std::string& name : rs.groups) {
            if (!contains(names, name)) {
                problems.push_back(head + "unknown residue group " + repr(name) + ".");
            }
        }
        if (rs.size <= 0 || !(0 < rs.opacity && rs.opacity <= 1)) {
            problems.push_back(head + "size must be positive and opacity in (0, 1].");
        }
    }
    if (!oneOf(v.renderer, {"auto", "software", "vtk"})) {
        problems.push_back("Unknown renderer " + repr(v.renderer) + " (auto, software or vtk).");
    }
    if (!(0.05 <= v.zoom && v.zoom <= 50)) {
        problems.push_back("Structure zoom must be between 0.05 and 50.");
    }
    const BasePairParams& bp = a.basepairs;
    if (!oneOf(bp.engine, {"builtin", "barnaba"})) {
        problems.push_back("Unknown base-pair engine " + repr(bp.engine) +
                           " (builtin or barnaba).");
    }
    if (bp.stride < 1) {
        problems.push_back("Base-pair stride must be at least 1.");
    }
    if (!(2.0 <= bp.cutoff && bp.cutoff <= 5.0)) {
        problems.push_back("Base-pair H-bond cutoff must be between 2.0 and 5.0 Å.");
    }
    if (!oneOf(v.rep_rest, kRepRest)) {
        problems.push_back("Unknown setting for residues without a representation " +
                           repr(v.rep_rest) + " (trace or hide).");
    }
    if (!oneOf(v.backbone_atom, {"P", "C4'"})) {
        problems.push_back("The backbone atom of ribbons and cartoons must be P or C4'.");
    }
    if (v.proximity_cutoff > a.proximity_search_radius) {
        problems.push_back("Proximity cutoff exceeds the proximity search radius "
                           "(raise analyses.proximity_search_radius).");
    }
    for (std::string& problem : densityProblems(project)) {
        problems.push_back(std::move(problem));
    }
    if (project.output.video.fps <= 0) {
        problems.push_back("Video FPS must be positive.");
    }
    return problems;
}

std::vector<std::string> layoutWarnings(const Project& project) {
    const auto& panels = project.layout.panels;
    const auto count = std::count_if(panels.begin(), panels.end(), [](const PanelSpec& panel) {
        return panel.type == "structure";
    });
    if (count > 1) {
        return {std::to_string(count) + " structure panels in the layout; usually one is intended "
                                        "(check each slot's panel type)"};
    }
    return {};
}

}  // namespace mdbridge::core
