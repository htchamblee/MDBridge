// The parameter model (mdatk/config.py of MDATK v0.0.4.1): plain data in and out, path
// resolution, folders and checks. Expected results come from MDATK itself
// (tests/parity/make_model_cases.py -> data/project-model, python_cases.py -> data/python; see
// data/README.md); the model parts of MDATK's test_config_selections.py, test_projectfile.py and
// test_density.py are ported at the end.

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#if !defined(_WIN32)
#include <pwd.h>
#endif

#include <gtest/gtest.h>

#include "mdbridge/core/ntpath.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/project.hpp"
#include "mdbridge/core/project_file.hpp"
#include "mdbridge/core/python.hpp"
#include "test_support.hpp"

using namespace mdbridge::core;
using namespace mdbridge::test;

namespace {

const fs::path kModel = kData / "project-model";

const Array& items(const Value& value, std::string_view key) {
    return *value.find(key)->get<Array>();
}

std::vector<std::string> strings(const Value& value) {
    std::vector<std::string> out;
    if (const Array* array = value.get<Array>()) {
        for (const Value& item : *array) {
            out.push_back(item.stringOr());
        }
    }
    return out;
}

// The entry at `path` (keys and list indices).
const Value& at(const Value& root, const Array& path) {
    const Value* node = &root;
    for (const Value& key : path) {
        if (const auto* index = key.get<std::int64_t>()) {
            node = &(*node->get<Array>())[static_cast<std::size_t>(*index)];
        } else {
            node = node->find(*key.get<std::string>());
        }
        if (!node) {
            throw std::runtime_error("no such entry");
        }
    }
    return *node;
}

// The input that sets `path` (an index stands for a one-item list of objects).
Value nest(const Array& path, Value value) {
    for (auto it = path.rbegin(); it != path.rend(); ++it) {
        if (it->get<std::int64_t>()) {
            value = Array{std::move(value)};
        } else {
            value = Table{{*it->get<std::string>(), std::move(value)}};
        }
    }
    return value;
}

std::string describe(const Value& c) {
    if (const Value* name = c.find("name")) {
        return name->stringOr();
    }
    return pythonRepr(*c.find("path")) + " = " + pythonRepr(tagged(*c.find("value")));
}

// A path resolved on Linux (MDATK's expectation) as MDBridge's Windows rules write it there: each
// ';'-separated piece normalised in generic form (posixpath keeps an absolute path as given, the
// Windows rules normalise it: "/abs/../x" is "/x", "a\b" is "a/b"). Unchanged on POSIX.
std::string platformResolved(const std::string& path) {
#if defined(_WIN32)
    std::string out;
    std::size_t start = 0;
    for (;;) {
        const std::size_t end = path.find(';', start);
        const std::string piece = path.substr(start, end == std::string::npos ? std::string::npos : end - start);
        if (start) {
            out += ';';
        }
        out += piece.empty() ? piece : ntpath::normpath(piece);
        if (end == std::string::npos) {
            return out;
        }
        start = end + 1;
    }
#else
    return path;
#endif
}

// Runs one from_dict case: the result's repr at result_path, or the error; and the warnings.
void checkFromDict(const Value& c, const Value& input, const Array& resultPath) {
#if defined(_WIN32)
    if (c.find("name")->stringOr() == "base dir, absolute and dotted paths") {
        return;  // posixpath keeps "/x/../t"; the Windows rules normalise it (test_ntpath.cpp)
    }
#endif
    const Value* base = c.find("base_dir");
    const std::string baseDir = base && base->isString() ? base->stringOr() : std::string();
    std::vector<std::string> warnings;
    if (const Value* error = c.find("error")) {
        try {
            projectFromPlain(input, baseDir, &warnings);
            ADD_FAILURE() << describe(c) << ": no error, expected " << error->stringOr();
        } catch (const ProjectError& e) {
            EXPECT_EQ(e.what(), error->stringOr()) << describe(c);
        }
    } else {
        try {
            const Project project = projectFromPlain(input, baseDir, &warnings);
            EXPECT_EQ(pythonRepr(at(toPlain(project), resultPath)), c.find("result")->stringOr())
                << describe(c);
        } catch (const std::exception& e) {
            ADD_FAILURE() << describe(c) << ": " << e.what();
        }
    }
    const Value* expected = c.find("warnings");
    EXPECT_EQ(warnings, expected ? strings(*expected) : std::vector<std::string>()) << describe(c);
}

// The current directory for the life of the object.
class ScopedCwd {
public:
    explicit ScopedCwd(const fs::path& path) : m_old(fs::current_path()) { fs::current_path(path); }
    ~ScopedCwd() { fs::current_path(m_old); }
    ScopedCwd(const ScopedCwd&) = delete;
    ScopedCwd& operator=(const ScopedCwd&) = delete;

private:
    fs::path m_old;
};

#if defined(_WIN32)
constexpr const char* kHomeVariable = "USERPROFILE";
#else
constexpr const char* kHomeVariable = "HOME";
#endif

class ScopedHome {
public:
    explicit ScopedHome(const std::string& value) {
        if (const char* old = std::getenv(kHomeVariable)) {
            m_old = old;
        }
        set(value);
    }
    ~ScopedHome() { set(m_old); }
    ScopedHome(const ScopedHome&) = delete;
    ScopedHome& operator=(const ScopedHome&) = delete;

private:
    static void set(const std::optional<std::string>& value) {
#if defined(_WIN32)
        _putenv_s(kHomeVariable, value ? value->c_str() : "");
#else
        if (value) {
            ::setenv(kHomeVariable, value->c_str(), 1);
        } else {
            ::unsetenv(kHomeVariable);
        }
#endif
    }

    std::optional<std::string> m_old;
};

std::string errorOf(const Value& data, std::string_view baseDir = {}) {
    try {
        std::vector<std::string> warnings;
        projectFromPlain(data, baseDir, &warnings);
    } catch (const ProjectError& e) {
        return e.what();
    }
    return "(no error)";
}

const std::vector<double> kIdentity{1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0};

}  // namespace

// ---- vectors of make_model_cases.py ----------------------------------------------------------

TEST(ProjectModel, DefaultProject) {
    const Value expected = readJson(kModel / "default_project.json");
    EXPECT_EQ(pythonRepr(toPlain(Project{})), expected.find("repr")->stringOr());
    EXPECT_EQ(projectFromPlain(Table{}), Project{});
    EXPECT_EQ(projectFromPlain(toPlain(Project{})), Project{});
}

TEST(ProjectModel, FromDictCases) {
    const Value vectors = readJson(kModel / "from_dict_cases.json");
    for (const Value& c : items(vectors, "cases")) {
        checkFromDict(c, tagged(*c.find("input")), *c.find("result_path")->get<Array>());
    }
}

TEST(ProjectModel, FromDictEveryField) {
    // every field of the model given a value of every TOML type (and None, NaN, infinity ...)
    const Value vectors = readJson(kModel / "from_dict_cases.json");
    const Array& cases = items(vectors, "leaf_cases");
    EXPECT_GT(cases.size(), 5000U);
    for (const Value& c : cases) {
        const Array& path = *c.find("path")->get<Array>();
        checkFromDict(c, nest(path, tagged(*c.find("value"))), path);
    }
}

TEST(ProjectModel, FromDictIntegralFloats) {
    const Value vectors = readJson(kModel / "from_dict_cases.json");
    for (const Value& c : items(vectors, "float_int_cases")) {
        const Array& path = *c.find("path")->get<Array>();
        checkFromDict(c, nest(path, tagged(*c.find("value"))), path);
    }
}

TEST(ProjectModel, FromDictBeyondInt64IsAnError) {
    // deviation: Python keeps the big integer; MDBridge's fields are 64-bit
    const Value vectors = readJson(kModel / "from_dict_cases.json");
    for (const Value& v : items(vectors, "beyond_int64")) {
        const Value value = tagged(v);
        EXPECT_EQ(errorOf(Table{{"inputs", Table{{"start", value}}}}),
                  "project.inputs.start: expected an integer, got " + pythonRepr(value));
        EXPECT_EQ(errorOf(Table{{"inputs", Table{{"stop", value}}}}),
                  "project.inputs.stop: expected an integer, got " + pythonRepr(value));
    }
    EXPECT_EQ(errorOf(Table{{"output", Table{{"video", Table{{"fps", 9223372036854775808.0}}}}}}),
              "project.output.video.fps: expected an integer, got 9.223372036854776e+18");
    const Project p = projectFromPlain(Table{{"inputs", Table{{"start", -9223372036854775808.0}}}});
    EXPECT_EQ(p.inputs.start, std::numeric_limits<std::int64_t>::min());
}

TEST(ProjectModel, ValidateCases) {
    const Value vectors = readJson(kModel / "validate_cases.json");
    TempDir dir;
    for (const std::string& name : strings(*vectors.find("files"))) {
        writeFile(dir / name, "");
    }
    for (const std::string& name : strings(*vectors.find("folders"))) {
        fs::create_directory(dir / name);
    }
    const ScopedCwd cwd(dir.path());
    for (const Value& c : items(vectors, "cases")) {
        const Project p = projectFromPlain(tagged(*c.find("input")));
        const std::string name = c.find("name")->stringOr();
        EXPECT_EQ(validate(p), strings(*c.find("problems"))) << name;
        EXPECT_EQ(densityProblems(p), strings(*c.find("density_problems"))) << name;
        EXPECT_EQ(densityActive(p), *c.find("density_active")->get<bool>()) << name;
        EXPECT_EQ(layoutWarnings(p), strings(*c.find("layout_warnings"))) << name;
        EXPECT_EQ(groupNames(p), strings(*c.find("group_names"))) << name;
    }
}

TEST(ProjectModel, ConfigDocOfEveryPath) {
    const Value vectors = readJson(kModel / "config_doc.json");
    const Array& paths = items(vectors, "paths");
    EXPECT_EQ(paths.size(), 240U);
    for (const Value& c : paths) {
        const std::string path = c.find("path")->stringOr();
        const auto comment = modelComment(path);
        ASSERT_TRUE(comment.has_value()) << path;
        EXPECT_EQ(*comment, c.find("doc")->stringOr()) << path;
    }
    for (const char* other :
         {"", "format", "format_version", "mdatk_version", "saved", "log", "log_file", "project",
          "inputs.nope", "inputs.topology.x", "layout.panels.options.mode", "groups.0.name",
          "groups[0].name", "workflow.nodes.settings.path", "workflow.nodes.files.structure",
          "analyses.contacts.", ".notes", "Notes"}) {
        EXPECT_FALSE(modelComment(other).has_value()) << other;
    }
}

// ---- vectors of python_cases.py --------------------------------------------------------------

TEST(ProjectModel, ResolvePathsCases) {
    const Value vectors = readJson(kData / "python" / "resolve_paths_cases.json");
    const ScopedHome home(vectors.find("home")->stringOr());
    const std::string base = vectors.find("base_dir")->stringOr();
    const Value before = *vectors.find("project_before");
    Project p = projectFromPlain(before);
    EXPECT_EQ(pythonRepr(toPlain(p)), pythonRepr(before));
    resolvePaths(p, base);

    // the changes, as paths like "inputs.trajectories[1]" or "workflow.nodes[0].files.structure"
    Value expected = before;
    for (const Value& change : items(vectors, "changes")) {
        const std::string path = change.find("path")->stringOr();
        Value* node = &expected;
        std::size_t i = 0;
        while (i < path.size()) {
            if (path[i] == '.') {
                ++i;
            } else if (path[i] == '[') {
                const std::size_t close = path.find(']', i);
                node = &(*node->get<Array>())[std::stoul(path.substr(i + 1, close - i - 1))];
                i = close + 1;
            } else {
                const std::size_t end = path.find_first_of(".[", i);
                node = node->find(path.substr(i, end == std::string::npos ? end : end - i));
                ASSERT_NE(node, nullptr) << path;
                i = end == std::string::npos ? path.size() : end;
            }
        }
        EXPECT_EQ(*node, change.find("before")->stringOr()) << path;
        *node = platformResolved(change.find("after")->stringOr());
    }
#if defined(_WIN32)
    // "/abs//x/../c.nc", absolute, is not among the changes (posixpath keeps it as it is); the
    // Windows rules normalise it too
    for (Value& trajectory : *expected.find("inputs")->find("trajectories")->get<Array>()) {
        trajectory = platformResolved(trajectory.stringOr());
    }
#endif
    EXPECT_EQ(pythonRepr(toPlain(p)), pythonRepr(expected));

    for (const Value& c : items(vectors, "cpptraj_path")) {
        Project q;
        q.cleaning.cpptraj_path = c.find("before")->stringOr();
        resolvePaths(q, base);
        EXPECT_EQ(q.cleaning.cpptraj_path, platformResolved(c.find("after")->stringOr()));
    }
#if !defined(_WIN32)
    const passwd* root = ::getpwnam("root");
    const bool rootIsRoot = root && root->pw_dir && std::string(root->pw_dir) == "/root";
#else
    const bool rootIsRoot = false;
#endif
    for (const Value& c : items(vectors, "fix")) {
        const std::string path = c.find("before")->stringOr();
        if (path.rfind("~root", 0) == 0 && !rootIsRoot) {
            continue;
        }
#if defined(_WIN32)
        if (path.rfind("~", 0) == 0 && path.rfind("~/", 0) != 0 && path != "~") {
            continue;  // "~user": no reference behaviour on Windows
        }
#endif
        EXPECT_EQ(resolvePath(path, base), platformResolved(c.find("after")->stringOr())) << path;
    }
    for (const Value& c : items(vectors, "fix_relative_base")) {
        EXPECT_EQ(resolvePath(c.find("before")->stringOr(), c.find("base_dir")->stringOr()),
                  platformResolved(c.find("after")->stringOr()));
    }
}

TEST(ProjectModel, DirsCases) {
    const Value vectors = readJson(kData / "python" / "dirs_cases.json");
    const std::string cwd = currentDirectory();
    const auto expected = [&](const Value& c, std::string_view key) {
        std::string s = c.find(key)->stringOr();
        if (s.rfind("{cwd}", 0) == 0) {
            s.replace(0, 5, cwd);
        }
#if defined(_WIN32)
        s = ntpath::abspath(s, cwd);  // the Windows rules: on the current drive, generic form
#endif
        return s;
    };
    for (const Value& c : items(vectors, "dirs")) {
        Project p;
        p.output.out_dir = c.find("output.out_dir")->stringOr();
        p.output.cache_dir = c.find("output.cache_dir")->stringOr();
        p.workflow.out_dir = c.find("workflow.out_dir")->stringOr();
        p.workflow.project_name = c.find("workflow.project_name")->stringOr();
        EXPECT_EQ(outDir(p), expected(c, "out_dir"));
        EXPECT_EQ(cacheDir(p), expected(c, "cache_dir"));
        EXPECT_EQ(workflowDir(p), expected(c, "workflow_dir"));
    }
    for (const Value& c : items(vectors, "group")) {
        Project p;
        const std::vector<std::string> names = strings(*c.find("groups"));
        if (names != groupNames(p)) {
            p.groups.clear();
            for (const std::string& name : names) {
                p.groups.push_back(ResidueGroup{name});
            }
        }
        const std::string name = c.find("name")->stringOr();
        if (const Value* found = c.find("found")) {
            const ResidueGroup& g = group(p, name);
            EXPECT_EQ(&g, &p.groups[static_cast<std::size_t>(
                              std::find(names.begin(), names.end(), name) - names.begin())]);
            EXPECT_EQ(pythonRepr(Table{{"name", g.name},
                                       {"residues", g.residues},
                                       {"selection", g.selection},
                                       {"color", g.color}}),
                      pythonRepr(*found));
        } else {
            try {
                group(p, name);
                ADD_FAILURE() << name;
            } catch (const std::out_of_range& e) {
                EXPECT_EQ(e.what(), c.find("error")->find("args0")->stringOr());
            }
        }
    }
}

// ---- behaviour without a Python reference ----------------------------------------------------

TEST(ProjectModel, NodeFileThatIsNotAPath) {
    // MDATK crashes (AttributeError); MDBridge names the entry
    const Value node = Table{{"tool", "input-file"}, {"files", Table{{"a", "x"}, {"b", 5}}}};
    EXPECT_EQ(errorOf(Table{{"workflow", Table{{"nodes", Array{Table{}, node}}}}}, "/base"),
              "project.workflow.nodes[1].files.b: expected a file path, got 5");
    EXPECT_EQ(errorOf(Table{{"workflow", Table{{"scratch", Array{node}}}}}, "/base"),
              "project.workflow.scratch[0].files.b: expected a file path, got 5");
    // without a base folder nothing is resolved, and free-form values are kept
    const Project p = projectFromPlain(Table{{"workflow", Table{{"nodes", Array{node}}}}});
    EXPECT_EQ(p.workflow.nodes[0].files, (Table{{"a", "x"}, {"b", 5}}));

    PlanNode bad;
    bad.files = Table{{"ok", "y"}, {"none", Value()}};
    EXPECT_THROW(resolveNodePaths(bad, "/b", "here"), ProjectError);
    EXPECT_EQ(bad.files, (Table{{"ok", "y"}, {"none", Value()}}));  // unchanged
    try {
        resolveNodePaths(bad, "/b", "here");
    } catch (const ProjectError& e) {
        EXPECT_STREQ(e.what(), "here.files.none: expected a file path, got None");
    }
}

TEST(ProjectModel, ResolveNodePaths) {
    PlanNode node;
    node.tool = "input-file";
    node.files = Table{{"s", " a.pdb ;　; /c.pdb;"}, {"e", ""}};
    node.settings = Table{{"kind", "pdb"}, {"path", "m/x.pdb"}, {"other", "o.pdb"}};
    resolveNodePaths(node, "/b", "n");
    EXPECT_EQ(node.files, (Table{{"s", "/b/a.pdb;/c.pdb"}, {"e", ""}}));
    EXPECT_EQ(node.settings, (Table{{"kind", "pdb"}, {"path", "/b/m/x.pdb"}, {"other", "o.pdb"}}));
    node.tool = "note";
    node.settings = Table{{"path", "x"}};
    resolveNodePaths(node, "/b", "n");
    EXPECT_EQ(node.settings, (Table{{"path", "x"}}));
}

TEST(ProjectModel, UnknownKeysAreLoggedWithoutAList) {
    const LogCapture log;
    const Project p = projectFromPlain(Table{{"zzz", 1}, {"inputs", Table{{"q", 2}}}});
    EXPECT_EQ(p, Project{});
    EXPECT_EQ(log.lines(), (std::vector<std::string>{
                               "WARNING: project: ignoring unknown key 'zzz'",
                               "WARNING: project.inputs: ignoring unknown key 'q'"}));
}

TEST(ProjectModel, ToPlainDropsNullsOfFreeFormTables) {
    Project p;
    const Table options{{"a", Value()},
                        {"b", Table{{"c", Value()}, {"d", 1}}},
                        {"l", Array{Value(), Table{{"e", Value()}}}}};
    p.layout.panels = {PanelSpec{"structure", {}, "", options}};
    p.workflow.nodes = {PlanNode{}};
    p.workflow.nodes[0].settings = Table{{"x", Value()}, {"y", 2}};
    p.workflow.nodes[0].files = Table{{"f", Value()}};
    p.inputs.stop = 3;
    const Value plain = toPlain(p);
    const Table kept{{"b", Table{{"d", 1}}}, {"l", Array{Value(), Table{}}}};
    EXPECT_EQ(*plain.path({"layout", "panels"}),
              (Array{Table{
                  {"type", "structure"}, {"series", Array{}}, {"title", ""}, {"options", kept}}}));
    const Value& node = (*plain.path({"workflow", "nodes"})->get<Array>())[0];
    EXPECT_EQ(*node.find("settings"), (Table{{"y", 2}}));
    EXPECT_EQ(*node.find("files"), Table{});
    EXPECT_EQ(*plain.path({"inputs", "stop"}), Value(std::int64_t{3}));
    EXPECT_TRUE(plain.path({"inputs", "dt_ps"})->isNull());
    EXPECT_TRUE(plain.path({"analyses", "energy", "target_value"})->isNull());
}

TEST(ProjectModel, NonTableInput) {
    EXPECT_EQ(errorOf(Array{}), "project: expected an object, got list");
    EXPECT_EQ(errorOf(Value()), "project: expected an object, got NoneType");
}

TEST(ProjectModel, RigidTransformProblem) {
    EXPECT_EQ(rigidTransformProblem(kIdentity), std::nullopt);
    EXPECT_EQ(rigidTransformProblem(std::vector<double>{1, 2, 3}),
              "a rigid transform needs 12 numbers (3x3 rotation + translation)");
    EXPECT_EQ(rigidTransformProblem(std::vector<double>{}),
              "a rigid transform needs 12 numbers (3x3 rotation + translation)");
    std::vector<double> t = kIdentity;
    t[11] = std::numeric_limits<double>::quiet_NaN();
    EXPECT_EQ(rigidTransformProblem(t), "the transform contains non-finite numbers");
    t = kIdentity;
    t[0] = 1e300;
    EXPECT_EQ(rigidTransformProblem(t), "the rotation part of the transform is not a rotation");
    t = kIdentity;
    t[4] = -1;
    EXPECT_EQ(rigidTransformProblem(t), "the rotation part of the transform is a reflection");
    // rotations about every axis, rounded to 6 digits as stored, are accepted
    for (int k = 0; k < 64; ++k) {
        const double a = 0.3 * k, b = -0.17 * k, c = 0.05 * k * k;
        const double ca = std::cos(a), sa = std::sin(a), cb = std::cos(b), sb = std::sin(b),
                     cc = std::cos(c), sc = std::sin(c);
        // Rz(c) Ry(b) Rx(a)
        std::vector<double> r{cc * cb, cc * sb * sa - sc * ca, cc * sb * ca + sc * sa,
                              sc * cb, sc * sb * sa + cc * ca, sc * sb * ca - cc * sa,
                              -sb,     cb * sa,                cb * ca,
                              1.5,     -2,                     1e6};
        for (double& x : r) {
            x = std::round(x * 1e6) / 1e6;
        }
        EXPECT_EQ(rigidTransformProblem(r), std::nullopt) << k;
        for (int i = 0; i < 9; ++i) {
            r[static_cast<std::size_t>(i)] *= i < 3 ? -1 : 1;  // a reflected first row
        }
        EXPECT_EQ(rigidTransformProblem(r), "the rotation part of the transform is a reflection")
            << k;
    }
}

// ---- MDATK's tests (model parts) ---------------------------------------------------------------

// test_config_selections.py::test_project_roundtrip, without the file
TEST(ProjectModel, ProjectRoundTrip) {
    Project p;
    p.inputs.topology = "/tmp/x/sys.parm7";
    p.inputs.trajectories = {"/tmp/x/a.nc", "/tmp/x/b.nc"};
    p.output.out_dir = "/tmp/x/out";
    p.analyses.contacts.cutoff = 5.25;
    p.layout.panels.push_back(
        PanelSpec{"contact_map", {"contacts/map"}, "Map", Table{{"mode", "frequency"}}});
    const Project q = projectFromPlain(toPlain(p), "/elsewhere");
    EXPECT_EQ(toPlain(q), toPlain(p));
    EXPECT_EQ(q, p);
}

// test_config_selections.py::test_relative_paths_and_aliases, from the data of that file
TEST(ProjectModel, RelativePathsAndAliases) {
    const Value data = parseToml("not_a_key = 1\n[input]\ntopology = \"top.parm7\"\n"
                                 "trajectory = \"traj.nc\"\n[output]\nout_dir = \"out\"\n",
                                 "p.mdatk.toml");
    const LogCapture log;
    const Project p = projectFromPlain(data, "/work/proj");
    EXPECT_EQ(p.inputs.topology, "/work/proj/top.parm7");
    EXPECT_EQ(p.inputs.trajectories, std::vector<std::string>{"/work/proj/traj.nc"});
    EXPECT_EQ(p.output.out_dir, "/work/proj/out");
    ASSERT_EQ(log.lines().size(), 1U);
    EXPECT_NE(log.lines()[0].find("not_a_key"), std::string::npos);
}

// test_config_selections.py::test_type_errors_are_reported
TEST(ProjectModel, TypeErrorsAreReported) {
    EXPECT_NE(errorOf(Table{{"analyses", Table{{"contacts", Table{{"cutoff", "far"}}}}}})
                  .find("analyses.contacts.cutoff"),
              std::string::npos);
    EXPECT_NE(errorOf(Table{{"analyses", Table{{"rmsd", Table{{"enabled", 1}}}}}}).find("enabled"),
              std::string::npos);
}

// test_config_selections.py::test_params_log_can_be_replayed: the "parameters" table loads
TEST(ProjectModel, ParamsLogParametersLoad) {
    Project p;
    p.analyses.hbonds.angle_cutoff = 150.0;
    EXPECT_EQ(projectFromPlain(toPlain(p)).analyses.hbonds.angle_cutoff, 150.0);
}

// test_config_selections.py::test_validate_flags_problems
TEST(ProjectModel, ValidateFlagsProblems) {
    Project p;
    const auto any = [](const std::vector<std::string>& problems, std::string_view word) {
        return std::any_of(problems.begin(), problems.end(),
                           [&](const std::string& m) { return m.find(word) != std::string::npos; });
    };
    EXPECT_TRUE(any(validate(p), "topology"));
    EXPECT_EQ(validate(p), (std::vector<std::string>{"No topology file selected.",
                                                     "No trajectory file selected."}));
    p.groups.push_back(ResidueGroup{"Bad name", "1-3"});
    p.visual.glow_mode = "sparkle";
    const std::vector<std::string> problems = validate(p);
    EXPECT_TRUE(any(problems, "spaces"));
    EXPECT_TRUE(any(problems, "glow"));
}

// test_projectfile.py::test_optional_values_and_round_trip, without the file
TEST(ProjectModel, OptionalValuesRoundTrip) {
    Project p;
    p.inputs.stop = std::nullopt;
    p.inputs.dt_ps = 2.5;
    p.analyses.energy.target_value = std::nullopt;
    p.layout.panels.push_back(
        PanelSpec{"contact_map", {"contacts/map"}, "Map", Table{{"transpose", true}}});
    const Project q = projectFromPlain(toPlain(p));
    EXPECT_FALSE(q.inputs.stop.has_value());
    EXPECT_EQ(q.inputs.dt_ps, 2.5);
    EXPECT_EQ(toPlain(projectFromPlain(toPlain(q))), toPlain(q));
}

// test_density.py::test_config_density_validation
TEST(ProjectModel, ConfigDensityValidation) {
    TempDir dir;
    Project p;
    EXPECT_TRUE(densityProblems(p).empty());
    EXPECT_FALSE(densityActive(p));
    p.density.map_path = pathToUtf8(dir / "missing.mrc");
    p.density.fit_frame = "first";
    p.density.transform = {1.0, 2.0};
    std::string probs;
    for (const std::string& m : densityProblems(p)) {
        probs += m + " | ";
    }
    EXPECT_NE(probs.find("not found"), std::string::npos);
    EXPECT_NE(probs.find("'ref' or a frame number"), std::string::npos);
    EXPECT_NE(probs.find("12 numbers"), std::string::npos);
    EXPECT_EQ(projectFromPlain(toPlain(p)).density, p.density);
}

// test_density.py: the angular step check and the placement checks of the pipeline tests
TEST(ProjectModel, DensityAngularStepAndPlacement) {
    TempDir dir;
    writeFile(dir / "map.mrc", "");
    Project p;
    p.density.map_path = pathToUtf8(dir / "map.mrc");
    p.density.angular_step = 90.0;
    const auto any = [](const std::vector<std::string>& problems, std::string_view word) {
        return std::any_of(problems.begin(), problems.end(),
                           [&](const std::string& m) { return m.find(word) != std::string::npos; });
    };
    EXPECT_TRUE(any(densityProblems(p), "angular step"));
    p.density.angular_step = 30.0;
    EXPECT_TRUE(densityProblems(p).empty());
    p.density.show = false;
    p.density.map_path = pathToUtf8(dir / "gone.mrc");
    EXPECT_TRUE(densityProblems(p).empty());  // a hidden missing map does not block runs
    p.density.show = true;
    p.density.map_path = pathToUtf8(dir / "map.mrc");
    p.density.transform = kIdentity;
    p.density.transform[9] = std::numeric_limits<double>::infinity();
    EXPECT_TRUE(any(densityProblems(p), "invalid"));
}
