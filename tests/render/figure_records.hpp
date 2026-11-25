#pragma once
// Helpers shared by the panel tests (wave 2): the figure cases of tests/data/render/figures
// (tests/parity/make_figure_cases.py) and their inputs, read back as MDBridge's objects:
//
//   untag(v)            the generators' tagged JSON ({"$float": "nan"}) as plain values
//   loadRenderInput     tests/data/render/inputs/<name>: dataset.json + <result>.npz -> core::Dataset
//                       (loadRenderInputFolder: any such folder)
//   caseProject         a case record's "project" (Project.to_dict(), paths {data}) -> core::Project
//   buildCase           build_figure of a case: its input, project, size, dpi, watermark
//
// and the comparison of a drawn MDBridge figure (figure_dump.hpp) with MDATK's record
// (figure_compare.hpp).

#include <cmath>
#include <filesystem>
#include <string>
#include <vector>

#include "mdbridge/analysis/cache.hpp"
#include "mdbridge/core/dataset.hpp"
#include "mdbridge/core/project.hpp"
#include "mdbridge/core/value.hpp"
#include "mdbridge/render/panels/build.hpp"
#include "mdbridge/render/text.hpp"
#include "test_support.hpp"

namespace mdbridge::test {

inline const std::filesystem::path kRenderData = kData / "render";
inline const std::filesystem::path kFigureCases = kRenderData / "figures";

inline const render::FontSet& figureFonts() {
    static const render::FontSet set = render::FontSet::load(
        [](std::string_view f) { return readFile(std::filesystem::path(MDBRIDGE_FONT_DIR) / std::string(f)); });
    return set;
}

// {"$float": "nan"} -> NaN, "{data}" in strings -> the test data folder
inline core::Value untag(const core::Value& v) {
    if (const core::Table* t = v.get<core::Table>()) {
        if (t->size() == 1 && (*t)[0].first == "$float") {
            const std::string s = (*t)[0].second.stringOr();
            if (s == "nan") return std::nan("");
            if (s == "inf") return INFINITY;
            if (s == "-inf") return -INFINITY;
            return std::stod(s);
        }
        core::Table out;
        for (const auto& [k, x] : *t) out.emplace_back(k, untag(x));
        return out;
    }
    if (const core::Array* a = v.get<core::Array>()) {
        core::Array out;
        for (const core::Value& x : *a) out.push_back(untag(x));
        return out;
    }
    if (const std::string* s = v.get<std::string>()) {
        std::string out = *s;
        const std::string data = kData.generic_string();
        for (std::size_t at = out.find("{data}"); at != std::string::npos; at = out.find("{data}", at + data.size())) {
            out.replace(at, 6, data);
        }
        return out;
    }
    return v;
}

// a dataset folder written by render_common.save_dataset (dataset.json + <result>.npz)
inline core::Dataset loadRenderInputFolder(const std::filesystem::path& folder) {
    const core::Value doc = untag(readJson(folder / "dataset.json"));
    std::vector<std::int64_t> frames;
    for (const core::Value& f : *doc.find("frames")->get<core::Array>()) frames.push_back(*f.get<std::int64_t>());
    std::vector<double> time;
    for (const core::Value& t : *doc.find("time_ps")->get<core::Array>()) {
        time.push_back(t.get<double>() ? *t.get<double>() : static_cast<double>(*t.get<std::int64_t>()));
    }
    core::Dataset ds(std::move(frames), std::move(time), {}, doc.find("time_unit")->stringOr());
    for (const core::Value& r : *doc.find("results")->get<core::Array>()) {
        core::AnalysisResult res;
        res.name = r.find("name")->stringOr();
        res.key = r.find("key")->stringOr();
        res.params = *r.find("params")->get<core::Table>();
        res.meta = *r.find("meta")->get<core::Table>();
        core::NamedArrays arrays = analysis::readNpz(folder / (res.name + ".npz"));
        for (const core::Value& n : *r.find("arrays")->get<core::Array>()) {
            for (auto& [k, a] : arrays) {
                if (k == n.stringOr()) res.arrays.emplace_back(k, a);
            }
        }
        ds.add(std::move(res));
    }
    return ds;
}

inline core::Dataset loadRenderInput(const std::string& name) { return loadRenderInputFolder(kRenderData / "inputs" / name); }

inline core::Project caseProject(const core::Value& record) {
    return core::projectFromPlain(untag(*record.find("project")));
}

inline std::vector<double> numbersOf(const core::Value& v) {
    std::vector<double> out;
    for (const core::Value& x : *v.get<core::Array>()) {
        if (const double* d = x.get<double>()) out.push_back(*d);
        else if (const std::int64_t* i = x.get<std::int64_t>()) out.push_back(static_cast<double>(*i));
        else out.push_back(std::nan(""));
    }
    return out;
}

// build_figure(dataset, project, w, h, dpi, watermark) of a case record
inline render::panels::BuiltFigure buildCase(const core::Value& record, const core::Dataset& dataset,
                                             const core::Project& project) {
    const std::vector<double> size = numbersOf(*record.find("size"));
    const core::Value* dpi = record.find("dpi");
    const double d = dpi->get<double>() ? *dpi->get<double>() : static_cast<double>(*dpi->get<std::int64_t>());
    return render::panels::buildFigure(dataset, project, figureFonts(), static_cast<int>(size[0]),
                                       static_cast<int>(size[1]), d, record.find("watermark")->stringOr());
}

}  // namespace mdbridge::test
