// Port of MDATK's mdatk/render/coloring.py and residues.py (the owner's code).

#include "mdbridge/render/shading.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <mutex>
#include <set>
#include <stdexcept>

#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/log.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/core/selections.hpp"
#include "mdbridge/core/text.hpp"
#include "mdbridge/io/pdb.hpp"
#include "mdbridge/render/colormap.hpp"

namespace mdbridge::render {

namespace {

using core::AnalysisResult;
using core::NumArray;

// The values of an integer or bool array as int64.
std::vector<std::int64_t> ints(const NumArray& a) {
    switch (a.dtype()) {
        case NumArray::DType::Bool: {
            auto v = a.values<std::uint8_t>();
            return {v.begin(), v.end()};
        }
        case NumArray::DType::Int8: {
            auto v = a.values<std::int8_t>();
            return {v.begin(), v.end()};
        }
        case NumArray::DType::Int32: {
            auto v = a.values<std::int32_t>();
            return {v.begin(), v.end()};
        }
        case NumArray::DType::Int64: {
            auto v = a.values<std::int64_t>();
            return {v.begin(), v.end()};
        }
        default:
            break;
    }
    std::vector<std::int64_t> out;
    for (double d : a.toDouble()) {
        out.push_back(static_cast<std::int64_t>(d));
    }
    return out;
}

std::vector<bool> bools(const NumArray& a) {
    std::vector<bool> out;
    if (a.dtype() == NumArray::DType::Float32 || a.dtype() == NumArray::DType::Float64) {
        for (double d : a.toDouble()) {
            out.push_back(d != 0.0);
        }
        return out;
    }
    for (std::int64_t v : ints(a)) {
        out.push_back(v != 0);
    }
    return out;
}

const NumArray& required(const AnalysisResult& r, std::string_view name) {
    const NumArray* a = r.array(name);
    if (a == nullptr) {
        throw std::out_of_range(std::string(name));
    }
    return *a;
}

std::string repr(std::string_view s) { return core::pythonRepr(core::Value(std::string(s))); }

}  // namespace

std::vector<float> FrameMatrix::column(std::size_t frame) const {
    std::vector<float> out(rows);
    for (std::size_t r = 0; r < rows; ++r) {
        out[r] = at(r, frame);
    }
    return out;
}

FrameMatrix glowIntensity(std::span<const double> minDistance, std::size_t rows, std::size_t frames,
                          const std::vector<bool>& member, double cutoff, std::string_view mode,
                          std::int64_t fadeFrames) {
    if (minDistance.size() != rows * frames || member.size() != rows) {
        throw std::invalid_argument("glowIntensity: min_distance must be rows x frames and member one per row");
    }
    FrameMatrix out{rows, frames, std::vector<float>(rows * frames, 0.0f)};
    auto near = [&](std::size_t r, std::size_t f) { return minDistance[r * frames + f] <= cutoff && !member[r]; };
    if (mode == "none") {
        return out;
    }
    if (mode == "current") {
        for (std::size_t r = 0; r < rows; ++r) {
            for (std::size_t f = 0; f < frames; ++f) {
                out.values[r * frames + f] = near(r, f) ? 1.0f : 0.0f;
            }
        }
        return out;
    }
    if (mode == "ever") {
        for (std::size_t r = 0; r < rows; ++r) {
            bool any = false;
            for (std::size_t f = 0; f < frames && !any; ++f) {
                any = near(r, f);
            }
            std::fill_n(out.values.begin() + static_cast<std::ptrdiff_t>(r * frames), frames, any ? 1.0f : 0.0f);
        }
        return out;
    }
    if (mode == "fade") {
        for (std::size_t r = 0; r < rows; ++r) {
            double since = std::numeric_limits<double>::infinity();
            for (std::size_t f = 0; f < frames; ++f) {
                since = near(r, f) ? 0.0 : since + 1.0;
                double v = 0.0;
                if (fadeFrames > 0) {
                    v = std::clamp(1.0 - since / (static_cast<double>(fadeFrames) + 1.0), 0.0, 1.0);
                } else {
                    v = near(r, f) ? 1.0 : 0.0;
                }
                out.values[r * frames + f] = static_cast<float>(v);
            }
        }
        return out;
    }
    throw core::ValueError("unknown glow mode " + repr(mode) + " (choose from none, ever, fade, current)");
}

std::vector<Rgba> staticColors(std::string_view mode, std::span<const std::int64_t> chain,
                               std::span<const std::int64_t> groupIndex, const std::vector<std::string>& groupColors,
                               std::string_view baseColor, std::string_view chainPalette) {
    const std::size_t n = chain.size();
    std::vector<Rgba> out(n, toRgba(baseColor));
    if (mode == "chain") {
        const std::set<std::int64_t> chains(chain.begin(), chain.end());  // np.unique: sorted
        const std::vector<Rgba> cols = palette(chainPalette, static_cast<int>(chains.size()));
        std::size_t k = 0;
        for (std::int64_t c : chains) {
            for (std::size_t i = 0; i < n; ++i) {
                if (chain[i] == c) {
                    out[i] = cols[k];
                }
            }
            ++k;
        }
    } else if (mode == "groups" || mode == "proximity") {
        if (groupIndex.size() != n) {
            throw std::invalid_argument("staticColors: group_index must have one entry per residue");
        }
        for (std::size_t gi = 0; gi < groupColors.size(); ++gi) {
            const Rgba c = toRgba(groupColors[gi]);
            for (std::size_t i = 0; i < n; ++i) {
                if (groupIndex[i] == static_cast<std::int64_t>(gi)) {
                    out[i] = c;
                }
            }
        }
    }
    return out;
}

std::vector<Rgba> frameColors(const std::vector<Rgba>& colors, std::string_view mode, std::string_view baseColor,
                              std::string_view proximityColor, std::span<const float> intensity,
                              const std::vector<bool>& member) {
    if (mode != "proximity") {
        return colors;
    }
    const Rgba base = toRgba(baseColor);
    const Rgba hi = toRgba(proximityColor);
    if (intensity.size() != colors.size() || member.size() != colors.size()) {
        throw std::invalid_argument("frameColors: one intensity and member flag per residue");
    }
    std::vector<Rgba> out(colors.size());
    for (std::size_t k = 0; k < out.size(); ++k) {
        if (member[k]) {
            out[k] = colors[k];
            continue;
        }
        const double i = static_cast<double>(intensity[k]);
        out[k] = {base.r * (1.0 - i) + hi.r * i, base.g * (1.0 - i) + hi.g * i, base.b * (1.0 - i) + hi.b * i,
                  base.a * (1.0 - i) + hi.a * i};
    }
    return out;
}

namespace {
std::mutex& warnedMutex() {
    static std::mutex m;
    return m;
}
std::set<std::string>& warned() {
    static std::set<std::string> s;
    return s;
}
}  // namespace

void warnOnce(const std::string& message) {
    {
        std::lock_guard lock(warnedMutex());
        if (!warned().insert(message).second) {
            return;
        }
    }
    core::logWarning(message);
}

void resetWarnOnce() {
    std::lock_guard lock(warnedMutex());
    warned().clear();
}

std::string oneLetter(std::string_view resname) {
    if (const auto info = io::parseResidueName(resname)) {
        return std::string(1, info->base);
    }
    static const std::map<std::string, char, std::less<>> amino{
        {"ALA", 'A'}, {"ARG", 'R'}, {"ASN", 'N'}, {"ASP", 'D'}, {"CYS", 'C'}, {"GLN", 'Q'}, {"GLU", 'E'},
        {"GLY", 'G'}, {"HIS", 'H'}, {"ILE", 'I'}, {"LEU", 'L'}, {"LYS", 'K'}, {"MET", 'M'}, {"PHE", 'F'},
        {"PRO", 'P'}, {"SER", 'S'}, {"THR", 'T'}, {"TRP", 'W'}, {"TYR", 'Y'}, {"VAL", 'V'}, {"HID", 'H'},
        {"HIE", 'H'}, {"HIP", 'H'}, {"CYX", 'C'}};
    // strip().upper(): only ASCII names can match the table (a non-ASCII letter whose upper case is
    // ASCII, U+0131 dotless i -> I and U+017F long s -> S, is mapped too, as str.upper maps it)
    std::string key;
    const std::string stripped = core::pythonStrip(resname);
    for (std::size_t k = 0; k < stripped.size(); ++k) {
        const auto c = static_cast<unsigned char>(stripped[k]);
        if (c == 0xC4 && k + 1 < stripped.size() && static_cast<unsigned char>(stripped[k + 1]) == 0xB1) {
            key += 'I';
            ++k;
        } else if (c == 0xC5 && k + 1 < stripped.size() && static_cast<unsigned char>(stripped[k + 1]) == 0xBF) {
            key += 'S';
            ++k;
        } else if (c >= 'a' && c <= 'z') {
            key += static_cast<char>(c - 'a' + 'A');
        } else {
            key += static_cast<char>(c);
        }
    }
    const auto it = amino.find(key);
    return std::string(1, it != amino.end() ? it->second : 'X');
}

std::vector<std::string> chainNames(const AnalysisResult& trace) {
    if (const NumArray* names = trace.array("chain_name")) {
        if (names->dtype() == NumArray::DType::String) {
            auto v = names->values<std::string>();
            if (std::any_of(v.begin(), v.end(), [](const std::string& s) { return !s.empty(); })) {
                return {v.begin(), v.end()};
            }
        }
    }
    std::vector<std::string> out;
    for (std::int64_t c : ints(required(trace, "chain"))) {
        const std::int64_t m = ((c % 26) + 26) % 26;  // Python's modulo
        out.emplace_back(1, static_cast<char>('A' + m));
    }
    return out;
}

std::vector<std::string> residueNames(const AnalysisResult& trace) {
    if (const NumArray* names = trace.array("resname")) {
        if (names->dtype() == NumArray::DType::String) {
            auto v = names->values<std::string>();
            return {v.begin(), v.end()};
        }
        std::vector<std::string> out;
        for (std::int64_t i : ints(*names)) {
            out.push_back(std::to_string(i));
        }
        return out;
    }
    const std::size_t n = required(trace, "chain").size();
    std::vector<std::string> labels;
    for (const auto& [k, v] : trace.meta) {
        if (k == "residue_labels") {
            if (const core::Array* a = v.get<core::Array>()) {
                for (const core::Value& e : *a) {
                    labels.push_back(e.get<std::string>() != nullptr ? *e.get<std::string>() : core::pythonStr(e));
                }
            }
        }
    }
    if (!labels.empty() && labels.size() == n) {  // "G12" / "A:G12"
        std::vector<std::string> out;
        for (const std::string& lab : labels) {
            const std::size_t colon = lab.rfind(':');
            std::string last = colon == std::string::npos ? lab : lab.substr(colon + 1);
            while (!last.empty() && ((last.back() >= '0' && last.back() <= '9') || last.back() == '-')) {
                last.pop_back();
            }
            out.push_back(last.empty() ? "X" : last);
        }
        return out;
    }
    return std::vector<std::string>(n, "X");
}

std::vector<std::int64_t> residueRows(const AnalysisResult& trace, std::string_view spec) {
    const std::size_t n = required(trace, "chain").size();
    std::vector<std::int64_t> out;
    if (core::pythonStrip(spec).empty()) {
        for (std::size_t k = 0; k < n; ++k) {
            out.push_back(static_cast<std::int64_t>(k));
        }
        return out;
    }
    const std::vector<std::int64_t> resid = ints(required(trace, "resid"));
    const std::vector<std::string> chains = chainNames(trace);
    std::vector<bool> keep(n, false);
    for (const core::ResidueRange& r : core::parseResidueRanges(spec)) {
        for (std::size_t k = 0; k < n; ++k) {
            if (resid[k] >= r.first && resid[k] <= r.last && (!r.chain || chains[k] == *r.chain)) {
                keep[k] = true;
            }
        }
    }
    for (std::size_t k = 0; k < n; ++k) {
        if (keep[k]) {
            out.push_back(static_cast<std::int64_t>(k));
        }
    }
    if (out.empty()) {
        throw core::SelectionError("no residue of the structure view matches " + repr(spec));
    }
    return out;
}

Shading residueShading(const core::Project& project, const core::Dataset& dataset, const AnalysisResult& trace,
                       std::optional<bool> needProximity) {
    const core::VisualParams& v = project.visual;
    const std::vector<std::int64_t> chain = ints(required(trace, "chain"));
    const std::size_t nRes = chain.size();
    const std::size_t nFrames = dataset.frameCount() != 0 ? dataset.frameCount() : 1;
    std::vector<std::string> groupColors;
    for (const core::ResidueGroup& g : project.groups) {
        groupColors.push_back(g.color);
    }
    Shading out;
    out.staticColors = staticColors(v.color_mode, chain, ints(required(trace, "group_index")), groupColors,
                                    v.base_color, v.chain_palette);
    out.member.assign(nRes, false);
    out.glow = {nRes, nFrames, std::vector<float>(nRes * nFrames, 0.0f)};
    out.colorIntensity = out.glow;
    const bool need = needProximity.value_or(v.color_mode == "proximity" || v.glow_mode != "none");
    if (need) {
        const AnalysisResult* prox = dataset.result("proximity");
        if (prox == nullptr) {
            throw std::out_of_range("proximity");
        }
        const NumArray& d = required(*prox, "min_distance");
        std::vector<double> dist = d.toDouble();
        std::size_t rows = d.shape().empty() ? 0 : d.shape()[0];
        const std::size_t frames = d.shape().size() > 1 ? d.shape()[1] : 0;
        std::vector<bool> member = bools(required(*prox, "member"));
        const std::vector<std::int64_t> proxIndex = ints(required(*prox, "resindex"));
        const std::vector<std::int64_t> traceIndex = ints(required(trace, "resindex"));
        if (proxIndex != traceIndex) {
            std::map<std::int64_t, std::size_t> order;
            for (std::size_t k = 0; k < proxIndex.size(); ++k) {
                order[proxIndex[k]] = k;  // the last of repeated indices, as a dict comprehension
            }
            std::vector<double> remapped;
            std::vector<bool> remember;
            for (std::int64_t ix : traceIndex) {
                const auto it = order.find(ix);
                if (it == order.end()) {
                    remapped.insert(remapped.end(), frames, std::numeric_limits<double>::infinity());
                    remember.push_back(false);
                } else {
                    remapped.insert(remapped.end(), dist.begin() + static_cast<std::ptrdiff_t>(it->second * frames),
                                    dist.begin() + static_cast<std::ptrdiff_t>((it->second + 1) * frames));
                    remember.push_back(member[it->second]);
                }
            }
            dist = std::move(remapped);
            member = std::move(remember);
            rows = traceIndex.size();
        }
        out.member = member;
        out.glow = glowIntensity(dist, rows, frames, member, v.proximity_cutoff, v.glow_mode, v.fade_frames);
        const std::string mode = v.glow_mode != "none" ? v.glow_mode : "current";
        out.colorIntensity = glowIntensity(dist, rows, frames, member, v.proximity_cutoff, mode, v.fade_frames);
    }
    const Rgba g = toRgba(v.glow_color);
    out.glowRgb = {g.r, g.g, g.b};
    return out;
}

std::vector<std::int64_t> tagRows(const AnalysisResult& trace, const core::VisualParams& visual) {
    const std::string spec = core::pythonStrip(visual.label_residues);
    const NumArray* residArray = trace.array("resid");
    if (residArray == nullptr || !(visual.label_every > 0 || !spec.empty())) {
        return {};
    }
    const std::vector<std::int64_t> resid = ints(*residArray);
    std::vector<bool> wanted(resid.size(), false);
    if (visual.label_every > 0) {
        for (std::size_t k = 0; k < resid.size(); ++k) {
            wanted[k] = wanted[k] || resid[k] % visual.label_every == 0;
        }
    }
    if (!spec.empty()) {
        try {
            const std::vector<core::ResidueRange> ranges = core::parseResidueRanges(visual.label_residues);
            std::vector<std::string> chains;
            for (const core::ResidueRange& r : ranges) {
                if (r.chain && chains.empty() && trace.array("chain") != nullptr) {
                    chains = chainNames(trace);
                }
                for (std::size_t k = 0; k < resid.size(); ++k) {
                    // RN-24: a chain-qualified range tags that chain only (MDATK ignored the chain)
                    const bool chainOk = !r.chain || chains.empty() || chains[k] == *r.chain;
                    if (resid[k] >= r.first && resid[k] <= r.last && chainOk) {
                        wanted[k] = true;
                    }
                }
            }
        } catch (const core::ValueError& e) {
            // RN-23: MDATK ignored a spec it could not read without a word
            warnOnce("visual.label_residues: " + std::string(e.what()) + "; no residue-number tags from it");
        }
    }
    std::vector<std::int64_t> out;
    for (std::size_t k = 0; k < wanted.size(); ++k) {
        if (wanted[k]) {
            out.push_back(static_cast<std::int64_t>(k));
        }
    }
    return out;
}

}  // namespace mdbridge::render
