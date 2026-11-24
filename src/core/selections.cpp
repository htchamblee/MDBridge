#include "mdbridge/core/selections.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>

#include "mdbridge/core/atom_selection.hpp"
#include "mdbridge/core/log.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/core/text.hpp"

namespace mdbridge::core {

namespace {

std::string repr(std::string_view text) { return pythonRepr(Value(std::string(text))); }

bool isAsciiLetter(char32_t c) { return (c >= U'a' && c <= U'z') || (c >= U'A' && c <= U'Z'); }
bool isAsciiDigit(char32_t c) { return c >= U'0' && c <= U'9'; }

char asciiUpper(char c) { return c >= 'a' && c <= 'z' ? static_cast<char>(c - 'a' + 'A') : c; }
char asciiLower(char c) { return c >= 'A' && c <= 'Z' ? static_cast<char>(c - 'A' + 'a') : c; }

std::string upper(std::string text) {
    std::transform(text.begin(), text.end(), text.begin(), asciiUpper);
    return text;
}

std::string lower(std::string text) {
    std::transform(text.begin(), text.end(), text.begin(), asciiLower);
    return text;
}

// re.split(r"[,;\s]+", text) without the empty pieces
std::vector<std::string> splitItems(std::string_view text) {
    std::vector<std::string> out;
    std::u32string current;
    for (const char32_t c : decodeUtf8(text)) {
        if (c == U',' || c == U';' || isPythonSpace(c)) {
            if (!current.empty()) {
                out.push_back(encodeUtf8(current));
                current.clear();
            }
        } else {
            current += c;
        }
    }
    if (!current.empty()) {
        out.push_back(encodeUtf8(current));
    }
    return out;
}

// str(int(digits)) of a run of \d characters
std::string decimalText(std::u32string_view digits) {
    std::string out;
    for (const char32_t c : digits) {
        // Unicode decimal digits come in runs of ten starting at a zero
        char32_t zero = c;
        while (zero > 0 && isPythonDigit(zero - 1)) {
            --zero;
        }
        out += static_cast<char>('0' + (c - zero) % 10);
    }
    const std::size_t nonZero = out.find_first_not_of('0');
    return nonZero == std::string::npos ? "0" : out.substr(nonZero);
}

// a < b for non-negative decimal texts without leading zeros
bool decimalLess(const std::string& a, const std::string& b) {
    return a.size() != b.size() ? a.size() < b.size() : a < b;
}

struct RangeText {
    std::optional<std::string> chain;
    std::string first, last;  // decimal, first <= last
};

std::size_t digitRun(const std::u32string& s, std::size_t at) {
    std::size_t end = at;
    while (end < s.size() && isPythonDigit(s[end])) {
        ++end;
    }
    return end;
}

// parse_residue_ranges with the numbers as Python ints (decimal texts)
std::vector<RangeText> parseRangeTexts(std::string_view text) {
    std::vector<RangeText> items;
    const std::string stripped = pythonStrip(text);
    if (stripped.empty()) {
        return items;
    }
    // re.sub(r"\s*([-:])\s*", r"\1", ...)
    std::u32string normalised;
    bool skipSpace = false;
    for (const char32_t c : decodeUtf8(stripped)) {
        if (c == U'-' || c == U':') {
            while (!normalised.empty() && isPythonSpace(normalised.back())) {
                normalised.pop_back();
            }
            normalised += c;
            skipSpace = true;
        } else if (skipSpace && isPythonSpace(c)) {
            continue;
        } else {
            skipSpace = false;
            normalised += c;
        }
    }
    for (const std::string& raw : splitItems(encodeUtf8(normalised))) {
        // ^(?:(?P<chain>[A-Za-z][A-Za-z0-9]*):)?(?P<a>\d+)(?:[-:](?P<b>\d+))?$
        const std::u32string s = decodeUtf8(raw);
        RangeText item;
        std::size_t pos = 0;
        if (!s.empty() && isAsciiLetter(s[0])) {
            std::size_t end = 1;
            while (end < s.size() && (isAsciiLetter(s[end]) || isAsciiDigit(s[end]))) {
                ++end;
            }
            if (end < s.size() && s[end] == U':') {
                item.chain = encodeUtf8(s.substr(0, end));
                pos = end + 1;
            }
        }
        const std::size_t aEnd = digitRun(s, pos);
        bool ok = aEnd > pos;
        std::size_t bEnd = aEnd;
        if (ok && aEnd < s.size()) {
            ok = (s[aEnd] == U'-' || s[aEnd] == U':');
            bEnd = digitRun(s, aEnd + 1);
            ok = ok && bEnd > aEnd + 1 && bEnd == s.size();
        }
        if (!ok) {
            throw SelectionError("invalid residue range item " + repr(raw));
        }
        item.first = decimalText(std::u32string_view(s).substr(pos, aEnd - pos));
        item.last = bEnd > aEnd
                        ? decimalText(std::u32string_view(s).substr(aEnd + 1, bEnd - aEnd - 1))
                        : item.first;
        if (decimalLess(item.last, item.first)) {
            std::swap(item.first, item.last);
        }
        items.push_back(std::move(item));
    }
    return items;
}

const std::map<std::string, std::string, std::less<>>& ionElements() {
    static const std::map<std::string, std::string, std::less<>> table = {
        {"MG", "MG"},  {"MG2", "MG"}, {"NA", "NA"},  {"NA+", "NA"}, {"SOD", "NA"},
        {"K", "K"},    {"K+", "K"},   {"POT", "K"},  {"CL", "CL"},  {"CL-", "CL"},
        {"CLA", "CL"}, {"CA", "CA"},  {"CAL", "CA"}, {"ZN", "ZN"},  {"MN", "MN"},
        {"FE", "FE"},  {"LI", "LI"},  {"RB", "RB"},  {"CS", "CS"},  {"BR", "BR"},
    };
    return table;
}

std::string rstripChars(std::string text, std::string_view chars) {
    while (!text.empty() && chars.find(text.back()) != std::string_view::npos) {
        text.pop_back();
    }
    return text;
}

bool isAlpha(char32_t c) {
    if (c < 0x80) {
        return isAsciiLetter(c);
    }
    return isPythonWordCharacter(c) && !isPythonDigit(c);
}

// "%.3f"
std::string fixed3(double x) {
    char buffer[400];
    std::snprintf(buffer, sizeof buffer, "%.3f", x);
    return buffer;
}

// A frame of parse_frame_spec: an int64, or a Python int beyond 64 bits (always outside)
struct FrameValue {
    std::int64_t value = 0;
    std::optional<std::string> big;  // str() of the Python int
};

// |big| - small for a decimal text big > small >= 0
std::string decimalMinus(std::string big, std::int64_t small) {
    std::string s = std::to_string(small);
    int borrow = 0;
    for (std::size_t i = 0; i < big.size(); ++i) {
        const std::size_t at = big.size() - 1 - i;
        int d = big[at] - '0' - borrow - (i < s.size() ? s[s.size() - 1 - i] - '0' : 0);
        borrow = d < 0 ? 1 : 0;
        big[at] = static_cast<char>('0' + (d + 10) % 10);
    }
    const std::size_t nonZero = big.find_first_not_of('0');
    return nonZero == std::string::npos ? "0" : big.substr(nonZero);
}

std::string listRepr(const std::vector<FrameValue>& values) {
    std::string out = "[";
    for (std::size_t i = 0; i < values.size(); ++i) {
        out += (i ? ", " : "") + (values[i].big ? *values[i].big : std::to_string(values[i].value));
    }
    return out + "]";
}

std::optional<std::int64_t> toInt64(const std::string& decimal) {
    try {
        return pythonInt(decimal);
    } catch (const std::overflow_error&) {
        return std::nullopt;
    }
}

}  // namespace

std::vector<ResidueRange> parseResidueRanges(std::string_view text) {
    std::vector<ResidueRange> out;
    for (const RangeText& item : parseRangeTexts(text)) {
        const auto first = toInt64(item.first), last = toInt64(item.last);
        if (!first || !last) {
            throw SelectionError("residue number too large in " + repr(text));
        }
        out.push_back({item.chain, *first, *last});
    }
    return out;
}

std::string rangesToSelection(std::string_view text, const Topology* topology) {
    const std::vector<RangeText> items = parseRangeTexts(text);
    if (items.empty()) {
        throw SelectionError("empty residue range");
    }
    const std::string chainAttr = topology && topology->chainIDs ? "chainID" : "segid";
    std::string out;
    for (const RangeText& item : items) {
        const std::string resid = item.first == item.last ? "resid " + item.first
                                                          : "resid " + item.first + ":" + item.last;
        out += out.empty() ? "" : " or ";
        out += item.chain ? "(" + chainAttr + " " + *item.chain + " and " + resid + ")"
                          : "(" + resid + ")";
    }
    return out;
}

std::string groupSelection(const ResidueGroup& group, const Topology* topology) {
    if (const std::string selection = pythonStrip(group.selection); !selection.empty()) {
        return "(" + selection + ")";
    }
    if (!pythonStrip(group.residues).empty()) {
        return "(" + rangesToSelection(group.residues, topology) + ")";
    }
    throw SelectionError("group " + repr(group.name) + " has neither residues nor a selection");
}

namespace {

bool isGroupNameCharacter(char c) {
    return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9') || c == '_' ||
           c == '.' || c == '+' || c == '-' || c == '\'';
}

std::string expand(std::string_view text, const std::vector<ResidueGroup>& groups,
                   const Topology* topology, int depth) {
    if (depth > 8) {
        throw SelectionError("group definitions are nested too deeply (circular reference?)");
    }
    std::string out;
    std::size_t pos = 0;
    for (;;) {
        std::size_t at = text.find("group:", pos);
        std::size_t end = at;
        while (at != std::string_view::npos) {
            end = at + 6;
            while (end < text.size() && isGroupNameCharacter(text[end])) {
                ++end;
            }
            if (end > at + 6) {
                break;
            }
            at = text.find("group:", at + 1);
        }
        if (at == std::string_view::npos) {
            out += text.substr(pos);
            return out;
        }
        out += text.substr(pos, at - pos);
        const std::string_view name = text.substr(at + 6, end - at - 6);
        // {g.name: g}: the last group of a name wins, names listed in first-seen order
        const ResidueGroup* found = nullptr;
        std::vector<std::string> names;
        for (const ResidueGroup& g : groups) {
            if (std::find(names.begin(), names.end(), g.name) == names.end()) {
                names.push_back(g.name);
            }
            if (g.name == name) {
                found = &g;
            }
        }
        if (!found) {
            std::string defined;
            for (const std::string& n : names) {
                defined += (defined.empty() ? "" : ", ") + n;
            }
            throw SelectionError("unknown residue group " + repr(name) +
                                 " (defined: " + (defined.empty() ? "none" : defined) + ")");
        }
        out += expand(groupSelection(*found, topology), groups, topology, depth + 1);
        pos = end;
    }
}

}  // namespace

std::string expandGroups(std::string_view text, const std::vector<ResidueGroup>& groups,
                         const Topology* topology) {
    return expand(text, groups, topology, 0);
}

std::string resolveSelection(std::string_view text, const std::vector<ResidueGroup>& groups,
                             const Topology* topology) {
    std::string resolved = pythonStrip(text);
    if (resolved.empty()) {
        throw SelectionError("empty selection");
    }
    if (resolved.starts_with("sel:")) {
        resolved = pythonStrip(std::string_view(resolved).substr(4));
    } else if (std::any_of(groups.begin(), groups.end(),
                           [&](const ResidueGroup& g) { return g.name == resolved; })) {
        resolved = "group:" + resolved;
    }
    return expandGroups(resolved, groups, topology);
}

AtomGroup select(const std::shared_ptr<const Topology>& topology, const Frame* current,
                 std::string_view text, const std::vector<ResidueGroup>& groups, bool allowEmpty) {
    const std::string resolved = resolveSelection(text, groups, topology.get());
    std::vector<std::int32_t> atoms;
    try {
        atoms = selectAtoms(*topology, current, resolved);
    } catch (const std::exception& e) {
        throw SelectionError("invalid selection " + repr(text) + " -> " + repr(resolved) + ": " +
                             e.what());
    }
    if (!allowEmpty && atoms.empty()) {
        throw SelectionError("selection " + repr(text) + " matched no atoms");
    }
    return AtomGroup(topology, std::move(atoms));
}

AtomGroup select(const System& system, std::string_view text,
                 const std::vector<ResidueGroup>& groups, bool allowEmpty) {
    return select(system.topology, system.currentFrame(), text, groups, allowEmpty);
}

std::string guessElement(std::string_view name, std::string_view resname, bool singleAtomResidue) {
    const std::string nm = upper(pythonStrip(name));
    const std::string rn = upper(pythonStrip(resname));
    if (singleAtomResidue) {
        for (const std::string& key :
             {nm, rn, rstripChars(nm, "+-0123456789"), rstripChars(rn, "+-0123456789")}) {
            if (const auto it = ionElements().find(key); it != ionElements().end()) {
                return it->second;
            }
        }
    }
    const std::size_t start = nm.find_first_not_of("0123456789");
    if (start == std::string::npos) {
        return "";
    }
    for (const char32_t c : decodeUtf8(std::string_view(nm).substr(start))) {
        if (isAlpha(c)) {
            std::string out;
            appendUtf8(out, c);
            return out;
        }
    }
    return "";
}

std::vector<std::string> atomElements(const AtomGroup& atoms) {
    const Topology& t = atoms.topology();
    std::vector<std::string> out;
    out.reserve(atoms.size());
    if (t.elements && !atoms.empty()) {
        bool all = true;
        for (const std::int32_t a : atoms.indices()) {
            std::string e = pythonStrip(t.strings[(*t.elements)[static_cast<std::size_t>(a)]]);
            if (e.empty()) {
                all = false;
                break;
            }
            out.push_back(upper(std::move(e)));
        }
        if (all) {
            return out;
        }
        out.clear();
    }
    for (const std::int32_t a : atoms.indices()) {
        const std::int32_t r = t.atomResidue[static_cast<std::size_t>(a)];
        out.push_back(guessElement(t.name(a), t.resname(r), t.residueAtoms(r).size() == 1));
    }
    return out;
}

std::vector<bool> heavyMask(const AtomGroup& atoms) {
    std::vector<bool> out;
    for (const std::string& e : atomElements(atoms)) {
        out.push_back(e != "H");
    }
    return out;
}

AtomGroup heavy(const AtomGroup& atoms) { return atoms.masked(heavyMask(atoms)); }

AtomGroup elementFilter(const AtomGroup& atoms, const std::vector<std::string>& elements) {
    std::set<std::string> wanted;
    for (const std::string& e : elements) {
        if (const std::string s = pythonStrip(e); !s.empty()) {
            wanted.insert(upper(s));
        }
    }
    std::vector<bool> mask;
    for (const std::string& e : atomElements(atoms)) {
        mask.push_back(wanted.count(e) > 0);
    }
    return atoms.masked(mask);
}

std::vector<std::string> splitCsv(std::string_view text) {
    std::vector<std::string> out;
    std::size_t start = 0;
    for (;;) {
        const std::size_t comma = text.find(',', start);
        if (std::string piece = pythonStrip(text.substr(start, comma - start)); !piece.empty()) {
            out.push_back(std::move(piece));
        }
        if (comma == std::string_view::npos) {
            return out;
        }
        start = comma + 1;
    }
}

std::vector<std::string> residueLabels(const Topology& topology,
                                       std::span<const std::int32_t> residues) {
    std::set<std::int32_t> segments;
    for (const std::int32_t r : residues) {
        segments.insert(topology.residueSegment[static_cast<std::size_t>(r)]);
    }
    std::set<std::string> segids;
    for (const std::int32_t s : segments) {
        segids.insert(topology.segid(s));
    }
    const bool multi = segids.size() > 1;
    std::vector<std::string> out;
    for (const std::int32_t r : residues) {
        std::string label =
            topology.resname(r) + std::to_string(topology.resids[static_cast<std::size_t>(r)]);
        if (multi) {
            label =
                topology.segid(topology.residueSegment[static_cast<std::size_t>(r)]) + ":" + label;
        }
        out.push_back(std::move(label));
    }
    return out;
}

std::vector<std::int64_t> parseFrameSpec(std::string_view text, std::int64_t frameCount) {
    std::vector<FrameValue> frames;
    if (pythonStrip(text).empty()) {
        return {};
    }
    auto add = [&](std::int64_t v) { frames.push_back({v, std::nullopt}); };
    for (const std::string& raw : splitItems(text)) {
        const std::string token = lower(raw);
        if (token == "last" || token == "end") {
            add(frameCount - 1);
            continue;
        }
        if (token == "first" || token == "start") {
            add(0);
            continue;
        }
        // (\d+)-(\d+|last)(?::(\d+))?
        const std::u32string s = decodeUtf8(token);
        const std::size_t aEnd = digitRun(s, 0);
        if (aEnd > 0 && aEnd < s.size() && s[aEnd] == U'-') {
            std::size_t bEnd = digitRun(s, aEnd + 1);
            bool isLast = false;
            if (bEnd == aEnd + 1 && s.compare(aEnd + 1, 4, U"last") == 0) {
                isLast = true;
                bEnd = aEnd + 5;
            }
            std::optional<std::string> step;
            bool ok = bEnd > aEnd + 1;
            if (ok && bEnd < s.size()) {
                const std::size_t stepEnd = digitRun(s, bEnd + 1);
                ok = s[bEnd] == U':' && stepEnd > bEnd + 1 && stepEnd == s.size();
                if (ok) {
                    step = decimalText(std::u32string_view(s).substr(bEnd + 1, stepEnd - bEnd - 1));
                }
            }
            if (ok) {
                const auto a = toInt64(decimalText(std::u32string_view(s).substr(0, aEnd)));
                const auto b = isLast ? std::optional<std::int64_t>(frameCount - 1)
                                      : toInt64(decimalText(std::u32string_view(s).substr(
                                            aEnd + 1, bEnd - aEnd - 1)));
                const auto stride = step ? toInt64(*step) : std::optional<std::int64_t>(1);
                if (stride && *stride < 1) {
                    throw SelectionError("invalid frame step in " + repr(raw));
                }
                if (!b && a) {
                    throw SelectionError("frame range too large in " + repr(raw));
                }
                if (a && b && !stride && *a <= *b) {
                    add(*a);  // a step beyond 64 bits: the first frame only
                }
                if (a && b && stride && *a <= *b) {
                    if ((*b - *a) / *stride > 100'000'000) {
                        throw SelectionError("frame range too large in " + repr(raw));
                    }
                    for (std::int64_t f = *a; f <= *b; f += *stride) {
                        add(f);
                        if (*b - f < *stride) {
                            break;
                        }
                    }
                }
                continue;
            }
        }
        std::int64_t v = 0;
        try {
            v = pythonInt(token);
        } catch (const ValueError&) {
            throw SelectionError("invalid frame specification " + repr(raw));
        } catch (const std::overflow_error&) {
            // a Python int beyond 64 bits: outside the trajectory either way
            std::string digits = pythonStrip(token);
            const bool negative = digits[0] == '-';
            std::u32string body = decodeUtf8(digits);
            body.erase(std::remove(body.begin(), body.end(), U'_'), body.end());
            if (body[0] == U'-' || body[0] == U'+') {
                body.erase(0, 1);
            }
            const std::string magnitude = decimalText(body);
            frames.push_back({0, negative ? "-" + decimalMinus(magnitude, frameCount) : magnitude});
            continue;
        }
        add(v < 0 ? v + frameCount : v);
    }
    std::vector<FrameValue> bad;
    std::set<std::int64_t> good;
    for (const FrameValue& f : frames) {
        if (f.big || f.value < 0 || f.value >= frameCount) {
            bad.push_back(f);
        } else {
            good.insert(f.value);
        }
    }
    if (!bad.empty()) {
        throw SelectionError("frame(s) " + listRepr(bad) + " outside trajectory (0-" +
                             std::to_string(frameCount - 1) + ")");
    }
    return {good.begin(), good.end()};
}

std::vector<double> parseTimeSpec(std::string_view text, std::string_view defaultUnit) {
    std::vector<double> times;
    if (pythonStrip(text).empty()) {
        return times;
    }
    if (defaultUnit == "frame") {
        throw SelectionError(
            "time specification needs a time unit when the display unit is 'frame'");
    }
    for (const std::string& raw : splitItems(text)) {
        // ([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s*([a-zA-Z]*) as a full match
        const std::u32string s = decodeUtf8(raw);
        std::size_t i = 0;
        if (i < s.size() && (s[i] == U'-' || s[i] == U'+')) {
            ++i;
        }
        const std::size_t intEnd = digitRun(s, i);
        std::size_t end = 0;
        if (intEnd < s.size() && s[intEnd] == U'.' && digitRun(s, intEnd + 1) > intEnd + 1) {
            end = digitRun(s, intEnd + 1);
        } else if (intEnd > i) {
            end = intEnd;
        }
        bool ok = end > 0;
        if (ok && end < s.size() && (s[end] == U'e' || s[end] == U'E')) {
            std::size_t j = end + 1;
            if (j < s.size() && (s[j] == U'-' || s[j] == U'+')) {
                ++j;
            }
            if (digitRun(s, j) > j) {
                end = digitRun(s, j);
            }
        }
        std::size_t unitStart = end;
        while (ok && unitStart < s.size() && isPythonSpace(s[unitStart])) {
            ++unitStart;
        }
        for (std::size_t k = unitStart; ok && k < s.size(); ++k) {
            ok = isAsciiLetter(s[k]);
        }
        if (!ok) {
            throw SelectionError("invalid time " + repr(raw));
        }
        const std::string number = encodeUtf8(std::u32string_view(s).substr(0, end));
        const std::string unitText = encodeUtf8(std::u32string_view(s).substr(unitStart));
        const std::string unit = lower(unitText.empty() ? std::string(defaultUnit) : unitText);
        if (unit != "ps" && unit != "ns" && unit != "us" && unit != "fs") {
            throw SelectionError("unknown time unit " + repr(unit) + " in " + repr(raw));
        }
        times.push_back(pythonFloat(number) * timeUnitFactor(unit));
    }
    return times;
}

double timeUnitFactor(std::string_view unit) {
    if (unit == "ps") {
        return 1.0;
    }
    if (unit == "ns") {
        return 1000.0;
    }
    if (unit == "us") {
        return 1.0e6;
    }
    if (unit == "fs") {
        return 1.0e-3;
    }
    throw ValueError("unknown time unit " + repr(unit));
}

AtomGroup resolveExportSelection(const System& system, std::string_view text,
                                 const std::vector<ResidueGroup>& groups) {
    std::string selection = pythonStrip(text);
    if (selection.empty()) {
        return system.atoms();
    }
    const std::vector<std::string> names = splitCsv(selection);
    const bool allGroups =
        !names.empty() && std::all_of(names.begin(), names.end(), [&](const std::string& n) {
            return std::any_of(groups.begin(), groups.end(),
                               [&](const ResidueGroup& g) { return g.name == n; });
        });
    if (allGroups) {
        selection.clear();
        for (const std::string& n : names) {
            selection += (selection.empty() ? "group:" : " or group:") + n;
        }
    }
    return select(system, selection, groups);
}

std::vector<std::int64_t> framesFromSpecs(Trajectory& trajectory, std::string_view framesText,
                                          std::string_view timesText, std::string_view timeUnit,
                                          std::optional<double> dtPs) {
    const std::int64_t n = trajectory.frameCount();
    const std::vector<std::int64_t> parsed = parseFrameSpec(framesText, n);
    std::set<std::int64_t> frames(parsed.begin(), parsed.end());
    const std::vector<double> times =
        parseTimeSpec(timesText, timeUnit != "frame" ? timeUnit : "ps");
    if (!times.empty()) {
        std::vector<double> all(static_cast<std::size_t>(n));
        if (dtPs && *dtPs != 0.0) {
            for (std::int64_t i = 0; i < n; ++i) {
                all[static_cast<std::size_t>(i)] = static_cast<double>(i) * *dtPs;
            }
        } else {
            for (std::int64_t i = 0; i < n; ++i) {
                all[static_cast<std::size_t>(i)] = trajectory.seek(i).time;
            }
            if (n > 0) {
                trajectory.seek(0);  // iterating an MDAnalysis trajectory rewinds it
            }
        }
        if (all.empty()) {
            throw ValueError("attempt to get argmin of an empty sequence");
        }
        double tolerance = 1.0;
        if (n > 1) {
            std::vector<double> steps;
            for (std::size_t i = 1; i < all.size(); ++i) {
                steps.push_back(all[i] - all[i - 1]);
            }
            std::sort(steps.begin(), steps.end());
            const std::size_t mid = steps.size() / 2;
            const double median =
                steps.size() % 2 ? steps[mid] : (steps[mid - 1] + steps[mid]) / 2.0;
            tolerance = 0.5 * median;
        }
        tolerance = tolerance > 1e-6 ? tolerance : 1e-6;
        for (const double t : times) {
            std::size_t k = 0;
            for (std::size_t i = 1; i < all.size(); ++i) {
                if (std::fabs(all[i] - t) < std::fabs(all[k] - t)) {
                    k = i;
                }
            }
            if (std::fabs(all[k] - t) > tolerance) {
                logWarning("time " + fixed3(t) + " ps is outside the trajectory; using frame " +
                           std::to_string(k) + " (" + fixed3(all[k]) + " ps)");
            }
            frames.insert(static_cast<std::int64_t>(k));
        }
    }
    return {frames.begin(), frames.end()};
}

}  // namespace mdbridge::core
