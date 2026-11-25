#include "mdbridge/core/suggest.hpp"

#include <array>
#include <cstddef>
#include <utility>
#include <vector>

#include "mdbridge/core/text.hpp"

namespace mdbridge::core {

namespace {

struct Block {
    std::size_t i = 0, j = 0, size = 0;
};

// The longest common block of a[alo, ahi) and b[blo, bhi): of the longest, the one that starts
// earliest in a, then earliest in b.
Block longestBlock(std::u32string_view a, std::u32string_view b, std::size_t alo, std::size_t ahi,
                   std::size_t blo, std::size_t bhi) {
    Block best{alo, blo, 0};
    // run[j + 1]: length of the common block ending at (i - 1, j) / (i, j)
    std::vector<std::size_t> previous(bhi - blo + 1, 0), current(bhi - blo + 1, 0);
    for (std::size_t i = alo; i < ahi; ++i) {
        for (std::size_t j = blo; j < bhi; ++j) {
            const std::size_t k = a[i] == b[j] ? previous[j - blo] + 1 : 0;
            current[j - blo + 1] = k;
            if (k > best.size) {
                best = {i + 1 - k, j + 1 - k, k};
            }
        }
        std::swap(previous, current);
    }
    return best;
}

std::size_t matchingCharacters(std::u32string_view a, std::u32string_view b) {
    std::size_t total = 0;
    std::vector<std::array<std::size_t, 4>> queue{{0, a.size(), 0, b.size()}};
    while (!queue.empty()) {
        const auto [alo, ahi, blo, bhi] = queue.back();
        queue.pop_back();
        const Block m = longestBlock(a, b, alo, ahi, blo, bhi);
        if (m.size == 0) {
            continue;
        }
        total += m.size;
        if (alo < m.i && blo < m.j) {
            queue.push_back({alo, m.i, blo, m.j});
        }
        if (m.i + m.size < ahi && m.j + m.size < bhi) {
            queue.push_back({m.i + m.size, ahi, m.j + m.size, bhi});
        }
    }
    return total;
}

}  // namespace

double similarity(std::u32string_view a, std::u32string_view b) {
    const std::size_t length = a.size() + b.size();
    if (length == 0) {
        return 1.0;
    }
    return 2.0 * static_cast<double>(matchingCharacters(a, b)) / static_cast<double>(length);
}

std::optional<std::string> suggestName(std::string_view name, const std::vector<std::string>& candidates) {
    const std::string lowerName = pythonLower(name);
    for (const std::string& c : candidates) {
        if (!c.empty() && pythonLower(c) == lowerName) {
            return c;
        }
    }
    const std::u32string word = decodeUtf8(lowerName);
    if (word.size() >= 200) {
        return std::nullopt;
    }
    // lower-case spelling -> the last candidate of that spelling
    std::vector<std::pair<std::u32string, const std::string*>> lowered;
    for (const std::string& c : candidates) {
        if (c.empty()) {
            continue;
        }
        std::u32string key = decodeUtf8(pythonLower(c));
        bool seen = false;
        for (auto& [k, original] : lowered) {
            if (k == key) {
                original = &c;
                seen = true;
            }
        }
        if (!seen) {
            lowered.emplace_back(std::move(key), &c);
        }
    }
    const std::pair<std::u32string, const std::string*>* best = nullptr;
    double bestScore = 0.0;
    for (const auto& entry : lowered) {
        const double score = similarity(entry.first, word);
        if (score < 0.6) {
            continue;
        }
        if (!best || score > bestScore || (score == bestScore && entry.first > best->first)) {
            best = &entry;
            bestScore = score;
        }
    }
    if (!best) {
        return std::nullopt;
    }
    return *best->second;
}

}  // namespace mdbridge::core
