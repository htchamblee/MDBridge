#pragma once
// "Did you mean ...?" for names (MDATK 0.1.x mdatk/selections.py suggest): the candidate a name
// was probably meant to be, by case or by a typo.

#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace mdbridge::core {

// Similarity of two strings (code points) in [0, 1]: 2 M / (len(a) + len(b)), M the characters of
// the matching blocks found by taking the longest common block (the earliest in `a`, then in `b`,
// of the longest) and recursing on both sides of it (Ratcliff / Obershelp; the value of Python's
// difflib.SequenceMatcher(None, a, b).ratio() while `b` is shorter than 200 characters). 1 for two
// empty strings.
double similarity(std::u32string_view a, std::u32string_view b);

// The candidate `name` was probably meant to be: the first candidate equal to it ignoring case,
// else the candidate whose lower-case spelling is most similar to name's (similarity at least 0.6;
// on a tie the one whose lower-case spelling sorts last; of candidates with the same lower-case
// spelling the last). Empty candidates are skipped. Nothing for a name of 200 characters or more.
std::optional<std::string> suggestName(std::string_view name, const std::vector<std::string>& candidates);

}  // namespace mdbridge::core
