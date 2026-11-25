// Citations (mdatk/citations.py of MDATK v0.0.4.1): the tables and every function, against what
// MDATK computes with MDBridge's changes CIT-01 / CIT-02 (tests/parity/improvements_citations.py;
// tests/parity/make_citation_cases.py -> data/citations; see data/README.md).

#include <algorithm>
#include <optional>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "mdbridge/core/citations.hpp"
#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/project.hpp"
#include "test_support.hpp"

using namespace mdbridge::core;
using namespace mdbridge::test;
namespace cit = mdbridge::core::citations;

namespace {

const Value& cases() {
    static const Value data = readJson(kData / "citations" / "citation_cases.json");
    return data;
}

const Array& items(const Value& value, std::string_view key) {
    return *value.find(key)->get<Array>();
}

std::vector<std::string> strings(const Value& value) {
    std::vector<std::string> out;
    for (const Value& item : *value.get<Array>()) {
        out.push_back(item.stringOr());
    }
    return out;
}

std::string joined(const std::vector<std::string>& items) {
    std::string out;
    for (const std::string& item : items) {
        out += (out.empty() ? "[" : ", ") + item;
    }
    return out + "]";
}

}  // namespace

TEST(Citations, ReferenceTable) {
    const Array& expected = items(cases(), "references");
    const auto table = cit::referenceTable();
    ASSERT_EQ(table.size(), expected.size());
    for (std::size_t i = 0; i < table.size(); ++i) {
        const Value& e = expected[i];
        EXPECT_EQ(table[i].key, text(e, "key"));
        EXPECT_EQ(table[i].text, text(e, "text")) << text(e, "key");
        EXPECT_EQ(table[i].bibtex, text(e, "bibtex")) << text(e, "key");
        EXPECT_EQ(table[i].doi, text(e, "doi")) << text(e, "key");
        EXPECT_EQ(table[i].url, text(e, "url")) << text(e, "key");
        EXPECT_EQ(table[i].verified, text(e, "verified")) << text(e, "key");
        EXPECT_EQ(cit::findReference(table[i].key), &table[i]);
    }
    EXPECT_EQ(cit::findReference("nope"), nullptr);
    EXPECT_EQ(cit::findReference(""), nullptr);
}

TEST(Citations, ResourceTable) {
    const Array& expected = items(cases(), "resources");
    const auto table = cit::resourceTable();
    ASSERT_EQ(table.size(), expected.size());
    for (std::size_t i = 0; i < table.size(); ++i) {
        const Value& e = expected[i];
        EXPECT_EQ(table[i].id, text(e, "id"));
        EXPECT_EQ(table[i].name, text(e, "name")) << text(e, "id");
        EXPECT_EQ(table[i].note, text(e, "note")) << text(e, "id");
        EXPECT_EQ(std::vector<std::string>(table[i].keys.begin(), table[i].keys.end()),
                  strings(*e.find("keys")))
            << text(e, "id");
        EXPECT_EQ(cit::findResource(table[i].id), &table[i]);
    }
    EXPECT_EQ(cit::findResource("nope"), nullptr);
}

TEST(Citations, UniqueAndKeys) {
    for (const Value& c : items(cases(), "unique")) {
        const auto input = strings(*c.find("items"));
        EXPECT_EQ(cit::unique(input), strings(*c.find("out"))) << joined(input);
    }
    for (const Value& c : items(cases(), "keys_for")) {
        const auto input = strings(*c.find("ids"));
        EXPECT_EQ(cit::keysFor(input), strings(*c.find("out"))) << joined(input);
    }
}

TEST(Citations, ReferencesBibtexText) {
    for (const Value& c : items(cases(), "references_of")) {
        const auto keys = strings(*c.find("keys"));
        std::vector<std::string> got;
        for (const cit::Reference& r : cit::references(keys)) {
            got.emplace_back(r.key);
        }
        EXPECT_EQ(got, strings(*c.find("out"))) << joined(keys);
    }
    for (const Value& c : items(cases(), "bibtex")) {
        const auto keys = strings(*c.find("keys"));
        EXPECT_EQ(cit::bibtex(keys), text(c, "out")) << joined(keys);
    }
    for (const Value& c : items(cases(), "text_lines")) {
        const auto keys = strings(*c.find("keys"));
        EXPECT_EQ(cit::textLines(keys), strings(*c.find("out"))) << joined(keys);
    }
}

TEST(Citations, Grouped) {
    for (const Value& c : items(cases(), "grouped")) {
        const auto ids = strings(*c.find("ids"));
        const auto got = cit::grouped(ids);
        const Array& expected = items(c, "out");
        ASSERT_EQ(got.size(), expected.size()) << joined(ids);
        for (std::size_t i = 0; i < got.size(); ++i) {
            const Array& e = *expected[i].get<Array>();
            EXPECT_EQ(got[i].first.id, e[0].stringOr()) << joined(ids);
            std::vector<std::string> keys;
            for (const cit::Reference& r : got[i].second) {
                keys.emplace_back(r.key);
            }
            EXPECT_EQ(keys, strings(e[1])) << joined(ids);
        }
    }
}

TEST(Citations, CommentLines) {
    for (const Value& c : items(cases(), "comment_lines")) {
        const auto ids = strings(*c.find("ids"));
        const std::string prefix = text(c, "prefix");
        const int width = static_cast<int>(*c.find("width")->get<std::int64_t>());
        std::optional<std::string> heading;
        if (const std::string* h = c.find("heading")->get<std::string>()) {
            heading = *h;
        }
        const std::string what = joined(ids) + " '" + prefix + "' " + std::to_string(width);
        try {
            const auto lines = cit::commentLines(
                ids, prefix, width,
                heading ? std::optional<std::string_view>(*heading) : std::nullopt);
            if (const Value* error = c.find("error")) {
                ADD_FAILURE() << what << ": no error, expected " << error->stringOr();
            } else {
                EXPECT_EQ(lines, strings(*c.find("lines"))) << what;
            }
        } catch (const ValueError& e) {
            EXPECT_EQ(e.what(), text(c, "error")) << what;
        }
    }
    // the defaults
    const std::vector<std::string> ids = {"ffmpeg"};
    EXPECT_EQ(cit::commentLines(ids),
              (std::vector<std::string>{"# Please cite (FFmpeg):",
                                        "#   - FFmpeg Developers. FFmpeg. https://ffmpeg.org"}));
}

TEST(Citations, PipelineResources) {
    for (const Value& c : items(cases(), "pipeline")) {
        const Project project = projectFromPlain(*c.find("project"));
        EXPECT_EQ(cit::pipelineResources(project), strings(*c.find("resources")))
            << toJson(*c.find("project"), -1);
    }
}

TEST(Citations, ProjectResources) {
    for (const Value& c : items(cases(), "project")) {
        const Project project = projectFromPlain(*c.find("project"));
        std::vector<std::vector<std::string>> steps;
        for (const Value& step : items(c, "generating")) {
            steps.push_back(strings(step));
        }
        EXPECT_EQ(cit::projectResources(project, steps), strings(*c.find("resources")))
            << toJson(c, -1);
    }
}

// CIT-01 / CIT-02 (docs/DIFFERENCES.md "Citations"): the software MDATK ran on is not cited and is
// not in the table; the RMSD's superposition methods are cited by their originators, only when used.
TEST(Citations, MethodsByTheirOriginators) {
    for (const char* id : {"mdanalysis", "numpy", "scipy", "matplotlib"}) {
        EXPECT_EQ(cit::findResource(id), nullptr) << id;
    }
    for (const char* key : {"michaudagrawal2011", "gowers2016", "harris2020", "virtanen2020", "hunter2007"}) {
        EXPECT_EQ(cit::findReference(key), nullptr) << key;
    }
    using Ids = std::vector<std::string>;
    EXPECT_EQ(cit::keysFor(Ids{"qcp"}), (Ids{"theobald2005", "liu2010"}));
    EXPECT_EQ(cit::keysFor(Ids{"kabsch"}), (Ids{"kabsch1976", "kabsch1978"}));
    // MDATK 0.1.1's resource of its generic Lennard-Jones table resolves (MDBridge's typed layer's
    // references), but is never a project resource: MDBridge cites the energy parameters per run
    EXPECT_EQ(cit::keysFor(Ids{"generic_lj"}), (Ids{"cornell1995", "wang2000", "jorgensen1983", "joung2008", "li2013"}));
    for (const char* key : {"theobald2005", "liu2010", "kabsch1976", "kabsch1978"}) {
        const cit::Reference* r = cit::findReference(key);
        ASSERT_NE(r, nullptr) << key;
        EXPECT_FALSE(r->doi.empty()) << key;
        EXPECT_FALSE(r->verified.empty()) << key;
    }

    Project p;  // defaults: RMSD on, two targets without a fit selection, superposition, video
    EXPECT_EQ(cit::pipelineResources(p), (Ids{"qcp", "ffmpeg"}));
    p.output.video.enabled = false;
    p.analyses.rmsd.targets[1].fit_selection = "name P";
    EXPECT_EQ(cit::pipelineResources(p), (Ids{"qcp", "kabsch"}));
    p.analyses.rmsd.targets[0].fit_selection = " \t";
    EXPECT_EQ(cit::pipelineResources(p), (Ids{"qcp", "kabsch"}));
    p.analyses.rmsd.targets[0].fit_selection = "name C4'";
    EXPECT_EQ(cit::pipelineResources(p), (Ids{"kabsch"}));
    p.analyses.rmsd.superposition = false;  // plain RMSD: no superposition method
    EXPECT_EQ(cit::pipelineResources(p), Ids{});
    p.analyses.rmsd.superposition = true;
    p.analyses.rmsd.enabled = false;
    EXPECT_EQ(cit::pipelineResources(p), Ids{});
    p.analyses.rmsd.enabled = true;
    p.analyses.rmsd.targets.clear();
    EXPECT_EQ(cit::pipelineResources(p), Ids{});

    // CIT-02: Barnaba is never run, so never cited
    Project b;
    b.analyses.basepairs.engine = "barnaba";
    PanelSpec secondary;
    secondary.type = "secondary";
    b.layout.panels = {secondary};
    const Ids ids = cit::pipelineResources(b);
    EXPECT_EQ(std::find(ids.begin(), ids.end(), "barnaba"), ids.end());
}
