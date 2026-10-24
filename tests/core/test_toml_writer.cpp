// The commented TOML writer against the vectors MDATK v0.0.4.1 produced
// (tests/parity/make_writer_vectors.py -> data/project-file/vectors).

#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "mdbridge/core/text.hpp"
#include "mdbridge/core/toml_writer.hpp"
#include "test_support.hpp"

using namespace mdbridge::core;
using namespace mdbridge::test;

namespace {

const fs::path kVectors = kData / "project-file" / "vectors";

std::vector<std::string> strings(const Value& array) {
    std::vector<std::string> out;
    for (const Value& item : *array.get<Array>()) {
        out.push_back(*item.get<std::string>());
    }
    return out;
}

}  // namespace

TEST(TomlWriterVectors, Textwrap) {
    const Value cases = readJson(kVectors / "textwrap_vectors.json");
    ASSERT_GT(cases.size(), 1000u);
    for (const Value& c : *cases.get<Array>()) {
        const int width = static_cast<int>(*c.find("width")->get<std::int64_t>());
        EXPECT_EQ(pythonTextwrap(text(c, "text"), width), strings(*c.find("lines")))
            << "text " << text(c, "text") << " width " << width;
    }
}

TEST(TomlWriterVectors, Comment) {
    const Value cases = readJson(kVectors / "comment_vectors.json");
    for (const Value& c : *cases.get<Array>()) {
        TomlWriter writer;
        const std::string* commentText = c.find("text")->get<std::string>();  // null or 0: no comment
        writer.comment(commentText ? std::optional<std::string>(*commentText) : std::nullopt, text(c, "indent"));
        EXPECT_EQ(writer.lines(), strings(*c.find("lines"))) << (commentText ? *commentText : "(none)");
    }
}

TEST(TomlWriterVectors, Values) {
    const Value cases = readJson(kVectors / "value_vectors.json");
    for (const Value& c : *cases.find("values")->get<Array>()) {
        EXPECT_EQ(tomlValue(tagged(*c.find("input"))), text(c, "output")) << toJson(*c.find("input"), -1);
    }
}

TEST(TomlWriterVectors, Keys) {
    const Value cases = readJson(kVectors / "key_vectors.json");
    for (const Value& c : *cases.get<Array>()) {
        const std::string key = *tagged(*c.find("key")).get<std::string>();
        if (!key.empty() && key.back() == '\n' && text(c, "output") == key) {
            // MDATK leaves "abc\n" bare (Python's "$" matches before a final newline), which no TOML
            // reader accepts; MDBridge quotes it (an accepted deviation)
            EXPECT_EQ(tomlKey(key), tomlQuote(key));
            continue;
        }
        EXPECT_EQ(tomlKey(key), text(c, "output")) << toJson(*c.find("key"), -1);
    }
}

TEST(TomlWriterVectors, Dumps) {
    const Value cases = readJson(kVectors / "dumps_vectors.json");
    for (const Value& c : *cases.get<Array>()) {
        const Value docs = *c.find("docs");
        auto doc = [&docs](const std::vector<std::string>& path, const std::vector<const Value*>&) {
            std::string key;
            for (const std::string& part : path) {
                key += (key.empty() ? "" : "/") + part;
            }
            const Value* text = docs.find(key);
            return text ? std::optional<std::string>(*text->get<std::string>()) : std::nullopt;
        };
        EXPECT_EQ(dumpToml(tagged(*c.find("data")), doc, strings(*c.find("header"))), text(c, "output"))
            << text(c, "name");
    }
}
