// mmCIF documents (the mmCIF part of mdatk/dictionaries.py) against what MDATK gives for the
// cases of tests/parity/make_cif_cases.py (data/io/cif, see data/README.md): tokens, parsing,
// writing, lookups, cif_value, is_dictionary_file and find_dictionary; MDATK's own parser test
// (tests/test_dictionaries.py) is ported at the end. core::pythonLower, which the lookups rely on,
// is checked here too.

#include <chrono>
#include <optional>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/hash.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/text.hpp"
#include "mdbridge/io/cif.hpp"
#include "test_support.hpp"

using namespace mdbridge;
using namespace mdbridge::test;
using core::Array;
using core::Table;
using core::Value;
using io::CifBlock;
using io::CifDocument;
using io::CifLoop;
using io::CifToken;
using io::DictionaryError;

namespace {

const Value& cases() {
    static const Value data = readJson(kData / "io" / "cif" / "cif_cases.json");
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

Value stringArray(const std::vector<std::string>& values) {
    Array out;
    for (const std::string& v : values) {
        out.emplace_back(v);
    }
    return out;
}

Value plain(const CifDocument& doc) {
    Array blocks;
    for (const CifBlock& b : doc.blocks) {
        Array loops;
        for (const CifLoop& item : b.items) {
            Array rows;
            for (const auto& row : item.rows) {
                rows.push_back(stringArray(row));
            }
            loops.emplace_back(Table{{"tags", stringArray(item.tags)},
                                     {"rows", std::move(rows)},
                                     {"pairs", item.pairs}});
        }
        blocks.emplace_back(Table{{"name", b.name}, {"items", std::move(loops)}});
    }
    return Table{{"header", stringArray(doc.header)}, {"blocks", std::move(blocks)}};
}

CifDocument fromPlain(const Value& data) {
    CifDocument doc;
    doc.header = strings(*data.find("header"));
    for (const Value& b : items(data, "blocks")) {
        CifBlock block;
        block.name = text(b, "name");
        for (const Value& i : items(b, "items")) {
            CifLoop loop;
            loop.tags = strings(*i.find("tags"));
            for (const Value& row : items(i, "rows")) {
                loop.rows.push_back(strings(row));
            }
            loop.pairs = *i.find("pairs")->get<bool>();
            block.items.push_back(std::move(loop));
        }
        doc.blocks.push_back(std::move(block));
    }
    return doc;
}

std::string caseText(const Value& c) {
    if (const Value* hex = c.find("hex")) {
        return fromHex(hex->stringOr());
    }
    return text(c, "text");
}

std::string_view kindName(CifToken::Kind kind) {
    switch (kind) {
    case CifToken::Kind::Value:
        return "value";
    case CifToken::Kind::Comment:
        return "comment";
    case CifToken::Kind::Data:
        return "data";
    case CifToken::Kind::Loop:
        return "loop";
    case CifToken::Kind::Other:
        return "other";
    case CifToken::Kind::Tag:
        return "tag";
    }
    return "?";
}

// What reparsing gave: {"ok": doc} or {"error": message}.
void checkAttempt(const Value& expected, const std::string& written, const std::string& what) {
    try {
        const CifDocument doc = io::parseCif(written);
        ASSERT_NE(expected.find("ok"), nullptr) << what << ": expected " << text(expected, "error");
        EXPECT_EQ(plain(doc), *expected.find("ok")) << what;
    } catch (const DictionaryError& e) {
        EXPECT_EQ(e.what(), text(expected, "error")) << what;
    }
}

template <class T, class U>
std::ptrdiff_t indexOf(const std::vector<T>& list, const U* found) {
    for (std::size_t k = 0; k < list.size(); ++k) {
        if (&list[k] == found) {
            return static_cast<std::ptrdiff_t>(k);
        }
    }
    return -1;
}

std::size_t number(const Value& value, std::string_view key) {
    return static_cast<std::size_t>(*value.find(key)->get<std::int64_t>());
}

void checkLookups(const CifDocument& doc, const Value& c, const std::string& what) {
    for (const Value& l : items(c, "lookups")) {
        const auto index = static_cast<std::ptrdiff_t>(*l.find("index")->get<std::int64_t>());
        if (const Value* name = l.find("block")) {
            EXPECT_EQ(indexOf(doc.blocks, doc.block(name->stringOr())), index) << what << " block "
                                                                              << name->stringOr();
            continue;
        }
        const CifBlock& block = doc.blocks[number(l, "in")];
        if (const Value* category = l.find("get")) {
            EXPECT_EQ(indexOf(block.items, block.get(category->stringOr())), index)
                << what << " get " << category->stringOr();
            continue;
        }
        const CifLoop& item = block.items[number(l, "item")];
        const std::string query = text(l, "col");
        EXPECT_EQ(item.col(query), index) << what << " col " << query;
        if (const Value* column = l.find("column")) {
            EXPECT_EQ(item.column(query), strings(*column)) << what << " column " << query;
        } else {
            EXPECT_THROW(item.column(query), std::out_of_range) << what << " column " << query;
        }
    }
}

}  // namespace

TEST(Cif, Tokens) {
    for (const Value& c : items(cases(), "parse")) {
        const std::string what = text(c, "name");
        try {
            const std::vector<CifToken> tokens = io::cifTokens(caseText(c));
            ASSERT_NE(c.find("tokens"), nullptr)
                << what << ": expected " << text(c, "tokens_error");
            const Array& expected = items(c, "tokens");
            ASSERT_EQ(tokens.size(), expected.size()) << what;
            for (std::size_t k = 0; k < tokens.size(); ++k) {
                const Array& e = *expected[k].get<Array>();
                EXPECT_EQ(kindName(tokens[k].kind), e[0].stringOr()) << what << " token " << k;
                EXPECT_EQ(tokens[k].text, e[1].stringOr()) << what << " token " << k;
            }
        } catch (const DictionaryError& e) {
            EXPECT_EQ(e.what(), text(c, "tokens_error")) << what;
        }
    }
}

TEST(Cif, ParseWriteAndLookups) {
    for (const Value& c : items(cases(), "parse")) {
        const std::string what = text(c, "name");
        try {
            const CifDocument doc = io::parseCif(caseText(c));
            ASSERT_NE(c.find("doc"), nullptr) << what << ": expected " << text(c, "error");
            EXPECT_EQ(plain(doc), *c.find("doc")) << what;
            const std::string written = io::writeCif(doc);
            EXPECT_EQ(written, text(c, "written")) << what;
            checkAttempt(*c.find("reparsed"), written, what + " (reparsed)");
            checkLookups(doc, c, what);
        } catch (const DictionaryError& e) {
            EXPECT_EQ(e.what(), text(c, "error")) << what;
            EXPECT_STREQ(e.pythonName(), "DictionaryError");
        }
    }
}

TEST(Cif, Write) {
    int k = 0;
    for (const Value& c : items(cases(), "write")) {
        const std::string what = "write case " + std::to_string(k++);
        const CifDocument doc = fromPlain(*c.find("doc"));
        const std::string written = io::writeCif(doc);
        EXPECT_EQ(written, text(c, "written")) << what;
        checkAttempt(*c.find("reparsed"), written, what);
    }
}

TEST(Cif, Value) {
    for (const Value& c : items(cases(), "cif_value")) {
        const Array& pair = *c.get<Array>();
        EXPECT_EQ(io::cifValue(pair[0].stringOr()), pair[1].stringOr())
            << core::pythonStringRepr(pair[0].stringOr());
    }
}

TEST(Cif, PythonLower) {
    const Value& lower = *cases().find("lower");
    for (const Value& c : items(lower, "strings")) {
        const Array& pair = *c.get<Array>();
        EXPECT_EQ(core::pythonLower(pair[0].stringOr()), pair[1].stringOr())
            << core::pythonStringRepr(pair[0].stringOr());
    }
    std::string single;
    std::string context;
    const std::string sigma = "\xce\xa3";
    for (char32_t code = 0; code < 0x110000; ++code) {
        if (code >= 0xD800 && code <= 0xDFFF) {
            continue;
        }
        std::string c;
        core::appendUtf8(c, code);
        single += (code ? "\n" : "") + core::pythonLower(c);
        context += (code ? "\n" : "") + core::pythonLower(c + sigma) + "|" +
                   core::pythonLower("A" + sigma + c) + "|" + core::pythonLower("A" + c + sigma);
    }
    EXPECT_EQ(core::sha256Hex(single), text(lower, "single_sha256"));
    EXPECT_EQ(core::sha256Hex(context), text(lower, "context_sha256"));
}

TEST(Cif, IsDictionaryFile) {
    TempDir dir;
    int k = 0;
    for (const Value& c : items(cases(), "dictionary")) {
        const bool expected = *c.find("result")->get<bool>();
        const fs::path path = dir / ("d" + std::to_string(k++) + ".cif");
        if (c.find("folder")) {
            fs::create_directory(path);
        } else if (const Value* parts = c.find("parts")) {
            std::string bytes;
            for (const Value& part : *parts->get<Array>()) {
                const Array& p = *part.get<Array>();
                const std::string unit = fromHex(p[0].stringOr());
                const auto count = *p[1].get<std::int64_t>();
                for (std::int64_t n = 0; n < count; ++n) {
                    bytes += unit;
                }
            }
            writeFile(path, bytes);
        }
        EXPECT_EQ(io::isDictionaryFile(core::pathToUtf8(path)), expected) << "case " << k - 1;
    }
}

#if !defined(_WIN32)
// The scenarios use POSIX file names ("*" in a folder name); parity is defined on Linux.
TEST(Cif, FindDictionary) {
    const fs::path oldCwd = fs::current_path();
    for (const Value& scenario : items(cases(), "find")) {
        TempDir dir;
        const std::string root = core::pathToUtf8(dir.path());
        for (const Value& d : items(scenario, "dirs")) {
            fs::create_directories(dir.path() / core::pathFromUtf8(d.stringOr()));
        }
        const auto now = fs::file_time_type::clock::now();
        for (const Value& f : items(scenario, "files")) {
            const Array& file = *f.get<Array>();
            const fs::path path = dir.path() / core::pathFromUtf8(file[0].stringOr());
            writeFile(path, file[1].stringOr());
            fs::last_write_time(path, now - std::chrono::seconds(*file[2].get<std::int64_t>()));
        }
        const auto withRoot = [&](const Value& list) {
            std::vector<std::string> out;
            for (std::string s : strings(list)) {
                for (auto at = s.find("{root}"); at != std::string::npos; at = s.find("{root}")) {
                    s.replace(at, 6, root);
                }
                out.push_back(s);
            }
            return out;
        };
        const auto withoutRoot = [&](std::string s) {
            for (std::size_t at = s.find(root); at != std::string::npos; at = s.find(root)) {
                s.replace(at, root.size(), "{root}");
            }
            return s;
        };
        fs::current_path(dir.path());
        for (const Value& call : items(scenario, "calls")) {
            const std::vector<std::string> sources = withRoot(*call.find("sources"));
            const Value* exclude = call.find("exclude");
            const std::vector<std::string> skip =
                exclude ? withRoot(*exclude) : std::vector<std::string>();
            const std::string code = text(call, "code");
            const std::string what = toJson(call, -1);
            try {
                const std::string found = io::findDictionary(sources, code, skip);
                EXPECT_EQ(withoutRoot(found), text(call, "result")) << what;
            } catch (const DictionaryError& e) {
                EXPECT_EQ(withoutRoot(e.what()), text(call, "error")) << what;
            } catch (const core::NamedError& e) {
                EXPECT_EQ(e.pythonName() + ": " + e.what(), text(call, "python_error")) << what;
            }
        }
        fs::current_path(oldCwd);
    }
}
#endif

// MDATK's test_cif_parser_quotes_text_fields_pairs_and_round_trip
TEST(Cif, MdatkParserRoundTrip) {
    const std::string text = R"(# header comment
data_comp_list
_chem_comp.id   XYZ
_chem_comp.name 'name with "quotes" inside'
_chem_comp.group non-polymer
data_comp_XYZ
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
XYZ "O5'" O
XYZ 'C1' C # trailing comment
XYZ N1'x N
loop_
_x.text
;
line one
line two
;
)";
    const CifDocument doc = io::parseCif(text);
    EXPECT_EQ(doc.header, std::vector<std::string>{"# header comment"});
    const CifLoop* clist = doc.block("comp_list")->get("_chem_comp");
    ASSERT_NE(clist, nullptr);
    EXPECT_TRUE(clist->pairs);
    EXPECT_EQ(clist->column("name"), std::vector<std::string>{"name with \"quotes\" inside"});
    const CifLoop* atoms = doc.block("comp_XYZ")->get("_chem_comp_atom");
    EXPECT_EQ(atoms->column("atom_id"), (std::vector<std::string>{"O5'", "C1", "N1'x"}));
    EXPECT_EQ(doc.block("comp_XYZ")->get("_x")->rows,
              std::vector<std::vector<std::string>>{{"line one\nline two"}});
    const CifDocument again = io::parseCif(io::writeCif(doc));
    EXPECT_EQ(again.block("comp_XYZ")->get("_chem_comp_atom")->rows, atoms->rows);
    EXPECT_EQ(again.block("comp_list")->get("_chem_comp")->rows, clist->rows);
    EXPECT_EQ(again.block("comp_XYZ")->get("_x")->rows,
              std::vector<std::vector<std::string>>{{"line one\nline two"}});
    EXPECT_NE(io::writeCif(doc).find("\"O5'\""), std::string::npos);
    EXPECT_THROW(io::parseCif("data_a\nloop_\n_a.b\n_a.c\n1 2 3\n"), DictionaryError);
    EXPECT_THROW(io::parseCif("data_a\n_a.b 'open\n"), DictionaryError);
}

#if !defined(_WIN32)
// MDATK's test_find_dictionary (the file contents only need the two loop markers)
TEST(Cif, MdatkFindDictionary) {
    TempDir dir;
    const fs::path work = dir / "work_0A5";
    fs::create_directories(work / "0A5_TMP");
    const std::string dict =
        "data_x\nloop_\n_chem_comp_atom.atom_id\nA\nloop_\n_chem_comp_bond.atom_id_1\nA\n";
    writeFile(work / "0A5_TMP" / "inter.cif", dict);
    writeFile(work / "notes.cif", "data_x\n_a.b 1\n");
    writeFile(work / "0A5.pdb", "END\n");
    writeFile(work / "0A5.cif", dict);
    const std::string w = core::pathToUtf8(work);
    const std::string expected = core::pathToUtf8(work / "0A5.cif");
    EXPECT_EQ(io::findDictionary(std::vector<std::string>{w}, "0A5"), expected);
    EXPECT_EQ(io::findDictionary(std::vector<std::string>{w + "/*.cif"}), expected);
    EXPECT_THROW(io::findDictionary(std::vector<std::string>{w}, "",
                                    std::vector<std::string>{expected}),
                 DictionaryError);
    EXPECT_THROW(io::findDictionary(std::vector<std::string>{core::pathToUtf8(dir / "missing")}),
                 DictionaryError);
}
#endif
