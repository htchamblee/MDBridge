// MDBridge's Windows rules for project paths (ntpath.hpp). MDATK never ran on Windows, so there are
// no reference vectors: the expectations below are the rules written out. Pure string functions,
// so these run on every platform.

#include <optional>
#include <ostream>
#include <string>

#include <gtest/gtest.h>

#include "mdbridge/core/ntpath.hpp"

using namespace mdbridge::core;

namespace {

struct Split {
    std::string drive, root, tail;
    bool operator==(const Split&) const = default;
};

Split split(std::string_view path) {
    const ntpath::Parts parts = ntpath::splitroot(path);
    return {std::string(parts.drive), std::string(parts.root), std::string(parts.tail)};
}

std::ostream& operator<<(std::ostream& out, const Split& s) {
    return out << "('" << s.drive << "', '" << s.root << "', '" << s.tail << "')";
}

const std::optional<std::string> kHome = std::string("C:\\Users\\Ann");

std::string resolve(std::string_view path, std::string_view base) { return ntpath::resolve(path, base, kHome); }

}  // namespace

TEST(NtPath, SplitRoot) {
    EXPECT_EQ(split("C:\\a\\b"), (Split{"C:", "\\", "a\\b"}));
    EXPECT_EQ(split("C:/a/b"), (Split{"C:", "/", "a/b"}));
    EXPECT_EQ(split("c:a"), (Split{"c:", "", "a"}));
    EXPECT_EQ(split("C:"), (Split{"C:", "", ""}));
    EXPECT_EQ(split("\\a"), (Split{"", "\\", "a"}));
    EXPECT_EQ(split("/a/b"), (Split{"", "/", "a/b"}));
    EXPECT_EQ(split("\\\\srv\\share\\a"), (Split{"\\\\srv\\share", "\\", "a"}));
    EXPECT_EQ(split("//srv/share/a/b"), (Split{"//srv/share", "/", "a/b"}));
    EXPECT_EQ(split("\\/srv\\share/a"), (Split{"\\/srv\\share", "/", "a"}));
    EXPECT_EQ(split("//srv/share"), (Split{"//srv/share", "", ""}));
    EXPECT_EQ(split("//srv"), (Split{"//srv", "", ""}));
    EXPECT_EQ(split("\\\\?\\UNC\\srv\\share\\a"), (Split{"\\\\?\\UNC\\srv\\share", "\\", "a"}));
    EXPECT_EQ(split("a/b"), (Split{"", "", "a/b"}));
    EXPECT_EQ(split("1:a"), (Split{"", "", "1:a"}));  // only an ASCII letter makes a drive
    EXPECT_EQ(split(""), (Split{"", "", ""}));
}

TEST(NtPath, IsAbs) {
    for (const char* path : {"C:\\x", "C:/x", "c:/", "\\x", "/x", "\\\\srv\\share\\x", "//srv/share/x",
                             "//srv/share", "\\\\?\\C:\\x", "//./pipe/x"}) {
        EXPECT_TRUE(ntpath::isabs(path)) << path;
    }
    for (const char* path : {"", "x", "x\\y", "C:", "C:x", ".\\x", "..", "~"}) {
        EXPECT_FALSE(ntpath::isabs(path)) << path;
    }
}

TEST(NtPath, Normpath) {
    EXPECT_EQ(ntpath::normpath("C:\\a\\.\\b\\..\\c\\"), "C:/a/c");
    EXPECT_EQ(ntpath::normpath("C:/a//b/"), "C:/a/b");
    EXPECT_EQ(ntpath::normpath("C:\\..\\..\\x"), "C:/x");  // never above the drive's root
    EXPECT_EQ(ntpath::normpath("C:\\"), "C:/");
    EXPECT_EQ(ntpath::normpath("C:/a/.."), "C:/");
    EXPECT_EQ(ntpath::normpath("c:\\Data\\X.pdb"), "c:/Data/X.pdb");  // case kept
    EXPECT_EQ(ntpath::normpath("\\\\srv\\share\\a\\..\\..\\b"), "//srv/share/b");  // nor the share
    EXPECT_EQ(ntpath::normpath("//srv/share/"), "//srv/share/");
    EXPECT_EQ(ntpath::normpath("//srv/share"), "//srv/share");
    EXPECT_EQ(ntpath::normpath("\\a\\..\\..\\b"), "/b");
    EXPECT_EQ(ntpath::normpath("a\\..\\..\\b"), "../b");  // relative: kept
    EXPECT_EQ(ntpath::normpath("C:a\\..\\..\\b"), "C:../b");
    EXPECT_EQ(ntpath::normpath("a/./b\\c"), "a/b/c");
    EXPECT_EQ(ntpath::normpath(""), ".");
    EXPECT_EQ(ntpath::normpath("a\\.."), ".");
    EXPECT_EQ(ntpath::normpath("..."), "...");  // not a parent reference
    EXPECT_EQ(ntpath::normpath("\\\\?\\C:\\a\\..\\b"), "\\\\?\\C:\\a\\..\\b");  // device paths as given
    EXPECT_EQ(ntpath::normpath("//./pipe/x"), "//./pipe/x");
}

TEST(NtPath, Join) {
    EXPECT_EQ(ntpath::join("C:\\a", "b"), "C:\\a\\b");
    EXPECT_EQ(ntpath::join("C:/a/", "b"), "C:/a/b");
    EXPECT_EQ(ntpath::join("C:/a", "D:\\b"), "D:\\b");
    EXPECT_EQ(ntpath::join("C:/a", "\\b"), "C:\\b");
    EXPECT_EQ(ntpath::join("//srv/share/a", "/b"), "//srv/share/b");
    EXPECT_EQ(ntpath::join("C:/a", "C:b"), "C:/a\\b");
    EXPECT_EQ(ntpath::join("C:/a", "c:b"), "c:/a\\b");
    EXPECT_EQ(ntpath::join("C:/a", "D:b"), "D:b");
    EXPECT_EQ(ntpath::join("//srv/share", "a"), "//srv/share\\a");
    EXPECT_EQ(ntpath::join("", "a"), "a");
    EXPECT_EQ(ntpath::join("C:", "a"), "C:a");
}

TEST(NtPath, Dirname) {
    EXPECT_EQ(ntpath::dirname("C:/a/b/p.mdatk.toml"), "C:/a/b");
    EXPECT_EQ(ntpath::dirname("C:\\a\\p.mdatk.toml"), "C:\\a");
    EXPECT_EQ(ntpath::dirname("C:/p.mdatk.toml"), "C:/");
    EXPECT_EQ(ntpath::dirname("C:/a//b"), "C:/a");
    EXPECT_EQ(ntpath::dirname("//srv/share/p.toml"), "//srv/share/");
    EXPECT_EQ(ntpath::dirname("//srv/share/d/p.toml"), "//srv/share/d");
    EXPECT_EQ(ntpath::dirname("p.toml"), "");
}

TEST(NtPath, ProfileHome) {
    const std::optional<std::string> none;
    EXPECT_EQ(ntpath::profileHome(std::string("C:\\Users\\Ann"), std::string("D:"), std::string("\\x")),
              "C:\\Users\\Ann");
    EXPECT_EQ(ntpath::profileHome(none, std::string("D:"), std::string("\\Users\\Ann")), "D:\\Users\\Ann");
    EXPECT_EQ(ntpath::profileHome(none, none, std::string("\\Users\\Ann")), "\\Users\\Ann");
    EXPECT_EQ(ntpath::profileHome(none, std::string("D:"), none), std::nullopt);
}

TEST(NtPath, Expanduser) {
    EXPECT_EQ(ntpath::expanduser("~", kHome), "C:\\Users\\Ann");
    EXPECT_EQ(ntpath::expanduser("~\\x", kHome), "C:\\Users\\Ann\\x");
    EXPECT_EQ(ntpath::expanduser("~/x", kHome), "C:\\Users\\Ann/x");
    EXPECT_EQ(ntpath::expanduser("~ann\\x", kHome), "~ann\\x");  // "~user" stays
    EXPECT_EQ(ntpath::expanduser("x\\~", kHome), "x\\~");
    EXPECT_EQ(ntpath::expanduser("~\\x", std::nullopt), "~\\x");
}

// resolvePath on Windows: the project's folder "C:/proj/p"
TEST(NtPath, ResolveAgainstTheProjectFolder) {
    const std::string base = "C:/proj/p";
    // relative, either separator
    EXPECT_EQ(resolve("top.parm7", base), "C:/proj/p/top.parm7");
    EXPECT_EQ(resolve("sub\\traj.nc", base), "C:/proj/p/sub/traj.nc");
    EXPECT_EQ(resolve("sub/x\\y.nc", base), "C:/proj/p/sub/x/y.nc");
    EXPECT_EQ(resolve(".\\a\\.\\b\\", base), "C:/proj/p/a/b");
    EXPECT_EQ(resolve("..\\data\\a.pdb", base), "C:/proj/data/a.pdb");
    EXPECT_EQ(resolve("..\\..\\..\\..\\a.pdb", base), "C:/a.pdb");  // never above the drive's root
    EXPECT_EQ(resolve(".", base), "C:/proj/p");
    EXPECT_EQ(resolve("", base), "");
    // drive paths: absolute, normalised, generic form, letter as given
    EXPECT_EQ(resolve("D:\\data\\a.pdb", base), "D:/data/a.pdb");
    EXPECT_EQ(resolve("d:/data/../a.pdb", base), "d:/a.pdb");
    EXPECT_EQ(resolve("C:\\", base), "C:/");
    // UNC
    EXPECT_EQ(resolve("\\\\srv\\share\\md\\a.nc", base), "//srv/share/md/a.nc");
    EXPECT_EQ(resolve("//srv/share/md/../b.nc", base), "//srv/share/b.nc");
    EXPECT_EQ(resolve("\\\\srv\\share\\..\\x", base), "//srv/share/x");
    // root-relative: the drive of the base
    EXPECT_EQ(resolve("\\data\\a.pdb", base), "C:/data/a.pdb");
    EXPECT_EQ(resolve("/data/a.pdb", base), "C:/data/a.pdb");
    EXPECT_EQ(resolve("/data/a.pdb", "//srv/share/proj"), "//srv/share/data/a.pdb");
    // drive-relative: the base when on that drive (letters without case), else the drive's root
    EXPECT_EQ(resolve("C:sub\\a.pdb", base), "C:/proj/p/sub/a.pdb");
    EXPECT_EQ(resolve("c:sub\\a.pdb", base), "c:/proj/p/sub/a.pdb");
    EXPECT_EQ(resolve("C:..\\a.pdb", base), "C:/proj/a.pdb");
    EXPECT_EQ(resolve("D:sub\\a.pdb", base), "D:/sub/a.pdb");
    EXPECT_EQ(resolve("D:..\\a.pdb", base), "D:/a.pdb");
    EXPECT_EQ(resolve("D:", base), "D:/");
    EXPECT_EQ(resolve("C:a.pdb", "//srv/share/proj"), "C:/a.pdb");
    // "~": the profile folder; "~user" is a relative name
    EXPECT_EQ(resolve("~", base), "C:/Users/Ann");
    EXPECT_EQ(resolve("~\\md\\..\\a.nc", base), "C:/Users/Ann/a.nc");
    EXPECT_EQ(resolve("~/a.nc", base), "C:/Users/Ann/a.nc");
    EXPECT_EQ(resolve("~ann\\a.nc", base), "C:/proj/p/~ann/a.nc");
    EXPECT_EQ(ntpath::resolve("~\\a.nc", base, std::nullopt), "C:/proj/p/~/a.nc");
    // case and other text kept as stored
    EXPECT_EQ(resolve("Data\\A.PDB", base), "C:/proj/p/Data/A.PDB");
    EXPECT_EQ(resolve(" a ", base), "C:/proj/p/ a ");
    EXPECT_EQ(resolve("caf\xC3\xA9.pdb", base), "C:/proj/p/caf\xC3\xA9.pdb");
    // device paths as given
    EXPECT_EQ(resolve("\\\\?\\C:\\long\\a.pdb", base), "\\\\?\\C:\\long\\a.pdb");
}

TEST(NtPath, AbspathAgainstTheCurrentDirectory) {
    EXPECT_EQ(ntpath::abspath("MDATK_output", "C:\\Users\\Ann\\work"), "C:/Users/Ann/work/MDATK_output");
    EXPECT_EQ(ntpath::abspath("C:\\Users\\Ann\\work\\p.mdatk.toml", "D:\\x"), "C:/Users/Ann/work/p.mdatk.toml");
    EXPECT_EQ(ntpath::abspath("/abs/out/", "C:\\Users\\Ann"), "C:/abs/out");
    EXPECT_EQ(ntpath::abspath("rel//out", "C:\\Users\\Ann\\work/x"), "C:/Users/Ann/work/x/rel/out");
    EXPECT_EQ(ntpath::abspath("~", "C:\\w"), "C:/w/~");  // no "~" expansion here
    // relative base (not used by the project code): stays relative
    EXPECT_EQ(ntpath::abspath("..\\a", "rel\\sub"), "rel/a");
}

TEST(NtPath, Generic) {
    EXPECT_EQ(ntpath::generic("C:\\a/b\\c"), "C:/a/b/c");
    EXPECT_EQ(ntpath::generic("\\\\srv\\share"), "//srv/share");
}
