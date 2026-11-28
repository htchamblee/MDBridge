// mdbridge_lua_lint FILE_OR_FOLDER...: the determinism lint (lint.hpp) over Lua files; folders are
// searched recursively for *.lua. Prints each finding as "file:line:column: rule: message" and
// exits with 0 if there is none, 1 if there are findings or a file cannot be read, 2 on bad usage.

#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <system_error>
#include <vector>

#include "mdbridge/workflow/lua/lint.hpp"

namespace fs = std::filesystem;

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: mdbridge_lua_lint FILE_OR_FOLDER...\n");
        return 2;
    }
    std::vector<fs::path> files;
    for (int i = 1; i < argc; ++i) {
        const fs::path path(argv[i]);
        std::error_code error;
        if (fs::is_directory(path, error)) {
            for (fs::recursive_directory_iterator it(path, error), end; !error && it != end; it.increment(error)) {
                if (it->is_regular_file(error) && it->path().extension() == ".lua") {
                    files.push_back(it->path());
                }
            }
        } else {
            files.push_back(path);
        }
        if (error) {
            std::fprintf(stderr, "error: %s: %s\n", path.string().c_str(), error.message().c_str());
            return 1;
        }
    }
    std::sort(files.begin(), files.end());
    std::size_t total = 0;
    bool failed = false;
    for (const fs::path& file : files) {
        std::ifstream in(file, std::ios::binary);
        if (!in) {
            std::fprintf(stderr, "error: cannot read %s\n", file.string().c_str());
            failed = true;
            continue;
        }
        std::ostringstream text;
        text << in.rdbuf();
        for (const auto& finding : mdbridge::workflow::lua::lint(text.str())) {
            std::printf("%s\n", mdbridge::workflow::lua::formatFinding(file.generic_string(), finding).c_str());
            ++total;
        }
    }
    std::printf("%zu Lua file(s), %zu finding(s)\n", files.size(), total);
    return failed || total > 0 ? 1 : 0;
}
