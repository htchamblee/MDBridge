#include "mdbridge/core/project_file.hpp"

#include <algorithm>
#include <array>
#include <fstream>
#include <sstream>
#include <system_error>
#include <tuple>

#include <toml++/toml.hpp>

#include "mdbridge/core/json.hpp"
#include "mdbridge/core/paths.hpp"

namespace mdbridge::core {

namespace {

namespace fs = std::filesystem;

// entries of a project file that describe the file, not the project (projectfile.read pops them)
constexpr std::array<std::string_view, 6> kMetaKeys = {"format",   "format_version", "mdatk_version",
                                                       "saved",    "log",            "log_file"};

std::string display(const fs::path& path) { return pathToUtf8(path); }

std::string readText(const fs::path& path) {
    std::error_code ec;
    if (!fs::is_regular_file(path, ec)) {
        throw ProjectError(display(path) + (fs::exists(path, ec) ? ": not a file" : ": no such file"));
    }
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw ProjectError(display(path) + ": cannot open the file");
    }
    std::ostringstream text;
    text << in.rdbuf();
    if (in.bad()) {
        throw ProjectError(display(path) + ": cannot read the file");
    }
    return std::move(text).str();
}

bool endsWithIgnoringCase(std::string_view text, std::string_view suffix) {
    if (text.size() < suffix.size()) {
        return false;
    }
    return std::equal(suffix.begin(), suffix.end(), text.end() - static_cast<std::ptrdiff_t>(suffix.size()),
                      [](char a, char b) {
                          auto lower = [](char c) { return c >= 'A' && c <= 'Z' ? static_cast<char>(c - 'A' + 'a') : c; };
                          return lower(a) == lower(b);
                      });
}

// projectfile._reject_json
void rejectJson(const fs::path& path) {
    if (endsWithIgnoringCase(display(path), ".json")) {
        throw ProjectError(display(path) + ": JSON presets are no longer supported; MDBridge uses TOML "
                                           "project files (" + std::string(kProjectSuffix) + ")");
    }
}

DateTime convertDate(const toml::date& d, DateTime::Kind kind) {
    DateTime out;
    out.kind = kind;
    out.year = d.year;
    out.month = d.month;
    out.day = d.day;
    return out;
}

void setTime(DateTime& out, const toml::time& t) {
    out.hour = t.hour;
    out.minute = t.minute;
    out.second = t.second;
    out.microsecond = static_cast<int>(t.nanosecond / 1000);  // tomllib truncates to microseconds
}

Value convert(const toml::node& node);

Value convertTable(const toml::table& table) {
    // toml++ keeps keys sorted; Python's tomllib keeps them in the order of the text
    std::vector<std::pair<const toml::key*, const toml::node*>> entries;
    entries.reserve(table.size());
    for (auto&& [key, node] : table) {
        entries.emplace_back(&key, &node);
    }
    std::stable_sort(entries.begin(), entries.end(), [](const auto& a, const auto& b) {
        const auto& pa = a.first->source().begin;
        const auto& pb = b.first->source().begin;
        return std::tie(pa.line, pa.column) < std::tie(pb.line, pb.column);
    });
    Table out;
    out.reserve(entries.size());
    for (const auto& [key, node] : entries) {
        out.emplace_back(std::string(key->str()), convert(*node));
    }
    return out;
}

Value convert(const toml::node& node) {
    switch (node.type()) {
    case toml::node_type::table:
        return convertTable(*node.as_table());
    case toml::node_type::array: {
        Array out;
        for (const toml::node& element : *node.as_array()) {
            out.push_back(convert(element));
        }
        return out;
    }
    case toml::node_type::string:
        return node.as_string()->get();
    case toml::node_type::integer:
        return node.as_integer()->get();
    case toml::node_type::floating_point:
        return node.as_floating_point()->get();
    case toml::node_type::boolean:
        return node.as_boolean()->get();
    case toml::node_type::date:
        return convertDate(node.as_date()->get(), DateTime::Kind::LocalDate);
    case toml::node_type::time: {
        DateTime out;
        out.kind = DateTime::Kind::LocalTime;
        setTime(out, node.as_time()->get());
        return out;
    }
    case toml::node_type::date_time: {
        const toml::date_time& dt = node.as_date_time()->get();
        DateTime out = convertDate(dt.date, dt.offset ? DateTime::Kind::OffsetDateTime : DateTime::Kind::LocalDateTime);
        setTime(out, dt.time);
        if (dt.offset) {
            out.offsetMinutes = dt.offset->minutes;
        }
        return out;
    }
    case toml::node_type::none:
        break;
    }
    return {};
}

// The first `count` characters (code points) of UTF-8 text, like Python's fh.read(count).
std::string_view firstCharacters(std::string_view text, std::size_t count) {
    std::size_t seen = 0;
    for (std::size_t i = 0; i < text.size(); ++i) {
        if ((static_cast<unsigned char>(text[i]) & 0xC0) != 0x80 && seen++ == count) {
            return text.substr(0, i);
        }
    }
    return text;
}

}  // namespace

Value parseToml(std::string_view text, std::string_view source) {
    try {
        return convertTable(toml::parse(text, source));
    } catch (const toml::parse_error& error) {
        const auto& where = error.source().begin;
        throw ProjectError(std::string(source) + ": not valid TOML: " + std::string(error.description()) +
                           " (at line " + std::to_string(where.line) + ", column " +
                           std::to_string(where.column) + ")");
    }
}

ProjectFile readProject(const fs::path& path) {
    rejectJson(path);
    Value data = parseToml(readText(path), display(path));

    ProjectFile file;
    file.path = path;
    Value fmt;
    if (const Value* f = data.find("format")) {
        fmt = *f;
    }
    if (fmt == Value(std::string(kParamsLogFormat))) {  // replay a previous run
        file.kind = FileKind::ParamsLog;
        Table header;
        Value parameters = Table{};
        for (auto& [key, value] : *data.get<Table>()) {
            if (key == "parameters") {
                if (value.isTable()) {
                    parameters = std::move(value);
                }
            } else {
                header.emplace_back(key, std::move(value));
            }
        }
        file.header = std::move(header);
        data = std::move(parameters);
        fmt = std::string(kProjectFormat);
    }
    if (!fmt.isNull() && fmt != Value(std::string(kProjectFormat))) {
        const std::string* text = fmt.get<std::string>();
        throw ProjectError(display(path) + ": not an MDATK project (format=" +
                           (text ? pythonStringRepr(*text) : toJson(fmt, -1)) + ")");
    }
    Table& header = *file.header.get<Table>();
    Table settings;
    for (auto& [key, value] : *data.get<Table>()) {
        const bool meta = std::find(kMetaKeys.begin(), kMetaKeys.end(), key) != kMetaKeys.end();
        if (!meta) {
            settings.emplace_back(key, std::move(value));
        } else if (file.kind == FileKind::Project && key != "log") {
            header.emplace_back(key, std::move(value));
        }
    }
    file.settings = std::move(settings);
    return file;
}

bool isParamsLog(const fs::path& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        return false;
    }
    std::string head(4 * 4096, '\0');  // 4096 characters are at most 4 * 4096 bytes of UTF-8
    in.read(head.data(), static_cast<std::streamsize>(head.size()));
    head.resize(static_cast<std::size_t>(in.gcount()));
    const std::string marker = "format = \"" + std::string(kParamsLogFormat) + "\"";
    return firstCharacters(head, 4096).find(marker) != std::string_view::npos;
}

std::string logFileName(const fs::path& project) {
    std::error_code ec;
    if (!fs::is_regular_file(project, ec)) {
        return {};
    }
    std::ifstream in(project, std::ios::binary);
    std::string line;
    for (int n = 0; n < 60 && std::getline(in, line); ++n) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (line.rfind("log_file = ", 0) == 0) {
            try {
                const Value entry = parseToml(line, display(project));
                const Value* name = entry.find("log_file");
                return name ? name->stringOr() : std::string();
            } catch (const ProjectError&) {
                return {};
            }
        }
    }
    return {};
}

fs::path logPath(const fs::path& project) {
    const std::string name = logFileName(project);
    if (name.empty()) {
        return project;
    }
    std::error_code ec;
    fs::path absolute = fs::absolute(project, ec);
    return (ec ? project : absolute).parent_path() / pathFromUtf8(name);
}

Array readLog(const fs::path& project, std::string* warning) {
    const fs::path target = logPath(project);
    try {
        Value data = parseToml(readText(target), display(target));
        if (Value* log = data.find("log"); log && log->isArray()) {
            return std::move(*log->get<Array>());
        }
        return {};
    } catch (const ProjectError& error) {
        if (warning) {
            *warning = "could not read the step log " + display(target) + ": " + error.what();
        }
        return {};
    }
}

}  // namespace mdbridge::core
