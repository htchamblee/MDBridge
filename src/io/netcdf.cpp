#include "mdbridge/io/netcdf.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <map>
#include <stdexcept>
#include <utility>

#include "binary_file.hpp"
#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/json.hpp"
#include "python_text.hpp"
#include "readers.hpp"
#include "text_file.hpp"

namespace mdbridge::io {

namespace {

using core::ValueError;
using detail::ReaderError;

// ---------------------------------------------------------------- netCDF classic (CDF-1 / CDF-2)

enum NcType : std::uint32_t {
    NcByte = 1,
    NcChar = 2,
    NcShort = 3,
    NcInt = 4,
    NcFloat = 5,
    NcDouble = 6
};

constexpr std::uint32_t kDimensionTag = 0x0A;
constexpr std::uint32_t kVariableTag = 0x0B;
constexpr std::uint32_t kAttributeTag = 0x0C;

std::size_t typeSize(std::uint32_t type) {
    switch (type) {
    case NcByte:
    case NcChar:
        return 1;
    case NcShort:
        return 2;
    case NcInt:
    case NcFloat:
        return 4;
    case NcDouble:
        return 8;
    default:
        throw ReaderError("KeyError", "unknown netCDF type " + std::to_string(type));
    }
}

// One value of type `type` at `p` (big-endian), exactly as a double.
double loadValue(const char* p, std::uint32_t type) {
    switch (type) {
    case NcByte:
    case NcChar:
        return static_cast<double>(static_cast<signed char>(*p));
    case NcShort:
        return static_cast<double>(static_cast<std::int16_t>(detail::loadU16(p, true)));
    case NcInt:
        return static_cast<double>(detail::loadI32(p, true));
    case NcFloat:
        return static_cast<double>(detail::loadF32(p, true));
    default:
        return detail::loadF64(p, true);
    }
}

struct Attribute {
    std::uint32_t type = NcChar;
    std::size_t count = 0;
    std::string bytes;  // the values, big-endian, without padding

    bool isText() const { return type == NcChar; }
    // scipy's value of a char attribute: the bytes without trailing NULs
    std::string text() const {
        std::string out = bytes;
        while (!out.empty() && out.back() == '\0') {
            out.pop_back();
        }
        return out;
    }
    double number(std::size_t i) const {
        return loadValue(bytes.data() + i * typeSize(type), type);
    }
};

using Attributes = std::vector<std::pair<std::string, Attribute>>;

const Attribute* findAttribute(const Attributes& attributes, std::string_view name) {
    for (const auto& [key, value] : attributes) {
        if (key == name) {
            return &value;
        }
    }
    return nullptr;
}

struct Dimension {
    std::string name;
    std::uint64_t length = 0;  // 0: the record dimension
};

struct Variable {
    std::string name;
    std::vector<std::size_t> dims;
    Attributes attributes;
    std::uint32_t type = NcFloat;
    std::uint64_t begin = 0;
    std::uint64_t vsize = 0;  // the header's value
    bool record = false;
    std::vector<std::uint64_t> shape;  // record dimension: the record count
    std::uint64_t itemBytes = 0;       // bytes of one entry of the first axis
    std::uint64_t recordOffset = 0;    // record variables: offset inside a record
};

struct Header {
    int version = 0;
    std::uint64_t records = 0;
    std::vector<Dimension> dims;
    Attributes attributes;
    std::vector<Variable> variables;
    std::uint64_t recordBegin = 0;
    std::uint64_t recordStride = 0;

    const Variable* variable(std::string_view name) const {
        for (const Variable& v : variables) {
            if (v.name == name) {
                return &v;
            }
        }
        return nullptr;
    }
    const Dimension* dimension(std::string_view name) const {
        for (const Dimension& d : dims) {
            if (d.name == name) {
                return &d;
            }
        }
        return nullptr;
    }
};

// The header ends before its last field.
struct TruncatedHeader {};

// The header, read front to back as scipy.io.netcdf_file reads it.
class HeaderReader {
public:
    explicit HeaderReader(detail::BinaryFile& file) : m_file(file) {}

    Header read() {
        Header header;
        const std::string magic = take(4);
        if (magic.substr(0, 3) != "CDF") {
            throw std::invalid_argument("not netcdf");
        }
        header.version = static_cast<signed char>(magic[3]);
        header.records = static_cast<std::uint32_t>(integer());
        // dimensions
        expectTag(kDimensionTag);
        const std::int32_t dimCount = integer();
        for (std::int32_t i = 0; i < dimCount; ++i) {
            Dimension d;
            d.name = name();
            d.length = static_cast<std::uint32_t>(integer());
            header.dims.push_back(std::move(d));
        }
        header.attributes = attributes();
        expectTag(kVariableTag);
        const std::int32_t varCount = integer();
        std::uint64_t recordBytes = 0;
        std::size_t recordVariables = 0;
        for (std::int32_t i = 0; i < varCount; ++i) {
            Variable v;
            v.name = name();
            const std::int32_t rank = integer();
            for (std::int32_t k = 0; k < rank; ++k) {
                const std::int32_t id = integer();
                if (id < 0 || static_cast<std::size_t>(id) >= header.dims.size()) {
                    throw ValueError("variable '" + v.name + "' refers to dimension " +
                                     std::to_string(id) + " of " +
                                     std::to_string(header.dims.size()));
                }
                v.dims.push_back(static_cast<std::size_t>(id));
                v.shape.push_back(header.dims[static_cast<std::size_t>(id)].length);
            }
            v.attributes = attributes();
            v.type = static_cast<std::uint32_t>(integer());
            const std::size_t size = typeSize(v.type);
            v.vsize = static_cast<std::uint32_t>(integer());
            v.begin = header.version == 1 ? static_cast<std::uint32_t>(integer())
                                          : static_cast<std::uint64_t>(integer64());
            v.record = !v.shape.empty() && v.shape[0] == 0;
            std::uint64_t entry = size;
            for (std::size_t k = 1; k < v.shape.size(); ++k) {
                entry *= v.shape[k];
            }
            v.itemBytes = v.shape.empty() ? size : entry;
            if (v.record) {
                v.shape[0] = header.records;
                if (recordVariables == 0) {
                    header.recordBegin = v.begin;
                }
                v.recordOffset = recordBytes;
                recordBytes += entry;
                if (v.type == NcByte || v.type == NcChar || v.type == NcShort) {
                    recordBytes += (4 - entry % 4) % 4;
                }
                ++recordVariables;
            }
            header.variables.push_back(std::move(v));
        }
        if (recordVariables == 1) {
            for (const Variable& v : header.variables) {
                if (v.record) {
                    recordBytes = v.itemBytes;
                }
            }
        }
        header.recordStride = recordBytes;
        return header;
    }

private:
    std::string take(std::size_t n) {
        std::string out = m_file.read(m_pos, n);
        if (out.size() != n) {
            throw TruncatedHeader();
        }
        m_pos += n;
        return out;
    }
    std::int32_t integer() { return detail::loadI32(take(4).data(), true); }
    std::int64_t integer64() { return detail::loadI64(take(8).data(), true); }
    void skipPadding(std::size_t n) { take((4 - n % 4) % 4); }
    void expectTag(std::uint32_t tag) {
        const auto value = static_cast<std::uint32_t>(integer());
        if (value != 0 && value != tag) {
            throw ValueError("Unexpected header.");
        }
    }
    std::string name() {
        const std::int32_t length = integer();
        if (length < 0) {
            throw ValueError("negative length in netCDF header");
        }
        std::string out = take(static_cast<std::size_t>(length));
        skipPadding(out.size());
        return out;
    }
    Attributes attributes() {
        expectTag(kAttributeTag);
        const std::int32_t count = integer();
        Attributes out;
        for (std::int32_t i = 0; i < count; ++i) {
            std::string key = name();
            Attribute a;
            a.type = static_cast<std::uint32_t>(integer());
            const std::size_t size = typeSize(a.type);
            a.count = static_cast<std::uint32_t>(integer());
            a.bytes = take(a.count * size);
            skipPadding(a.bytes.size());
            out.emplace_back(std::move(key), std::move(a));
        }
        return out;
    }

    detail::BinaryFile& m_file;
    std::uint64_t m_pos = 0;
};

// ---------------------------------------------------------------- Amber conventions (reader)

[[noreturn]] void typeError(const std::string& message) {
    throw ReaderError("TypeError", message);
}

// An attribute's text as str: bytes.decode() (UTF-8, strict)
std::string decodedText(const Attribute& attribute) {
    const std::string bytes = attribute.text();
    if (const auto error = detail::utf8Error(bytes)) {
        throw detail::UnicodeDecodeError(*error);
    }
    return bytes;
}

// repr() of a bytes value
std::string bytesRepr(const std::string& bytes) {
    const bool doubleQuote =
        bytes.find('\'') != std::string::npos && bytes.find('"') == std::string::npos;
    const char quote = doubleQuote ? '"' : '\'';
    static const char hex[] = "0123456789abcdef";
    std::string out = "b";
    out += quote;
    for (const char ch : bytes) {
        const auto c = static_cast<unsigned char>(ch);
        if (ch == '\\' || ch == quote) {
            out += '\\';
            out += ch;
        } else if (ch == '\t') {
            out += "\\t";
        } else if (ch == '\n') {
            out += "\\n";
        } else if (ch == '\r') {
            out += "\\r";
        } else if (c < 0x20 || c >= 0x7F) {
            out += "\\x";
            out += hex[c >> 4];
            out += hex[c & 15];
        } else {
            out += ch;
        }
    }
    return out + quote;
}

// str() of a numeric attribute value: a scalar for one element, else NumPy's array text
// (fixed notation, elements padded to a common width).
std::string numericText(const Attribute& attribute) {
    const bool integral = attribute.type != NcFloat && attribute.type != NcDouble;
    std::vector<std::string> items;
    for (std::size_t i = 0; i < attribute.count; ++i) {
        const double x = attribute.number(i);
        if (integral) {
            items.push_back(std::to_string(static_cast<std::int64_t>(x)));
        } else {
            std::string text = core::pythonFloatRepr(x);
            if (text.size() > 2 && text.compare(text.size() - 2, 2, ".0") == 0) {
                text.pop_back();
            }
            items.push_back(std::move(text));
        }
    }
    if (attribute.count == 1) {
        std::string text = items.front();
        if (!integral && text.back() == '.') {
            text += '0';
        }
        return text;
    }
    std::size_t head = 0;
    std::size_t tail = 0;
    for (const std::string& item : items) {
        const std::size_t dot = std::min(item.find('.'), item.size());
        head = std::max(head, dot);
        tail = std::max(tail, item.size() - dot);
    }
    std::string out = "[";
    for (std::size_t i = 0; i < items.size(); ++i) {
        const std::string& item = items[i];
        const std::size_t dot = std::min(item.find('.'), item.size());
        out += (i ? " " : "") + std::string(head - dot, ' ') + item +
               std::string(tail - (item.size() - dot), ' ');
    }
    return out + "]";
}

std::string valueText(const Attribute& attribute) {
    return attribute.isText() ? bytesRepr(attribute.text()) : numericText(attribute);
}

void checkUnits(const std::string& file, const Variable& variable, const std::string& expected) {
    const Attribute* units = findAttribute(variable.attributes, "units");
    if (units == nullptr) {
        throw ReaderError("OSError", file + ": variable '" + variable.name +
                                         "' has no 'units' attribute (expected '" + expected +
                                         "')");
    }
    if (!units->isText()) {
        throw ReaderError("OSError", file + ": the 'units' attribute of variable '" +
                                         variable.name + "' is not text (expected '" + expected +
                                         "', found " + valueText(*units) + ")");
    }
    const std::string found = decodedText(*units);
    if (found != expected) {
        throw ReaderError("NotImplementedError",
                          "NETCDFReader currently assumes that the trajectory was written in "
                          "units of " +
                              found + " instead of " + expected);
    }
}

// Whether AMBER is one of the comma-separated pieces or one of the whitespace-separated words.
bool hasAmberToken(const std::string& conventions) {
    std::size_t start = 0;
    for (std::size_t comma = 0; comma != std::string::npos; start = comma + 1) {
        comma = conventions.find(',', start);
        const std::size_t length = comma == std::string::npos ? std::string::npos : comma - start;
        if (conventions.compare(start, length, "AMBER") == 0) {
            return true;
        }
    }
    for (const detail::U32& word : detail::split(detail::toU32(conventions))) {
        if (word == U"AMBER") {
            return true;
        }
    }
    return false;
}

// A variable's scale factor: applied unless it is within 1e-9 (relative) of 1.
struct Scale {
    double factor = 1.0;
    bool single = false;  // a float32 factor
    bool applies = false;
};

Scale scaleOf(const Attribute& attribute) {
    Scale s;
    s.factor = attribute.number(0);
    s.single = attribute.type == NcFloat;
    const bool close = std::isfinite(s.factor) &&
                       std::fabs(s.factor - 1.0) <= 1e-9 * std::max(std::fabs(s.factor), 1.0);
    s.applies = !close;
    return s;
}

// A value of `type` times the scale factor in NumPy's promoted type, exactly as a double: float32
// for float32 / int8 / int16 data with a float32 factor, else float64.
double scaled(double raw, std::uint32_t type, const Scale& scale) {
    if (!scale.applies) {
        return raw;
    }
    if (scale.single && (type == NcFloat || type == NcShort || type == NcByte || type == NcChar)) {
        const float product = static_cast<float>(raw) * static_cast<float>(scale.factor);
        return static_cast<double>(product);
    }
    return raw * scale.factor;
}

constexpr std::array<std::string_view, 6> kScaledVariables = {
    "time", "cell_lengths", "cell_angles", "coordinates", "velocities", "forces"};

class NetcdfTrajectory : public detail::FileTrajectory {
public:
    NetcdfTrajectory(const std::filesystem::path& path, std::optional<double> dt,
                     std::size_t expectedAtoms)
        : m_file(path) {
        const std::string name = detail::pathText(path);
        readContainer(name);
        checkConventions(name);
        // CDF-1 (32-bit offsets) is read like CDF-2 (RD-N08).
        // TODO(LD-03): for version 1, warn "<file>: NetCDF file uses 32-bit offsets (CDF-1); the
        // Amber convention asks for 64-bit offsets (CDF-2)"
        const Dimension* spatial = m_header.dimension("spatial");
        if (spatial == nullptr) {
            throw ValueError("NCDF trajectory does not contain spatial dimension");
        }
        if (spatial->length != 3) {
            typeError("Incorrect spatial value for NCDF trajectory file");
        }
        // (a missing program / programVersion attribute is only a warning)
        const Dimension* atom = m_header.dimension("atom");
        if (atom == nullptr) {
            throw ValueError("NCDF trajectory " + name + " does not contain atom information");
        }
        m_atoms = static_cast<std::size_t>(atom->length);
        if (expectedAtoms != 0 && expectedAtoms != m_atoms) {
            throw ValueError("The topology and NCDF trajectory files don't have the same number "
                             "of atoms!\nTopology number of atoms " +
                             std::to_string(expectedAtoms) + "\nTrajectory: " + name +
                             " Number of atoms " + std::to_string(m_atoms));
        }
        const Dimension* frame = m_header.dimension("frame");
        const Variable* coordinates = m_header.variable("coordinates");
        if (frame == nullptr || (frame->length == 0 && coordinates == nullptr)) {
            throw ValueError("NCDF trajectory " + name + " does not contain frame information");
        }
        if (frame->length != 0) {
            m_frameCount = static_cast<std::int64_t>(frame->length);
        } else {
            m_frameCount =
                coordinates->shape.empty() ? 0 : static_cast<std::int64_t>(coordinates->shape[0]);
        }
        // (the title attribute is MDAnalysis' remarks, which MDATK does not use)
        m_time = m_header.variable("time");
        if (m_time != nullptr) {
            checkUnits(name, *m_time, "picosecond");
        }
        if (coordinates == nullptr) {
            throw ReaderError("KeyError", "'coordinates'");
        }
        m_coordinates = coordinates;
        checkUnits(name, *m_coordinates, "angstrom");
        readScaleFactors();
        if (const Variable* velocities = m_header.variable("velocities")) {
            checkUnits(name, *velocities, "angstrom/picosecond");
        }
        if (const Variable* forces = m_header.variable("forces")) {
            checkUnits(name, *forces, "kilocalorie/mole/angstrom");
        }
        m_lengths = m_header.variable("cell_lengths");
        if (m_lengths != nullptr) {
            checkUnits(name, *m_lengths, "angstrom");
            m_angles = m_header.variable("cell_angles");
            if (m_angles == nullptr) {
                throw ReaderError("KeyError", "'cell_angles'");
            }
            checkUnits(name, *m_angles, "degree");
        }
        checkCellShape(name, m_lengths, "cell_spatial");
        checkCellShape(name, m_angles, "cell_angular");
        // TODO(LD-03): without a time variable and a given dt, warn once per file (RD-N10)
        // "<file>: no time information; dt set to 1.0 ps"
        m_dt = dt ? *dt : fileDt();
        if (m_frameCount == 0) {
            throw ReaderError("OSError", name + ": trajectory has no frames");
        }
        seek(0);
    }

    bool storedTime(std::int64_t) const override { return m_time != nullptr; }

protected:
    void readFrame(std::int64_t index, core::Frame& frame) override {
        const auto k = static_cast<std::uint64_t>(index);
        const std::vector<double> xyz = row(*m_coordinates, k);
        if (xyz.size() != 3 * m_atoms) {
            throw ValueError("could not broadcast input array from shape (" +
                             std::to_string(xyz.size()) + ",) into shape (" +
                             std::to_string(m_atoms) + ",3)");
        }
        frame.positions.resize(xyz.size());
        for (std::size_t i = 0; i < xyz.size(); ++i) {
            frame.positions[i] = static_cast<float>(xyz[i]);
        }
        if (m_time != nullptr) {
            const std::vector<double> t = row(*m_time, k);
            if (t.size() != 1) {
                typeError("only length-1 arrays can be converted to Python scalars");
            }
            frame.time = t.front();
        } else {
            frame.time = static_cast<double>(index) * m_dt;
        }
        if (m_lengths != nullptr) {
            const std::array<double, 3> lengths = triple(*m_lengths, k);
            const std::array<double, 3> angles = triple(*m_angles, k);
            std::array<float, 6> box{};
            for (std::size_t i = 0; i < 3; ++i) {
                box[i] = static_cast<float>(lengths[i]);
                box[i + 3] = static_cast<float>(angles[i]);
            }
            if (box[0] != 0.0f || box[1] != 0.0f || box[2] != 0.0f) {  // (0, 0, 0): no cell
                frame.box = box;
            }
        }
    }

private:
    void readContainer(const std::string& name) {
        const auto fail = [&](const std::string& message) {
            throw ReaderError("OSError", name + ": " + message);
        };
        const std::uint64_t size = m_file.size();
        if (size == 0) {
            fail("file is empty");
        }
        const std::string magic = m_file.read(0, 4);
        if (magic.substr(0, 3) != "CDF") {
            typeError("Error: " + name + " is not a valid NetCDF 3 file");
        }
        if (magic.size() < 4) {
            fail("truncated header");
        }
        if (magic[3] == 5) {
            fail("unsupported NetCDF variant (CDF-5)");
        }
        if (magic[3] != 1 && magic[3] != 2) {
            fail("unsupported NetCDF variant (version byte " +
                 std::to_string(static_cast<int>(static_cast<unsigned char>(magic[3]))) + ")");
        }
        try {
            m_header = HeaderReader(m_file).read();
        } catch (const TruncatedHeader&) {
            fail("truncated header");
        } catch (const ValueError& error) {
            fail(std::string("corrupt NetCDF header (") + error.what() + ")");
        }
        if (m_header.records == 0xFFFFFFFFu) {
            // streaming (RD-N05): as many records as the file holds
            const std::uint64_t stride = m_header.recordStride;
            m_header.records = stride == 0 || size < m_header.recordBegin
                                   ? 0
                                   : (size - m_header.recordBegin) / stride;
            for (Variable& v : m_header.variables) {
                if (v.record) {
                    v.shape[0] = m_header.records;
                }
            }
        }
        const std::uint64_t needed = m_header.recordBegin + m_header.records * m_header.recordStride;
        if (m_header.records > 0 && needed > size) {
            const std::uint64_t complete = size < m_header.recordBegin
                                               ? 0
                                               : (size - m_header.recordBegin) / m_header.recordStride;
            fail("file ends inside frame " + std::to_string(complete) + " (expected " +
                 std::to_string(needed) + " bytes, found " + std::to_string(size) + ")");
        }
        for (const Variable& v : m_header.variables) {
            if (!v.record) {
                const std::uint64_t entries = v.shape.empty() ? 1 : v.shape[0];
                const std::uint64_t end = v.begin + entries * v.itemBytes;
                if (end > size) {
                    fail("file ends inside variable '" + v.name + "' (expected " +
                         std::to_string(end) + " bytes, found " + std::to_string(size) + ")");
                }
            }
        }
    }

    // A cell variable is (frame, <dimension>) or a constant (<dimension>) (RD-N07).
    void checkCellShape(const std::string& name, const Variable* v, std::string_view dimension) const {
        if (v == nullptr) {
            return;
        }
        const auto isDim = [&](std::size_t id, std::string_view which) {
            return m_header.dims[id].name == which;
        };
        const bool perFrame =
            v->dims.size() == 2 && isDim(v->dims[0], "frame") && isDim(v->dims[1], dimension);
        const bool constant = v->dims.size() == 1 && isDim(v->dims[0], dimension);
        if (!perFrame && !constant) {
            std::string shape;
            for (const std::size_t id : v->dims) {
                shape += (shape.empty() ? "" : ", ") + m_header.dims[id].name;
            }
            throw ReaderError("OSError", name + ": variable '" + v->name + "' has dimensions (" +
                                             shape + "); expected (frame, " +
                                             std::string(dimension) + ") or (" +
                                             std::string(dimension) + ")");
        }
    }

    void checkConventions(const std::string& name) {
        const Attribute* conventions = findAttribute(m_header.attributes, "Conventions");
        if (conventions == nullptr || !conventions->isText()) {
            throw ValueError("NCDF trajectory " + name + " is missing Conventions");
        }
        if (!hasAmberToken(decodedText(*conventions))) {
            typeError("NCDF trajectory " + name +
                      " does not conform to AMBER specifications, "
                      "http://ambermd.org/netcdf/nctraj.xhtml ('AMBER' must be one of the token "
                      "in attribute Conventions)");
        }
        const Attribute* version = findAttribute(m_header.attributes, "ConventionVersion");
        if (version == nullptr || !version->isText()) {
            throw ValueError("NCDF trajectory " + name + " is missing ConventionVersion");
        }
        // a version other than 1.0 is only a warning, but the text must decode
        decodedText(*version);
    }

    void readScaleFactors() {
        for (const Variable& v : m_header.variables) {
            const Attribute* factor = findAttribute(v.attributes, "scale_factor");
            if (factor == nullptr) {
                continue;
            }
            if (std::find(kScaledVariables.begin(), kScaledVariables.end(), v.name) ==
                kScaledVariables.end()) {
                throw ReaderError("NotImplementedError",
                                  "scale_factors for variable " + v.name + " are not implemented");
            }
            if ((factor->type != NcFloat && factor->type != NcDouble) || factor->count != 1) {
                typeError(valueText(*factor) + " is not a float");
            }
            m_scales[v.name] = scaleOf(*factor);
        }
    }

    Scale scaleFor(const Variable& v) const {
        const auto it = m_scales.find(v.name);
        return it == m_scales.end() ? Scale{} : it->second;
    }

    std::uint64_t entryOffset(const Variable& v, std::uint64_t index) const {
        return v.begin + index * (v.record ? m_header.recordStride : v.itemBytes);
    }

    // Entry `index` along the first axis of `v`, scaled, each value exactly as a double.
    std::vector<double> row(const Variable& v, std::uint64_t index) {
        if (v.shape.empty()) {
            throw ReaderError("IndexError", "invalid index to scalar variable.");
        }
        if (index >= v.shape[0]) {
            throw ReaderError("IndexError", "index " + std::to_string(index) +
                                                " is out of bounds for axis 0 with size " +
                                                std::to_string(v.shape[0]));
        }
        const std::string bytes =
            m_file.read(entryOffset(v, index), static_cast<std::size_t>(v.itemBytes));
        const std::size_t size = typeSize(v.type);
        const Scale scale = scaleFor(v);
        std::vector<double> out(bytes.size() / size);
        for (std::size_t i = 0; i < out.size(); ++i) {
            out[i] = scaled(loadValue(bytes.data() + i * size, v.type), v.type, scale);
        }
        return out;
    }

    // The three values of a cell variable for frame `index`; a variable without the frame
    // dimension holds the cell of every frame.
    std::array<double, 3> triple(const Variable& v, std::uint64_t index) {
        std::vector<double> values;
        if (v.dims.size() == 1) {  // one cell for every frame
            for (std::uint64_t k = 0; k < 3; ++k) {
                values.push_back(row(v, k).front());
            }
        } else {
            values = row(v, index);
        }
        if (values.size() != 3) {
            throw ValueError("could not broadcast input array from shape (" +
                             std::to_string(values.size()) + ",) into shape (3,)");
        }
        return {values[0], values[1], values[2]};
    }

    // time[1] - time[0] of the scaled times (as the frames' times), in their type; 1.0 without a
    // time variable or with one frame.
    double fileDt() {
        if (m_time == nullptr || m_frameCount < 2 || m_time->shape.empty() ||
            m_time->shape[0] < 2 || m_time->itemBytes != typeSize(m_time->type)) {
            return 1.0;
        }
        const std::size_t size = typeSize(m_time->type);
        const double t0 = loadValue(m_file.read(entryOffset(*m_time, 0), size).data(), m_time->type);
        const double t1 = loadValue(m_file.read(entryOffset(*m_time, 1), size).data(), m_time->type);
        const Scale scale = scaleFor(*m_time);
        if (scale.applies) {  // RD-N06
            const double s0 = scaled(t0, m_time->type, scale);
            const double s1 = scaled(t1, m_time->type, scale);
            const bool single = scale.single && m_time->type != NcInt && m_time->type != NcDouble;
            return single ? static_cast<double>(static_cast<float>(s1) - static_cast<float>(s0))
                          : s1 - s0;
        }
        if (m_time->type == NcFloat) {
            return static_cast<double>(static_cast<float>(t1) - static_cast<float>(t0));
        }
        if (m_time->type == NcDouble) {
            return t1 - t0;
        }
        // integers: exact in a double, the difference wraps as NumPy's does
        const auto difference = static_cast<std::int64_t>(t1) - static_cast<std::int64_t>(t0);
        switch (m_time->type) {
        case NcInt:
            return static_cast<double>(static_cast<std::int32_t>(difference));
        case NcShort:
            return static_cast<double>(static_cast<std::int16_t>(difference));
        default:
            return static_cast<double>(static_cast<std::int8_t>(difference));
        }
    }

    detail::BinaryFile m_file;
    Header m_header;
    std::size_t m_atoms = 0;
    const Variable* m_coordinates = nullptr;
    const Variable* m_time = nullptr;
    const Variable* m_lengths = nullptr;
    const Variable* m_angles = nullptr;
    std::map<std::string, Scale> m_scales;
};

// ---------------------------------------------------------------- NCDFWriter (scipy backend)

void putName(std::string& out, std::string_view name) {
    detail::storeU32BE(out, static_cast<std::uint32_t>(name.size()));
    out += name;
    out.append((4 - name.size() % 4) % 4, '\0');
}

void putTextAttribute(std::string& out, std::string_view name, std::string_view value) {
    putName(out, name);
    detail::storeU32BE(out, NcChar);
    putName(out, value);
}

}  // namespace

std::unique_ptr<core::Trajectory> openNetcdfTrajectory(const std::filesystem::path& path,
                                                       std::optional<double> dt,
                                                       std::size_t expectedAtoms) {
    return std::make_unique<NetcdfTrajectory>(path, dt, expectedAtoms);
}

void writeNetcdfTrajectory(const std::filesystem::path& path,
                           const std::vector<core::Frame>& frames, std::size_t atoms,
                           const NetcdfWriteOptions& options) {
    if (frames.empty()) {
        return;  // NCDFWriter creates its file with the first frame
    }
    for (const core::Frame& frame : frames) {
        if (frame.positions.size() != atoms * 3) {
            throw ValueError("NCDFWriter: Timestep does not have the correct number of atoms");
        }
    }
    const bool periodic = frames.front().box.has_value();
    if (periodic) {
        for (const core::Frame& frame : frames) {
            if (!frame.box) {
                throw ValueError("NCDFWriter: a frame without a box in a trajectory with boxes");
            }
        }
    }

    struct Var {
        std::string_view name;
        std::vector<std::uint32_t> dims;
        std::string_view units;
        std::uint32_t type;
        std::uint32_t vsize;
    };
    // dimensions: frame, atom, spatial, cell_spatial, cell_angular, label
    std::vector<Var> fixed;
    if (periodic) {
        fixed.push_back({"cell_angular", {4, 5}, "", NcChar, 16});
    }
    fixed.push_back({"spatial", {2}, "", NcChar, 4});
    if (periodic) {
        fixed.push_back({"cell_spatial", {3}, "", NcChar, 4});
    }
    std::vector<Var> records;
    const auto coordinateBytes = static_cast<std::uint32_t>(atoms * 12);
    records.push_back({"coordinates", {0, 1, 2}, "angstrom", NcFloat, coordinateBytes});
    records.push_back({"time", {0}, "picosecond", NcFloat, 4});
    if (periodic) {
        records.push_back({"cell_lengths", {0, 3}, "angstrom", NcDouble, 24});
        records.push_back({"cell_angles", {0, 4}, "degree", NcDouble, 24});
    }

    std::string header = "CDF";
    header += '\x02';
    detail::storeU32BE(header, static_cast<std::uint32_t>(frames.size()));
    detail::storeU32BE(header, kDimensionTag);
    detail::storeU32BE(header, 6);
    const std::pair<std::string_view, std::uint32_t> dims[] = {
        {"frame", 0},        {"atom", static_cast<std::uint32_t>(atoms)},
        {"spatial", 3},      {"cell_spatial", 3},
        {"cell_angular", 3}, {"label", 5}};
    for (const auto& [name, length] : dims) {
        putName(header, name);
        detail::storeU32BE(header, length);
    }
    detail::storeU32BE(header, kAttributeTag);
    detail::storeU32BE(header, 5);
    putTextAttribute(header, "program", options.program);
    putTextAttribute(header, "programVersion", options.programVersion);
    putTextAttribute(header, "Conventions", "AMBER");
    putTextAttribute(header, "ConventionVersion", "1.0");
    putTextAttribute(header, "application", options.application);
    detail::storeU32BE(header, kVariableTag);
    detail::storeU32BE(header, static_cast<std::uint32_t>(fixed.size() + records.size()));
    std::vector<std::size_t> beginAt;
    for (const auto* group : {&fixed, &records}) {
        for (const Var& v : *group) {
            putName(header, v.name);
            detail::storeU32BE(header, static_cast<std::uint32_t>(v.dims.size()));
            for (const std::uint32_t d : v.dims) {
                detail::storeU32BE(header, d);
            }
            if (v.units.empty()) {
                header.append(8, '\0');
            } else {
                detail::storeU32BE(header, kAttributeTag);
                detail::storeU32BE(header, 1);
                putTextAttribute(header, "units", v.units);
            }
            detail::storeU32BE(header, v.type);
            detail::storeU32BE(header, v.vsize);
            beginAt.push_back(header.size());
            header.append(8, '\0');
        }
    }

    std::string data;
    std::size_t slot = 0;
    const auto setBegin = [&](std::size_t at, std::uint64_t begin) {
        std::string bytes;
        detail::storeU64BE(bytes, begin);
        header.replace(at, 8, bytes);
    };
    for (const Var& v : fixed) {
        setBegin(beginAt[slot++], header.size() + data.size());
        if (v.name == "cell_angular") {
            data += std::string("alphabeta gamma", 15);
            data += '\0';
        } else if (v.name == "spatial") {
            data += std::string("xyz\0", 4);
        } else {
            data += std::string("abc\0", 4);
        }
    }
    std::uint64_t offset = header.size() + data.size();
    for (const Var& v : records) {
        setBegin(beginAt[slot++], offset);
        offset += v.vsize;
    }
    for (const core::Frame& frame : frames) {
        for (const float x : frame.positions) {
            detail::storeF32BE(data, x);
        }
        detail::storeF32BE(data, static_cast<float>(frame.time));
        if (periodic) {
            for (const float x : *frame.box) {
                detail::storeF64BE(data, static_cast<double>(x));
            }
        }
    }
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    if (!out) {
        throw std::filesystem::filesystem_error("cannot write file", path,
                                                std::make_error_code(std::errc::permission_denied));
    }
    out.write(header.data(), static_cast<std::streamsize>(header.size()));
    out.write(data.data(), static_cast<std::streamsize>(data.size()));
    if (!out) {
        throw std::filesystem::filesystem_error("cannot write file", path,
                                                std::make_error_code(std::errc::io_error));
    }
}

}  // namespace mdbridge::io
