#include "mdbridge/io/netcdf.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <stdexcept>
#include <utility>

#include "binary_file.hpp"
#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/json.hpp"
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
                    throw ReaderError("IndexError", "list index out of range");
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
            throw ReaderError("IndexError", "index 0 is out of bounds for axis 0 with size 0");
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

// ---------------------------------------------------------------- Amber conventions (NCDFReader)

// A scale_factor attribute: numpy multiplies in float32 when both the data type and the factor
// allow it (byte / short / float data, float factor), else in float64.
struct Scale {
    double value = 1.0;
    bool isFloat32 = false;
    bool active = false;
};

bool productIsFloat32(std::uint32_t dataType, const Scale& scale) {
    return scale.isFloat32 &&
           (dataType == NcByte || dataType == NcChar || dataType == NcShort || dataType == NcFloat);
}

double scaled(double raw, std::uint32_t dataType, const Scale& scale) {
    if (!scale.active) {
        return raw;
    }
    if (productIsFloat32(dataType, scale)) {
        return static_cast<double>(static_cast<float>(raw) * static_cast<float>(scale.value));
    }
    return raw * scale.value;
}

std::string shapeText(const std::vector<std::uint64_t>& shape) {
    std::string out = "(";
    for (std::size_t i = 0; i < shape.size(); ++i) {
        out += (i ? "," : "") + std::to_string(shape[i]);
    }
    return out + (shape.size() == 1 ? ",)" : ")");
}

// numpy's check of `target[:] = source` (broadcasting).
void checkBroadcast(const std::vector<std::uint64_t>& source,
                    const std::vector<std::uint64_t>& target) {
    bool ok = source.size() <= target.size() ||
              std::all_of(source.begin(), source.end() - static_cast<std::ptrdiff_t>(target.size()),
                          [](std::uint64_t d) { return d == 1; });
    for (std::size_t i = 0; ok && i < std::min(source.size(), target.size()); ++i) {
        const std::uint64_t s = source[source.size() - 1 - i];
        ok = s == 1 || s == target[target.size() - 1 - i];
    }
    if (!ok) {
        throw ValueError("could not broadcast input array from shape " + shapeText(source) +
                         " into shape " + shapeText(target));
    }
}

std::string attributeStr(const Attribute& a) {
    if (a.isText()) {
        return core::pythonStringRepr(a.text()).insert(0, "b");
    }
    const bool integral = a.type != NcFloat && a.type != NcDouble;
    const auto one = [&](std::size_t i) {
        const double v = a.number(i);
        if (integral) {
            return std::to_string(static_cast<long long>(v));
        }
        std::string text = core::pythonFloatRepr(v);
        if (text.size() > 2 && text.compare(text.size() - 2, 2, ".0") == 0) {
            text.pop_back();
        }
        return text;
    };
    if (a.count == 1) {
        return one(0);
    }
    std::string out = "[";
    for (std::size_t i = 0; i < a.count; ++i) {
        out += (i ? " " : "") + one(i);
    }
    return out + "]";
}

std::vector<std::string> splitOn(const std::string& text, char sep) {
    std::vector<std::string> out(1);
    for (const char c : text) {
        if (c == sep) {
            out.emplace_back();
        } else {
            out.back() += c;
        }
    }
    return out;
}

std::vector<std::string> splitWhitespace(const std::string& text) {
    std::vector<std::string> out;
    std::string word;
    for (const char c : text) {
        const bool space = c == ' ' || (c >= '\t' && c <= '\r') || (c >= '\x1c' && c <= '\x1f');
        if (space) {
            if (!word.empty()) {
                out.push_back(std::move(word));
                word.clear();
            }
        } else {
            word += c;
        }
    }
    if (!word.empty()) {
        out.push_back(std::move(word));
    }
    return out;
}

class NetcdfTrajectory : public detail::FileTrajectory {
public:
    NetcdfTrajectory(const std::filesystem::path& path, std::optional<double> dt,
                     std::size_t expectedAtoms)
        : m_file(path), m_path(detail::pathText(path)) {
        if (m_file.size() == 0) {
            throw ValueError("cannot mmap an empty file");
        }
        try {
            m_header = HeaderReader(m_file).read();
        } catch (const std::invalid_argument&) {
            throw ReaderError("TypeError", "Error: " + m_path + " is not a valid NetCDF 3 file");
        }
        checkData();
        checkConventions();
        const Dimension* atom = m_header.dimension("atom");
        if (atom == nullptr) {
            throw ValueError("NCDF trajectory " + m_path + " does not contain atom information");
        }
        m_atoms = atom->length;
        if (expectedAtoms != 0 && expectedAtoms != m_atoms) {
            throw ValueError("Supplied n_atoms (" + std::to_string(expectedAtoms) +
                             ") != natom from ncdf (" + std::to_string(m_atoms) +
                             "). Note: n_atoms can be None and then the ncdf value is used!");
        }
        const Dimension* frame = m_header.dimension("frame");
        const Variable* coordinates = m_header.variable("coordinates");
        if (frame == nullptr || (frame->length == 0 && coordinates == nullptr)) {
            throw ValueError("NCDF trajectory " + m_path + " does not contain frame information");
        }
        m_frameCount =
            static_cast<std::int64_t>(frame->length != 0 ? frame->length : coordinates->shape[0]);
        m_time = m_header.variable("time");
        if (m_time != nullptr) {
            verifyUnits(*m_time, "picosecond");
        }
        if (coordinates == nullptr) {
            throw ReaderError("KeyError", "'coordinates'");
        }
        m_coordinates = coordinates;
        verifyUnits(*m_coordinates, "angstrom");
        readScaleFactors();
        if (const Variable* v = m_header.variable("velocities")) {
            verifyUnits(*v, "angstrom/picosecond");
        }
        if (const Variable* f = m_header.variable("forces")) {
            verifyUnits(*f, "kilocalorie/mole/angstrom");
        }
        m_lengths = m_header.variable("cell_lengths");
        if (m_lengths != nullptr) {
            verifyUnits(*m_lengths, "angstrom");
            m_angles = m_header.variable("cell_angles");
            if (m_angles == nullptr) {
                throw ReaderError("KeyError", "'cell_angles'");
            }
            verifyUnits(*m_angles, "degree");
        }
        if (dt) {
            m_dt = *dt;
        } else if (m_time != nullptr && m_time->shape.size() >= 1 && m_time->shape[0] >= 2) {
            m_dt = timeStep();
        }
        if (m_frameCount == 0) {
            throw ReaderError("IndexError", "frame index must be 0 <= frame < 0");
        }
        seek(0);
    }

    bool storedTime(std::int64_t) const override { return m_time != nullptr; }

protected:
    void readFrame(std::int64_t index, core::Frame& frame) override {
        const auto k = static_cast<std::uint64_t>(index);
        if (k >= m_coordinates->shape[0]) {
            throw ReaderError("IndexError", "index " + std::to_string(k) +
                                                " is out of bounds for axis 0 with size " +
                                                std::to_string(m_coordinates->shape[0]));
        }
        std::vector<std::uint64_t> shape(m_coordinates->shape.begin() + 1,
                                         m_coordinates->shape.end());
        checkBroadcast(shape, {m_atoms, 3});
        const std::vector<double> raw = values(*m_coordinates, k);
        const Scale& scale = m_scales["coordinates"];
        frame.positions.assign(m_atoms * 3, 0.0f);
        const bool f32 = productIsFloat32(m_coordinates->type, scale);
        for (std::size_t i = 0; i < frame.positions.size(); ++i) {
            const double x = raw[broadcastIndex(shape, {m_atoms, 3}, i)];
            if (!scale.active) {
                frame.positions[i] = static_cast<float>(x);
            } else if (f32) {
                frame.positions[i] = static_cast<float>(x) * static_cast<float>(scale.value);
            } else {
                frame.positions[i] = static_cast<float>(x * scale.value);
            }
        }
        if (m_time != nullptr) {
            frame.time = scaled(values(*m_time, k).at(0), m_time->type, m_scales["time"]);
        } else {
            frame.time = static_cast<double>(index) * m_dt;
        }
        if (m_lengths != nullptr) {
            std::array<float, 6> box{};
            const Variable* parts[2] = {m_lengths, m_angles};
            for (int p = 0; p < 2; ++p) {
                const Variable& v = *parts[p];
                std::vector<std::uint64_t> partShape(v.shape.begin() + 1, v.shape.end());
                checkBroadcast(partShape, {3});
                const std::vector<double> raw3 = values(v, k);
                for (std::size_t i = 0; i < 3; ++i) {
                    box[p * 3 + i] = static_cast<float>(
                        scaled(raw3[broadcastIndex(partShape, {3}, i)], v.type, m_scales[v.name]));
                }
            }
            if (box[0] == 0.0f && box[1] == 0.0f && box[2] == 0.0f) {
                throw ReaderError("TypeError", "'NoneType' object is not subscriptable");
            }
            frame.box = box;
        }
    }

private:
    void checkData() {
        const std::uint64_t size = m_file.size();
        for (const Variable& v : m_header.variables) {
            if (v.record) {
                continue;
            }
            std::uint64_t bytes = typeSize(v.type);
            for (const std::uint64_t d : v.shape) {
                bytes *= d;
            }
            const std::uint64_t available = v.begin >= size ? 0 : std::min(bytes, size - v.begin);
            if (available != bytes) {
                throw ValueError("cannot reshape array of size " +
                                 std::to_string(available / typeSize(v.type)) + " into shape " +
                                 shapeText(v.shape));
            }
        }
        if (m_header.recordStride == 0) {
            return;
        }
        std::uint64_t recordSize = 0;
        for (const Variable& v : m_header.variables) {
            if (v.record) {
                recordSize += v.vsize;
            }
        }
        const std::uint64_t wanted = m_header.records * recordSize;
        const std::uint64_t available =
            m_header.recordBegin >= size ? 0 : std::min(wanted, size - m_header.recordBegin);
        if (available % m_header.recordStride != 0) {
            throw ValueError(
                "When changing to a larger dtype, its size must be a divisor of the total size in "
                "bytes of the last axis of the array.");
        }
        if (available / m_header.recordStride < m_header.records) {
            throw ValueError("cannot reshape array of size " +
                             std::to_string(available / m_header.recordStride) + " into shape (" +
                             std::to_string(m_header.records) + ",)");
        }
    }

    void checkConventions() {
        const Attribute* conventions = findAttribute(m_header.attributes, "Conventions");
        if (conventions == nullptr || !conventions->isText()) {
            throw ValueError("NCDF trajectory " + m_path + " is missing Conventions");
        }
        const std::string text = conventions->text();
        const auto commas = splitOn(text, ',');
        const auto words = splitWhitespace(text);
        if (std::find(commas.begin(), commas.end(), "AMBER") == commas.end() &&
            std::find(words.begin(), words.end(), "AMBER") == words.end()) {
            throw ReaderError("TypeError",
                              "NCDF trajectory " + m_path +
                                  " does not conform to AMBER specifications, "
                                  "http://ambermd.org/netcdf/nctraj.xhtml ('AMBER' must be one "
                                  "of the token in attribute Conventions)");
        }
        const Attribute* version = findAttribute(m_header.attributes, "ConventionVersion");
        if (version == nullptr || !version->isText()) {
            throw ValueError("NCDF trajectory " + m_path + " is missing ConventionVersion");
        }
        if (m_header.version != 2) {
            throw ReaderError("TypeError",
                              "NCDF trajectory " + m_path +
                                  " does not conform to AMBER specifications, as detailed in "
                                  "https://ambermd.org/netcdf/nctraj.xhtml (NetCDF file does not "
                                  "use 64 bit offsets [version_byte = 2])");
        }
        const Dimension* spatial = m_header.dimension("spatial");
        if (spatial == nullptr) {
            throw ValueError("NCDF trajectory does not contain spatial dimension");
        }
        if (spatial->length != 3) {
            throw ReaderError("TypeError", "Incorrect spatial value for NCDF trajectory file");
        }
    }

    static void verifyUnits(const Variable& v, const std::string& expected) {
        const Attribute* units = findAttribute(v.attributes, "units");
        if (units == nullptr) {
            throw ReaderError("AttributeError",
                              "'netcdf_variable' object has no attribute 'units'");
        }
        if (!units->isText()) {
            throw ReaderError("AttributeError", "'numpy.ndarray' object has no attribute 'decode'");
        }
        if (units->text() != expected) {
            throw ReaderError(
                "NotImplementedError",
                "NETCDFReader currently assumes that the trajectory was written in units of " +
                    units->text() + " instead of " + expected);
        }
    }

    void readScaleFactors() {
        static const std::array<std::string_view, 6> known = {
            "time", "cell_lengths", "cell_angles", "coordinates", "velocities", "forces"};
        for (const Variable& v : m_header.variables) {
            const Attribute* factor = findAttribute(v.attributes, "scale_factor");
            if (factor == nullptr) {
                continue;
            }
            if (std::find(known.begin(), known.end(), v.name) == known.end()) {
                throw ReaderError("NotImplementedError",
                                  "scale_factors for variable " + v.name + " are not implemented");
            }
            if ((factor->type != NcFloat && factor->type != NcDouble) || factor->count != 1) {
                throw ReaderError("TypeError", attributeStr(*factor) + " is not a float");
            }
            Scale scale;
            scale.value = factor->number(0);
            scale.isFloat32 = factor->type == NcFloat;
            // math.isclose(scale_factor, 1): such factors are not applied
            scale.active =
                !(std::fabs(scale.value - 1.0) <= 1e-9 * std::max(std::fabs(scale.value), 1.0));
            m_scales[v.name] = scale;
        }
    }

    // time[1] - time[0] in the variable's own type (no scale factor)
    double timeStep() {
        const double t1 = values(*m_time, 1).at(0);
        const double t0 = values(*m_time, 0).at(0);
        if (m_time->type == NcFloat) {
            return static_cast<double>(static_cast<float>(t1) - static_cast<float>(t0));
        }
        if (m_time->type == NcDouble) {
            return t1 - t0;
        }
        return static_cast<double>(static_cast<std::int64_t>(t1) - static_cast<std::int64_t>(t0));
    }

    // The values of entry k of the first axis.
    std::vector<double> values(const Variable& v, std::uint64_t k) {
        const std::size_t size = typeSize(v.type);
        const std::uint64_t offset =
            v.record ? m_header.recordBegin + k * m_header.recordStride + v.recordOffset
                     : v.begin + k * v.itemBytes;
        const std::string bytes = m_file.read(offset, static_cast<std::size_t>(v.itemBytes));
        if (bytes.size() != v.itemBytes) {
            throw ValueError("NCDF trajectory " + m_path + " is truncated");
        }
        std::vector<double> out(bytes.size() / size);
        for (std::size_t i = 0; i < out.size(); ++i) {
            out[i] = loadValue(bytes.data() + i * size, v.type);
        }
        return out;
    }

    static std::size_t broadcastIndex(const std::vector<std::uint64_t>& source,
                                      const std::vector<std::uint64_t>& target, std::size_t flat) {
        std::size_t index = 0;
        std::size_t stride = 1;
        for (std::size_t i = 0; i < source.size() && i < target.size(); ++i) {
            const std::uint64_t t = target[target.size() - 1 - i];
            const std::uint64_t s = source[source.size() - 1 - i];
            const std::size_t coordinate = flat % t;
            flat /= t;
            if (s != 1) {
                index += coordinate * stride;
            }
            stride *= s;
        }
        return index;
    }

    detail::BinaryFile m_file;
    std::string m_path;
    Header m_header;
    std::uint64_t m_atoms = 0;
    const Variable* m_coordinates = nullptr;
    const Variable* m_time = nullptr;
    const Variable* m_lengths = nullptr;
    const Variable* m_angles = nullptr;
    std::map<std::string, Scale, std::less<>> m_scales;
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
