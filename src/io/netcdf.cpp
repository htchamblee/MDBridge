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

// ---------------------------------------------------------------- Amber conventions (reader)
// CLEAN ROOM: reimplement from scratchpad/cleanroom/netcdf-amber.md

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
    (void)path; (void)dt; (void)expectedAtoms;
    throw std::logic_error("clean room: NetCDF reader not implemented");
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
