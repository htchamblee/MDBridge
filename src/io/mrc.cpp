#include "mdbridge/io/mrc.hpp"

#include <zlib.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <deque>
#include <fstream>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <system_error>
#include <tuple>
#include <utility>

#include "binary_file.hpp"
#include "mdbridge/core/detail/numpy_sum.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/python.hpp"
#include "python_text.hpp"

// Compiled without floating-point contraction (see CMakeLists.txt): the statistics and simulated
// maps round every operation as numpy does, in numpy's order. The sums are numpy's pairwise
// summation (mdbridge/core/detail/numpy_sum.hpp, a port of numpy's code under its licence).

namespace mdbridge::io {

namespace fs = std::filesystem;

namespace {

constexpr std::size_t kHeaderBytes = 1024;
constexpr std::int32_t kImodStamp = 1146047817;
constexpr std::size_t kChunk = 1 << 16;

[[noreturn]] void fail(const std::string& path, const std::string& message) {
    throw MapError(path + ": " + message);
}

std::ifstream openRead(const fs::path& path) {
    std::error_code ec;
    if (fs::is_directory(path, ec)) {
        throw fs::filesystem_error("cannot read file", path,
                                   std::make_error_code(std::errc::is_a_directory));
    }
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw fs::filesystem_error("cannot read file", path,
                                   std::make_error_code(std::errc::no_such_file_or_directory));
    }
    return in;
}

// Python's repr() of a bytes object.
std::string bytesRepr(std::string_view bytes) {
    const bool doubleQuotes =
        bytes.find('\'') != std::string_view::npos && bytes.find('"') == std::string_view::npos;
    const char quote = doubleQuotes ? '"' : '\'';
    std::string out = "b";
    out += quote;
    for (const char ch : bytes) {
        const auto c = static_cast<unsigned char>(ch);
        if (c == '\\' || c == static_cast<unsigned char>(quote)) {
            out += '\\';
            out += ch;
        } else if (c == '\t') {
            out += "\\t";
        } else if (c == '\n') {
            out += "\\n";
        } else if (c == '\r') {
            out += "\\r";
        } else if (c < 0x20 || c >= 0x7f) {
            static const char* digits = "0123456789abcdef";
            out += "\\x";
            out += digits[c >> 4];
            out += digits[c & 15];
        } else {
            out += ch;
        }
    }
    out += quote;
    return out;
}

// ---- the bytes of a map file, plain or gzip ---------------------------------------------------

class Source {
public:
    virtual ~Source() = default;
    // Up to n bytes; fewer only at the end of the data.
    virtual std::size_t read(char* out, std::size_t n) = 0;
    // The bytes left, when known without decompressing.
    virtual std::optional<std::uint64_t> remaining() const = 0;

    std::uint64_t skip(std::uint64_t n) {
        char buffer[4096];
        std::uint64_t done = 0;
        while (done < n) {
            const std::size_t want =
                static_cast<std::size_t>(std::min<std::uint64_t>(n - done, sizeof buffer));
            const std::size_t got = read(buffer, want);
            done += got;
            if (got < want) {
                break;
            }
        }
        return done;
    }
};

class PlainSource : public Source {
public:
    PlainSource(std::ifstream in, std::uint64_t size) : m_in(std::move(in)), m_size(size) {}

    std::size_t read(char* out, std::size_t n) override {
        m_in.read(out, static_cast<std::streamsize>(n));
        const auto got = static_cast<std::size_t>(m_in.gcount());
        m_pos += got;
        return got;
    }

    std::optional<std::uint64_t> remaining() const override {
        return m_size > m_pos ? m_size - m_pos : 0;
    }

private:
    std::ifstream m_in;
    std::uint64_t m_size;
    std::uint64_t m_pos = 0;
};

// Python's gzip reader: concatenated members, zero padding between them, each member decompressed
// only as far as the data is needed.
class GzipSource : public Source {
public:
    GzipSource(std::ifstream in, std::string path)
            : m_in(std::move(in)), m_path(std::move(path)), m_buffer(kChunk) {
        if (inflateInit2(&m_z, -MAX_WBITS) != Z_OK) {
            throw std::bad_alloc();
        }
    }
    ~GzipSource() override { inflateEnd(&m_z); }
    GzipSource(const GzipSource&) = delete;
    GzipSource& operator=(const GzipSource&) = delete;

    std::optional<std::uint64_t> remaining() const override { return std::nullopt; }

    std::size_t read(char* out, std::size_t n) override {
        std::size_t done = 0;
        while (done < n && m_state != State::End) {
            if (m_state == State::Trailer) {
                finishMember();
                continue;
            }
            if (m_state == State::Header) {
                if (!readHeader()) {
                    m_state = State::End;
                    break;
                }
                inflateReset(&m_z);
                m_crc = crc32(0, nullptr, 0);
                m_size = 0;
                m_state = State::Body;
            }
            if (!fill()) {  // the stream ends inside a member
                m_state = State::End;
                break;
            }
            m_z.next_in = m_buffer.data() + m_pos;
            m_z.avail_in = static_cast<uInt>(m_end - m_pos);
            auto* target = reinterpret_cast<Bytef*>(out + done);
            m_z.next_out = target;
            m_z.avail_out = static_cast<uInt>(std::min<std::size_t>(n - done, 1u << 30));
            const int rc = inflate(&m_z, Z_NO_FLUSH);
            const auto produced = static_cast<std::size_t>(m_z.next_out - target);
            m_pos = m_end - m_z.avail_in;
            m_crc = crc32(m_crc, target, static_cast<uInt>(produced));
            m_size += produced;
            done += produced;
            if (rc == Z_STREAM_END) {
                m_state = State::Trailer;
            } else if (rc == Z_DATA_ERROR) {
                fail(m_path, std::string("Error -3 while decompressing data: ") +
                                 (m_z.msg ? m_z.msg : "invalid input data"));
            } else if (rc == Z_MEM_ERROR) {
                throw std::bad_alloc();
            } else if (rc != Z_OK && rc != Z_BUF_ERROR) {
                fail(m_path, "Error " + std::to_string(rc) + " while decompressing data");
            }
        }
        return done;
    }

private:
    enum class State { Header, Body, Trailer, End };

    bool fill() {
        if (m_pos < m_end) {
            return true;
        }
        m_in.read(reinterpret_cast<char*>(m_buffer.data()), static_cast<std::streamsize>(kChunk));
        m_pos = 0;
        m_end = static_cast<std::size_t>(m_in.gcount());
        return m_end > 0;
    }

    int byte() {
        if (!fill()) {
            return -1;
        }
        return m_buffer[m_pos++];
    }

    bool take(std::size_t n, unsigned char* out = nullptr) {
        for (std::size_t k = 0; k < n; ++k) {
            const int b = byte();
            if (b < 0) {
                return false;
            }
            if (out) {
                out[k] = static_cast<unsigned char>(b);
            }
        }
        return true;
    }

    // gzip._read_gzip_header; false at the end of the file (or of a truncated header)
    bool readHeader() {
        const int b0 = byte();
        if (b0 < 0) {
            return false;
        }
        const int b1 = byte();
        if (b0 != 0x1f || b1 != 0x8b) {
            std::string magic(1, static_cast<char>(b0));
            if (b1 >= 0) {
                magic += static_cast<char>(b1);
            }
            fail(m_path, "Not a gzipped file (" + bytesRepr(magic) + ")");
        }
        unsigned char head[8];
        if (!take(8, head)) {
            return false;
        }
        if (head[0] != 8) {
            fail(m_path, "Unknown compression method");
        }
        const int flags = head[1];
        if (flags & 4) {  // FEXTRA
            unsigned char size[2];
            if (!take(2, size) || !take(static_cast<std::size_t>(size[0] | size[1] << 8))) {
                return false;
            }
        }
        for (const int flag : {8, 16}) {  // FNAME, FCOMMENT: zero-terminated
            if (flags & flag) {
                for (int b = byte(); b > 0; b = byte()) {}
            }
        }
        if ((flags & 2) && !take(2)) {  // FHCRC
            return false;
        }
        return true;
    }

    // gzip._GzipReader._read_eof: CRC and size, then zero padding
    void finishMember() {
        unsigned char trailer[8];
        if (!take(8, trailer)) {
            m_state = State::End;
            return;
        }
        const std::uint32_t crc = detail::loadU32(trailer, false);
        const std::uint32_t size = detail::loadU32(trailer + 4, false);
        if (crc != static_cast<std::uint32_t>(m_crc)) {
            const auto hex = [](std::uint32_t v) {
                static const char* digits = "0123456789abcdef";
                std::string out;
                do {
                    out.insert(out.begin(), digits[v & 15]);
                    v >>= 4;
                } while (v);
                return "0x" + out;
            };
            fail(m_path,
                 "CRC check failed " + hex(crc) + " != " + hex(static_cast<std::uint32_t>(m_crc)));
        }
        if (size != static_cast<std::uint32_t>(m_size & 0xFFFFFFFFu)) {
            fail(m_path, "Incorrect length of data produced");
        }
        for (;;) {
            const int b = byte();
            if (b != 0) {
                if (b > 0) {
                    --m_pos;  // the byte came from the buffer
                }
                break;
            }
        }
        m_state = State::Header;
    }

    std::ifstream m_in;
    std::string m_path;
    std::vector<unsigned char> m_buffer;
    std::size_t m_pos = 0, m_end = 0;
    z_stream m_z{};
    State m_state = State::Header;
    uLong m_crc = 0;
    std::uint64_t m_size = 0;
};

std::size_t itemSize(int mode) {
    switch (mode) {
    case 0:
        return 1;
    case 1:
    case 6:
    case 12:
        return 2;
    default:
        return 4;
    }
}

bool supportedMode(std::int32_t mode) {
    return mode == 0 || mode == 1 || mode == 2 || mode == 6 || mode == 12;
}

float halfToFloat(std::uint16_t h) {
    const std::uint32_t sign = static_cast<std::uint32_t>(h & 0x8000u) << 16;
    const std::uint32_t exponent = h >> 10 & 0x1Fu;
    std::uint32_t mantissa = h & 0x3FFu;
    std::uint32_t bits;
    if (exponent == 0x1F) {
        bits = sign | 0x7F800000u | mantissa << 13;
    } else if (exponent != 0) {
        bits = sign | (exponent + 112) << 23 | mantissa << 13;
    } else if (mantissa == 0) {
        bits = sign;
    } else {  // subnormal: normalise
        std::uint32_t e = 113;
        while (!(mantissa & 0x400u)) {
            mantissa <<= 1;
            --e;
        }
        bits = sign | e << 23 | (mantissa & 0x3FFu) << 13;
    }
    float value;
    std::memcpy(&value, &bits, 4);
    return value;
}

std::uint64_t saturatingMul(std::uint64_t a, std::uint64_t b) {
    if (a != 0 && b > std::numeric_limits<std::uint64_t>::max() / a) {
        return std::numeric_limits<std::uint64_t>::max();
    }
    return a * b;
}

std::string absolutePath(const fs::path& path) {
#ifdef _WIN32
    return core::pathToUtf8(fs::absolute(path).lexically_normal());
#else
    return core::posixpath::abspath(core::pathToUtf8(path), core::currentDirectory());
#endif
}

// np.add.reduce of a contiguous float32 array: the identity 0 plus numpy's pairwise sum
template <class F>
float numpySum(const F& value, std::size_t n) {
    return core::detail::numpySum<float>(value, n);
}

void checkShape(const DensityMap& map) {
    const auto& n = map.shapeXyz;
    if (n[0] < 0 || n[1] < 0 || n[2] < 0 ||
        saturatingMul(
            saturatingMul(static_cast<std::uint64_t>(n[0]), static_cast<std::uint64_t>(n[1])),
            static_cast<std::uint64_t>(n[2])) != map.data.size()) {
        throw std::invalid_argument("density map: the data do not match the grid shape");
    }
}

}  // namespace

// ---- DensityMap ------------------------------------------------------------------------------

float DensityMap::at(std::int64_t x, std::int64_t y, std::int64_t z) const {
    if (x < 0 || y < 0 || z < 0 || x >= shapeXyz[0] || y >= shapeXyz[1] || z >= shapeXyz[2]) {
        throw std::out_of_range("density map: voxel index out of range");
    }
    return data[static_cast<std::size_t>((z * shapeXyz[1] + y) * shapeXyz[0] + x)];
}

core::Vec3 DensityMap::extent() const {
    core::Vec3 out;
    for (std::size_t k = 0; k < 3; ++k) {
        out[k] = origin[k] + static_cast<double>(shapeXyz[k] - 1) * voxel[k];
    }
    return out;
}

MapStats DensityMap::stats() const {
    if (data.empty()) {
        throw core::ValueError(
            "zero-size array to reduction operation minimum which has no identity");
    }
    float lo = data[0], hi = data[0];
    for (const float x : data) {
        if (std::isnan(lo)) {
            break;
        }
        if (std::isnan(x)) {
            lo = hi = x;
        } else {
            lo = x < lo ? x : lo;
            hi = x > hi ? x : hi;
        }
    }
    const std::size_t n = data.size();
    const double count = static_cast<double>(n);
    const float* d = data.data();
    const float sum = numpySum([d](std::size_t i) { return d[i]; }, n);
    const auto mean = static_cast<float>(static_cast<double>(sum) / count);
    const float squares = numpySum(
        [d, mean](std::size_t i) {
            const float x = d[i] - mean;
            return x * x;
        },
        n);
    const auto variance = static_cast<float>(static_cast<double>(squares) / count);
    return {lo, hi, mean, std::sqrt(variance)};
}

core::Vec3 DensityMap::indexToXyz(const core::Vec3& index) const {
    return {origin[0] + index[0] * voxel[0], origin[1] + index[1] * voxel[1],
            origin[2] + index[2] * voxel[2]};
}

core::Vec3 DensityMap::xyzToIndex(const core::Vec3& xyz) const {
    return {(xyz[0] - origin[0]) / voxel[0], (xyz[1] - origin[1]) / voxel[1],
            (xyz[2] - origin[2]) / voxel[2]};
}

// ---- reading ---------------------------------------------------------------------------------

DensityMap readMap(const fs::path& path) {
    const std::string name = core::pathToUtf8(path);
    std::ifstream in = openRead(path);
    char magic[2] = {};
    in.read(magic, 2);
    const bool gzip = in.gcount() == 2 && static_cast<unsigned char>(magic[0]) == 0x1f &&
                      static_cast<unsigned char>(magic[1]) == 0x8b;
    in.clear();
    in.seekg(0, std::ios::end);
    const auto size = static_cast<std::uint64_t>(in.tellg());
    in.seekg(0);
    std::unique_ptr<Source> source;
    if (gzip) {
        source = std::make_unique<GzipSource>(std::move(in), name);
    } else {
        source = std::make_unique<PlainSource>(std::move(in), size);
    }

    std::string raw(kHeaderBytes, '\0');
    if (source->read(raw.data(), kHeaderBytes) < kHeaderBytes) {
        fail(name, "file too short for an MRC header");
    }
    const auto byteAt = [&raw](std::size_t k) {
        return static_cast<unsigned char>(raw[k]);
    };
    bool big = false;
    if (byteAt(212) == 0x44) {
        big = false;
    } else if (byteAt(212) == 0x11) {
        big = true;
    } else {  // no usable machine stamp: the byte order whose MODE and MAPC look valid
        bool found = false;
        for (const bool candidate : {false, true}) {
            const std::int32_t mode = detail::loadI32(raw.data() + 12, candidate);
            std::array<std::int32_t, 3> crs{};
            for (std::size_t k = 0; k < 3; ++k) {
                crs[k] = detail::loadI32(raw.data() + 64 + 4 * k, candidate);
            }
            std::sort(crs.begin(), crs.end());
            if (supportedMode(mode) && (crs == std::array<std::int32_t, 3>{1, 2, 3} ||
                                        crs == std::array<std::int32_t, 3>{0, 0, 0})) {
                big = candidate;
                found = true;
                break;
            }
        }
        if (!found) {
            throw MapError("cannot determine the byte order of the map header");
        }
    }
    const auto word = [&raw, big](std::size_t k) {
        return detail::loadI32(raw.data() + 4 * k, big);
    };
    const auto real = [&raw, big](std::size_t k) {
        return static_cast<double>(detail::loadF32(raw.data() + 4 * k, big));
    };
    const std::int32_t nc = word(0), nr = word(1), ns = word(2), mode = word(3);
    const core::Vec3 nstart{double(word(4)), double(word(5)), double(word(6))};
    core::Vec3 mxyz{double(word(7)), double(word(8)), double(word(9))};
    const core::Vec3 cell{real(10), real(11), real(12)};
    const core::Vec3 angles{real(13), real(14), real(15)};
    std::array<std::int32_t, 3> mapcrs{word(16), word(17), word(18)};
    const std::int32_t nsymbt = word(23);
    const core::Vec3 originField{real(49), real(50), real(51)};

    if (!supportedMode(mode)) {
        fail(name, "unsupported MRC mode " + std::to_string(mode) + " (supported: 0, 1, 2, 6, 12)");
    }
    if (std::min({nc, nr, ns}) <= 0) {
        fail(name, "invalid dimensions " + std::to_string(nc) + " x " + std::to_string(nr) + " x " +
                       std::to_string(ns));
    }
    std::array<std::int32_t, 3> sorted = mapcrs;
    std::sort(sorted.begin(), sorted.end());
    if (sorted != std::array<std::int32_t, 3>{1, 2, 3}) {
        if (mapcrs == std::array<std::int32_t, 3>{0, 0, 0}) {  // empty field: default order
            mapcrs = {1, 2, 3};
        } else {
            fail(name, "invalid axis order MAPC/MAPR/MAPS = [" + std::to_string(mapcrs[0]) + ", " +
                           std::to_string(mapcrs[1]) + ", " + std::to_string(mapcrs[2]) + "]");
        }
    }
    const bool skewed = std::any_of(angles.begin(), angles.end(),
                                    [](double a) { return std::abs(a - 90.0) > 0.01; });
    const bool positive = std::all_of(angles.begin(), angles.end(), [](double a) { return a > 0; });
    if (skewed && positive) {
        fail(name, "non-orthogonal cell angles [" + core::pythonFloatRepr(angles[0]) + ", " +
                       core::pythonFloatRepr(angles[1]) + ", " + core::pythonFloatRepr(angles[2]) +
                       "] are not supported");
    }
    const bool unsignedBytes = mode == 0 && (raw.compare(208, 4, "MAP ") != 0 ||
                                             (word(38) == kImodStamp && !(word(39) & 1)));
    const std::size_t item = itemSize(mode);
    const std::uint64_t count =
        saturatingMul(saturatingMul(static_cast<std::uint64_t>(nc), static_cast<std::uint64_t>(nr)),
                      static_cast<std::uint64_t>(ns));
    const std::uint64_t bytes = saturatingMul(count, item);
    if (nsymbt > 0) {
        source->skip(static_cast<std::uint64_t>(nsymbt));
    }
    const std::optional<std::uint64_t> left = source->remaining();
    if ((left && *left < bytes) || bytes > std::numeric_limits<std::size_t>::max() / 2) {
        fail(name, "file is shorter than its header says");
    }
    std::string buffer;
    if (left) {
        buffer.resize(static_cast<std::size_t>(bytes));
        buffer.resize(source->read(buffer.data(), buffer.size()));
    } else {  // compressed: grow with the data, never by the header's word alone
        while (buffer.size() < bytes) {
            const std::size_t have = buffer.size();
            const std::size_t want = static_cast<std::size_t>(
                std::min<std::uint64_t>(bytes - have, std::max<std::size_t>(have, kChunk)));
            buffer.resize(have + want);
            const std::size_t got = source->read(buffer.data() + have, want);
            buffer.resize(have + got);
            if (got < want) {
                break;
            }
        }
    }
    if (buffer.size() < bytes) {
        fail(name, "file is shorter than its header says");
    }

    // file axis k (column, row, section) holds Cartesian axis mapcrs[k] - 1
    const std::array<std::size_t, 3> axisOf{static_cast<std::size_t>(mapcrs[0] - 1),
                                            static_cast<std::size_t>(mapcrs[1] - 1),
                                            static_cast<std::size_t>(mapcrs[2] - 1)};
    const std::array<std::int64_t, 3> nFile{nc, nr, ns};
    DensityMap map;
    core::Vec3 nCart{}, startCart{};
    for (std::size_t k = 0; k < 3; ++k) {
        map.shapeXyz[axisOf[k]] = nFile[k];
        nCart[axisOf[k]] = static_cast<double>(nFile[k]);
        startCart[axisOf[k]] = nstart[k];
    }
    const std::array<std::size_t, 3> cartStride{
        1, static_cast<std::size_t>(map.shapeXyz[0]),
        static_cast<std::size_t>(map.shapeXyz[0] * map.shapeXyz[1])};
    const std::size_t sc = cartStride[axisOf[0]], sr = cartStride[axisOf[1]],
                      ss = cartStride[axisOf[2]];
    map.data.resize(static_cast<std::size_t>(count));
    const char* p = buffer.data();
    for (std::int64_t s = 0; s < ns; ++s) {
        for (std::int64_t r = 0; r < nr; ++r) {
            std::size_t out = static_cast<std::size_t>(s) * ss + static_cast<std::size_t>(r) * sr;
            for (std::int64_t c = 0; c < nc; ++c, out += sc, p += item) {
                float value;
                switch (mode) {
                case 0:
                    value = unsignedBytes ? float(static_cast<unsigned char>(*p))
                                          : float(static_cast<signed char>(*p));
                    break;
                case 1:
                    value = float(static_cast<std::int16_t>(detail::loadU16(p, big)));
                    break;
                case 6:
                    value = float(detail::loadU16(p, big));
                    break;
                case 12:
                    value = halfToFloat(detail::loadU16(p, big));
                    break;
                default:
                    value = detail::loadF32(p, big);
                }
                map.data[out] = value;
            }
        }
    }
    for (std::size_t k = 0; k < 3; ++k) {
        mxyz[k] = mxyz[k] > 0 ? mxyz[k] : nCart[k];
    }
    if (std::any_of(cell.begin(), cell.end(), [](double c) { return c <= 0; })) {
        map.voxel = {1, 1, 1};
    } else {
        for (std::size_t k = 0; k < 3; ++k) {
            map.voxel[k] = cell[k] / mxyz[k];
        }
    }
    const bool hasOrigin =
        std::any_of(originField.begin(), originField.end(), [](double o) { return o != 0; }) &&
        std::all_of(originField.begin(), originField.end(),
                    [](double o) { return std::isfinite(o); });
    for (std::size_t k = 0; k < 3; ++k) {
        map.origin[k] = hasOrigin ? originField[k] : startCart[k] * map.voxel[k];
    }
    map.path = name;
    return map;
}

// ---- writing ---------------------------------------------------------------------------------

void writeMap(const fs::path& path, const DensityMap& map) {
    checkShape(map);
    const MapStats st = map.stats();  // numpy raises for empty maps before the file is opened
    std::array<std::uint32_t, 256> words{};
    const auto setInt = [&words](std::size_t k, std::int64_t v) {
        words[k] = static_cast<std::uint32_t>(v);
    };
    const auto setFloat = [&words](std::size_t k, double v) {
        const auto f = static_cast<float>(v);
        std::memcpy(&words[k], &f, 4);
    };
    const auto& n = map.shapeXyz;
    for (std::size_t k = 0; k < 3; ++k) {
        setInt(k, n[k]);
        setInt(7 + k, n[k]);
    }
    setInt(3, 2);
    // np.allclose(start, np.round(start), atol=1e-4)
    core::Vec3 start{}, rounded{};
    bool whole = true;
    for (std::size_t k = 0; k < 3; ++k) {
        start[k] = map.origin[k] / map.voxel[k];
        rounded[k] = std::nearbyint(start[k]);
        const bool close = (std::abs(start[k] - rounded[k]) <= 1e-4 + 1e-5 * std::abs(rounded[k]) &&
                            std::isfinite(rounded[k])) ||
                           start[k] == rounded[k];
        whole = whole && close;
    }
    if (whole) {
        for (std::size_t k = 0; k < 3; ++k) {
            // numpy's float64 -> int32 cast: INT32_MIN outside the range
            const double r = rounded[k];
            setInt(4 + k, r >= -2147483648.0 && r <= 2147483647.0 ? static_cast<std::int64_t>(r)
                                                                  : std::int64_t{-2147483647 - 1});
        }
    }
    for (std::size_t k = 0; k < 3; ++k) {
        setFloat(10 + k, static_cast<double>(n[k]) * map.voxel[k]);
        setFloat(13 + k, 90.0);
        setInt(16 + k, static_cast<std::int64_t>(k + 1));
        setFloat(49 + k, map.origin[k]);
    }
    setFloat(19, st.min);
    setFloat(20, st.max);
    setFloat(21, st.mean);
    setInt(22, 1);
    setInt(27, 20140);
    setFloat(54, st.rms);
    std::string raw;
    raw.reserve(kHeaderBytes + 4 * map.data.size());
    const auto put = [&raw](std::uint32_t v) {
        for (int shift = 0; shift < 32; shift += 8) {
            raw += static_cast<char>(v >> shift & 0xFF);
        }
    };
    for (const std::uint32_t w : words) {
        put(w);
    }
    raw.replace(208, 4, "MAP ");
    raw.replace(212, 4, std::string("\x44\x44\x00\x00", 4));
    for (const float x : map.data) {
        std::uint32_t bits;
        std::memcpy(&bits, &x, 4);
        put(bits);
    }
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    if (!out) {
        throw fs::filesystem_error("cannot write file", path,
                                   std::make_error_code(std::errc::permission_denied));
    }
    out.write(raw.data(), static_cast<std::streamsize>(raw.size()));
    if (!out) {
        throw fs::filesystem_error("cannot write file", path,
                                   std::make_error_code(std::errc::io_error));
    }
}

// ---- cache, level, info ----------------------------------------------------------------------

std::shared_ptr<const DensityMap> loadMapCached(const fs::path& path) {
    using Key = std::tuple<std::string, std::uintmax_t, std::int64_t>;
    static std::mutex mutex;
    static std::deque<std::pair<Key, std::shared_ptr<const DensityMap>>> cache;
    const std::uintmax_t size = fs::file_size(path);
    const auto mtime =
        static_cast<std::int64_t>(fs::last_write_time(path).time_since_epoch().count());
    const Key key{absolutePath(path), size, mtime};
    const auto find = [&key]() -> std::shared_ptr<const DensityMap> {
        for (const auto& [k, map] : cache) {
            if (k == key) {
                return map;
            }
        }
        return nullptr;
    };
    {
        const std::lock_guard<std::mutex> lock(mutex);
        if (auto hit = find()) {
            return hit;
        }
    }
    auto map = std::make_shared<const DensityMap>(readMap(path));
    const std::lock_guard<std::mutex> lock(mutex);
    if (auto hit = find()) {  // another thread read it meanwhile
        return hit;
    }
    if (cache.size() >= 2) {
        cache.pop_front();
    }
    cache.emplace_back(key, map);
    return map;
}

double mapLevel(const DensityMap& map, double sigmas) {
    const MapStats st = map.stats();
    return st.mean + sigmas * st.rms;
}

core::Value densityInfo(const fs::path& path) {
    const DensityMap map = readMap(path);
    const MapStats st = map.stats();
    const auto rounded = [](const core::Vec3& v, int digits) {
        return core::Array{core::pythonRound(v[0], digits), core::pythonRound(v[1], digits),
                           core::pythonRound(v[2], digits)};
    };
    return core::Table{
        {"map", absolutePath(path)},
        {"grid_xyz", core::Array{map.shapeXyz[0], map.shapeXyz[1], map.shapeXyz[2]}},
        {"voxel_A", rounded(map.voxel, 4)},
        {"origin_A", rounded(map.origin, 3)},
        {"extent_A", rounded(map.extent(), 3)},
        {"min", core::pythonRound(st.min, 6)},
        {"max", core::pythonRound(st.max, 6)},
        {"mean", core::pythonRound(st.mean, 6)},
        {"rms", core::pythonRound(st.rms, 6)},
        {"level_2sigma", core::pythonRound(mapLevel(map, 2.0), 6)},
    };
}

// ---- simulated maps --------------------------------------------------------------------------

std::vector<double> atomWeights(std::span<const std::string> elements) {
    static const std::pair<const char*, double> table[] = {
        {"H", 1},   {"C", 6},   {"N", 7},   {"O", 8},   {"F", 9},   {"NA", 11},
        {"MG", 12}, {"P", 15},  {"S", 16},  {"CL", 17}, {"K", 19},  {"CA", 20},
        {"MN", 25}, {"FE", 26}, {"CO", 27}, {"NI", 28}, {"CU", 29}, {"ZN", 30},
        {"SE", 34}, {"BR", 35}, {"RB", 37}, {"SR", 38}, {"CS", 55}, {"BA", 56}};
    std::vector<double> out;
    out.reserve(elements.size());
    for (const std::string& element : elements) {
        const std::string key =
            detail::toUtf8(detail::upper(detail::strip(detail::toU32(element))));
        double weight = 6;
        for (const auto& [symbol, z] : table) {
            if (key == symbol) {
                weight = z;
                break;
            }
        }
        out.push_back(weight);
    }
    return out;
}

DensityMap simulateMap(std::span<const core::Vec3> xyz, std::span<const double> weights,
                       double resolution, std::optional<double> spacing, double sigmaFactor,
                       double padSigmas) {
    if (xyz.empty()) {
        throw FitError("no atoms to simulate a map from");
    }
    if (resolution <= 0) {
        throw FitError("the resolution must be positive");
    }
    if (weights.size() != xyz.size()) {
        throw std::invalid_argument("simulateMap: one weight per atom is needed");
    }
    const double sigma = sigmaFactor * resolution;
    const double h = spacing && *spacing != 0 ? *spacing : resolution / 3.0;
    const double pad = padSigmas * sigma + h;
    core::Vec3 lo{}, hi{};
    for (std::size_t k = 0; k < 3; ++k) {  // numpy's min / max: NaN propagates
        double mn = xyz[0][k], mx = xyz[0][k];
        for (const core::Vec3& p : xyz) {
            if (std::isnan(mn)) {
                break;
            }
            if (std::isnan(p[k])) {
                mn = mx = p[k];
            } else {
                mn = p[k] < mn ? p[k] : mn;
                mx = p[k] > mx ? p[k] : mx;
            }
        }
        lo[k] = std::floor((mn - pad) / h) * h;  // on whole voxels (CCP4 NSTART)
        hi[k] = mx;
    }
    std::array<std::int64_t, 3> n{};
    std::uint64_t total = 1;
    for (std::size_t k = 0; k < 3; ++k) {
        const double cells = std::ceil((hi[k] + pad - lo[k]) / h);
        if (!(cells >= 0) || cells > 9.0e15) {
            throw core::ValueError("negative dimensions are not allowed");
        }
        n[k] = static_cast<std::int64_t>(cells) + 1;
        total = saturatingMul(total, static_cast<std::uint64_t>(n[k]));
    }
    if (total > std::numeric_limits<std::size_t>::max() / sizeof(double)) {
        throw std::bad_alloc();
    }
    std::vector<double> grid(static_cast<std::size_t>(total), 0.0);

    // splat_gaussians(grid, (xyz - lo) / h, weights, sigma / h, pad_sigmas)
    const double sigmaVox = sigma / h;
    const double cutoff = padSigmas;
    const auto rad =
        std::max<std::int64_t>(1, static_cast<std::int64_t>(std::ceil(cutoff * sigmaVox)));
    const auto width = static_cast<std::size_t>(2 * rad + 1);
    const std::size_t atoms = xyz.size();
    std::vector<std::array<std::int64_t, 3>> base(atoms);
    std::array<std::vector<double>, 3> factor;  // [axis][atom * width + offset]
    const double norm = 1.0 / (2.0 * sigmaVox * sigmaVox);
    for (std::size_t k = 0; k < 3; ++k) {
        factor[k].resize(atoms * width);
    }
    for (std::size_t a = 0; a < atoms; ++a) {
        for (std::size_t k = 0; k < 3; ++k) {
            const double ijk = (xyz[a][k] - lo[k]) / h;
            const double fl = std::floor(ijk);
            base[a][k] = static_cast<std::int64_t>(fl);
            const double frac = ijk - static_cast<double>(base[a][k]);
            for (std::size_t o = 0; o < width; ++o) {
                const double d = static_cast<double>(static_cast<std::int64_t>(o) - rad) - frac;
                factor[k][a * width + o] = std::exp(-(d * d) * norm);
            }
        }
    }
    const std::size_t nx = static_cast<std::size_t>(n[0]), ny = static_cast<std::size_t>(n[1]);
    for (std::size_t dzi = 0; dzi < width; ++dzi) {
        const std::int64_t dz = static_cast<std::int64_t>(dzi) - rad;
        for (std::size_t dyi = 0; dyi < width; ++dyi) {
            const std::int64_t dy = static_cast<std::int64_t>(dyi) - rad;
            for (std::size_t a = 0; a < atoms; ++a) {
                const std::int64_t z = base[a][2] + dz, y = base[a][1] + dy;
                if (z < 0 || z >= n[2] || y < 0 || y >= n[1]) {
                    continue;
                }
                const double amp =
                    weights[a] * factor[2][a * width + dzi] * factor[1][a * width + dyi];
                const std::size_t row =
                    (static_cast<std::size_t>(z) * ny + static_cast<std::size_t>(y)) * nx;
                for (std::size_t dxi = 0; dxi < width; ++dxi) {
                    const std::int64_t x = base[a][0] + static_cast<std::int64_t>(dxi) - rad;
                    if (x < 0 || x >= n[0]) {
                        continue;
                    }
                    grid[row + static_cast<std::size_t>(x)] += amp * factor[0][a * width + dxi];
                }
            }
        }
    }
    DensityMap map;
    map.shapeXyz = n;
    map.data.resize(grid.size());
    std::transform(grid.begin(), grid.end(), map.data.begin(),
                   [](double v) { return static_cast<float>(v); });
    map.origin = lo;
    map.voxel = {h, h, h};
    return map;
}

}  // namespace mdbridge::io
