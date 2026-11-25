// The analysis cache (mdatk/io/cache.py): .npz archives of the arrays, JSON descriptions, the
// manifest. See cache.hpp for the file layout.

#include "mdbridge/analysis/cache.hpp"

#include <zlib.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstring>
#include <fstream>
#include <limits>
#include <map>
#include <stdexcept>
#include <type_traits>

#include "mdbridge/core/clock.hpp"
#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/log.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/core/text.hpp"
#include "mdbridge/version.hpp"
#ifdef _WIN32
#include "mdbridge/core/ntpath.hpp"
#endif

namespace mdbridge::analysis {

namespace fs = std::filesystem;
using core::NumArray;
using core::Table;
using core::Value;

static_assert(std::endian::native == std::endian::little, "the cache files are written little-endian");

namespace {

[[noreturn]] void fail(const fs::path& path, const std::string& what) {
    throw std::runtime_error(core::pathToUtf8(path) + ": " + what);
}

// ---- little-endian fields ---------------------------------------------------------------------

void put16(std::string& out, std::uint32_t v) {
    out.push_back(static_cast<char>(v & 0xFF));
    out.push_back(static_cast<char>((v >> 8) & 0xFF));
}
void put32(std::string& out, std::uint32_t v) {
    put16(out, v & 0xFFFF);
    put16(out, v >> 16);
}
void put64(std::string& out, std::uint64_t v) {
    put32(out, static_cast<std::uint32_t>(v & 0xFFFFFFFFu));
    put32(out, static_cast<std::uint32_t>(v >> 32));
}
std::uint64_t getLe(const unsigned char* p, int bytes) {
    std::uint64_t v = 0;
    for (int n = bytes - 1; n >= 0; --n) {
        v = (v << 8) | p[n];
    }
    return v;
}

// ---- .npy -------------------------------------------------------------------------------------

std::string npyDescr(const NumArray& a) {
    switch (a.dtype()) {
    case NumArray::DType::Bool: return "|b1";
    case NumArray::DType::Int8: return "|i1";
    case NumArray::DType::Int32: return "<i4";
    case NumArray::DType::Int64: return "<i8";
    case NumArray::DType::Float32: return "<f4";
    case NumArray::DType::Float64: return "<f8";
    case NumArray::DType::String: return "<U" + std::to_string(a.stringWidth());
    }
    throw std::logic_error("npyDescr");
}

// Python's repr of the shape tuple: (3,), (2, 3), ()
std::string shapeRepr(const std::vector<std::size_t>& shape) {
    std::string s = "(";
    for (std::size_t n = 0; n < shape.size(); ++n) {
        s += (n ? ", " : "") + std::to_string(shape[n]);
    }
    return s + (shape.size() == 1 ? ",)" : ")");
}

// The header numpy's format.write_array writes (format 1.0, or 2.0 when it does not fit).
std::string npyHeader(const NumArray& a) {
    const std::string dict =
        "{'descr': '" + npyDescr(a) + "', 'fortran_order': False, 'shape': " + shapeRepr(a.shape()) + ", }";
    for (const std::uint32_t version : {1u, 2u}) {
        const std::size_t prefix = 6 + 2 + (version == 1 ? 2 : 4);
        const std::size_t total = (prefix + dict.size() + 1 + 63) / 64 * 64;
        const std::size_t length = total - prefix;
        if (version == 1 && length > 0xFFFF) {
            continue;
        }
        std::string out("\x93NUMPY", 6);
        out.push_back(static_cast<char>(version));
        out.push_back('\0');
        if (version == 1) {
            put16(out, static_cast<std::uint32_t>(length));
        } else {
            put32(out, static_cast<std::uint32_t>(length));
        }
        out += dict;
        out.append(length - dict.size() - 1, ' ');
        out.push_back('\n');
        return out;
    }
    throw std::logic_error("npyHeader");
}

// Emits the bytes of an array's data (after the header) in pieces.
template <class Sink>
void emitData(const NumArray& a, Sink&& sink) {
    const auto raw = [&](const auto span) {
        const auto* bytes = reinterpret_cast<const unsigned char*>(span.data());
        std::size_t left = span.size_bytes();
        constexpr std::size_t kChunk = std::size_t{1} << 20;
        while (left > 0) {
            const std::size_t n = std::min(left, kChunk);
            sink(bytes, n);
            bytes += n;
            left -= n;
        }
    };
    switch (a.dtype()) {
    case NumArray::DType::Bool: raw(a.values<std::uint8_t>()); return;
    case NumArray::DType::Int8: raw(a.values<std::int8_t>()); return;
    case NumArray::DType::Int32: raw(a.values<std::int32_t>()); return;
    case NumArray::DType::Int64: raw(a.values<std::int64_t>()); return;
    case NumArray::DType::Float32: raw(a.values<float>()); return;
    case NumArray::DType::Float64: raw(a.values<double>()); return;
    case NumArray::DType::String: {
        const std::size_t width = a.stringWidth();
        std::vector<unsigned char> buffer;
        for (const std::string& s : a.values<std::string>()) {
            const std::u32string codes = core::decodeUtf8(s);
            for (std::size_t n = 0; n < width; ++n) {
                const std::uint32_t c = n < codes.size() ? static_cast<std::uint32_t>(codes[n]) : 0;
                for (int b = 0; b < 4; ++b) {
                    buffer.push_back(static_cast<unsigned char>((c >> (8 * b)) & 0xFF));
                }
            }
            if (buffer.size() >= (std::size_t{1} << 20)) {
                sink(buffer.data(), buffer.size());
                buffer.clear();
            }
        }
        if (!buffer.empty()) {
            sink(buffer.data(), buffer.size());
        }
        return;
    }
    }
}

std::size_t dataBytes(const NumArray& a) {
    switch (a.dtype()) {
    case NumArray::DType::Bool:
    case NumArray::DType::Int8: return a.size();
    case NumArray::DType::Int32:
    case NumArray::DType::Float32: return a.size() * 4;
    case NumArray::DType::Int64:
    case NumArray::DType::Float64: return a.size() * 8;
    case NumArray::DType::String: return a.size() * a.stringWidth() * 4;
    }
    return 0;
}

// ---- writing a ZIP archive --------------------------------------------------------------------

constexpr std::uint32_t kMax32 = 0xFFFFFFFFu;
constexpr std::uint16_t kDosTime = 0;  // 00:00:00
constexpr std::uint16_t kDosDate = (0 << 9) | (1 << 5) | 1;  // 1980-01-01, as numpy writes members

struct CentralEntry {
    std::string name;
    std::uint32_t crc = 0;
    std::uint64_t compressed = 0;
    std::uint64_t size = 0;
    std::uint64_t offset = 0;
    bool zip64 = false;
};

bool isAscii(std::string_view s) {
    return std::all_of(s.begin(), s.end(), [](char c) { return static_cast<unsigned char>(c) < 0x80; });
}

class ZipWriter {
public:
    explicit ZipWriter(fs::path path) : m_path(std::move(path)), m_out(m_path, std::ios::binary | std::ios::trunc) {
        if (!m_out) {
            fail(m_path, "cannot write the file");
        }
    }

    void add(const std::string& member, const NumArray& array) {
        const std::string header = npyHeader(array);
        CentralEntry e;
        e.name = member;
        e.size = header.size() + dataBytes(array);
        e.offset = position();
        // deflate may grow incompressible data a little: ZIP64 well before the 32-bit limits
        e.zip64 = e.size >= 0xF0000000u || e.offset >= kMax32;
        const std::string local = localHeader(e);
        write(local);

        z_stream z{};
        if (deflateInit2(&z, Z_DEFAULT_COMPRESSION, Z_DEFLATED, -MAX_WBITS, 8, Z_DEFAULT_STRATEGY) != Z_OK) {
            fail(m_path, "zlib deflateInit2 failed");
        }
        std::vector<unsigned char> buffer(std::size_t{1} << 18);
        uLong crc = crc32(0L, Z_NULL, 0);
        std::uint64_t compressed = 0;
        const auto pump = [&](const unsigned char* data, std::size_t n, int flush) {
            z.next_in = const_cast<Bytef*>(data);
            z.avail_in = static_cast<uInt>(n);
            do {
                z.next_out = buffer.data();
                z.avail_out = static_cast<uInt>(buffer.size());
                const int rc = deflate(&z, flush);
                if (rc == Z_STREAM_ERROR) {
                    deflateEnd(&z);
                    fail(m_path, "zlib deflate failed");
                }
                const std::size_t produced = buffer.size() - z.avail_out;
                m_out.write(reinterpret_cast<const char*>(buffer.data()), static_cast<std::streamsize>(produced));
                compressed += produced;
            } while (z.avail_out == 0);
        };
        const auto sink = [&](const unsigned char* data, std::size_t n) {
            crc = crc32(crc, data, static_cast<uInt>(n));
            pump(data, n, Z_NO_FLUSH);
        };
        sink(reinterpret_cast<const unsigned char*>(header.data()), header.size());
        emitData(array, sink);
        pump(nullptr, 0, Z_FINISH);
        deflateEnd(&z);
        e.crc = static_cast<std::uint32_t>(crc);
        e.compressed = compressed;
        if (!e.zip64 && compressed >= kMax32) {
            fail(m_path, "member " + member + " too large");
        }
        // the sizes and CRC into the local header
        const std::uint64_t end = position();
        m_out.seekp(static_cast<std::streamoff>(e.offset));
        write(localHeader(e));
        m_out.seekp(static_cast<std::streamoff>(end));
        m_entries.push_back(std::move(e));
    }

    void finish() {
        const std::uint64_t cdOffset = position();
        for (const CentralEntry& e : m_entries) {
            std::string h;
            put32(h, 0x02014b50);
            const bool offset64 = e.offset >= kMax32;
            const bool z64 = e.zip64 || offset64;
            put16(h, z64 ? 45 : 20);  // made by (MS-DOS attributes)
            put16(h, z64 ? 45 : 20);  // needed
            put16(h, isAscii(e.name) ? 0 : 0x0800);
            put16(h, 8);
            put16(h, kDosTime);
            put16(h, kDosDate);
            put32(h, e.crc);
            std::string extra;
            if (e.zip64) {
                put32(h, kMax32);
                put32(h, kMax32);
                put64(extra, e.size);
                put64(extra, e.compressed);
            } else {
                put32(h, static_cast<std::uint32_t>(e.compressed));
                put32(h, static_cast<std::uint32_t>(e.size));
            }
            if (offset64) {
                put64(extra, e.offset);
            }
            if (!extra.empty()) {
                std::string field;
                put16(field, 0x0001);
                put16(field, static_cast<std::uint32_t>(extra.size()));
                extra = field + extra;
            }
            put16(h, static_cast<std::uint32_t>(e.name.size()));
            put16(h, static_cast<std::uint32_t>(extra.size()));
            put16(h, 0);  // comment
            put16(h, 0);  // disk
            put16(h, 0);  // internal attributes
            put32(h, 0x01800000u);  // external attributes: -rw------- (numpy's 0o600 << 16)
            put32(h, offset64 ? kMax32 : static_cast<std::uint32_t>(e.offset));
            h += e.name;
            h += extra;
            write(h);
        }
        const std::uint64_t cdSize = position() - cdOffset;
        const std::uint64_t count = m_entries.size();
        std::string tail;
        if (count >= 0xFFFF || cdSize >= kMax32 || cdOffset >= kMax32) {
            const std::uint64_t record = position();
            put32(tail, 0x06064b50);
            put64(tail, 44);
            put16(tail, 45);
            put16(tail, 45);
            put32(tail, 0);
            put32(tail, 0);
            put64(tail, count);
            put64(tail, count);
            put64(tail, cdSize);
            put64(tail, cdOffset);
            put32(tail, 0x07064b50);
            put32(tail, 0);
            put64(tail, record);
            put32(tail, 1);
        }
        put32(tail, 0x06054b50);
        put16(tail, 0);
        put16(tail, 0);
        put16(tail, static_cast<std::uint32_t>(std::min<std::uint64_t>(count, 0xFFFF)));
        put16(tail, static_cast<std::uint32_t>(std::min<std::uint64_t>(count, 0xFFFF)));
        put32(tail, static_cast<std::uint32_t>(std::min<std::uint64_t>(cdSize, kMax32)));
        put32(tail, static_cast<std::uint32_t>(std::min<std::uint64_t>(cdOffset, kMax32)));
        put16(tail, 0);
        write(tail);
        m_out.close();
        if (!m_out) {
            fail(m_path, "cannot write the file");
        }
    }

private:
    std::uint64_t position() { return static_cast<std::uint64_t>(m_out.tellp()); }

    void write(const std::string& bytes) {
        m_out.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
        if (!m_out) {
            fail(m_path, "cannot write the file");
        }
    }

    static std::string localHeader(const CentralEntry& e) {
        std::string h;
        put32(h, 0x04034b50);
        put16(h, e.zip64 ? 45 : 20);
        put16(h, isAscii(e.name) ? 0 : 0x0800);
        put16(h, 8);
        put16(h, kDosTime);
        put16(h, kDosDate);
        put32(h, e.crc);
        std::string extra;
        if (e.zip64) {
            put32(h, kMax32);
            put32(h, kMax32);
            put16(extra, 0x0001);
            put16(extra, 16);
            put64(extra, e.size);
            put64(extra, e.compressed);
        } else {
            put32(h, static_cast<std::uint32_t>(e.compressed));
            put32(h, static_cast<std::uint32_t>(e.size));
        }
        put16(h, static_cast<std::uint32_t>(e.name.size()));
        put16(h, static_cast<std::uint32_t>(extra.size()));
        return h + e.name + extra;
    }

    fs::path m_path;
    std::ofstream m_out;
    std::vector<CentralEntry> m_entries;
};

// ---- reading a ZIP archive --------------------------------------------------------------------

struct ZipMember {
    std::string name;
    std::uint16_t method = 0;
    std::uint32_t crc = 0;
    std::uint64_t compressed = 0;
    std::uint64_t size = 0;
    std::uint64_t offset = 0;
};

class ZipReader {
public:
    explicit ZipReader(fs::path path) : m_path(std::move(path)), m_in(m_path, std::ios::binary) {
        if (!m_in) {
            fail(m_path, "cannot open the file");
        }
        readDirectory();
    }

    const std::vector<ZipMember>& members() const { return m_members; }

    // Streams one member's uncompressed bytes.
    class Stream {
    public:
        Stream(ZipReader& zip, const ZipMember& m) : m_zip(zip), m_member(m) {
            unsigned char h[30];
            zip.readAt(m.offset, h, sizeof h);
            if (getLe(h, 4) != 0x04034b50) {
                zip.bad("bad local header of " + m.name);
            }
            m_pos = m.offset + 30 + getLe(h + 26, 2) + getLe(h + 28, 2);
            m_left = m.compressed;
            if (m.method == 8) {
                if (inflateInit2(&m_z, -MAX_WBITS) != Z_OK) {
                    zip.bad("zlib inflateInit2 failed");
                }
                m_inflating = true;
            } else if (m.method != 0) {
                zip.bad("unsupported compression of " + m.name);
            }
            m_crc = crc32(0L, Z_NULL, 0);
        }
        ~Stream() {
            if (m_inflating) {
                inflateEnd(&m_z);
            }
        }
        Stream(const Stream&) = delete;
        Stream& operator=(const Stream&) = delete;

        // Exactly n bytes, or an error.
        void read(void* dst, std::size_t n) {
            auto* out = static_cast<unsigned char*>(dst);
            std::size_t done = 0;
            while (done < n) {
                const std::size_t want = std::min<std::size_t>(n - done, std::numeric_limits<uInt>::max() / 2);
                std::size_t got = 0;
                if (!m_inflating) {
                    if (m_left < want) {
                        m_zip.bad("member " + m_member.name + " is truncated");
                    }
                    m_zip.readAt(m_pos, out + done, want);
                    m_pos += want;
                    m_left -= want;
                    got = want;
                } else {
                    if (m_z.avail_in == 0 && m_left > 0) {
                        const auto chunk =
                            static_cast<std::size_t>(std::min<std::uint64_t>(m_left, m_buffer.size()));
                        m_zip.readAt(m_pos, m_buffer.data(), chunk);
                        m_pos += chunk;
                        m_left -= chunk;
                        m_z.next_in = m_buffer.data();
                        m_z.avail_in = static_cast<uInt>(chunk);
                    }
                    m_z.next_out = out + done;
                    m_z.avail_out = static_cast<uInt>(want);
                    const int rc = inflate(&m_z, Z_NO_FLUSH);
                    got = want - m_z.avail_out;
                    if (rc == Z_STREAM_END && got < want && done + got < n) {
                        m_zip.bad("member " + m_member.name + " is truncated");
                    }
                    if (rc != Z_OK && rc != Z_STREAM_END && !(rc == Z_BUF_ERROR && got > 0)) {
                        m_zip.bad("member " + m_member.name + " is corrupt (zlib error " + std::to_string(rc) + ")");
                    }
                    if (got == 0 && m_z.avail_in == 0 && m_left == 0) {
                        m_zip.bad("member " + m_member.name + " is truncated");
                    }
                }
                m_crc = crc32(m_crc, out + done, static_cast<uInt>(got));
                m_total += got;
                done += got;
            }
        }

        // After the last read: every byte was read and the CRC matches.
        void finish() {
            if (m_total != m_member.size || static_cast<std::uint32_t>(m_crc) != m_member.crc) {
                m_zip.bad("member " + m_member.name + " is corrupt (size or CRC)");
            }
        }

    private:
        ZipReader& m_zip;
        const ZipMember& m_member;
        z_stream m_z{};
        bool m_inflating = false;
        std::uint64_t m_pos = 0;
        std::uint64_t m_left = 0;
        std::uint64_t m_total = 0;
        uLong m_crc = 0;
        std::vector<unsigned char> m_buffer = std::vector<unsigned char>(std::size_t{1} << 18);
    };

    [[noreturn]] void bad(const std::string& what) const { fail(m_path, what); }

    void readAt(std::uint64_t offset, void* dst, std::size_t n) {
        m_in.clear();
        m_in.seekg(static_cast<std::streamoff>(offset));
        m_in.read(static_cast<char*>(dst), static_cast<std::streamsize>(n));
        if (static_cast<std::size_t>(m_in.gcount()) != n) {
            bad("unexpected end of file");
        }
    }

private:
    void readDirectory() {
        m_in.seekg(0, std::ios::end);
        const auto fileSize = static_cast<std::uint64_t>(m_in.tellg());
        const std::uint64_t tailSize = std::min<std::uint64_t>(fileSize, 22 + 0xFFFF);
        std::vector<unsigned char> tail(static_cast<std::size_t>(tailSize));
        readAt(fileSize - tailSize, tail.data(), tail.size());
        std::ptrdiff_t eocd = -1;
        for (std::ptrdiff_t p = static_cast<std::ptrdiff_t>(tail.size()) - 22; p >= 0; --p) {
            if (getLe(&tail[static_cast<std::size_t>(p)], 4) == 0x06054b50) {
                eocd = p;
                break;
            }
        }
        if (eocd < 0) {
            bad("not a ZIP archive (.npz)");
        }
        const unsigned char* e = &tail[static_cast<std::size_t>(eocd)];
        std::uint64_t count = getLe(e + 10, 2);
        std::uint64_t cdSize = getLe(e + 12, 4);
        std::uint64_t cdOffset = getLe(e + 16, 4);
        const std::uint64_t eocdPos = fileSize - tailSize + static_cast<std::uint64_t>(eocd);
        if ((count == 0xFFFF || cdSize == kMax32 || cdOffset == kMax32) && eocdPos >= 20) {
            unsigned char loc[20];
            readAt(eocdPos - 20, loc, sizeof loc);
            if (getLe(loc, 4) == 0x07064b50) {
                unsigned char rec[56];
                readAt(getLe(loc + 8, 8), rec, sizeof rec);
                if (getLe(rec, 4) != 0x06064b50) {
                    bad("bad ZIP64 end record");
                }
                count = getLe(rec + 32, 8);
                cdSize = getLe(rec + 40, 8);
                cdOffset = getLe(rec + 48, 8);
            }
        }
        if (cdOffset + cdSize > fileSize) {
            bad("bad central directory");
        }
        std::vector<unsigned char> cd(static_cast<std::size_t>(cdSize));
        readAt(cdOffset, cd.data(), cd.size());
        std::size_t p = 0;
        for (std::uint64_t n = 0; n < count; ++n) {
            if (p + 46 > cd.size() || getLe(&cd[p], 4) != 0x02014b50) {
                bad("bad central directory");
            }
            const unsigned char* h = &cd[p];
            ZipMember m;
            m.method = static_cast<std::uint16_t>(getLe(h + 10, 2));
            m.crc = static_cast<std::uint32_t>(getLe(h + 16, 4));
            m.compressed = getLe(h + 20, 4);
            m.size = getLe(h + 24, 4);
            const std::size_t nameLen = getLe(h + 28, 2);
            const std::size_t extraLen = getLe(h + 30, 2);
            const std::size_t commentLen = getLe(h + 32, 2);
            m.offset = getLe(h + 42, 4);
            if (p + 46 + nameLen + extraLen + commentLen > cd.size()) {
                bad("bad central directory");
            }
            m.name.assign(reinterpret_cast<const char*>(h + 46), nameLen);
            // ZIP64 extra: the fields whose 32-bit value is 0xFFFFFFFF, in this order
            const unsigned char* x = h + 46 + nameLen;
            for (std::size_t q = 0; q + 4 <= extraLen;) {
                const std::uint64_t id = getLe(x + q, 2);
                const std::size_t len = getLe(x + q + 2, 2);
                if (q + 4 + len > extraLen) {
                    break;
                }
                if (id == 0x0001) {
                    std::size_t r = q + 4;
                    const auto take = [&](std::uint64_t& field) {
                        if (field == kMax32 && r + 8 <= q + 4 + len) {
                            field = getLe(x + r, 8);
                            r += 8;
                        }
                    };
                    take(m.size);
                    take(m.compressed);
                    take(m.offset);
                }
                q += 4 + len;
            }
            m_members.push_back(std::move(m));
            p += 46 + nameLen + extraLen + commentLen;
        }
    }

    fs::path m_path;
    std::ifstream m_in;
    std::vector<ZipMember> m_members;
};

// The value of 'key' in a .npy header dict (the text after "'key':", stripped).
std::string headerField(const std::string& header, const std::string& key) {
    const std::size_t at = header.find("'" + key + "'");
    if (at == std::string::npos) {
        return {};
    }
    std::size_t p = header.find(':', at);
    if (p == std::string::npos) {
        return {};
    }
    ++p;
    while (p < header.size() && header[p] == ' ') {
        ++p;
    }
    if (p < header.size() && header[p] == '\'') {
        const std::size_t end = header.find('\'', p + 1);
        return end == std::string::npos ? std::string() : header.substr(p + 1, end - p - 1);
    }
    if (p < header.size() && header[p] == '(') {
        const std::size_t end = header.find(')', p);
        return end == std::string::npos ? std::string() : header.substr(p, end - p + 1);
    }
    std::size_t end = p;
    while (end < header.size() && header[end] != ',' && header[end] != '}') {
        ++end;
    }
    return header.substr(p, end - p);
}

template <class T>
NumArray readValues(ZipReader::Stream& s, std::size_t count, std::vector<std::size_t> shape) {
    std::vector<T> values(count);
    if (count > 0) {
        s.read(values.data(), count * sizeof(T));
    }
    return NumArray(std::move(values), std::move(shape));
}

NumArray readNpy(ZipReader& zip, const ZipMember& m) {
    ZipReader::Stream s(zip, m);
    unsigned char magic[8];
    s.read(magic, sizeof magic);
    if (std::memcmp(magic, "\x93NUMPY", 6) != 0) {
        zip.bad("member " + m.name + " is not a .npy array");
    }
    std::size_t headerLen = 0;
    if (magic[6] == 1) {
        unsigned char b[2];
        s.read(b, 2);
        headerLen = getLe(b, 2);
    } else {
        unsigned char b[4];
        s.read(b, 4);
        headerLen = getLe(b, 4);
    }
    std::string header(headerLen, '\0');
    s.read(header.data(), headerLen);
    const std::string descr = headerField(header, "descr");
    const std::string fortran = headerField(header, "fortran_order");
    const std::string shapeText = headerField(header, "shape");
    std::vector<std::size_t> shape;
    std::size_t count = 1;
    for (std::size_t p = 1; p < shapeText.size();) {
        while (p < shapeText.size() && (shapeText[p] == ' ' || shapeText[p] == ',')) {
            ++p;
        }
        std::size_t end = p;
        while (end < shapeText.size() && shapeText[end] >= '0' && shapeText[end] <= '9') {
            ++end;
        }
        if (end == p) {
            break;
        }
        shape.push_back(static_cast<std::size_t>(std::stoull(shapeText.substr(p, end - p))));
        count *= shape.back();
        p = end;
    }
    if (shapeText.empty() || shapeText.front() != '(' || (fortran == "True" && shape.size() > 1)) {
        zip.bad("member " + m.name + ": unsupported .npy header " + header);
    }
    if (shape.empty()) {
        shape = {1};  // a 0-d array as one value
    }
    NumArray out;
    if (descr == "|b1") {
        std::vector<std::uint8_t> bytes(count);
        if (count > 0) {
            s.read(bytes.data(), count);
        }
        out = NumArray(std::vector<bool>(bytes.begin(), bytes.end()), shape);
    } else if (descr == "|i1") {
        out = readValues<std::int8_t>(s, count, shape);
    } else if (descr == "<i4") {
        out = readValues<std::int32_t>(s, count, shape);
    } else if (descr == "<i8") {
        out = readValues<std::int64_t>(s, count, shape);
    } else if (descr == "<f4") {
        out = readValues<float>(s, count, shape);
    } else if (descr == "<f8") {
        out = readValues<double>(s, count, shape);
    } else if (descr.size() > 2 && descr.compare(0, 2, "<U") == 0) {
        const std::size_t width = static_cast<std::size_t>(std::stoull(descr.substr(2)));
        std::vector<std::string> values(count);
        std::vector<unsigned char> buffer(width * 4);
        for (std::string& v : values) {
            if (width > 0) {
                s.read(buffer.data(), buffer.size());
            }
            std::u32string codes;
            for (std::size_t n = 0; n < width; ++n) {
                codes.push_back(static_cast<char32_t>(getLe(&buffer[n * 4], 4)));
            }
            while (!codes.empty() && codes.back() == 0) {
                codes.pop_back();
            }
            v = core::encodeUtf8(codes);
        }
        out = NumArray::strings(std::move(values), std::max<std::size_t>(width, 1), shape);
    } else {
        zip.bad("member " + m.name + ": unsupported dtype '" + descr + "'");
    }
    s.finish();
    return out;
}

// ---- JSON files -------------------------------------------------------------------------------

std::string readFile(const fs::path& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        fail(path, "cannot read the file");
    }
    std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    return text;
}

// write_json: indent 2, a final newline, <path>.tmp then rename.
void writeJson(const std::string& path, const Value& data) {
    const fs::path target = core::pathFromUtf8(path);
    const fs::path tmp = core::pathFromUtf8(path + ".tmp");
    {
        std::ofstream out(tmp, std::ios::binary | std::ios::trunc);
        const std::string text = core::toJson(data, 2) + "\n";
        out.write(text.data(), static_cast<std::streamsize>(text.size()));
        out.close();
        if (!out) {
            fail(tmp, "cannot write the file");
        }
    }
    fs::rename(tmp, target);
}

const Value* entryOf(const Value& table, std::string_view key) { return table.isTable() ? table.find(key) : nullptr; }

std::string absolute(const std::string& dir) {
#ifdef _WIN32
    return core::ntpath::abspath(dir, core::currentDirectory());
#else
    return core::posixpath::abspath(dir, core::currentDirectory());
#endif
}

// os.path.join: the folders are absolute in generic form ('/' separators) on every platform and the
// names relative, so '/' joins them.
std::string joinPath(const std::string& dir, std::string_view name) { return core::posixpath::join(dir, name); }

// float(np.median(np.diff(t)))
double medianStep(const std::vector<double>& t) {
    std::vector<double> d;
    for (std::size_t n = 1; n < t.size(); ++n) {
        d.push_back(t[n] - t[n - 1]);
    }
    if (std::any_of(d.begin(), d.end(), [](double x) { return std::isnan(x); })) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    std::sort(d.begin(), d.end());
    const std::size_t m = d.size() / 2;
    return d.size() % 2 ? d[m] : (d[m - 1] + d[m]) / 2.0;
}

}  // namespace

// ---- .npz ---------------------------------------------------------------------------------------

void writeNpz(const fs::path& path, const core::NamedArrays& arrays) {
    ZipWriter zip(path);
    for (const auto& [name, array] : arrays) {
        zip.add(name + ".npy", array);
    }
    zip.finish();
}

core::NamedArrays readNpz(const fs::path& path) {
    ZipReader zip(path);
    core::NamedArrays out;
    for (const ZipMember& m : zip.members()) {
        std::string name = m.name;
        if (name.size() >= 4 && name.compare(name.size() - 4, 4, ".npy") == 0) {
            name.resize(name.size() - 4);
        }
        out.emplace_back(std::move(name), readNpy(zip, m));
    }
    return out;
}

// ---- CacheStore ----------------------------------------------------------------------------------

CacheStore::CacheStore(std::string cacheDir) : m_dir(absolute(cacheDir)) {}

std::string CacheStore::manifestPath() const { return joinPath(m_dir, kManifestName); }

std::string CacheStore::entryFile(std::string_view name, std::string_view key) {
    return std::string(name) + "-" + std::string(key) + std::string(kEntryTag) + ".npz";
}

std::string CacheStore::entryBase(std::string_view name, std::string_view key) const {
    return joinPath(m_dir, std::string(name) + "-" + std::string(key) + std::string(kEntryTag));
}

void CacheStore::ensure() const { fs::create_directories(core::pathFromUtf8(m_dir)); }

bool CacheStore::has(std::string_view name, std::string_view key) const {
    const std::string base = entryBase(name, key);
    std::error_code ec;
    return fs::exists(core::pathFromUtf8(base + ".npz"), ec) && fs::exists(core::pathFromUtf8(base + ".json"), ec);
}

std::string CacheStore::saveResult(const core::AnalysisResult& result, std::span<const std::int64_t> frames,
                                   std::span<const double> timePs) const {
    ensure();
    const std::string base = entryBase(result.name, result.key);
    core::NamedArrays arrays = result.arrays;
    arrays.emplace_back("_frames", NumArray(std::vector<std::int64_t>(frames.begin(), frames.end())));
    arrays.emplace_back("_time_ps", NumArray(std::vector<double>(timePs.begin(), timePs.end())));
    const fs::path tmp = core::pathFromUtf8(base + ".tmp.npz");
    writeNpz(tmp, arrays);
    fs::rename(tmp, core::pathFromUtf8(base + ".npz"));
    std::vector<std::string> names;
    for (const auto& [name, array] : result.arrays) {
        names.push_back(name);
    }
    std::sort(names.begin(), names.end());
    core::Array sorted(names.begin(), names.end());
    writeJson(base + ".json", Table{{"project", kName},
                                    {"mdbridge_version", kVersion},
                                    {"created", core::clockNow().isoformat()},
                                    {"name", result.name},
                                    {"key", result.key},
                                    {"params", result.params},
                                    {"meta", result.meta},
                                    {"arrays", sorted},
                                    {"n_frames", static_cast<std::int64_t>(frames.size())}});
    return base + ".npz";
}

std::optional<CachedResult> CacheStore::loadResult(std::string_view name, std::string_view key) const {
    if (!has(name, key)) {
        return std::nullopt;
    }
    const std::string base = entryBase(name, key);
    try {
        const Value info = core::parseJson(readFile(core::pathFromUtf8(base + ".json")), base + ".json");
        core::NamedArrays arrays = readNpz(core::pathFromUtf8(base + ".npz"));
        CachedResult out;
        bool haveFrames = false;
        bool haveTimes = false;
        for (auto& [arrayName, array] : arrays) {
            if (arrayName == "_frames") {
                const auto v = array.values<std::int64_t>();
                out.frames.assign(v.begin(), v.end());
                haveFrames = true;
            } else if (arrayName == "_time_ps") {
                const auto v = array.values<double>();
                out.timePs.assign(v.begin(), v.end());
                haveTimes = true;
            } else {
                out.result.arrays.emplace_back(arrayName, std::move(array));
            }
        }
        if (!haveFrames || !haveTimes) {
            throw std::runtime_error("no _frames / _time_ps arrays");
        }
        const Value* infoName = entryOf(info, "name");
        const Value* infoKey = entryOf(info, "key");
        if (!infoName || !infoName->isString() || !infoKey || !infoKey->isString()) {
            throw std::runtime_error("no name / key in the description");
        }
        out.result.name = *infoName->get<std::string>();
        out.result.key = *infoKey->get<std::string>();
        if (const Value* params = entryOf(info, "params"); params && params->isTable()) {
            out.result.params = *params->get<Table>();
        }
        if (const Value* meta = entryOf(info, "meta"); meta && meta->isTable()) {
            out.result.meta = *meta->get<Table>();
        }
        return out;
    } catch (const std::exception& e) {  // a corrupt / partial cache entry
        core::logWarning("ignoring unreadable cache entry " + base + ": " + e.what());
        return std::nullopt;
    }
}

std::optional<Value> CacheStore::readManifest() const {
    const fs::path path = core::pathFromUtf8(manifestPath());
    std::error_code ec;
    if (!fs::exists(path, ec)) {
        return std::nullopt;
    }
    return core::parseJson(readFile(path), manifestPath());
}

std::string CacheStore::writeManifest(const Value& system, const Table& frameSelection, const core::Dataset& dataset,
                                      const core::Array& groups, const CacheEntries& entries,
                                      const Table& extra) const {
    ensure();
    std::optional<Value> old;
    try {
        old = readManifest();
    } catch (const std::exception& e) {
        core::logWarning("ignoring the unreadable cache manifest " + manifestPath() + ": " + e.what());
    }
    // all_entries by "<name>-<key>"
    std::map<std::string, std::pair<std::pair<std::string, std::string>, Value>> known;
    const auto remember = [&](const std::string& name, const std::string& key, Value entry) {
        known[name + "-" + key] = {{name, key}, std::move(entry)};
    };
    if (old) {
        if (const Value* all = entryOf(*old, "all_entries"); all && all->isArray()) {
            for (const Value& e : *all->get<core::Array>()) {
                const Value* n = entryOf(e, "name");
                const Value* k = entryOf(e, "key");
                if (n && n->isString() && k && k->isString()) {
                    remember(*n->get<std::string>(), *k->get<std::string>(), e);
                }
            }
        }
    }
    core::Array current;
    std::vector<std::string> seen;
    for (const auto& [name, key] : entries) {
        if (std::find(seen.begin(), seen.end(), name) != seen.end()) {
            continue;  // a dict: one entry per name
        }
        seen.push_back(name);
        const std::string& entryKey = key;
        const core::AnalysisResult* res = dataset.result(name);
        Table outputs;
        if (res) {
            const Value* described = nullptr;
            for (const auto& [k, v] : res->meta) {
                if (k == "outputs") {
                    described = &v;
                }
            }
            if (described && described->isTable()) {
                for (const auto& [kind, outs] : *described->get<Table>()) {
                    std::vector<std::string> names;
                    if (const Table* t = outs.get<Table>()) {
                        for (const auto& [out, d] : *t) {
                            names.push_back(out);
                        }
                    }
                    std::sort(names.begin(), names.end());
                    outputs.emplace_back(kind, core::Array(names.begin(), names.end()));
                }
            } else if (!described) {
                outputs = {{"series", core::Array{}}, {"matrices", core::Array{}}, {"pairmaps", core::Array{}}};
            }
        }
        const std::string file = entryFile(name, entryKey);
        current.push_back(Table{{"name", name},
                                {"key", entryKey},
                                {"file", file},
                                {"params", res ? Value(res->params) : Value(Table{})},
                                {"outputs", outputs}});
        remember(name, entryKey, Table{{"name", name}, {"key", entryKey}, {"file", file}});
    }
    std::vector<std::pair<std::pair<std::string, std::string>, Value>> all;
    for (auto& [id, entry] : known) {
        all.push_back(std::move(entry));
    }
    std::stable_sort(all.begin(), all.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
    core::Array allEntries;
    for (auto& [id, entry] : all) {
        allEntries.push_back(std::move(entry));
    }

    const std::size_t n = dataset.frameCount();
    Table frames = frameSelection;
    const auto set = [](Table& t, const std::string& key, Value v) {
        for (auto& [k, x] : t) {
            if (k == key) {
                x = std::move(v);
                return;
            }
        }
        t.emplace_back(key, std::move(v));
    };
    set(frames, "n_frames", static_cast<std::int64_t>(n));
    set(frames, "first_frame", n ? Value(dataset.frames().front()) : Value());
    set(frames, "last_frame", n ? Value(dataset.frames().back()) : Value());
    const std::vector<double>& t = dataset.timePs();
    Table timeAxis{{"unit", "ps"},
                   {"n", static_cast<std::int64_t>(t.size())},
                   {"start", t.empty() ? Value() : Value(t.front())},
                   {"end", t.empty() ? Value() : Value(t.back())},
                   {"dt", t.size() > 1 ? Value(medianStep(t)) : Value()}};
    Table manifest{{"project", kName},
                   {"mdbridge_version", kVersion},
                   {"updated", core::clockNow().isoformat()},
                   {"system", system},
                   {"frames", frames},
                   {"time_axis", timeAxis},
                   {"groups", groups},
                   {"available", dataset.describe()},
                   {"entries", current},
                   {"all_entries", allEntries}};
    for (const auto& [k, v] : extra) {
        set(manifest, k, v);
    }
    writeJson(manifestPath(), manifest);
    return manifestPath();
}

core::Dataset CacheStore::loadDatasetFromManifest(const std::string& timeUnit) const {
    const std::optional<Value> manifest = readManifest();
    if (!manifest) {
        throw core::NamedError("FileNotFoundError", "no " + std::string(kManifestName) + " in " + m_dir);
    }
    std::optional<core::Dataset> dataset;
    if (const Value* entries = entryOf(*manifest, "entries"); entries && entries->isArray()) {
        for (const Value& entry : *entries->get<core::Array>()) {
            const std::string name = entryOf(entry, "name") ? entryOf(entry, "name")->stringOr() : std::string();
            const std::string key = entryOf(entry, "key") ? entryOf(entry, "key")->stringOr() : std::string();
            std::optional<CachedResult> loaded = loadResult(name, key);
            if (!loaded) {
                const Value* file = entryOf(entry, "file");
                throw core::NamedError("FileNotFoundError", "cache entry " + (file ? core::pythonStr(*file) : name) +
                                                                " is missing");
            }
            if (!dataset) {
                dataset.emplace(loaded->frames, loaded->timePs, std::vector<core::AnalysisResult>{}, timeUnit);
            }
            dataset->add(std::move(loaded->result));
        }
    }
    if (!dataset) {
        throw core::NamedError("FileNotFoundError", "the cache manifest lists no analyses");
    }
    return std::move(*dataset);
}

}  // namespace mdbridge::analysis
