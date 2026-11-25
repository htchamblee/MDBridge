// GROMACS XTC and TRR trajectories as MDAnalysis 2.11.0's XTCReader / TRRReader read them.
//
// The XTC coordinate decompression (sizeOfInt, sizeOfInts, decodeBits, decodeInts, kMagicInts and
// the loop of decompressXtc) is ported from chemfiles (src/files/XDRFile.cpp, commit 920a1347),
// whose code derives from the xdrfile library. Both are used under their licences:
//
// chemfiles: Copyright (c) 2020, Guillaume Fraux and contributors (see the AUTHORS file)
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted
// provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this list of
//    conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice, this list of
//    conditions and the following disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its contributors may be used to
//    endorse or promote products derived from this software without specific prior written
//    permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
// IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
// FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// xdrfile: Copyright (c) 2009-2014, Erik Lindahl & David van der Spoel
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted
// provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this list of
//    conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice, this list of
//    conditions and the following disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
// IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
// FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#include <algorithm>
#include <array>
#include <climits>
#include <cstdint>
#include <stdexcept>
#include <utility>

#include "binary_file.hpp"
#include "mdbridge/core/errors.hpp"
#include "readers.hpp"
#include "text_file.hpp"

namespace mdbridge::io::detail {

namespace {

using core::ValueError;

// ---------------------------------------------------------------- xtc decompression (ported)

std::uint32_t sizeOfInt(std::uint32_t size) {
    std::uint32_t num = 1;
    std::uint32_t bits = 0;
    while (size >= num && bits < 4 * CHAR_BIT) {
        ++bits;
        num <<= 1;
    }
    return bits;
}

std::uint32_t sizeOfInts(const std::uint32_t sizes[3]) {
    std::uint32_t byteCount = 1;
    std::uint8_t bytes[32];
    bytes[0] = 1;
    std::uint32_t bits = 0;
    for (std::size_t i = 0; i < 3; ++i) {
        std::uint32_t tmp = 0;
        std::uint32_t k = 0;
        for (; k < byteCount; ++k) {
            tmp = bytes[k] * sizes[i] + tmp;
            bytes[k] = static_cast<std::uint8_t>(tmp & 0xff);
            tmp >>= CHAR_BIT;
        }
        while (tmp != 0) {
            bytes[k++] = static_cast<std::uint8_t>(tmp & 0xff);
            tmp >>= CHAR_BIT;
        }
        byteCount = k;
    }
    std::uint32_t num = 1;
    --byteCount;
    while (bytes[byteCount] >= num) {
        ++bits;
        num *= 2;
    }
    return bits + byteCount * CHAR_BIT;
}

struct BitReader {
    const std::string& buffer;
    std::size_t count = 0;
    std::size_t lastBits = 0;
    std::uint32_t lastByte = 0;

    std::uint8_t byte() {
        if (count >= buffer.size()) {
            throw ReaderError("OSError", "XTC read error = compression");
        }
        return static_cast<std::uint8_t>(buffer[count++]);
    }
};

// receivebits: the next `bits` (< 32) bits of the buffer as an integer
std::uint32_t decodeBits(BitReader& state, std::uint32_t bits) {
    const std::uint32_t mask = bits >= 32 ? 0xFFFFFFFFu : (std::uint32_t{1} << bits) - 1;
    std::size_t lastBits = state.lastBits;
    std::uint32_t lastByte = state.lastByte;
    std::uint32_t num = 0;
    while (bits >= CHAR_BIT) {
        lastByte = (lastByte << CHAR_BIT) | state.byte();
        num |= (lastByte >> lastBits) << (bits - CHAR_BIT);
        bits -= CHAR_BIT;
    }
    if (bits > 0) {
        if (lastBits < bits) {
            lastBits += CHAR_BIT;
            lastByte = (lastByte << CHAR_BIT) | state.byte();
        }
        lastBits -= bits;
        num |= (lastByte >> lastBits) & ((std::uint32_t{1} << bits) - 1);
    }
    num &= mask;
    state.lastBits = lastBits;
    state.lastByte = lastByte & 0xff;
    return num;
}

template <typename T>
void unpackFromInt(BitReader& state, std::uint32_t bits, const std::uint32_t sizes[3],
                   std::int32_t nums[3]) {
    T v = 0;
    std::size_t byteCount = 0;
    while (bits >= CHAR_BIT) {
        const T b = decodeBits(state, CHAR_BIT);
        v |= b << (CHAR_BIT * byteCount++);
        bits -= CHAR_BIT;
    }
    if (bits > 0) {
        const T b = decodeBits(state, bits);
        v |= b << (CHAR_BIT * byteCount);
    }
    const T sz = sizes[2];
    const T sy = sizes[1];
    const T szy = sz * sy;
    const T x1 = v / szy;
    const T q1 = v - x1 * szy;
    const T y1 = q1 / sz;
    const T z1 = q1 - y1 * sz;
    nums[0] = static_cast<std::int32_t>(x1);
    nums[1] = static_cast<std::int32_t>(y1);
    nums[2] = static_cast<std::int32_t>(z1);
}

// receiveints: three small integers packed with the given sizes in `bits` bits
void decodeInts(BitReader& state, std::uint32_t bits, const std::uint32_t sizes[3],
                std::int32_t nums[3]) {
    if (sizes[0] == 0 || sizes[1] == 0 || sizes[2] == 0) {
        throw ReaderError("OSError", "XTC read error = compression");
    }
    if (bits <= 32) {
        unpackFromInt<std::uint32_t>(state, bits, sizes, nums);
        return;
    }
    if (bits <= 64) {
        unpackFromInt<std::uint64_t>(state, bits, sizes, nums);
        return;
    }
    std::uint8_t bytes[32] = {};
    std::size_t byteCount = 0;
    while (bits >= CHAR_BIT) {
        bytes[byteCount++] = static_cast<std::uint8_t>(decodeBits(state, CHAR_BIT));
        bits -= CHAR_BIT;
    }
    if (bits > 0) {
        bytes[byteCount++] = static_cast<std::uint8_t>(decodeBits(state, bits));
    }
    for (std::size_t i = 2; i > 0; --i) {
        std::uint32_t num = 0;
        for (std::size_t j = 0; j < byteCount; ++j) {
            const std::size_t k = byteCount - 1 - j;
            num = (num << CHAR_BIT) | bytes[k];
            const std::uint32_t p = num / sizes[i];
            bytes[k] = static_cast<std::uint8_t>(p);
            num = num - p * sizes[i];
        }
        nums[i] = static_cast<std::int32_t>(num);
    }
    nums[0] =
        static_cast<std::int32_t>(std::uint32_t{bytes[0]} | std::uint32_t{bytes[1]} << 8 |
                                  std::uint32_t{bytes[2]} << 16 | std::uint32_t{bytes[3]} << 24);
}

constexpr std::array<int, 73> kMagicInts = {
    0,        0,        0,       0,       0,       0,       0,       0,       0,       8,
    10,       12,       16,      20,      25,      32,      40,      50,      64,      80,
    101,      128,      161,     203,     256,     322,     406,     512,     645,     812,
    1024,     1290,     1625,    2048,    2580,    3250,    4096,    5060,    6501,    8192,
    10321,    13003,    16384,   20642,   26007,   32768,   41285,   52015,   65536,   82570,
    104031,   131072,   165140,  208063,  262144,  330280,  416127,  524287,  660561,  832255,
    1048576,  1321122,  1664510, 2097152, 2642245, 3329021, 4194304, 5284491, 6658042, 8388607,
    10568983, 13316085, 16777216};
constexpr std::uint32_t kFirstIndex = 9;

// Writers can step the small-integer size one past the table (xdrfile then reads beyond its
// array); such a size is only an error when a run of small integers uses it.
int magicInt(std::uint32_t index) {
    return index < kMagicInts.size() ? kMagicInts[index] : 0;
}

struct Compressed {
    float precision = 0.0f;
    std::int32_t minInt[3] = {};
    std::int32_t maxInt[3] = {};
    std::uint32_t smallIndex = 0;
    std::string data;
};

// The coordinates (nm, float32) of `atoms` atoms: decoded integers times 1 / precision.
std::vector<float> decompressXtc(const Compressed& c, std::size_t atoms) {
    if (c.smallIndex >= kMagicInts.size()) {
        throw ReaderError("OSError", "XTC read error = compression");
    }
    std::uint32_t sizeInt[3];
    std::uint32_t bitSizeInt[3] = {0, 0, 0};
    for (int i = 0; i < 3; ++i) {
        sizeInt[i] = static_cast<std::uint32_t>(c.maxInt[i] - c.minInt[i]) + 1;
    }
    std::uint32_t bitSize = 0;
    if ((sizeInt[0] | sizeInt[1] | sizeInt[2]) > 0xffffff) {
        for (int i = 0; i < 3; ++i) {
            bitSizeInt[i] = sizeOfInt(sizeInt[i]);
        }
    } else {
        bitSize = sizeOfInts(sizeInt);
    }
    std::uint32_t smallIndex = c.smallIndex;
    int smaller = magicInt(std::max(kFirstIndex, smallIndex - 1)) / 2;
    int smallNum = magicInt(smallIndex) / 2;
    std::uint32_t sizeSmall[3];
    sizeSmall[0] = sizeSmall[1] = sizeSmall[2] = static_cast<std::uint32_t>(magicInt(smallIndex));

    const float inversePrecision = static_cast<float>(1.0 / static_cast<double>(c.precision));
    std::vector<float> out(3 * atoms);
    BitReader state{c.data};
    int run = 0;
    std::size_t written = 0;
    std::int32_t thisCoord[3];
    std::int32_t prevCoord[3];
    const auto put = [&](const std::int32_t xyz[3]) {
        if (written >= atoms) {
            throw ReaderError("OSError", "XTC read error = compression");
        }
        for (int d = 0; d < 3; ++d) {
            out[3 * written + d] = static_cast<float>(xyz[d]) * inversePrecision;
        }
        ++written;
    };
    for (std::size_t read = 0; read < atoms; ++read) {
        if (bitSize == 0) {
            for (int d = 0; d < 3; ++d) {
                thisCoord[d] = static_cast<std::int32_t>(decodeBits(state, bitSizeInt[d]));
            }
        } else {
            decodeInts(state, bitSize, sizeInt, thisCoord);
        }
        for (int d = 0; d < 3; ++d) {
            thisCoord[d] += c.minInt[d];
            prevCoord[d] = thisCoord[d];
        }
        int isSmaller = 0;
        if (decodeBits(state, 1) != 0) {
            run = static_cast<int>(decodeBits(state, 5));
            isSmaller = run % 3;
            run -= isSmaller;
            --isSmaller;
        }
        if (run > 0) {
            if (written * 3 + static_cast<std::size_t>(run) > out.size()) {
                throw ReaderError("OSError", "XTC read error = compression");
            }
            for (int k = 0; k < run; k += 3) {
                decodeInts(state, smallIndex, sizeSmall, thisCoord);
                ++read;
                for (int d = 0; d < 3; ++d) {
                    thisCoord[d] += prevCoord[d] - smallNum;
                }
                if (k == 0) {
                    // the first two atoms are swapped (better compression of water molecules)
                    for (int d = 0; d < 3; ++d) {
                        std::swap(thisCoord[d], prevCoord[d]);
                    }
                    put(prevCoord);
                } else {
                    for (int d = 0; d < 3; ++d) {
                        prevCoord[d] = thisCoord[d];
                    }
                }
                put(thisCoord);
            }
        } else {
            put(thisCoord);
        }
        if (isSmaller < 0) {
            --smallIndex;
            smallNum = smaller;
            smaller = smallIndex > kFirstIndex ? magicInt(smallIndex - 1) / 2 : 0;
        } else if (isSmaller > 0) {
            ++smallIndex;
            smaller = smallNum;
            smallNum = magicInt(smallIndex) / 2;
        }
        sizeSmall[0] = sizeSmall[1] = sizeSmall[2] =
            static_cast<std::uint32_t>(magicInt(smallIndex));
    }
    return out;
}

// ---------------------------------------------------------------- frames

constexpr std::int32_t kXtcMagic = 1995;
constexpr std::int32_t kXtcNewMagic = 2023;  // 64-bit byte count of the compressed data
constexpr std::size_t kXtcSmallHeader = 56;  // magic, natoms, step, time, box, natoms
constexpr std::size_t kXtcHeader = kXtcSmallHeader + 32;  // precision, minint, maxint, smallidx

std::uint64_t roundUp4(std::uint64_t n) {
    return (n + 3) & ~std::uint64_t{3};
}

// MDAnalysis: dimensions = triclinic_box(*box) of the nm vectors, lengths then * 10 in float32.
std::optional<std::array<float, 6>> boxOf(const float matrix[9]) {
    std::array<std::array<double, 3>, 3> v{};
    for (int i = 0; i < 9; ++i) {
        v[static_cast<std::size_t>(i / 3)][static_cast<std::size_t>(i % 3)] =
            static_cast<double>(matrix[i]);
    }
    std::array<float, 6> box = triclinicBox(v[0], v[1], v[2]);
    if (box[0] == 0.0f && box[1] == 0.0f && box[2] == 0.0f) {
        return std::nullopt;
    }
    for (int i = 0; i < 3; ++i) {
        box[static_cast<std::size_t>(i)] *= 10.0f;
    }
    return box;
}

void toAngstrom(std::vector<float>& positions) {
    for (float& x : positions) {
        x *= 10.0f;
    }
}

class XdrTrajectory : public FileTrajectory {
public:
    bool storedTime(std::int64_t) const override { return true; }

protected:
    XdrTrajectory(const std::filesystem::path& path, std::optional<double> dt)
        : m_file(path), m_name(pathText(path)), m_givenDt(dt) {}

    // An error found while opening the file, naming it (RD-D03).
    [[noreturn]] void fail(const std::string& message) const {
        throw ReaderError("OSError", m_name + ": " + message);
    }

    // Decodes frame `index`; readFrame adds the file and the frame to its errors (RD-D03).
    virtual void decodeFrame(std::int64_t index, core::Frame& frame) = 0;

    void readFrame(std::int64_t index, core::Frame& frame) override {
        try {
            decodeFrame(index, frame);
        } catch (const ReaderError& error) {
            throw ReaderError(error.pythonName(), m_name + " frame " + std::to_string(index) +
                                                      ": " + error.what());
        }
    }

    // dt: the given one, else time[1] - time[0]; 0 for one frame (even when given)
    void finish() {
        if (m_offsets.empty()) {
            fail("XDR read error = endoffile");
        }
        m_frameCount = static_cast<std::int64_t>(m_offsets.size());
        core::Frame first;
        readFrame(0, first);
        m_firstFileTime = first.time;
        if (m_frameCount > 1) {
            core::Frame second;
            readFrame(1, second);
            m_dt = m_givenDt ? *m_givenDt : second.time - first.time;
        } else {
            m_dt = 0.0;
        }
        seek(0);
    }

    void setTime(std::int64_t index, core::Frame& frame, float fileTime) const {
        frame.time =
            m_givenDt ? static_cast<double>(index) * *m_givenDt : static_cast<double>(fileTime);
    }

    BinaryFile m_file;
    std::string m_name;
    std::optional<double> m_givenDt;
    std::vector<std::uint64_t> m_offsets;
    std::size_t m_atoms = 0;
    double m_firstFileTime = 0.0;
};

class XtcTrajectory : public XdrTrajectory {
public:
    XtcTrajectory(const std::filesystem::path& path, std::optional<double> dt)
        : XdrTrajectory(path, dt) {
        const std::string head = m_file.read(0, 8);
        if (head.size() < 4) {
            fail("XDR read error = endoffile");
        }
        const std::int32_t magic = loadI32(head.data(), true);
        if (magic != kXtcMagic && magic != kXtcNewMagic) {
            fail("XDR read error = magic");
        }
        if (head.size() < 8) {
            fail("XDR read error = endoffile");
        }
        const std::int32_t atoms = loadI32(head.data() + 4, true);
        if (atoms < 0) {
            fail("XDR read error = header");
        }
        m_atoms = static_cast<std::size_t>(atoms);
        const std::uint64_t size = m_file.size();
        if (m_atoms <= 9) {
            const std::uint64_t frameBytes = kXtcSmallHeader + 12 * m_atoms;
            for (std::uint64_t k = 0; k < size / frameBytes; ++k) {
                m_offsets.push_back(k * frameBytes);
            }
        } else {
            const std::size_t countBytes = magic == kXtcNewMagic ? 8 : 4;
            std::uint64_t pos = 0;
            while (pos + kXtcHeader + countBytes <= size) {
                const std::string count = m_file.read(pos + kXtcHeader, countBytes);
                const std::uint64_t bytes =
                    countBytes == 8 ? loadU64(count.data(), true) : loadU32(count.data(), true);
                m_offsets.push_back(pos);
                pos += kXtcHeader + countBytes + roundUp4(bytes);
            }
        }
        finish();
    }

protected:
    void decodeFrame(std::int64_t index, core::Frame& frame) override {
        const std::uint64_t pos = m_offsets[static_cast<std::size_t>(index)];
        const std::string head = m_file.read(pos, kXtcSmallHeader);
        if (head.size() != kXtcSmallHeader) {
            throw ReaderError("OSError", "XTC read error = header");
        }
        const std::int32_t magic = loadI32(head.data(), true);
        if ((magic != kXtcMagic && magic != kXtcNewMagic) ||
            loadI32(head.data() + 4, true) != static_cast<std::int32_t>(m_atoms)) {
            throw ReaderError("OSError", "XTC read error = header");
        }
        const float time = loadF32(head.data() + 12, true);
        float matrix[9];
        for (int i = 0; i < 9; ++i) {
            matrix[i] = loadF32(head.data() + 16 + 4 * i, true);
        }
        std::vector<float> positions;
        if (m_atoms <= 9) {
            const std::string values = m_file.read(pos + kXtcSmallHeader, 12 * m_atoms);
            if (values.size() != 12 * m_atoms) {
                throw ReaderError("OSError", "XTC read error = float");
            }
            positions.resize(3 * m_atoms);
            for (std::size_t i = 0; i < positions.size(); ++i) {
                positions[i] = loadF32(values.data() + 4 * i, true);
            }
        } else {
            const std::size_t countBytes = magic == kXtcNewMagic ? 8 : 4;
            const std::string params = m_file.read(pos + kXtcSmallHeader, 32 + countBytes);
            if (params.size() != 32 + countBytes) {
                throw ReaderError("OSError", "XTC read error = compression");
            }
            Compressed c;
            c.precision = loadF32(params.data(), true);
            for (int i = 0; i < 3; ++i) {
                c.minInt[i] = loadI32(params.data() + 4 + 4 * i, true);
                c.maxInt[i] = loadI32(params.data() + 16 + 4 * i, true);
            }
            c.smallIndex = loadU32(params.data() + 28, true);
            const std::uint64_t bytes = countBytes == 8 ? loadU64(params.data() + 32, true)
                                                        : loadU32(params.data() + 32, true);
            c.data = m_file.read(pos + kXtcHeader + countBytes, static_cast<std::size_t>(bytes));
            if (c.data.size() != bytes) {
                throw ReaderError("OSError", "XTC read error = compression");
            }
            positions = decompressXtc(c, m_atoms);
        }
        toAngstrom(positions);
        frame.positions = std::move(positions);
        frame.box = boxOf(matrix);
        setTime(index, frame, time);
    }
};

class TrrTrajectory : public XdrTrajectory {
public:
    TrrTrajectory(const std::filesystem::path& path, std::optional<double> dt)
        : XdrTrajectory(path, dt) {
        std::uint64_t pos = 0;
        while (true) {
            const auto header = readHeader(pos);
            if (!header) {
                if (pos == 0) {
                    fail("XDR read error = integer");
                }
                break;
            }
            if (pos == 0) {
                m_atoms = header->atoms;
            }
            m_offsets.push_back(pos);
            pos += header->bytes + header->dataBytes();
        }
        finish();
    }

protected:
    void decodeFrame(std::int64_t index, core::Frame& frame) override {
        const std::uint64_t pos = m_offsets[static_cast<std::size_t>(index)];
        const auto header = readHeader(pos);
        if (!header || header->atoms != m_atoms) {
            throw ReaderError("OSError", "TRR read error = header");
        }
        const std::size_t real = header->real;
        std::uint64_t at = pos + header->bytes;
        const auto reals = [&](std::size_t count) {
            const std::string bytes = m_file.read(at, count * real);
            if (bytes.size() != count * real) {
                throw ReaderError("OSError", "TRR read error = float");
            }
            at += count * real;
            std::vector<float> out(count);
            for (std::size_t i = 0; i < count; ++i) {
                out[i] = real == 8 ? static_cast<float>(loadF64(bytes.data() + 8 * i, true))
                                   : loadF32(bytes.data() + 4 * i, true);
            }
            return out;
        };
        float matrix[9] = {};
        if (header->boxSize != 0) {
            const std::vector<float> box = reals(9);
            std::copy(box.begin(), box.end(), matrix);
        }
        at += header->virialSize + header->pressureSize;
        if (header->xSize == 0) {
            throw ReaderError("NoDataError", "This Timestep has no position information");
        }
        std::vector<float> positions = reals(3 * m_atoms);
        toAngstrom(positions);
        frame.positions = std::move(positions);
        frame.box = boxOf(matrix);
        setTime(index, frame, header->time);
    }

private:
    struct Header {
        std::uint64_t bytes = 0;  // of the header
        std::size_t real = 4;
        std::uint64_t boxSize = 0;
        std::uint64_t virialSize = 0;
        std::uint64_t pressureSize = 0;
        std::uint64_t xSize = 0;
        std::uint64_t vSize = 0;
        std::uint64_t fSize = 0;
        std::size_t atoms = 0;
        float time = 0.0f;

        std::uint64_t dataBytes() const {
            return boxSize + virialSize + pressureSize + xSize + vSize + fSize;
        }
    };

    // magic, version string, 13 sizes and counts, time and lambda (float or double)
    std::optional<Header> readHeader(std::uint64_t pos) {
        const std::string start = m_file.read(pos, 12);
        if (start.size() != 12) {
            return std::nullopt;
        }
        const std::uint64_t stringBytes = roundUp4(loadU32(start.data() + 8, true));
        const std::uint64_t intsAt = pos + 12 + stringBytes;
        const std::string ints = m_file.read(intsAt, 52);
        if (ints.size() != 52) {
            return std::nullopt;
        }
        const auto value = [&](int i) { return loadI32(ints.data() + 4 * i, true); };
        Header h;
        h.boxSize = static_cast<std::uint32_t>(value(2));
        h.virialSize = static_cast<std::uint32_t>(value(3));
        h.pressureSize = static_cast<std::uint32_t>(value(4));
        h.xSize = static_cast<std::uint32_t>(value(7));
        h.vSize = static_cast<std::uint32_t>(value(8));
        h.fSize = static_cast<std::uint32_t>(value(9));
        const std::int32_t atoms = value(10);
        if (atoms < 0) {
            return std::nullopt;
        }
        h.atoms = static_cast<std::size_t>(atoms);
        if (h.boxSize != 0) {
            h.real = h.boxSize / 9;
        } else if (h.xSize != 0 && atoms > 0) {
            h.real = h.xSize / (3 * h.atoms);
        } else if (h.vSize != 0 && atoms > 0) {
            h.real = h.vSize / (3 * h.atoms);
        } else if (h.fSize != 0 && atoms > 0) {
            h.real = h.fSize / (3 * h.atoms);
        }
        if (h.real != 4 && h.real != 8) {
            return std::nullopt;
        }
        const std::string times = m_file.read(intsAt + 52, 2 * h.real);
        if (times.size() != 2 * h.real) {
            return std::nullopt;
        }
        h.time = h.real == 8 ? static_cast<float>(loadF64(times.data(), true))
                             : loadF32(times.data(), true);
        h.bytes = 12 + stringBytes + 52 + 2 * h.real;
        return h;
    }
};

}  // namespace

std::unique_ptr<core::Trajectory> openXtc(const std::filesystem::path& path,
                                          std::optional<double> dt, std::size_t) {
    return std::make_unique<XtcTrajectory>(path, dt);
}

std::unique_ptr<core::Trajectory> openTrr(const std::filesystem::path& path,
                                          std::optional<double> dt, std::size_t) {
    return std::make_unique<TrrTrajectory>(path, dt);
}

}  // namespace mdbridge::io::detail
