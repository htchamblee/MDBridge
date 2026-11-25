// numpy.random's SeedSequence, PCG64 and Generator.uniform / normal / integers (NumPy 2.5.5), bit
// for bit. Port of numpy/random/bit_generator.pyx, src/pcg64/pcg64.h + pcg64.c (PCG: Copyright 2014
// Melissa O'Neill, MIT), src/distributions/distributions.c (Copyright (c) 2005-2017 NumPy
// Developers, BSD-3-Clause; ziggurat derived from Julia, MIT) and _generator.pyx /
// _bounded_integers.pyx.in (Copyright (c) 2019 Kevin Sheppard, NCSA / BSD-3-Clause). Licence texts:
// packaging/licenses/numpy-random-LICENSE.txt. Changes: C++, scalar parameters, emulated 128-bit
// arithmetic everywhere.

#include "mdbridge/render/numpy_random.hpp"

#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "mdbridge/core/errors.hpp"

namespace mdbridge::render {

namespace {

#include "ziggurat_constants.inc"

// crc32: the reflected polynomial 0xEDB88320, as zlib computes it
const std::array<std::uint32_t, 256>& crcTable() {
    static const std::array<std::uint32_t, 256> table = [] {
        std::array<std::uint32_t, 256> t{};
        for (std::uint32_t n = 0; n < 256; ++n) {
            std::uint32_t c = n;
            for (int k = 0; k < 8; ++k) {
                c = (c & 1U) != 0 ? 0xEDB88320U ^ (c >> 1) : c >> 1;
            }
            t[n] = c;
        }
        return t;
    }();
    return table;
}

// SeedSequence constants (bit_generator.pyx)
constexpr std::uint32_t kInitA = 0x43b0d7e5U;
constexpr std::uint32_t kMultA = 0x931e8875U;
constexpr std::uint32_t kInitB = 0x8b51f9ddU;
constexpr std::uint32_t kMultB = 0x58f38dedU;
constexpr std::uint32_t kMixMultL = 0xca01f9ddU;
constexpr std::uint32_t kMixMultR = 0x4973f715U;
constexpr unsigned kXShift = 16;
constexpr std::size_t kPoolSize = 4;

std::uint32_t hashmix(std::uint32_t value, std::uint32_t& hashConst) {
    value ^= hashConst;
    hashConst *= kMultA;
    value *= hashConst;
    value ^= value >> kXShift;
    return value;
}

std::uint32_t mix(std::uint32_t x, std::uint32_t y) {
    std::uint32_t result = kMixMultL * x - kMixMultR * y;
    result ^= result >> kXShift;
    return result;
}

// PCG64's 128-bit multiplier
constexpr std::uint64_t kMultHigh = 2549297995355413924ULL;
constexpr std::uint64_t kMultLow = 4865540595714422341ULL;

// high and low 64 bits of x * y
void mult64(std::uint64_t x, std::uint64_t y, std::uint64_t& high, std::uint64_t& low) {
    low = x * y;
    const std::uint64_t x0 = x & 0xFFFFFFFFULL;
    const std::uint64_t x1 = x >> 32;
    const std::uint64_t y0 = y & 0xFFFFFFFFULL;
    const std::uint64_t y1 = y >> 32;
    const std::uint64_t w0 = x0 * y0;
    const std::uint64_t t = x1 * y0 + (w0 >> 32);
    std::uint64_t w1 = t & 0xFFFFFFFFULL;
    const std::uint64_t w2 = t >> 32;
    w1 += x0 * y1;
    high = x1 * y1 + w2 + (w1 >> 32);
}

std::uint64_t rotr64(std::uint64_t value, unsigned rot) { return (value >> rot) | (value << ((0U - rot) & 63U)); }

}  // namespace

std::uint32_t crc32(const void* data, std::size_t size, std::uint32_t crc) {
    const auto& table = crcTable();
    const auto* bytes = static_cast<const unsigned char*>(data);
    crc = ~crc;
    for (std::size_t k = 0; k < size; ++k) {
        crc = table[(crc ^ bytes[k]) & 0xFFU] ^ (crc >> 8);
    }
    return ~crc;
}

std::vector<std::uint32_t> seedSequenceState(std::uint64_t entropy, std::size_t nWords) {
    // _coerce_to_uint32_array(int): 32-bit words, lowest first; 0 is [0]
    std::vector<std::uint32_t> words;
    if (entropy == 0) {
        words.push_back(0);
    }
    for (std::uint64_t n = entropy; n > 0; n >>= 32) {
        words.push_back(static_cast<std::uint32_t>(n & 0xFFFFFFFFULL));
    }
    // mix_entropy
    std::array<std::uint32_t, kPoolSize> pool{};
    std::uint32_t hashConst = kInitA;
    for (std::size_t i = 0; i < kPoolSize; ++i) {
        pool[i] = hashmix(i < words.size() ? words[i] : 0U, hashConst);
    }
    for (std::size_t src = 0; src < kPoolSize; ++src) {
        for (std::size_t dst = 0; dst < kPoolSize; ++dst) {
            if (src != dst) {
                pool[dst] = mix(pool[dst], hashmix(pool[src], hashConst));
            }
        }
    }
    for (std::size_t src = kPoolSize; src < words.size(); ++src) {
        for (std::size_t dst = 0; dst < kPoolSize; ++dst) {
            pool[dst] = mix(pool[dst], hashmix(words[src], hashConst));
        }
    }
    // generate_state
    std::vector<std::uint32_t> state(nWords);
    std::uint32_t h = kInitB;
    for (std::size_t i = 0; i < nWords; ++i) {
        std::uint32_t v = pool[i % kPoolSize];
        v ^= h;
        h *= kMultB;
        v *= h;
        v ^= v >> kXShift;
        state[i] = v;
    }
    return state;
}

NumpyGenerator::NumpyGenerator(std::uint64_t seed) {
    // PCG64.__init__: generate_state(4, uint64) -> seed[0..1], inc[0..1] (little-endian word pairs)
    const std::vector<std::uint32_t> w = seedSequenceState(seed, 8);
    auto word = [&](std::size_t k) { return static_cast<std::uint64_t>(w[2 * k]) | (static_cast<std::uint64_t>(w[2 * k + 1]) << 32); };
    const std::uint64_t initStateHigh = word(0), initStateLow = word(1);
    const std::uint64_t initSeqHigh = word(2), initSeqLow = word(3);
    // pcg_setseq_128_srandom_r
    m_stateHigh = 0;
    m_stateLow = 0;
    m_incHigh = (initSeqHigh << 1U) | (initSeqLow >> 63U);
    m_incLow = (initSeqLow << 1U) | 1U;
    next64();
    const std::uint64_t low = m_stateLow + initStateLow;
    m_stateHigh = m_stateHigh + initStateHigh + (low < initStateLow ? 1U : 0U);
    m_stateLow = low;
    next64();
    m_hasUint32 = false;
    m_uinteger = 0;
}

std::uint64_t NumpyGenerator::next64() {
    // state = state * MULT + inc (mod 2^128), then the XSL-RR output of the new state
    std::uint64_t high = 0, low = 0;
    mult64(m_stateLow, kMultLow, high, low);
    high += m_stateHigh * kMultLow + m_stateLow * kMultHigh;
    const std::uint64_t sumLow = low + m_incLow;
    m_stateHigh = high + m_incHigh + (sumLow < m_incLow ? 1U : 0U);
    m_stateLow = sumLow;
    return rotr64(m_stateHigh ^ m_stateLow, static_cast<unsigned>(m_stateHigh >> 58U));
}

std::uint32_t NumpyGenerator::next32() {
    if (m_hasUint32) {
        m_hasUint32 = false;
        return m_uinteger;
    }
    const std::uint64_t next = next64();
    m_hasUint32 = true;
    m_uinteger = static_cast<std::uint32_t>(next >> 32);
    return static_cast<std::uint32_t>(next & 0xFFFFFFFFULL);
}

double NumpyGenerator::nextDouble() { return static_cast<double>(next64() >> 11) * (1.0 / 9007199254740992.0); }

double NumpyGenerator::uniform(double low, double high) {
    const double range = high - low;
    if (!std::isfinite(range)) {
        throw std::overflow_error("high - low range exceeds valid bounds");
    }
    return low + range * nextDouble();
}

std::vector<double> NumpyGenerator::uniform(double low, double high, std::size_t count) {
    std::vector<double> out;
    out.reserve(count);
    for (std::size_t k = 0; k < count; ++k) {
        out.push_back(uniform(low, high));
    }
    return out;
}

double NumpyGenerator::standardNormal() {
    for (;;) {
        std::uint64_t r = next64();
        const int idx = static_cast<int>(r & 0xFFU);
        r >>= 8;
        const bool negative = (r & 0x1U) != 0;
        const std::uint64_t rabs = (r >> 1) & 0x000fffffffffffffULL;
        double x = static_cast<double>(rabs) * wi_double[idx];
        if (negative) {
            x = -x;
        }
        if (rabs < ki_double[idx]) {
            return x;  // 99.3% of the time
        }
        if (idx == 0) {
            for (;;) {
                const double xx = -ziggurat_nor_inv_r * std::log1p(-nextDouble());
                const double yy = -std::log1p(-nextDouble());
                if (yy + yy > xx * xx) {
                    return ((rabs >> 8) & 0x1U) != 0 ? -(ziggurat_nor_r + xx) : ziggurat_nor_r + xx;
                }
            }
        }
        if (((fi_double[idx - 1] - fi_double[idx]) * nextDouble() + fi_double[idx]) < std::exp(-0.5 * x * x)) {
            return x;
        }
    }
}

double NumpyGenerator::normal(double loc, double scale) {
    if (!std::isnan(scale) && std::signbit(scale)) {
        throw core::ValueError("scale < 0");
    }
    return loc + scale * standardNormal();
}

std::vector<double> NumpyGenerator::normal(double loc, double scale, std::size_t count) {
    std::vector<double> out;
    out.reserve(count);
    for (std::size_t k = 0; k < count; ++k) {
        out.push_back(normal(loc, scale));
    }
    return out;
}

std::uint64_t NumpyGenerator::boundedUint64(std::uint64_t off, std::uint64_t rng) {
    // random_bounded_uint64_fill (use_masked = false, one value)
    if (rng == 0) {
        return off;
    }
    if (rng <= 0xFFFFFFFFULL) {
        if (rng == 0xFFFFFFFFULL) {
            return off + next32();
        }
        // buffered_bounded_lemire_uint32
        const auto rngExcl = static_cast<std::uint32_t>(rng + 1);
        std::uint64_t m = static_cast<std::uint64_t>(next32()) * rngExcl;
        auto leftover = static_cast<std::uint32_t>(m & 0xFFFFFFFFULL);
        if (leftover < rngExcl) {
            const std::uint32_t threshold = (std::numeric_limits<std::uint32_t>::max() - static_cast<std::uint32_t>(rng)) % rngExcl;
            while (leftover < threshold) {
                m = static_cast<std::uint64_t>(next32()) * rngExcl;
                leftover = static_cast<std::uint32_t>(m & 0xFFFFFFFFULL);
            }
        }
        return off + (m >> 32);
    }
    if (rng == std::numeric_limits<std::uint64_t>::max()) {
        return off + next64();
    }
    // bounded_lemire_uint64
    const std::uint64_t rngExcl = rng + 1;
    std::uint64_t x = next64();
    std::uint64_t leftover = x * rngExcl;
    if (leftover < rngExcl) {
        const std::uint64_t threshold = (std::numeric_limits<std::uint64_t>::max() - rng) % rngExcl;
        while (leftover < threshold) {
            x = next64();
            leftover = x * rngExcl;
        }
    }
    std::uint64_t high = 0, low = 0;
    mult64(x, rngExcl, high, low);
    return off + high;
}

std::int64_t NumpyGenerator::integers(std::int64_t low, std::int64_t high, bool endpoint) {
    // _rand_int64: closed interval [low, high - 1] (or [low, high] with endpoint)
    const bool invalid = endpoint ? low > high : low >= high;
    if (invalid) {
        if (low == 0) {
            throw core::ValueError(endpoint ? "high < 0" : "high <= 0");
        }
        throw core::ValueError(endpoint ? "low > high" : "low >= high");
    }
    const std::uint64_t top = static_cast<std::uint64_t>(high) - (endpoint ? 0U : 1U);
    const std::uint64_t rng = top - static_cast<std::uint64_t>(low);
    return static_cast<std::int64_t>(boundedUint64(static_cast<std::uint64_t>(low), rng));
}

std::vector<std::int64_t> NumpyGenerator::integerArray(std::int64_t low, std::int64_t high, std::size_t count, bool endpoint) {
    // the fill function shares one Lemire state across the values; with one value per call the
    // stream is the same (the uint32 buffer lives in the bit generator)
    std::vector<std::int64_t> out;
    out.reserve(count);
    for (std::size_t k = 0; k < count; ++k) {
        out.push_back(integers(low, high, endpoint));
    }
    return out;
}

}  // namespace mdbridge::render
