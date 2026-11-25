#pragma once
// numpy.random.default_rng(seed) as MDATK's example data uses it (mdatk/render/example.py): the
// PCG64 bit generator seeded through SeedSequence, and Generator.uniform / normal / integers, giving
// numpy's numbers bit for bit (the GUI's example preview looks exactly like MDATK's). Qt-free.
//
// Port of numpy.random (NumPy 2.5.5): SeedSequence (bit_generator.pyx), PCG64 (src/pcg64/pcg64.h,
// pcg64.c; PCG by Melissa O'Neill, MIT), random_standard_normal (ziggurat, derived from Julia, MIT),
// random_uniform, random_bounded_uint64_fill and the Lemire bounded integers
// (src/distributions/distributions.c), Generator.uniform / normal / integers (_generator.pyx,
// _bounded_integers.pyx.in). Copyright (c) 2005-2017 NumPy Developers; Copyright (c) 2019 Kevin
// Sheppard; BSD-3-Clause / NCSA. Licence texts: packaging/licenses/numpy-random-LICENSE.txt.
// Changes from numpy: C++ API for scalar parameters only; 128-bit arithmetic emulated with two
// 64-bit words (numpy's PCG_EMULATED_128BIT_MATH path) on every compiler.

#include <cstdint>
#include <vector>

namespace mdbridge::render {

// zlib.crc32(data) (CRC-32/ISO-HDLC, the seed of MDATK's example series: crc32 of the UTF-8 key).
std::uint32_t crc32(const void* data, std::size_t size, std::uint32_t crc = 0);

// SeedSequence(entropy).generate_state(nWords, uint32) for a non-negative integer entropy (split
// into 32-bit words, lowest first; 0 is one word) with the default pool size 4.
std::vector<std::uint32_t> seedSequenceState(std::uint64_t entropy, std::size_t nWords);

// numpy.random.Generator(PCG64(seed)) = numpy.random.default_rng(seed).
class NumpyGenerator {
public:
    explicit NumpyGenerator(std::uint64_t seed);

    // PCG64.random_raw(): the next 64-bit output (state step, then the XSL-RR output).
    std::uint64_t next64();
    // next_uint32: the low half of a 64-bit output, the high half kept for the next call.
    std::uint32_t next32();
    // next_double: (next64() >> 11) * 2^-53.
    double nextDouble();

    // Generator.random(): nextDouble().
    double random() { return nextDouble(); }
    // Generator.uniform(low, high): low + (high - low) * nextDouble().
    double uniform(double low, double high);
    std::vector<double> uniform(double low, double high, std::size_t count);
    // Generator.standard_normal() (ziggurat) and normal(loc, scale) = loc + scale * z.
    double standardNormal();
    double normal(double loc, double scale);
    std::vector<double> normal(double loc, double scale, std::size_t count);
    // Generator.integers(low, high) (dtype int64, high exclusive; endpoint: inclusive). Throws
    // core::ValueError("low >= high") (numpy's message; with endpoint "low > high").
    std::int64_t integers(std::int64_t low, std::int64_t high, bool endpoint = false);
    std::vector<std::int64_t> integerArray(std::int64_t low, std::int64_t high, std::size_t count, bool endpoint = false);

private:
    std::uint64_t boundedUint64(std::uint64_t off, std::uint64_t rng);

    std::uint64_t m_stateHigh = 0, m_stateLow = 0;
    std::uint64_t m_incHigh = 0, m_incLow = 0;
    bool m_hasUint32 = false;
    std::uint32_t m_uinteger = 0;
};

}  // namespace mdbridge::render
