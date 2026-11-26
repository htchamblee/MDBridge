#pragma once
// Small helpers shared by the numerics sources (not installed): BLAS-like sums in index order and a
// deterministic parallel loop.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <span>
#include <thread>
#include <vector>

namespace mdbridge::numerics::detail {

// np.dot of two vectors, summed in index order (BLAS' ddot may block the sum differently: the last
// bits can differ from numpy's)
inline double dot(std::span<const double> a, std::span<const double> b) {
    double s = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        s += a[i] * b[i];
    }
    return s;
}

// np.linalg.norm of a vector: sqrt(dot(x, x))
inline double norm2(std::span<const double> a) {
    return std::sqrt(dot(a, a));
}

// fn(i) for i in [0, count) on up to `threads` threads; every index is handled exactly once, so a
// result written per index does not depend on the thread count. The first exception is rethrown.
template <class F>
void parallelFor(std::size_t count, std::size_t threads, F&& fn) {
    threads = std::max<std::size_t>(1, std::min(threads, count));
    if (threads <= 1) {
        for (std::size_t i = 0; i < count; ++i) {
            fn(i);
        }
        return;
    }
    std::vector<std::exception_ptr> errors(threads);
    std::vector<std::thread> pool;
    pool.reserve(threads);
    for (std::size_t t = 0; t < threads; ++t) {
        pool.emplace_back([&, t]() {
            try {
                for (std::size_t i = t; i < count; i += threads) {
                    fn(i);
                }
            } catch (...) {
                errors[t] = std::current_exception();
            }
        });
    }
    for (std::thread& th : pool) {
        th.join();
    }
    for (const std::exception_ptr& e : errors) {
        if (e) {
            std::rethrow_exception(e);
        }
    }
}

}  // namespace mdbridge::numerics::detail
