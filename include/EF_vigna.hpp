/*
ELIAS-FANO IMPLEMENTATION USING VIGNA'S SUX LIBRARY
ground truth for benchmarking my custom implementations
*/
#pragma once
#include <span>
#include "concept_type.hpp"
#include "sux/bits/EliasFano.hpp"

class VignaEFSet final {
public:
    VignaEFSet(std::span<const u64> sorted_vals, u64 universe)
        : ef_(std::vector<uint64_t>(sorted_vals.begin(), sorted_vals.end()), universe), 
          m(sorted_vals.size()), 
          n(universe) {}

    u64 size()     const { return m; }
    u64 universe() const { return n; }

    size_t bytes() const { return ef_.bitCount() / 8; }

    std::optional<u64> access(u64 k) const {
        if (k >= m) return std::nullopt;
        return ef_.select(k);
    }

    std::optional<u64> predecessor(u64 x) const {
        if (m == 0) return std::nullopt;
        u64 r = ef_.rank(std::min(x + 1, n));
        if (r == 0) return std::nullopt;
        return ef_.select(r - 1);
    }

    std::optional<u64> successor(u64 x) const {
        if (x >= n) return std::nullopt;
        u64 r = ef_.rank(x);
        if (r >= m) return std::nullopt;
        return ef_.select(r);
    }

    bool contains(u64 x) const {
        if (x >= n) return false;
        u64 r = ef_.rank(x);
        return r < m && ef_.select(r) == x;
    }

private:
    mutable sux::bits::EliasFano<> ef_;
    u64 m, n;
};
