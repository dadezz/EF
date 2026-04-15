#pragma once
#include <concepts>
#include <optional>
#include <cstdint>

template <typename T>
concept EF = requires (const T t, uint64_t k, uint64_t x) {

    /// @brief size of the set
    /// @return the number of elements in the set
    { t.size() } -> std::convertible_to<uint64_t>;

    /// @brief universe size
    /// @return largest possible element + 1, i.e., the smallest integer u such that all elements in the set are < u
    { t.universe() } -> std::convertible_to<uint64_t>;

    /// @brief memory usage
    /// @return the number of bytes used by the data structure
    { t.bytes() } -> std::convertible_to<size_t>;

    /// @brief k-th smallest element (0-indexed), k must be < size()
    /// @param k index of the element to access
    /// @return the k-th smallest element wrapped in std::optional, or std::nullopt if k is out of bounds
    { t.access(k) } -> std::convertible_to<std::optional<uint64_t>>;

    /// @brief largest element <= x (predecessor query)
    /// @param x the value to find the predecessor of
    /// @return the largest element <= x wrapped in std::optional, or std::nullopt if no such element exists
    { t.predecessor(x) } -> std::convertible_to<std::optional<uint64_t>>;

    /// @brief smallest element => x (successor query)
    /// @param x the value to find the successor of
    /// @return the smallest element => x wrapped in std::optional, or std::nullopt if no such element exists
    { t.successor(x) } -> std::convertible_to<std::optional<uint64_t>>;

    /// @brief check if x is in the set
    /// @param x the value to check for membership
    /// @return true if x is in the set, false otherwise
    { t.contains(x) } -> std::same_as<bool>;
};

