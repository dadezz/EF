#include <algorithm>
#include <cassert>
#include <iostream>
#include <random>
#include <set>
#include <span>
#include <concepts>
#include <chrono>
#include <vector>

#include "concept_type.hpp"
#include "EF_bucket.hpp"
#include "EF_bucket2.hpp"
using u64 = uint64_t;

std::mt19937_64 rng(42);  // global random generator

template<typename T>
requires EF<T>
void checkCorrectness(std::span<const u64> sorted_vals, 
    u64 universe, const T& t, std::string_view label) {

    assert(t.size() == sorted_vals.size());

    // ################ ACCESS TESTS ################
    for (size_t k = 0; k < sorted_vals.size(); k++) {
        assert(t.access(k) == sorted_vals[k]);
    }

    // ######## PREDECESSOR & SUCCESSOR TESTS ########

    constexpr u64 EXHAUSTIVE_LIMIT = 100'000;
    constexpr int SPOT_CHECKS = 200'000;

    if (universe <= EXHAUSTIVE_LIMIT) {

        for (u64 x = 0; x < universe; x++) {
            auto it = std::upper_bound(sorted_vals.begin(), sorted_vals.end(), x);
            auto expected = (it != sorted_vals.begin()) ? std::optional(*std::prev(it)) : std::nullopt;
            assert(t.predecessor(x) == expected);

            auto sit = std::lower_bound(sorted_vals.begin(), sorted_vals.end(), x);
            auto expected_successor = (sit != sorted_vals.end()) ? std::optional(*sit) : std::nullopt;
            assert(t.successor(x) == expected_successor);
        }
    } else {
        // spot-check
        for (int i = 0; i < SPOT_CHECKS; i++) {
            uint64_t x = rng() % universe;
            
            auto it = std::upper_bound(sorted_vals.begin(), sorted_vals.end(), x);
            auto expected = (it != sorted_vals.begin()) ? std::optional(*std::prev(it)) : std::nullopt;
            assert(t.predecessor(x) == expected);

            auto sit = std::lower_bound(sorted_vals.begin(), sorted_vals.end(), x);
            auto expected_successor = (sit != sorted_vals.end()) ? std::optional(*sit) : std::nullopt;
            assert(t.successor(x) == expected_successor);
        }
    }
    std::cout << "OK: " << label << " (m=" << sorted_vals.size()
              << ", n=" << universe << ", " << t.bytes() << " bytes)\n";
}

enum class TestType {
    ACCESS, 
    PREDECESSOR, 
    SUCCESSOR,
    CONTAINS
};

template<typename T>
requires EF<T>
void benchmark(const T& t, TestType type, std::string_view label) {
    constexpr int NUM_QUERIES = 100000000;
    uint64_t n = t.universe();
    uint64_t size = t.size();

    auto tests = [&] (std::vector<uint64_t> const& queries) -> u64{
        u64 result = 0; // dummy variable to prevent optimization
        switch (type) {
            case TestType::ACCESS:
                for (uint64_t x : queries)
                    result ^= t.access(x).value_or(0);
            break;
            case TestType::PREDECESSOR:
                for (uint64_t x : queries)
                    result ^= t.predecessor(x).value_or(0);
            break;
            case TestType::SUCCESSOR:
                for (uint64_t x : queries)
                    result ^= t.successor(x).value_or(0);
            break;
            case TestType::CONTAINS:
                for (uint64_t x : queries)
                    result ^= t.contains(x);
            break;
        }
        return result;
    };

    // cold start
    std::vector<uint64_t> queries(1000);
    for (int i = 0; i < 1000; i++) {
        if (type == TestType::ACCESS) {
            queries[i] = rng() % size; // Access richiede indici < size
        } else {
            queries[i] = rng() % n;    // Gli altri richiedono valori < universe
        }
    }
    auto dummy = tests(queries);


    // generator of random queries (outside timed section to avoid bias due to RNG)
    std::vector<uint64_t> queries2(NUM_QUERIES);
    for (int i = 0; i < NUM_QUERIES; i++) {
        if (type == TestType::ACCESS) {
            queries2[i] = rng() % size; // Access richiede indici < size
        } else {
            queries2[i] = rng() % n;    // Gli altri richiedono valori < universe
        }
    }

    auto start = std::chrono::high_resolution_clock::now();
    dummy ^= tests(queries2);
    auto end = std::chrono::high_resolution_clock::now();
    
    std::chrono::duration<double, std::nano> elapsed = end - start;
    double ns_per_query = elapsed.count() / NUM_QUERIES;
    
    const char* type_str = "";
    if (type == TestType::ACCESS) type_str = "Access";
    if (type == TestType::PREDECESSOR) type_str = "Predecessor";
    if (type == TestType::SUCCESSOR) type_str = "Successor";
    if (type == TestType::CONTAINS) type_str = "Contains";
    
    std::cout << "Benchmark " << label << "\n"
        << "\ttype:" << type_str 
        << "\telapsed time: " << elapsed.count() 
        << "\tns per query: " << ns_per_query 
        << " (dummy=" << dummy << ")"<<std::endl;
}

template<typename T>
requires EF<T>
void run_test_suite(const std::vector<u64>& sorted_vals, u64 universe, std::string_view label) {
    T ef(sorted_vals, universe);
    checkCorrectness(sorted_vals, universe, ef, label);
    benchmark(ef, TestType::ACCESS, label);
    benchmark(ef, TestType::PREDECESSOR, label);
    benchmark(ef, TestType::SUCCESSOR, label);
    benchmark(ef, TestType::CONTAINS, label);
    std::cout << std::string(60, '-') << "\n";
}

int main() {

    // 1) Random uniform
    {
        uint64_t n = 1'000'000, m = 5000;
        std::set<uint64_t> dedup;
        while (dedup.size() < m) dedup.insert(rng() % n);
        std::vector<uint64_t> sorted(dedup.begin(), dedup.end());
        run_test_suite<BucketEFSet>(sorted, n, "1. Random uniform on EF1");
        run_test_suite<BucketEFSet2>(sorted, n, "1. Random uniform on EF2");
    }

    // 2) Dense cluster — most buckets empty
    {
        uint64_t n = 1'000'000, m = 1000;
        std::set<uint64_t> dedup;
        uint64_t base = 500'000;
        while (dedup.size() < m) dedup.insert(base + (rng() % 2000));
        std::vector<uint64_t> sorted(dedup.begin(), dedup.end());
        run_test_suite<BucketEFSet>(sorted, n, "2. Dense cluster on EF1");
        run_test_suite<BucketEFSet2>(sorted, n, "2. Dense cluster on EF2");
    }

    // 3) Extreme sparsity: few elements, huge universe
    {
        uint64_t n = 1ULL << 40, m = 10;
        std::set<uint64_t> dedup;
        while (dedup.size() < m) dedup.insert(rng() % n);
        std::vector<uint64_t> sorted(dedup.begin(), dedup.end());
        run_test_suite<BucketEFSet>(sorted, n, "3. Extreme sparsity (m=10, n=2^40) on EF1");
        run_test_suite<BucketEFSet2>(sorted, n, "3. Extreme sparsity (m=10, n=2^40) on EF2");
    }

    // 4) Small edge cases
    {
        run_test_suite<BucketEFSet>({0}, 1, "4a. Single {0} n=1 on EF1");
        run_test_suite<BucketEFSet2>({0}, 1, "4a. Single {0} n=1 on EF2");
        run_test_suite<BucketEFSet>({0}, 100, "4b. Single {0} n=100 on EF1");
        run_test_suite<BucketEFSet2>({0}, 100, "4b. Single {0} n=100 on EF2");
        run_test_suite<BucketEFSet>({99}, 100, "4c. Single {99} n=100 on EF1");
        run_test_suite<BucketEFSet2>({99}, 100, "4c. Single {99} n=100 on EF2");
        run_test_suite<BucketEFSet>({0, 99}, 100, "4d. Two {0,99} n=100 on EF1");
        run_test_suite<BucketEFSet2>({0, 99}, 100, "4d. Two {0,99} n=100 on EF2");
    }

    // 5) All elements at start of universe
    {
        std::vector<uint64_t> sorted;
        for (uint64_t i = 0; i < 500; i++) sorted.push_back(i);
        run_test_suite<BucketEFSet>(sorted, 1'000'000, "5. elements at start: first 500 on EF1");
        run_test_suite<BucketEFSet2>(sorted, 1'000'000, "5. elements at start: first 500 on EF2");
    }

    // 6) All elements at end of universe
    {
        std::vector<uint64_t> sorted;
        for (uint64_t i = 0; i < 500; i++) sorted.push_back(999'500 + i);
        run_test_suite<BucketEFSet>(sorted, 1'000'000, "6. elements at end: last 500 on EF1");
        run_test_suite<BucketEFSet2>(sorted, 1'000'000, "6. elements at end: last 500 on EF2");
    }

    // 7) Large random
    {
        uint64_t n = 1ULL << 32, m = 100'000;
        std::set<uint64_t> dedup;
        while (dedup.size() < m) dedup.insert(rng() % n);
        std::vector<uint64_t> sorted(dedup.begin(), dedup.end());
        run_test_suite<BucketEFSet>(sorted, n, "7. large random (n=2^32) on EF1");
        run_test_suite<BucketEFSet2>(sorted, n, "7. large random (n=2^32) on EF2");
    }

    // 8) Evenly spaced elements (stride)
    {
        uint64_t n = 1'000'000;
        std::vector<uint64_t> sorted;
        for (uint64_t i = 0; i < n; i += 100) sorted.push_back(i);
        run_test_suite<BucketEFSet>(sorted, n, "8. Evenly spaced (stride=100) on EF1");
        run_test_suite<BucketEFSet2>(sorted, n, "8. Evenly spaced (stride=100) on EF2");
    }

    // 9) Two disjoint dense clusters separated by huge gap
    {
        uint64_t n = 1'000'000;
        std::set<uint64_t> dedup;
        for(int i=0; i<500; ++i) dedup.insert(rng() % 1000); 
        for(int i=0; i<500; ++i) dedup.insert(999'000 + (rng() % 1000));
        std::vector<uint64_t> sorted(dedup.begin(), dedup.end());
        run_test_suite<BucketEFSet>(sorted, n, "9. Two disjoint clusters on EF1");
        run_test_suite<BucketEFSet2>(sorted, n, "9. Two disjoint clusters on EF2");
    }
    
    // 10) Exponential gaps
    {
        uint64_t n = 1ULL << 60;
        std::vector<uint64_t> sorted;
        for (uint64_t i = 1; i < 60; i++) sorted.push_back(1ULL << i);
        run_test_suite<BucketEFSet>(sorted, n, "10. Exponential gaps (powers of 2) on EF1");
        run_test_suite<BucketEFSet2>(sorted, n, "10. Exponential gaps (powers of 2) on EF2");
    }

    std::cout << "\nAll tests passed.\n";
    return 0;
}
