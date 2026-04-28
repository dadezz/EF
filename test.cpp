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
#include "EF_bucket2v2.hpp"
#include "EF_bucketv3.hpp"
#include "EF_vigna.hpp"
using u64 = uint64_t;

std::mt19937_64 rng(42);  // global random generator

template<typename T>
requires EF<T>
void checkCorrectness(std::span<const u64> sorted_vals, const T& t, std::string_view label) {
    std::cout << "Checking correctness of " << label << "...\n";
    assert(t.size() == sorted_vals.size());

    // ################ ACCESS TESTS ################
    for (size_t k = 0; k < sorted_vals.size(); k++) {
        assert(t.access(k) == sorted_vals[k]);
    }

    // ######## PREDECESSOR & SUCCESSOR TESTS ########

    constexpr u64 EXHAUSTIVE_LIMIT = 100'000;
    constexpr int SPOT_CHECKS = 200'000;
    u64 universe = sorted_vals.empty() ? 0 : sorted_vals.back() + 1;
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
u64 run_queries(const T& t, TestType type, const std::vector<u64>& queries) {
    u64 result = 0;
    switch (type) {
        case TestType::ACCESS:
            for (u64 x : queries) result ^= t.access(x).value_or(0);
            break;
        case TestType::PREDECESSOR:
            for (u64 x : queries) result ^= t.predecessor(x).value_or(0);
            break;
        case TestType::SUCCESSOR:
            for (u64 x : queries) result ^= t.successor(x).value_or(0);
            break;
        case TestType::CONTAINS:
            for (u64 x : queries) result ^= t.contains(x);
            break;
    }
    return result;
}

template<typename T>
requires EF<T>
void benchmark(const T& t, TestType type, std::string_view label,
               std::optional<std::pair<u64,u64>> query_range = std::nullopt) {
    constexpr int NUM_QUERIES = 100;

    u64 qlo, qhi;
    if (type == TestType::ACCESS) {
        qlo = 0;
        qhi = t.size();
    } else if (query_range) {
        qlo = query_range->first;
        qhi = query_range->second;
    } else {
        qlo = 0;
        qhi = t.universe();
    }
    u64 span = qhi - qlo;

    // cold start
    std::vector<u64> queries(1000);
    for (int i = 0; i < 1000; i++) queries[i] = qlo + rng() % span;
    u64 dummy = run_queries(t, type, queries);

    // timed queries
    std::vector<u64> queries2(NUM_QUERIES);
    for (int i = 0; i < NUM_QUERIES; i++) queries2[i] = qlo + rng() % span;

    auto start = std::chrono::high_resolution_clock::now();
    dummy ^= run_queries(t, type, queries2);
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
        << " (dummy=" << dummy << ")" << std::endl;
}

template<typename T>
requires EF<T>
void run_test_suite(const std::vector<u64>& sorted_vals, std::string_view label,
                    std::optional<std::pair<u64,u64>> query_range = std::nullopt) {
    T ef(sorted_vals);
    checkCorrectness(sorted_vals, ef, label);
    std::cout << "Running benchmarks for " << label << "...\n";
    benchmark(ef, TestType::ACCESS, label);
    benchmark(ef, TestType::PREDECESSOR, label, query_range);
    benchmark(ef, TestType::SUCCESSOR, label, query_range);
    benchmark(ef, TestType::CONTAINS, label, query_range);
    std::cout << std::string(60, '-') << "\n";
}

void test_every_method(const std::vector<u64>& sorted_vals, std::string_view label,
                    std::optional<std::pair<u64,u64>> query_range = std::nullopt) {
    std::string label_1 = std::string(label) + " on EF1";
    std::string label_2 = std::string(label) + " on EF2";
    std::string label_2v2 = std::string(label) + " on EF2v2";
    std::string label_2v3 = std::string(label) + " on EF2v3";
    std::string label_vigna = std::string(label) + " on Vigna's EF";
    
    run_test_suite<BucketEFSet>(sorted_vals, label_1, query_range);
    run_test_suite<BucketEFSet2>(sorted_vals, label_2, query_range);
    run_test_suite<BucketEFSet2v2>(sorted_vals, label_2v2, query_range);
    run_test_suite<BucketEFSet2v3>(sorted_vals, label_2v3, query_range);
    run_test_suite<VignaEFSet>(sorted_vals, label_vigna, query_range);
}
int main() {
    // 1) Random uniform
    {
        uint64_t n = 1'000'000, m = 5000;
        std::set<uint64_t> dedup;
        while (dedup.size() < m) dedup.insert(rng() % n);
        std::vector<uint64_t> sorted(dedup.begin(), dedup.end());
        test_every_method(sorted, "1. Random uniform");
        
        dedup.insert(n-1); 
        sorted.assign(dedup.begin(), dedup.end());
        test_every_method(sorted, "1b. Random uniform with n-1");
    }
    // 2) Dense cluster — most buckets empty
    {
        uint64_t n = 1'000'000, m = 1000;
        std::set<uint64_t> dedup;
        uint64_t base = 500'000;
        while (dedup.size() < m) dedup.insert(base + (rng() % 2000));
        std::vector<uint64_t> sorted(dedup.begin(), dedup.end());
        test_every_method(sorted, "2. Dense cluster");

        dedup.insert(n-1);
        sorted.assign(dedup.begin(), dedup.end());
        test_every_method(sorted, "2b. Dense cluster with n-1");
    }
    // 3) Extreme sparsity: few elements, huge universe
    {
        uint64_t n = 1ULL << 40, m = 10;
        std::set<uint64_t> dedup;
        while (dedup.size() < m) dedup.insert(rng() % n);
        std::vector<uint64_t> sorted(dedup.begin(), dedup.end());
        test_every_method(sorted, "3. Extreme sparsity (m=10, n=2^40), tight universe");

        dedup.insert(n-1);
        sorted.assign(dedup.begin(), dedup.end());
        test_every_method(sorted, "3b. Extreme sparsity (m=10, n=2^40) with n-1");
    }
    // 4) Small edge cases
    {
        test_every_method({0}, "4a. Single {0} n=1 ");
        test_every_method({99}, "4b. Single {99} n=100");
        test_every_method({0, 99}, "4c. Two {0,99} n=100");
    }
    // 5) All elements at start of universe
    {
        std::vector<uint64_t> sorted;
        for (uint64_t i = 0; i < 500; i++) sorted.push_back(i);
        test_every_method(sorted, "5. n=m, elements at start: first 500");
    }
    // 6) All elements at end of universe
    {
        std::vector<uint64_t> sorted;
        for (uint64_t i = 0; i < 500; i++) sorted.push_back(999'500 + i);
        test_every_method(sorted, "6. n=1000000, elements at end: last 500");
    }
    // 7) Large random
    {
        uint64_t n = 1ULL << 32, m = 100'000;
        std::set<uint64_t> dedup;
        while (dedup.size() < m) dedup.insert(rng() % n);
        std::vector<uint64_t> sorted(dedup.begin(), dedup.end());
        test_every_method(sorted, "7. large random (n=2^32) tight universe");

        dedup.insert(n-1);
        sorted.assign(dedup.begin(), dedup.end());
        test_every_method(sorted, "7b. large random (n=2^32) with n-1");
    }
    // 8) Evenly spaced elements (stride)
    {
        uint64_t n = 1'000'000;
        std::vector<uint64_t> sorted;
        for (uint64_t i = 0; i < n; i += 100) sorted.push_back(i);
        test_every_method(sorted, "8. Evenly spaced (stride=100)");
    }
    // 9) Two disjoint dense clusters separated by huge gap
    {
        uint64_t n = 1'000'000;
        std::set<uint64_t> dedup;
        for(int i=0; i<500; ++i) dedup.insert(rng() % 1000); 
        for(int i=0; i<500; ++i) dedup.insert(999'000 + (rng() % 1000));
        std::vector<uint64_t> sorted(dedup.begin(), dedup.end());
        test_every_method(sorted, "9. Two disjoint clusters, tight universe");

        dedup.insert(n-1);
        sorted.assign(dedup.begin(), dedup.end());
        test_every_method(sorted, "9b. Two disjoint clusters, with n-1");
    }
    // 10) Exponential gaps
    {
        uint64_t n = 1ULL << 60;
        std::vector<uint64_t> sorted;
        for (uint64_t i = 1; i < 60; i++) sorted.push_back(1ULL << i);
        test_every_method(sorted, "10. Exponential gaps (powers of 2), tight universe");
        if (!sorted.back() == n - 1) {
            sorted.push_back(n - 1);
            test_every_method(sorted, "10b. Exponential gaps (powers of 2) with n-1");
        }
    }
    // 11) m = 10M elements
    {
        uint64_t n = 1ULL << 63;
        uint64_t m = 10000000;
        std::set<uint64_t> dedup;
        while (dedup.size() < m) dedup.insert(rng() % n);
        std::vector<uint64_t> sorted(dedup.begin(), dedup.end());
        test_every_method(sorted, "11. m = 10M elements tight universe");

        dedup.insert(n-1);
        sorted.assign(dedup.begin(), dedup.end());
        test_every_method(sorted, "11b. m = 10M elements with n-1");
    }
    // 12) 50M Random Walk (90% small gaps, 10% massive gaps)
    {
        uint64_t n = 1ULL << 63;
        uint64_t m = 50'000'000;
        std::vector<uint64_t> sorted;
        sorted.reserve(m);
        
        uint64_t current = 0;
        for (uint64_t i = 0; i < m; i++) {
            sorted.push_back(current);
            if (rng() % 100 < 90) {
                current += (rng() % 100) + 1;           // 90% chance: gap of 1-100
            } else {
                current += (rng() % (1ULL << 30)) + 1;  // 10% chance: massive gap
            }
            if (current >= n) break; // Prevent overflow
        }
        
        // Ensure unique and strictly increasing
        sorted.erase(std::unique(sorted.begin(), sorted.end()), sorted.end());
        
        test_every_method(sorted, "12. 50M Random Walk tight universe");

        if (sorted.back() != n - 1) {
            sorted.push_back(n - 1);
            test_every_method(sorted, "12b. 50M Random Walk with n-1");
        }
    }
    // 13) 100M Contiguous Bursts
    {
        uint64_t n = 1ULL << 63;
        uint64_t num_clusters = 100'000;
        uint64_t burst_size = 1000;
        
        std::vector<uint64_t> sorted;
        sorted.reserve(num_clusters * burst_size);
        
        for (uint64_t i = 0; i < num_clusters; i++) {
            uint64_t start_val = rng() % (n - burst_size);
            for (uint64_t j = 0; j < burst_size; j++) {
                sorted.push_back(start_val + j);
            }
        }
        
        std::sort(sorted.begin(), sorted.end());
        sorted.erase(std::unique(sorted.begin(), sorted.end()), sorted.end());
        
        test_every_method(sorted, "13. 100M Contiguous Bursts, tight universe");

        if (sorted.back() != n - 1) {
            sorted.push_back(n - 1);
            test_every_method(sorted, "13b. 100M Contiguous Bursts with n-1");
        }
    }
    // 14) 50M Geometric (Power-Law) Gaps
    {
        uint64_t n = 1ULL << 63;
        uint64_t m = 50'000'000;
        std::vector<uint64_t> sorted;
        sorted.reserve(m);
        
        uint64_t current = 0;
        for (uint64_t i = 0; i < m; i++) {
            sorted.push_back(current);
            // Gap magnitude is 2^X, where X is random between 0 and 35
            uint64_t bit_width = rng() % 36; 
            current += (1ULL << bit_width) + (rng() % 10); // Base jump + noise
            if (current >= n) break;
        }
        
        sorted.erase(std::unique(sorted.begin(), sorted.end()), sorted.end());
        
        test_every_method(sorted, "14. 50M Geometric Gaps tight universe");

        if (sorted.back() != n - 1) {
            sorted.push_back(n - 1);
            test_every_method(sorted, "14b. 50M Geometric Gaps with n-1");
        }
    }
    // 15) 10M Adversarial Single-Prefix Cluster
    {
        uint64_t n = 1ULL << 63;
        uint64_t m = 10'000'000;
        std::vector<uint64_t> sorted;
        sorted.reserve(m);
        
        // Fix the top 40 bits to a random static value
        uint64_t fixed_prefix = (rng() % (1ULL << 40)) << 23; 
        
        for (uint64_t i = 0; i < m; i++) {
            // Only the bottom 23 bits vary
            sorted.push_back(fixed_prefix | (rng() % (1ULL << 23)));
        }
        
        std::sort(sorted.begin(), sorted.end());
        sorted.erase(std::unique(sorted.begin(), sorted.end()), sorted.end());
        
        test_every_method(sorted, "15. 10M Adversarial Single-Prefix tight universe");

        if (sorted.back() != n - 1) {
            sorted.push_back(n - 1);
            test_every_method(sorted, "15b. 10M Adversarial Single-Prefix with n-1");
        }
    }
    // 16) 100M True Uniform Random
    {
        uint64_t n = 1ULL << 63;
        uint64_t m = 100'000'000;
        std::vector<uint64_t> sorted(m);
        
        for (auto& val : sorted) {
            val = rng() % n;
        }
        
        std::sort(sorted.begin(), sorted.end());
        sorted.erase(std::unique(sorted.begin(), sorted.end()), sorted.end());
        
        test_every_method(sorted, "16. 100M Uniform Random tight universe");

        if (sorted.back() != n - 1) {
            sorted.push_back(n - 1);
            test_every_method(sorted, "16b. 100M Uniform Random with n-1");
        }
    }
    // 17) Uniform scan to map EF2 vs Vigna crossover
    {
        for (uint64_t m : {10000ULL, 50000ULL, 200000ULL, 500000ULL, 2000000ULL}) {
            uint64_t n = 1000 * m;  // density kept roughly constant
            std::set<uint64_t> dedup;
            while (dedup.size() < m) dedup.insert(rng() % n);
            std::vector<uint64_t> sorted(dedup.begin(), dedup.end());
            std::string label = "17. Uniform crossover m=" + std::to_string(m) + " tight universe";
            test_every_method(sorted, label);

            std::string label_b = "17b. Uniform crossover m=" + std::to_string(m) + " with n-1";
            if (sorted.back() != n - 1) {   
                sorted.push_back(n - 1);
                test_every_method(sorted, label_b);
            }
        }
    }

    // 18) Uniform sparse large m — stresses EF2 exp search in uniform regime
    {
        uint64_t m = 1'000'000, n = 1ULL << 50;  // n/m = 2^30, many empty buckets
        std::set<uint64_t> dedup;
        while (dedup.size() < m) dedup.insert(rng() % n);
        std::vector<uint64_t> sorted(dedup.begin(), dedup.end());
        test_every_method(sorted, "18. Uniform sparse m=1M tight universe");

        if (sorted.back() != n - 1) {
            sorted.push_back(n - 1);
            test_every_method(sorted, "18b. Uniform sparse m=1M with n-1");
        }
    }
    // 19) 50M Random Walk with queries restricted to [min, max] of data
    //     removes the "query in empty universe tail" artifact that inflated test 12
    {
        std::mt19937_64 local_rng(12345);
        uint64_t m = 50'000'000, n = 1ULL << 63;
        std::vector<uint64_t> sorted;
        sorted.reserve(m);
        uint64_t cur = 0;
        for (uint64_t i = 0; i < m; ++i) {
            cur += 1 + (local_rng() % 1000);
            sorted.push_back(cur);
        }
        auto range = std::make_pair(sorted.front(), sorted.back() + 1);
        test_every_method(sorted, "19. 50M Random Walk [bounded queries] tight universe", range);

        range = std::make_pair(sorted.front(), sorted.back() + 1);
        if (sorted.back() != n - 1) {
            sorted.push_back(n - 1);
            test_every_method(sorted, "19b. 50M Random Walk [bounded queries] with n-1", range);
        }
    }

    std::cout << "\nAll tests passed.\n";
    return 0;
}
