# Bucket Elias-Fano implementation.
This repo contains an academic research on alternative Elias-Fano implementations, using a bucket-based approach instead of the traditional rank/select approach.

The goal of this project is to provide a new Data Structure that is more efficient in terms of execution time, specifically for the *predecessor* queries, while maintaining a succint representation of the data.
The Elias-Fano implementation taken as a reference is the one implemented by Sebastiano Vigna in vigna/sux repository. 

Note that in the repo you'll find many different implementations, which are cronologically ordered. Each implementation is (or at least should be) an improvement of the previous one, based on the results of the tests and benchmarks.
The entire project is still WIP, so the code could be messy, and so could be the benchmarks.

# High level overview of the ideas
Input: a sorted list of _m_ unique integers in the range [0, n), where _n_ is the universe size. 

Query that we want to support:

1. access(i): returns the i-th integer in the list, aka random-access query.
1. predecessor(x): returns the largest integer in the list that is less than or equal to x.
1. successor(x): returns the smallest integer in the list that is greater than or equal to x.
1. contains(x): returns true if x is in the list, false otherwise.

High-level idea of the implementation:

the binary representation of numbers is split into high (most significant) and low (least significant) bits. Let's call the high bits "prefix" and the low bits "suffix".
The only purpose of the prefix is to decide which bucket the number belongs to. Central idea: If I'm looking at the bucket number 6, prefix is 6 by definition, so it doesn't need to be stored explicitly: We need just the suffixes, which are then stored compactly in a contiguous array.

There's a problem though: this array is just an indistinct bitvector, it doesn't contain any info about where each bucket starts and ends.
Every variant of the implementation in this repo is a different attempt to solve this problem, while keeping the representation as compact as possible and the queries as fast as possible.

# Data Structures investigated
each implementation will be called EFx, where x is an incremental number (there's a mismatch between the code and the readme, it'll be fixed). All the math briefly mentioned into the various sections, will have a formal discussion later.

## EF1
Single vector split in three logical sections: **S[]B[]F[]**

* **S = bucket starts**. S[i] = How many elements there are in the buckets before the i-th one. It contains therefore two informations: "how many numbers are there before this prefix?" and "what is the starting position of this bucket in the F vector?"
* **B = bitflag**. 1 if the bucket has at least one element, 0 if it's empty. Needed for the predecessor(x)/successor(x) queries, when the predecessor doesn't belong to the same bucket of x.
* **F = suffixes**. data.

### Problems encountered
B[] is a problem. If the data is dense, scanning B is really fast, since we can parallelize the scan over 128 bits using hardware-accelerated instructions.
Over sparse data, it's a disaster. For large datasets, there can be millions of buckets, leading to a size of many kilobytes for the bitflag and the cost is linear on the distance. Plus, every pred/succ query leads to a 3 cache misses by design

## EF2
Only 2 logical sections: **S[]F[]**.

We got rid of the bitmap. B[] is redundant. S[] contains in itself the information: bucket _i_ is empty if and only if S[i]==S[i+1].
So we can do an exponential search on S[]. Tradeoff: logaritmic cost on the distance. extremely faster over sparse data, slower on denser. 2 cache misses by design.

### Problems encountered
Well, it's a tradeoff. on clustered dense data, I measured some queries 3x slower than EF1 (90ns vs 30ns). Exponential search has some costs, but overall it's an improvement. The data structured could be further investigated, with tradeoffs over the memory

## EF3
Why *search* the predecessor cross-bucket if it can be cached?
Let's modify **S[]F[]** in **SP[]F[]**, where SP[] is the portion of the array containing the starting positions of the buckets (as always) and the cached largest element predecessor of the bucket, interleaved.
This waay, predecessor(x) cross-bucket is a single lookup. Is it still succint? log(n) bits per bucket, over m/logm buckets, so mlogn/logm = o(mlog(n/m)) bits, which is still succint.

### problems encountered
access, contains and successor never touch P[], so the interleaving is just streching the distance between S[i] and S[i+1], leading to more cache misses.

### for myself in the future
There's no need at all to store the predecessor entirely. what's missing is the information about the prefix. the suffix is already read into F[], so uoi can just store the prefix to retrieve the number, saving a lot of space

## EF4
Same as EF3, but without the interleaving. **S[]P[]F[]**. Successor -30%, contains -23%, predecessor +5%, same memory occupancy.

### problems ecountered
None. this is the best variant in the group variable-length-buckets. Still to be investigated the "eureka moment" of the previous "for myself in the future"

## EF5
Let's change completely the design. Let's think about how to obtain a single cache miss (on average) per query.
The underlying problem of the previous versions is that the buckets are always variable-length. In order to know where a bucket starts, you got to read S[], leading to at least 2 cache misses by design.
It would be better to leave this information implicit into a deterministic formula. Let's try to build fixed-length buckets. 
It brings A LOT of problems when there's not enough space in the bucket. I came up with a solution.
Layout of the single bucket:
```
┌────────┬──────────┬──────┬──────┬─────┬──────┐
│  size  │ max_prev │ slot │ slot │ ... │ slot │
└────────┴──────────┴──────┴──────┴─────┴──────┘
  log m     log n       b      b            b
            bit
C = log m + c·√(log m) slot
```
where:
* size: how many elements are there into the bucket
* max_prev: the cross-bucket predecessor, the same idea of P[] above
* slot: the suffixes

C = log m + c·√(log m) slot? on average, every bucket receives logm elements. It's a random variable, therefore, I added c·√(log m), squared rooted because it's the sd of a binomial distribution, and the probability of an overflow is limited to exp(−c²/2); the space is still o(m*b).

what about the overflows? there's an auxiliary data structure, at the moment leaved as a simple vector, it can be handled with an EF interface, with another API added (accessByPrefix), more on this later. Rather than a recursive EF5 structure, an EF4 auxiliary should be better in the degenerate case of all the m numbers in a single bucket.

There still one problem regarding the access query. how can it know what bucket look for? There should still be an **S[]** somehow, otherwise it should linearly scan every single bucket to retrieve the size.
a complete S[] vector would be expensive, so I opted for a sampled S[]: one sample every √B buckets, which is the best tradeoff, balancing sample array size against scan length.

### problems encountered
Still under investigation. Predecessor is truly 1 cache miss by design on average, but:
* memory occupancy exploded compared to the other solutions (there's still the variable of the overflow data structure, atm it inflates the total memory)
* access 5-15 times slower
* empty buckets are a waste of space, so is the extra-space c·√(log m), by design.
* not implemented yet: at the cost of losing the 1-cache miss assumption in the adversarial cases, i could store the max_prev directly in the empty slots of the previous bucket, if empty. otherwise (if not empty), search the predecessor directly on the previous bucket. It's expensive only when the predecessor has to be retrieved from the overflow DS.