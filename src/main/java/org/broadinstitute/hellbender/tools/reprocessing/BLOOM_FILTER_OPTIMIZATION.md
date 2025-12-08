# Bloom Filter Optimization for SplitReadsByRealignmentDifficulty

## Overview

This document describes the optimization of `SplitReadsByRealignmentDifficulty.java` from using exact hash sets to a probabilistic Bloom filter, achieving a **~65x memory reduction** while maintaining correctness for the use case.

## Problem Statement

The original implementation used two `HashedStringSet128` instances to track read names:
- `overlapReadNameSet`: Reads classified as "naive" (easy to realign)
- `noOverlapReadNameSet`: Reads classified as "uncertain" (difficult to realign)

For a typical whole-genome sequencing dataset with **1,162,717,900 read pairs** (~2.3 billion raw reads), this required approximately **20GB of heap memory**.

### Original Memory Analysis

| Component | Per Entry | Total Entries | Memory |
|-----------|-----------|---------------|--------|
| `Long2LongOpenHashMap` overhead | ~16 bytes | 1.16B | ~18.6 GB |
| Hash computation overhead | - | - | ~1-2 GB |
| **Total** | | | **~20 GB** |

## Solution: Bloom Filter with Minority Class Tracking

### Key Insight

Based on empirical data:
- **~85% of reads are naive** (overlapping uniquely mappable regions)
- **~15% of reads are uncertain** (require full realignment)

Rather than tracking both classes, we can:
1. **Only track the minority class** (uncertain reads) in a Bloom filter
2. **Default to the majority class** (naive) for anything not in the filter

### Why This Works

Bloom filters have an important property: **no false negatives**.

- If a read is **in the filter** → it's uncertain (or a false positive from naive)
- If a read is **not in the filter** → it's definitely naive

False positives (naive reads incorrectly classified as uncertain) are **safe** for this use case:
- They simply get extra processing (full realignment)
- No data is lost or corrupted
- The uncertain pipeline handles them correctly

### Memory Calculation

For 175 million uncertain read names with 0.1% false positive rate:

```
Bits per element = -ln(FPR) / (ln(2)²)
                 = -ln(0.001) / 0.4805
                 = 14.4 bits

Total bits = 175,000,000 × 14.4 = 2.52 billion bits
Total bytes = 2.52B / 8 = 315 MB
```

| Parameter | Value |
|-----------|-------|
| Expected uncertain reads | 175,000,000 |
| False positive rate | 0.001 (0.1%) |
| Bits per element | 14.4 |
| **Total memory** | **~314 MB** |

### Memory Comparison

| Implementation | Memory | Reduction |
|----------------|--------|-----------|
| Original (two hash sets) | ~20 GB | - |
| Bloom filter (uncertain only) | ~314 MB | **~65x** |

## Implementation Details

### New Command-Line Arguments

```
--expected_uncertain_reads (-EUR)
    Expected number of uncertain read names (approximately 15% of total read pairs).
    Used to size the Bloom filter. Overestimating is safe and uses slightly more memory;
    underestimating increases false positive rate but remains correct.
    Default: 250,000,000

--bloom_filter_fpr (-FPR)
    False positive rate for the Bloom filter (0.0 to 1.0). Lower values use more memory
    but result in fewer naive reads being incorrectly classified as uncertain.
    Default: 0.001 (0.1%)
```

### Algorithm Changes

#### Pass 1: Classification

**Before:**
```java
if (isReadOverlap && isPrimary) {
    overlapReadNameSet.add(read.getName());
} else {
    noOverlapReadNameSet.add(read.getName());
}
```

**After:**
```java
if (!(isReadOverlap && isPrimary)) {
    // Only track uncertain reads - naive reads require no storage
    uncertainBloomFilter.put(read.getName());
}
```

#### Pass 2: Writing

**Before:**
```java
if (overlapReadNameSet.contains(read.getName())) {
    outputWriter.addRead(read);  // naive
} else if (noOverlapReadNameSet.contains(read.getName())) {
    outputWriterUncertain.addRead(read);  // uncertain
} else {
    // error: unknown read
}
```

**After:**
```java
if (uncertainBloomFilter.mightContain(read.getName())) {
    outputWriterUncertain.addRead(read);  // uncertain (or false positive)
} else {
    outputWriter.addRead(read);  // definitely naive
}
```

### Removed Code

The entire `HashedStringSet128` inner class (~45 lines) was removed, including:
- `Long2LongOpenHashMap` usage
- Manual 128-bit MurmurHash computation
- Collision detection logic

## False Positive Impact Analysis

With 0.1% FPR and the expected read distribution:

| Metric | Value |
|--------|-------|
| Total read pairs | 1,162,717,900 |
| Naive reads (85%) | ~988,310,215 |
| Uncertain reads (15%) | ~174,407,685 |
| Expected false positives | ~988,310 |
| False positive rate (of naive) | 0.1% |

**Impact:** Approximately 988K naive reads will be incorrectly sent to the uncertain output. This is:
- **0.085% of total reads**
- **Safe:** They get extra processing but produce correct results
- **Acceptable trade-off** for 65x memory reduction

## Tuning Guidelines

### For Larger Datasets

If processing more reads, increase `--expected_uncertain_reads`:
```bash
# For 2 billion read pairs, expect ~300M uncertain
--expected_uncertain_reads 300000000
```

### For Lower Memory

Accept higher false positive rate:
```bash
# 1% FPR uses ~60% of the memory
--bloom_filter_fpr 0.01
```

### For Fewer False Positives

Use lower FPR (more memory):
```bash
# 0.01% FPR uses ~140% of the memory
--bloom_filter_fpr 0.0001
```

## Logging Output

The tool now logs:
```
Initialized Bloom filter: expectedInsertions=250000000, fpr=0.001, estimated size (MB)=XXX.X
Uncertain read names added to Bloom filter: XXXXXXXX
Naive reprocessing reads: XXXXXXXXX
Uncertain reprocessing reads: XXXXXXXXX
Estimated false positives (naive reads sent to uncertain due to Bloom filter FPR): ~XXXXXX
```

## Dependencies

Uses Guava's `BloomFilter` (already a GATK dependency):
```java
import com.google.common.hash.BloomFilter;
import com.google.common.hash.Funnels;
```

## References

- [Bloom Filter - Wikipedia](https://en.wikipedia.org/wiki/Bloom_filter)
- [Guava BloomFilter](https://guava.dev/releases/snapshot-jre/api/docs/com/google/common/hash/BloomFilter.html)
- Original implementation: `HashedStringSet128` inner class (removed)
