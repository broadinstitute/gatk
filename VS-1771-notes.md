# VS-1771 Investigation Notes
_Branch: `vs_1771_prepare_oh_so_many_intervals`_

## Goal

Make `GvsPrepareCallset` efficiently handle interval lists of any size
(WGS ~300, Exome ~143K merged, ClinVar ~100K merged, ACAF ~58M) by filtering
the `vet_%` and `ref_ranges_%` prepare data to only the requested intervals,
reducing cost and runtime for AoU Researcher Workbench interval-filtered extracts.

## Conclusion

**The goal is not achievable with the current BigQuery table structure.**
All filtering approaches explored either increase cost, increase runtime, or both.
See "Why there is no viable fix" below. This branch ships no production changes;
the Docker improvements developed alongside this work were extracted into
`vs_1771_leavings` for separate review.

---

## Approach 1: Temp table INNER JOIN

### Design

Three-way strategy in `get_location_filters()`:

1. Pad every interval's **start** by `interval_list_padding` bp to the left
   (default 1000 bp), merge overlapping results via
   `pybedtools.each(_pad).saveas().merge()`.
   Only left-side padding is needed: GVS filters on `location` (the leftmost
   VCF POS of a variant). A variant whose POS falls just left of an interval
   boundary needs left padding to be captured; there is no analogous case on
   the right.
2. Count merged intervals and bases covered on GVS chromosomes (1–22, X, Y).
3. Pick a strategy based on merged interval count and genome coverage:
   - **skip** (no filtering): coverage ≥ `skip_filter_coverage_threshold` (default 50%)
   - **inline SQL WHERE clause**: merged count ≤ 5,000
   - **temp table INNER JOIN**: everything else

Expected outcomes at default 1,000 bp padding:
| Interval list | Merged intervals | Coverage | Strategy   |
|---------------|-----------------|----------|------------|
| WGS           | ~291            | ~89%     | skip       |
| Exome / BGE   | ~143K           | ~18%     | temp_table |
| ClinVar       | ~100K           | ~35%     | temp_table |
| ACAF          | many            | ≥90%     | skip       |

The temp table JOIN used a `chrom_index` equality key to allow BQ to hash-partition
the join by chromosome before applying the BETWEEN range filter, avoiding a pure
O(n × intervals) nested loop. Both populate functions use `SELECT q_all.*` so
interval table columns are excluded from the output.

### Outcome

Integration testing confirmed correct results — prepare tables tied out to truth
data and were significantly smaller than baseline. But:

- **~3× longer runtime** for Exome and ClinVar paths.
- **Equal or higher cost**: `vet_%` / `ref_ranges_%` are partitioned by `sample_id`,
  so BQ cannot prune partitions based on location. The full sample data is always
  scanned; the JOIN adds overhead on top.

The `chrom_index` equality key helps prevent the worst-case nested loop but cannot
avoid the full scan — the BETWEEN predicate still applies to every row in each
chromosome's hash bucket.

---

## Approach 2: Chromosome-batched constant BETWEEN predicates

### Design

Instead of a JOIN against a temp table (variable bounds, no clustering benefit),
express each chromosome's intervals as **constant** `BETWEEN` predicates in the
WHERE clause. Constant bounds allow BigQuery's location clustering to skip blocks
whose `[min, max]` range falls entirely outside all intervals. One query per
chromosome, run in parallel.

BQ's 1 MiB query text limit caps usable batch size at ~8,000 intervals per query.
Chr1 (17,729 merged intervals) requires 3 sequential batches.

### Results

**chr22** (3,661 intervals, 100 samples, `vet_001`):
```
Bytes processed:  992 MB  (vs. 4,188 MB baseline — 76% reduction)
Server job time:  29.3s   (vs. 0.69s baseline)
  S00: Input   62,050,187 → 111 rows   read_ms=157   compute_ms=19,659
```

**chr1** (17,729 intervals in 3 batches of ≤8,000, 3,206 samples, `vet_001`):
```
batch 1/3 (8,000 intervals):  1,804,760,686 rows in →  83,901,828 out
                               28.9 GB  191s  read_ms=105  compute_ms=55,092
batch 2/3 (8,000 intervals):  2,283,382,330 rows in →  81,910,448 out
                               36.5 GB  229s  read_ms=44   compute_ms=68,715
batch 3/3 (1,729 intervals):  2,242,505,950 rows in →  20,533,112 out
                               35.9 GB   16s  read_ms=176  compute_ms=9,520
─────────────────────────────────────────────────────────────────────
Total chr1 filtered:          101.3 GB  437s sequential
Est. full-genome baseline:    ~134 GB
```

### Analysis

`read_ms` is negligible (44–176 ms) — clustering is reading data fast. But BQ
evaluates every OR condition against every row that survives block-skipping.
For chr1 that is 1.8–2.3 billion rows × 8,000 conditions = tens of trillions
of comparisons per batch. `compute_ms` of 55–69 seconds per batch dominates.

Chr22 got a 76% byte reduction because its 3,661 intervals are sparse on a small
chromosome — many blocks fall entirely outside all intervals and are skipped.
Chr1 got only ~25% byte reduction: 17,729 exome intervals are distributed
densely throughout chr1, so almost every block overlaps at least one interval
and must be read. Byte savings on the largest chromosomes (where they matter
most) are minimal.

437 seconds sequential for chr1 alone on one vet table. At AoU scale (~100 vet
tables), even with full parallelism across chromosomes and batches, total prepare
time would be many times longer than baseline.

---

## Why there is no viable fix

BigQuery allows only one partition dimension per table. The `vet_%` and
`ref_ranges_%` tables are partitioned by `sample_id`, which is the correct
choice for the primary access pattern (scan all locations for a given sample
set). Location-based partitioning would require abandoning `sample_id`
partitioning, breaking all other access patterns.

The tables are already clustered by location, providing block-level skipping for
constant-range WHERE filters. Our experiments show exactly what that buys —
76% byte reduction on a small sparse chromosome, ~25% on the densest chromosome —
while the per-row OR evaluation that dominates runtime is unaffected by clustering.
RANGE partitioning by location (a coarser form of the same idea) would not affect
per-row computation and is structurally incompatible with `sample_id` partitioning.

---

## Docker improvements (extracted to `vs_1771_leavings`)

Developed alongside this investigation and extracted into branch `vs_1771_leavings`
for separate review:

- **htslib / bcftools**: 1.22 → 1.23.1
- **bedtools 2.31.1** (new): built from source in `build_base.Dockerfile` with
  `CXXFLAGS="-g -Wall -O2 -std=c++11 -include cstdint"` to work around a missing
  `#include <cstdint>` in bedtools 2.31.1 (GCC ≥ 13). Flagged for reviewer consideration.
- **pyarrow 24.0.0** from PyPI musllinux wheel: eliminated the ~1 hour from-source
  Arrow build from `build_base.Dockerfile`.
