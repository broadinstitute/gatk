# VS-1771 Work-in-Progress Notes
_Branch: `vs_1771_prepare_oh_so_many_intervals`_

## What this ticket does

Makes `GvsPrepareCallset` automatically handle interval lists of any size
(WGS ~300, Exome ~300K, ClinVar ~1M, ACAF ~58M) by choosing a strategy
based on post-padding genome coverage, without any user intervention.

## Strategy logic (in `get_location_filters()`)

1. Pad every interval's **start** by `interval_list_padding` bp to the left
   (default 1000), merge overlapping results —
   `pybedtools.each(_pad).saveas().merge()`.
   Only left-side padding is needed: GVS filters on `location` (the leftmost
   VCF POS of a variant). A deletion or ref block whose POS falls just left of
   an interval boundary needs left padding to be captured; there is no analogous
   case on the right (a variant whose POS is right of the interval end starts
   beyond the interval regardless of its length).
2. Count merged intervals and bases covered on GVS chromosomes (1-22, X, Y).
3. Pick a strategy:
   - **skip** (return `""`): coverage ≥ `skip_filter_coverage_threshold` (default 50%)
   - **inline SQL**: merged count ≤ `INTERVAL_TEMP_TABLE_THRESHOLD` (5000)
   - **temp table INNER JOIN**: everything else

Expected outcomes for default 1000 bp padding:
| Interval list | Merged intervals | Coverage | Strategy   |
|---------------|-----------------|----------|------------|
| WGS           | ~291            | ~89%     | skip       |
| Exome / BGE   | ~143K           | ~18%     | temp_table |
| ClinVar       | ~100K           | ~35%     | temp_table |
| ACAF          | many, but >90%  | >90%     | skip       |

## Files changed

### `scripts/variantstore/scripts/create_ranges_cohort_extract_data_table.py`
- Added `HG38_GENOME_SIZE = 3_088_269_832` constant
- Rewrote `get_location_filters()` with three-way strategy
- Updated `make_extract_table()`: replaced `maximum_merge_distance` param with
  `interval_list_padding` and `skip_filter_threshold`
- Updated argparse: `--interval_list_padding` (default 1000),
  `--skip_filter_coverage_threshold` (default 0.5)
- **Temp table SQL**: uses `INNER JOIN ... ON location BETWEEN location_start AND location_end`
  (not `WHERE EXISTS`). BigQuery rejects EXISTS with non-equality correlated subqueries
  (translates to LEFT SEMI JOIN and requires an equality condition). Safe because
  intervals are always merged before loading, so a location matches at most one row.
- Both `populate_final_extract_table_with_ref` and `populate_final_extract_table_with_vet`
  use `SELECT q_all.*` (not `SELECT *`) so the interval table columns are excluded
  when the INNER JOIN is active.

### `scripts/variantstore/wdl/GvsPrepareRangesCallset.wdl`
- Removed the `PadIntervalList` GATK task call (padding now done in Python)
- Added `interval_list_padding = 1000` and `skip_filter_coverage_threshold = 0.5`
  to workflow inputs and plumbed them through to `PrepareRangesCallsetTask`
- Passes raw `interval_list` (not pre-padded) to the Python script

### `scripts/variantstore/wdl/GvsJointVariantCalling.wdl`
- **Bug fix**: was not passing `interval_list` to `GvsPrepareCallset` at all.
  Fixed by adding `interval_list = interval_list_to_use` to the call.
  Discovered during integration testing when prepare tables looked identical
  between branch and baseline.

### `merge_stats.py` (local analysis script, not deployed)
- Rewrote to use slop(padding) + merge instead of merge-distance sweep
- Sweeps over padding values (same default range as before)
- Filters coverage to GVS chromosomes only (matches production logic)
- Added `strategy` column showing what the production code would do

## Docker image work

### `scripts/variantstore/scripts/build_base.Dockerfile`
Changes vs. the prior build-base (`2025-06-08-alpine-build-base`):
- **Removed** entire Apache Arrow / pyarrow build from source (~1 hour of build time).
  Arrow 24+ ships `musllinux` wheels on PyPI so it can now be `pip install`ed.
- **htslib**: 1.22 → 1.23.1
- **bcftools**: 1.22 → 1.23.1
- **vcftools**: 0.1.17 (unchanged, already latest)
- **bedtools 2.31.1** (new — required by `pybedtools`):
  - Not available as an Alpine package; built from source
  - Uses `--strip-components=1` to avoid depending on tarball's internal dir name
  - Requires `CXXFLAGS="-g -Wall -O2 -std=c++11 -include cstdint"` to work around
    a missing `#include <cstdint>` in bedtools 2.31.1's `lineFileUtilities.h`
    (a known issue with GCC ≥ 13)
- Fixed `ENV PERL5LIB` to use modern `KEY=value` format

### `scripts/variantstore/scripts/Dockerfile`
- Added `COPY --from=build /bedtools /bedtools`
- Added `/bedtools/bin` to `PATH`
- Updated build-base tag to `2026-06-16-alpine-build-base`

### `scripts/variantstore/scripts/requirements.txt`
- Added `pyarrow==24.0.0` (replaces the from-source Arrow build)

## Outcome of integration testing

Both Docker images were built and pushed. Manual integration testing confirmed
that the prepare tables produced by this branch are **correct** — they tie out
to the truth data and are significantly smaller than the baseline tables for
interval-filtered extracts (Exome, ClinVar).

However, two performance problems were observed that block shipping this work:

### Slightly higher cost
BigQuery charges are based on **data scanned**, not data written. The `vet_%` and
`ref_ranges_%` tables are partitioned by **sample_id**, not by location. BigQuery
cannot prune partitions based on a location range filter, so the branch scans
exactly the same volume of vet/ref data as the baseline — plus the additional
intervals temp table. Cost is therefore equal-to-or-higher-than baseline.

### ~3x longer runtime for Exome and ClinVar
The `vet_%` and `ref_ranges_%` tables are clustered by location, which does help
for simple constant-range WHERE clause filters (e.g. `WHERE location BETWEEN
1000000001000 AND 1000000002000` — BigQuery uses block-level min/max metadata
to skip blocks). However, clustering **does not help** when the range bounds are
variable inputs from a JOIN partner. In our query:

```sql
INNER JOIN intervals ON chrom_index = ... AND location BETWEEN location_start AND location_end
```

`location_start` and `location_end` vary per interval row, so the optimizer
cannot precompute which blocks to skip — every vet/ref block must be read and
evaluated. This is distinct from BigQuery RANGE partitioning (a separate, newer
feature that would create discrete location-range partitions with guaranteed
pruning at query planning time, regardless of whether bounds are constant or
variable).

The `chrom_index` equality key added to the join does help somewhat by allowing
a hash join to be partitioned by chromosome (reducing each hash bucket to ~6K
intervals rather than 143K), but the BETWEEN evaluation within each chromosome
bucket is still a full scan of all vet/ref rows in that bucket. For an exome
this is still hundreds of millions of row evaluations, a cost that dwarfs the
baseline's simple scan-and-copy.

The inline SQL path (≤5000 intervals) would likely benefit from location
clustering because it emits literal constant BETWEEN predicates in the WHERE
clause. But 143K intervals (exome) or ~100K intervals (ClinVar) exceed the 5000
threshold at which that approach is feasible.

### Chromosome-batched constant-BETWEEN approach (also tried, also rejected)

A second approach was explored: instead of a temp table JOIN, express each
chromosome's intervals as constant `BETWEEN` predicates in the `WHERE` clause
and run one query per chromosome in parallel. Constant bounds allow BigQuery to
use location clustering for block-skipping, which the variable-bound JOIN cannot.

Smoke test on **chr22** (3,661 merged intervals, 100 samples, `vet_001`):

```
Bytes processed:  992 MB  (vs. 4,188 MB baseline — 76% reduction)
Server job time:  29.3s   (vs. 0.69s baseline — 42x slower)

Execution plan:
  S00: Input   62,050,187 → 111 rows   read_ms=157   compute_ms=19,659
```

Smoke test on **chr1** (17,729 merged intervals split into 3 batches of ≤8,000,
3,206 samples, `vet_001`):

```
batch 1/3 (8,000 intervals):  1,804,760,686 rows in →  83,901,828 out
                               28.9 GB  191s  read_ms=105  compute_ms=55,092
batch 2/3 (8,000 intervals):  2,283,382,330 rows in →  81,910,448 out
                               36.5 GB  229s  read_ms=44   compute_ms=68,715
batch 3/3 (1,729 intervals):  2,242,505,950 rows in →  20,533,112 out
                               35.9 GB   16s  read_ms=176  compute_ms=9,520
─────────────────────────────────────────────────────────────────────
Total chr1 filtered:          101.3 GB  437s sequential
Est. full-genome baseline:    ~134 GB   (scaling 100-sample/4.2 GB result)
```

**Why chr22 got 76% byte reduction but chr1 got only ~25%:**
Clustering block-skipping only helps when data blocks fall entirely outside all
interval conditions. Chr22's 3,661 intervals are sparse enough on a small
chromosome that many blocks can be skipped. Chr1's 17,729 exome intervals are
distributed densely throughout the chromosome — almost every block overlaps at
least one interval and must be read.

**The compute cost is catastrophic regardless:**
`read_ms` is negligible (44–176 ms) — clustering reads data fast. But BQ
evaluates every OR condition against every row that survives block-skipping:
55,092 / 68,715 ms of compute per batch on 1.8–2.3 billion rows each. BQ
parallelizes this internally but wall-clock is still 3–4 minutes per batch.

437 seconds sequential for chr1 alone on one vet table. With ~100 vet tables
and 23 chromosomes, even aggressive parallelism across chromosomes and batches
would multiply total prepare time many times over baseline. Rejected.

### Why there is no viable architectural fix

BigQuery allows only one partition dimension per table. The `vet_%` and
`ref_ranges_%` tables are partitioned by `sample_id`, which is the correct
choice for the primary access pattern (scan all locations for a given sample
set). Location-based partitioning would require abandoning `sample_id`
partitioning, breaking all other access patterns.

The `vet_%` / `ref_ranges_%` tables are already clustered by location, so
block-level skipping is already available. Our chr22 and chr1 smoke tests
showed exactly what clustering buys: a 76% byte reduction on chr22 (sparse
intervals, many blocks skippable) but only ~25% on chr1 (dense intervals,
most blocks contain at least one exome interval and cannot be skipped). More
critically, block-skipping is irrelevant to the per-row OR-condition evaluation
that dominates runtime. Even where clustering pruned bytes aggressively, the
compute cost was catastrophic. RANGE partitioning by location (a coarser form
of the same idea) would provide partition-level rather than block-level pruning,
but would not affect per-row computation — and is structurally incompatible with
`sample_id` partitioning for the same reason as above.

### Motivating use case
The original motivation was to reduce cost and runtime for AoU Researcher
Workbench interval-filtered extracts (e.g., a user requesting an exome-sized
extract from a WGS callset). Because we cannot make the prepare step cheaper,
the per-user extract savings are outweighed by the higher prepare cost and
runtime.

## What this branch does ship (Docker + local tooling)

Even though the interval-filtering approach is being set aside, the Docker image
work in this branch is independently useful and can be pulled into a separate PR:

- **htslib / bcftools**: 1.22 → 1.23.1
- **bedtools 2.31.1** (new): required by `pybedtools`; built from source in
  `build_base.Dockerfile` with `CXXFLAGS="-g -Wall -O2 -std=c++11 -include cstdint"`
  to work around a missing `#include <cstdint>` in bedtools 2.31.1 (GCC ≥ 13)
- **pyarrow 24.0.0** from PyPI musllinux wheel (eliminated the ~1 hour
  from-source Arrow build from `build_base.Dockerfile`)
- **`merge_stats.py` / `plot_merge_stats.py`**: local analysis scripts updated
  to model left-side-only padding + merge, matching the revised production logic

## Decision

The interval-filtering changes to `create_ranges_cohort_extract_data_table.py`,
`GvsPrepareRangesCallset.wdl`, and `GvsJointVariantCalling.wdl` will not be
merged. The Docker image updates will be extracted into a separate ticket/PR.
