# VS-1771 Work-in-Progress Notes
_Branch: `vs_1771_prepare_oh_so_many_intervals`_

## What this ticket does

Makes `GvsPrepareCallset` automatically handle interval lists of any size
(WGS ~300, Exome ~300K, ClinVar ~1M, ACAF ~58M) by choosing a strategy
based on post-padding genome coverage, without any user intervention.

## Strategy logic (in `get_location_filters()`)

1. Pad every interval by `interval_list_padding` bp per side (default 1000),
   merge overlapping results — using `pybedtools.each(_pad).saveas().merge()`.
   This is equivalent to GATK IntervalListTools `--PADDING`, NOT `merge -d 2000`
   (which only closes gaps, doesn't expand outer boundaries).
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

## Current status

Both Docker images are built and pushed. Short integration tests with WGS and
Exome/BGE data looked reasonable. A long-running integration test against truth
data is now running.

Next steps on return:
1. Check integration test results — verify prepare tables tie out to truth data
   for both WGS (skip path) and Exome/BGE (temp table path)
2. Open PR for VS-1771
