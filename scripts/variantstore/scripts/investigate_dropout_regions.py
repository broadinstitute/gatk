"""
investigate_dropout_regions.py
-------------------------------
Investigates the two known variant dropout regions in the AoU v9 VDS by
comparing variant_data and reference_data for the affected superpartitions.

The central question: if variant_data is absent for SP83/chr4 and SP64/chr19,
what does reference_data show?

Possible outcomes and their implications
-----------------------------------------
A) variant_data absent, reference_data ALSO absent for those samples
   → The entire genomic representation is missing. The Avro export simply
     omitted those samples' data for those regions entirely, leaving gaps in
     both the variant and reference layers.

B) variant_data absent, reference_data PRESENT with large reference blocks
   spanning the dropout window
   → The import process "filled in" those regions as all-reference when they
     should have had variant calls. This could happen if the variant Avro and
     reference Avro were exported separately and the variant shard covering
     those regions was silently dropped or empty, while the reference shard
     was intact.

C) variant_data absent, reference_data present but with blocks that stop
   BEFORE the dropout window and resume AFTER it (i.e., a gap in ref_ranges
   that coincides with the variant dropout)
   → Consistent with a complete data loss (both variant and ref) for that
     specific (superpartition, genomic region) combination, perhaps due to a
     scatter-boundary issue in the Avro export interval list.

Usage
-----
Run interactively in a Hail Jupyter notebook or as a standalone script on
a Dataproc cluster that has the VDS accessible.

  python investigate_dropout_regions.py \
      --vds-path     gs://.../aou_srwgs_short_variants_v9_r2_p1.vds \
      --samples-path gs://.../sample_superpartition_sampling.tsv
"""

import argparse
import hail as hl
import pandas as pd


# ---------------------------------------------------------------------------
# Known dropout regions (update as new regions are discovered)
# ---------------------------------------------------------------------------
DROPOUT_REGIONS = [
    dict(
        label       = "SP83 / chr4",
        interval    = "chr4:56580000-57030000",
        superpartition = 83,
    ),
    dict(
        label       = "SP64 / chr19",
        interval    = "chr19:40090000-40580000",
        superpartition = 64,
    ),
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_samples_for_sp(samples_path: str, superpartition: int) -> list[str]:
    """Return sample names (as str) for a single superpartition."""
    with hl.hadoop_open(samples_path, 'r') as f:
        df = pd.read_csv(f, sep='\t', dtype={'sample_name': str})
    return list(df.loc[df['superpartition'] == superpartition, 'sample_name'])


def investigate_window(vds, interval_str: str, samples: list[str], label: str):
    """
    Print a diagnostic summary for one (superpartition, genomic window) pair,
    comparing what variant_data and reference_data contain.
    """
    print(f"\n{'='*70}")
    print(f"  Region: {label}  —  {interval_str}")
    print(f"  Samples in SP: {len(samples)}")
    print(f"{'='*70}")

    interval  = hl.parse_locus_interval(interval_str, reference_genome='GRCh38')
    sp_set    = hl.literal(set(samples), hl.tset(hl.tstr))

    # Parse window boundaries as plain Python ints from the interval string
    # e.g. "chr4:56580000-57030000" -> start=56580000, end=57030000
    _contig, _range = interval_str.split(':')
    window_start, window_end = (int(x) for x in _range.split('-'))
    window_size = window_end - window_start

    # ------------------------------------------------------------------
    # 1. variant_data
    # ------------------------------------------------------------------
    vd = hl.filter_intervals(vds.variant_data, [interval])
    vd = vd.filter_cols(sp_set.contains(vd.s))

    vd_called = vd.filter_rows(hl.agg.any(hl.is_defined(vd.LGT)))
    n_vd_loci = vd_called.count_rows()

    vd_per_sample = (
        vd_called
        .annotate_cols(n_loci=hl.agg.count_where(hl.is_defined(vd_called.LGT)))
        .cols()
        .collect()
    )
    n_samples_with_any_call = sum(1 for r in vd_per_sample if r.n_loci > 0)
    total_variant_calls = sum(r.n_loci for r in vd_per_sample)

    # A healthy window this wide should have many hundreds to thousands of
    # variant calls across 10 samples.  Treat < 10 total as a de-facto dropout.
    VARIANT_DROPOUT_THRESHOLD = 10
    vd_is_dropout = total_variant_calls < VARIANT_DROPOUT_THRESHOLD

    print(f"\n  variant_data")
    print(f"    loci with ≥1 SP sample having defined LGT : {n_vd_loci:,}")
    print(f"    total calls across all SP samples          : {total_variant_calls:,}  "
          f"{'⚠ DROPOUT' if vd_is_dropout else '✓ present'}")
    print(f"    samples with ≥1 defined LGT in window     : {n_samples_with_any_call} / {len(samples)}")
    if vd_per_sample:
        counts = sorted([r.n_loci for r in vd_per_sample], reverse=True)
        print(f"    per-sample locus counts (sorted desc)     : {counts[:10]}{'...' if len(counts)>10 else ''}")

    # ------------------------------------------------------------------
    # 2. reference_data
    # ------------------------------------------------------------------
    rd = hl.filter_intervals(vds.reference_data, [interval])
    rd = rd.filter_cols(sp_set.contains(rd.s))

    rd_present = rd.filter_rows(hl.agg.any(hl.is_defined(rd.END)))
    n_rd_block_starts = rd_present.count_rows()

    # Per-sample: bases covered, and whether blocks reach both window edges.
    # Reaching the window edges tells us the region is continuously covered
    # by reference blocks (no internal gap), which is the hallmark of Outcome B.
    rd_entries = rd_present.entries()
    rd_entries = rd_entries.filter(hl.is_defined(rd_entries.END))

    rd_per_sample = (
        rd_entries
        .group_by(rd_entries.s)
        .aggregate(
            n_blocks=hl.agg.count(),
            bases_covered=hl.agg.sum(
                hl.min(rd_entries.END, window_end) -
                hl.max(rd_entries.locus.position, window_start) + 1
            ),
            max_end=hl.agg.max(rd_entries.END),
            min_start=hl.agg.min(rd_entries.locus.position),
        )
        .collect()
    )

    n_samples_with_ref  = sum(1 for r in rd_per_sample if r.n_blocks > 0)
    coverages           = sorted([r.bases_covered or 0 for r in rd_per_sample], reverse=True)
    # "Reaches start" = a block starts at or before window_start+1 (i.e. the
    #  block was already open when the window begins — normal for interior windows)
    n_reaching_start = sum(1 for r in rd_per_sample
                           if (r.min_start or window_end) <= window_start + 1)
    # "Reaches end" = a block whose END is at or past window_end-1
    n_reaching_end   = sum(1 for r in rd_per_sample
                           if (r.max_end   or 0)          >= window_end   - 1)

    print(f"\n  reference_data (blocks overlapping window)")
    print(f"    ref-block start rows with ≥1 SP sample     : {n_rd_block_starts:,}")
    print(f"    samples with ≥1 ref block in window        : {n_samples_with_ref} / {len(samples)}")
    print(f"    window size (bp)                           : {window_size:,}")
    print(f"    per-sample bases covered (sorted desc)     : {coverages}")
    print(f"    samples whose blocks reach window start    : {n_reaching_start} / {len(samples)}")
    print(f"    samples whose blocks reach window end      : {n_reaching_end} / {len(samples)}")

    # ------------------------------------------------------------------
    # 3. Interpretation
    # ------------------------------------------------------------------
    print(f"\n  Interpretation:")
    if not vd_is_dropout:
        print(f"    variant_data IS present ({total_variant_calls:,} calls) — no dropout detected.")
    elif n_samples_with_ref == 0:
        print(f"    BOTH variant_data and reference_data are absent.")
        print(f"    → Outcome A: complete data gap — the Avro export likely omitted")
        print(f"      this (superpartition, region) entirely.")
    elif n_reaching_start + n_reaching_end >= len(samples):
        # Blocks reach at least one edge for most samples → continuous coverage
        print(f"    variant_data ABSENT, reference_data continuously spans window.")
        print(f"    → Outcome B: the VDS treats this region as all-reference for these")
        print(f"      samples.  Variant rows were either absent from the Avro export or")
        print(f"      dropped during VDS creation (Hail combiner).")
    else:
        print(f"    variant_data ABSENT, reference blocks do NOT span both window edges.")
        print(f"    → Outcome C: the data gap has detectable boundaries within the window.")
        print(f"      The exact positions where ref blocks stop/start could pinpoint")
        print(f"      where the loss occurred.")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description='Investigate known variant dropout regions in a GVS VDS',
    )
    parser.add_argument('--vds-path',     type=str, required=True,
                        help='Path to the Hail VDS (GCS)')
    parser.add_argument('--samples-path', type=str, required=True,
                        help='GCS path to the superpartition samples TSV '
                             '(columns: sample_name, sample_id, superpartition)')
    parser.add_argument('--temp-path',    type=str, default=None,
                        help='Optional Hail temporary directory (GCS path)')
    args = parser.parse_args()

    init_kwargs = {}
    if args.temp_path:
        init_kwargs['tmp_dir'] = f'{args.temp_path}/hail_tmp_investigate'
    hl.init(**init_kwargs)

    print(f'Loading VDS from {args.vds_path} ...')
    vds = hl.vds.read_vds(args.vds_path)

    for region in DROPOUT_REGIONS:
        samples = load_samples_for_sp(args.samples_path, region['superpartition'])
        if not samples:
            print(f"\nWARNING: no samples found for superpartition {region['superpartition']} "
                  f"in {args.samples_path} — skipping {region['label']}")
            continue
        investigate_window(vds, region['interval'], samples, region['label'])

    print('\nDone.')
