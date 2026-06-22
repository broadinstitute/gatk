#!/usr/bin/env python3
"""
smoke_test_chr_clustering.py

Tests whether BigQuery location clustering reduces bytes scanned when a
chromosome's exome intervals are expressed as constant BETWEEN predicates
in the WHERE clause, vs. no location filter (baseline full scan).

For chromosomes with many intervals (e.g. chr1 with ~17K merged intervals),
the WHERE clause can exceed BQ's 1 MiB query limit. Use --max_intervals_per_batch
to split into multiple sequential queries; bytes and timing are summed.

Usage:
    python smoke_test_chr_clustering.py \\
        --project   broad-dsde-methods \\
        --dataset   vs_1234_my_dataset \\
        --vet_table vet_001 \\
        --sample_ids 100 200 300 \\           # OR --sample_ids_file ids.txt
        [--chrom chr1] \\
        [--max_intervals_per_batch 8000] \\
        [--padding 1000] \\
        [--interval_list bge_exome_calling_regions.v1.1.interval_list]

--sample_ids_file: plain text file, one integer sample_id per line.

Requires google-cloud-bigquery and GCP auth (gcloud auth application-default login).
"""

import argparse
import sys
import time

from google.cloud import bigquery


CHROM_MAP = {
    'chr1': 1, 'chr2': 2, 'chr3': 3, 'chr4': 4, 'chr5': 5,
    'chr6': 6, 'chr7': 7, 'chr8': 8, 'chr9': 9, 'chr10': 10,
    'chr11': 11, 'chr12': 12, 'chr13': 13, 'chr14': 14, 'chr15': 15,
    'chr16': 16, 'chr17': 17, 'chr18': 18, 'chr19': 19, 'chr20': 20,
    'chr21': 21, 'chr22': 22, 'chrX': 23, 'chrY': 24,
}

BQ_QUERY_LIMIT_BYTES = 1_048_576  # 1 MiB


def parse_interval_list(path, chrom, padding):
    """Read interval_list, return merged 0-based half-open intervals for `chrom`."""
    raw = []
    with open(path) as f:
        for line in f:
            if line.startswith('@'):
                continue
            parts = line.strip().split('\t')
            if parts[0] != chrom:
                continue
            # interval_list is 1-based closed; BED / pybedtools convention is 0-based half-open
            start = max(0, int(parts[1]) - 1 - padding)
            end = int(parts[2])
            raw.append((start, end))
    raw.sort()
    merged = []
    for s, e in raw:
        if merged and s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return merged


def encode_location(chrom_index, pos_0based):
    return chrom_index * 1_000_000_000_000 + pos_0based


def build_where_clause(intervals, chrom_index):
    """Build a SQL WHERE clause with constant BETWEEN predicates for `intervals`."""
    parts = [
        f"location BETWEEN {encode_location(chrom_index, s)} AND {encode_location(chrom_index, e)}"
        for s, e in intervals
    ]
    clause = "WHERE (" + " OR ".join(parts) + ")"
    clause_bytes = len(clause.encode())
    if clause_bytes > BQ_QUERY_LIMIT_BYTES - 512:
        raise ValueError(
            f"WHERE clause is {clause_bytes:,} bytes — exceeds BQ's 1 MiB query limit. "
            f"Reduce --max_intervals_per_batch (currently {len(intervals):,} intervals in this batch)."
        )
    return clause


def print_query_plan(job):
    stages = job.query_plan
    if not stages:
        print("    (no query plan — likely a cache hit)")
        return
    print(f"    {'Stage':<30s}  {'records_in':>12s}  {'records_out':>12s}  "
          f"{'read_ms':>8s}  {'compute_ms':>10s}  {'write_ms':>8s}")
    for st in stages:
        print(f"    {st.name:<30s}  {(st.records_read or 0):>12,}  {(st.records_written or 0):>12,}  "
              f"  {(st.read_ms_avg or 0):>6,.0f}  "
              f"  {(st.compute_ms_avg or 0):>8,.0f}  "
              f"  {(st.write_ms_avg or 0):>6,.0f}")


def _submit_and_report(client, sql, batch_label):
    """Submit one BQ query, wait for result, print stats. Returns (bytes, count, elapsed)."""
    print(f"  {batch_label}: query={len(sql.encode()):,} bytes  submitting...", flush=True)
    t0 = time.monotonic()
    job = client.query(sql)
    rows = list(job.result())
    elapsed = time.monotonic() - t0

    count = rows[0][0]
    server_s = None
    if job.started and job.ended:
        server_s = (job.ended - job.started).total_seconds()

    print(f"    rows={count:,}  bytes={job.total_bytes_processed:,}"
          f"  wall={elapsed:.2f}s"
          + (f"  server={server_s:.2f}s" if server_s is not None else ""))
    print_query_plan(job)
    return job.total_bytes_processed, count, elapsed


def run_baseline(client, project, dataset, vet_table, sample_ids):
    sample_stanza = ", ".join(str(s) for s in sample_ids)
    fq = f"`{project}.{dataset}.{vet_table}`"
    sql = f"SELECT COUNT(*) FROM {fq} WHERE sample_id IN ({sample_stanza})"
    print(f"\n[BASELINE (no location filter)]")
    return _submit_and_report(client, sql, "baseline")


def run_filtered(client, project, dataset, vet_table, sample_ids,
                 intervals, chrom_index, chrom, max_per_batch):
    """Run one or more batched queries covering all intervals; return summed totals."""
    batches = [intervals[i:i + max_per_batch]
               for i in range(0, len(intervals), max_per_batch)]
    n = len(batches)
    print(f"\n[FILTERED ({chrom} constant BETWEEN predicates)]")
    print(f"  {len(intervals):,} intervals → {n} batch(es) of ≤{max_per_batch:,}")

    sample_stanza = ", ".join(str(s) for s in sample_ids)
    fq = f"`{project}.{dataset}.{vet_table}`"

    total_bytes, total_count, total_time = 0, 0, 0.0
    for i, batch in enumerate(batches):
        try:
            where = build_where_clause(batch, chrom_index)
        except ValueError as e:
            sys.exit(str(e))
        location_cond = where[len("WHERE "):]
        sql = (f"SELECT COUNT(*) FROM {fq} "
               f"WHERE sample_id IN ({sample_stanza}) AND {location_cond}")
        b, c, t = _submit_and_report(client, sql, f"batch {i+1}/{n}")
        total_bytes += b
        total_count += c
        total_time += t

    if n > 1:
        print(f"  totals: rows={total_count:,}  bytes={total_bytes:,}  wall={total_time:.2f}s")
    return total_bytes, total_count, total_time


def load_sample_ids(args):
    if args.sample_ids_file:
        with open(args.sample_ids_file) as f:
            return [int(line.strip()) for line in f if line.strip()]
    return args.sample_ids


def ratio_str(a, b):
    return f"{a/b:>7.3f}x" if b else "  (cache)"


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--project",   required=True)
    p.add_argument("--dataset",   required=True)
    p.add_argument("--vet_table", required=True, help="e.g. vet_001")

    id_group = p.add_mutually_exclusive_group(required=True)
    id_group.add_argument("--sample_ids", nargs="+", type=int)
    id_group.add_argument("--sample_ids_file",
                          help="Plain text file, one integer sample_id per line")

    p.add_argument("--chrom", default="chr22")
    p.add_argument("--max_intervals_per_batch", type=int, default=8000,
                   help="Max intervals per WHERE clause batch (default: 8000, ~430 KB per query)")
    p.add_argument("--padding", type=int, default=1000,
                   help="Left-side padding in bp (default: 1000)")
    p.add_argument("--interval_list",
                   default="bge_exome_calling_regions.v1.1.interval_list")
    args = p.parse_args()

    if args.chrom not in CHROM_MAP:
        sys.exit(f"Unknown chromosome: {args.chrom}. Known: {sorted(CHROM_MAP)}")

    sample_ids = load_sample_ids(args)
    print(f"Loaded {len(sample_ids):,} sample_ids")

    chrom_index = CHROM_MAP[args.chrom]
    print(f"Parsing {args.chrom} from {args.interval_list} (padding={args.padding} bp)...")
    intervals = parse_interval_list(args.interval_list, args.chrom, args.padding)
    n_batches = -(-len(intervals) // args.max_intervals_per_batch)  # ceiling div
    print(f"  {len(intervals):,} merged intervals → {n_batches} batch(es) of ≤{args.max_intervals_per_batch:,}")

    client = bigquery.Client(project=args.project)

    bytes_b, count_b, time_b = run_baseline(
        client, args.project, args.dataset, args.vet_table, sample_ids)

    bytes_f, count_f, time_f = run_filtered(
        client, args.project, args.dataset, args.vet_table, sample_ids,
        intervals, chrom_index, args.chrom, args.max_intervals_per_batch)

    print("\n--- Summary ---")
    print(f"  {'':25s}  {'Baseline':>15s}  {'Filtered':>15s}  {'Ratio':>10s}")
    print(f"  {'Bytes processed':25s}  {bytes_b:>15,}  {bytes_f:>15,}  {ratio_str(bytes_f, bytes_b)}")
    print(f"  {'Rows returned':25s}  {count_b:>15,}  {count_f:>15,}  {ratio_str(count_f, count_b)}")
    print(f"  {'Wall-clock (s)':25s}  {time_b:>15.2f}  {time_f:>15.2f}  {ratio_str(time_f, time_b)}")
    if bytes_b > 0:
        print(f"\n  Byte reduction: {(1 - bytes_f/bytes_b)*100:.1f}%")


if __name__ == "__main__":
    main()
