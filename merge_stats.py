#!/usr/bin/env python3
"""
merge_stats.py — Compute interval count and % genome coverage after
applying left-side padding + merge, mirroring the processing done by
create_ranges_cohort_extract_data_table.py.

Only left-side (start) padding is applied; see that script for the rationale.

For each padding value tested, reports:
  padding              left-side padding in bp applied to every interval's start
  interval_count       merged interval count (GVS chromosomes only)
  bases_covered        total bp covered (GVS chromosomes only)
  pct_genome_covered   bases_covered / HG38_GENOME_SIZE * 100
  strategy             inline | temp_table | skip

Usage:
    python merge_stats.py INPUT [options]

Output: TSV to stdout (or --output FILE).
"""

import argparse
import subprocess
import sys
import os
import tempfile


# Total bases across GVS-encoded chromosomes 1-22, X (23), Y (24).
HG38_GENOME_SIZE = 3_088_269_832

# GVS chromosome name → chromosome index used in location encoding.
CHROM_MAP = {
    **{str(i): i for i in range(1, 23)},
    **{f"chr{i}": i for i in range(1, 23)},
    "X": 23, "chrX": 23,
    "Y": 24, "chrY": 24,
}

INTERVAL_TEMP_TABLE_THRESHOLD = 5000


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input", help="Input BED or Picard interval_list file")
    p.add_argument("-o", "--output", default="-",
                   help="Output TSV file (default: stdout)")

    pad_group = p.add_mutually_exclusive_group()
    pad_group.add_argument("--paddings", nargs="+", type=int, metavar="P",
                           help="Explicit per-side padding values to test (bp)")
    pad_group.add_argument("--range", nargs=3, type=int,
                           metavar=("MIN", "MAX", "STEP"),
                           help="Range of per-side padding values: min max step (bp)")

    p.add_argument("--skip-threshold", type=float, default=0.5,
                   help="Genome coverage fraction at which filtering is skipped "
                        "(default: 0.5)")
    p.add_argument("--bedtools", default="bedtools",
                   help="Path to bedtools executable (default: bedtools)")
    return p.parse_args()


def interval_list_to_bed_file(path, bed_fh):
    """Convert a Picard interval_list to BED (1-based closed → 0-based half-open)."""
    count = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            parts = line.split("\t")
            chrom, start1, end = parts[0], int(parts[1]), parts[2].rstrip()
            bed_fh.write(f"{chrom}\t{start1 - 1}\t{end}\n")
            count += 1
    return count


def prepare_bed_file(input_path):
    """Write input to a temp BED file. Returns (tmp_path, raw_interval_count)."""
    tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False)
    ext = os.path.splitext(input_path)[1].lower()
    try:
        if ext == ".interval_list":
            count = interval_list_to_bed_file(input_path, tmp)
        else:
            count = 0
            with open(input_path) as fh:
                for line in fh:
                    if not line.startswith("#"):
                        tmp.write(line)
                        count += 1
    finally:
        tmp.close()
    return tmp.name, count


def pad_and_merge(bed_path, padding, bedtools_bin):
    """
    Pad every interval's start by `padding` bp to the left (start clamped to 0),
    write to a temp file, then run `bedtools merge`.

    Only left-side padding is applied: GVS filters on `location` (the leftmost
    position of a variant), so only variants whose POS falls left of the interval
    boundary need to be captured by padding.

    Counts only intervals on GVS-encoded chromosomes (1-22, X, Y).
    Returns (interval_count, bases_covered).
    """
    # Pad in Python so we don't need a genome-sizes file for bedtools slop.
    padded_tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False)
    try:
        with open(bed_path) as fh:
            for line in fh:
                parts = line.rstrip().split("\t")
                if len(parts) < 3:
                    continue
                chrom = parts[0]
                start = max(0, int(parts[1]) - padding)
                end = int(parts[2])
                padded_tmp.write(f"{chrom}\t{start}\t{end}\n")
        padded_tmp.close()

        proc = subprocess.Popen(
            [bedtools_bin, "merge", "-i", padded_tmp.name],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        interval_count = 0
        bases_covered = 0
        for line in proc.stdout:
            parts = line.split("\t")
            if len(parts) >= 3 and parts[0] in CHROM_MAP:
                interval_count += 1
                bases_covered += int(parts[2]) - int(parts[1])
        proc.wait()
        if proc.returncode != 0:
            raise subprocess.CalledProcessError(proc.returncode, proc.args,
                                                stderr=proc.stderr.read())
    finally:
        os.unlink(padded_tmp.name)

    return interval_count, bases_covered


def strategy(interval_count, pct_covered, skip_threshold):
    if pct_covered / 100.0 >= skip_threshold:
        return "skip"
    if interval_count <= INTERVAL_TEMP_TABLE_THRESHOLD:
        return "inline"
    return "temp_table"


def main():
    args = parse_args()

    if args.paddings:
        paddings = sorted(args.paddings)
    elif args.range:
        mn, mx, step = args.range
        paddings = list(range(mn, mx + 1, step))
    else:
        # Default: 0, 100 bp steps to 1 kb, 1 kb steps to 10 kb, 10 kb steps to 100 kb
        paddings = (
            [0] +
            list(range(100, 1_000, 100)) +
            list(range(1_000, 10_001, 1_000)) +
            list(range(20_000, 100_001, 10_000))
        )

    print(f"# Input: {args.input}", file=sys.stderr)
    print("# Preparing BED file...", file=sys.stderr, flush=True)
    bed_path, raw_count = prepare_bed_file(args.input)

    print(f"# Raw interval count: {raw_count:,}", file=sys.stderr)
    print(f"# GVS genome size: {HG38_GENOME_SIZE:,} bp", file=sys.stderr)
    print(f"# Skip-filter threshold: {args.skip_threshold:.0%}", file=sys.stderr)
    print(f"# Inline SQL threshold: {INTERVAL_TEMP_TABLE_THRESHOLD:,} intervals",
          file=sys.stderr)
    print(f"# Testing {len(paddings)} padding value(s): "
          f"{paddings[0]} – {paddings[-1]} bp", file=sys.stderr)

    out = open(args.output, "w") if args.output != "-" else sys.stdout
    try:
        out.write("padding\tinterval_count\tbases_covered\tpct_genome_covered\tstrategy\n")
        for p in paddings:
            count, bases = pad_and_merge(bed_path, p, args.bedtools)
            pct = 100.0 * bases / HG38_GENOME_SIZE
            strat = strategy(count, pct, args.skip_threshold)
            out.write(f"{p}\t{count}\t{bases}\t{pct:.6f}\t{strat}\n")
            print(f"  padding={p:>7,}: {count:>8,} intervals, "
                  f"{bases:>14,} bp ({pct:.3f}%)  [{strat}]",
                  file=sys.stderr)
    finally:
        if args.output != "-":
            out.close()
        os.unlink(bed_path)


if __name__ == "__main__":
    main()
