#!/usr/bin/env python3
"""
plot_merge_stats.py — Plot interval count and % genome coverage vs. per-side
padding for one or more TSV files produced by merge_stats.py.

Usage:
    python plot_merge_stats.py [FILE ...] [-o OUTPUT]

If no files are given, all *_merge.tsv files in the current directory are used.
Labels are derived by stripping the '_merge.tsv' suffix from each filename.
"""

import argparse
import os
import sys
import glob

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("files", nargs="*",
                   help="TSV files to plot (default: all *_merge.tsv in cwd)")
    p.add_argument("-o", "--output", default="merge_stats_plot.png",
                   help="Output image file (default: merge_stats_plot.png)")
    p.add_argument("--dpi", type=int, default=150,
                   help="Output DPI (default: 150)")
    p.add_argument("--title", default="Interval padding + merge analysis",
                   help="Figure title")
    return p.parse_args()


def label_from_path(path):
    base = os.path.basename(path)
    for suffix in ("_merge.tsv", ".tsv"):
        if base.endswith(suffix):
            base = base[: -len(suffix)]
            break
    return base


def bp_formatter(x, _pos):
    """Format x-axis tick labels in human-readable bp/kb/Mb units."""
    if x == 0:
        return "0"
    if x < 1_000:
        return f"{int(x):,}"
    if x < 1_000_000:
        val = x / 1_000
        return f"{val:g}k"
    return f"{x / 1_000_000:g}M"


def main():
    args = parse_args()

    files = args.files or sorted(glob.glob("*_merge.tsv"))
    if not files:
        sys.exit("No *_merge.tsv files found. Pass file paths explicitly or "
                 "run from the directory containing them.")

    datasets = []
    for f in files:
        df = pd.read_csv(f, sep="\t")
        datasets.append((label_from_path(f), df))

    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(10, 8), sharex=True,
        gridspec_kw={"hspace": 0.08}
    )
    fig.suptitle(args.title, fontsize=13, y=0.98)

    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    for i, (label, df) in enumerate(datasets):
        color = colors[i % len(colors)]
        x = df["padding"]
        ax_top.plot(x, df["interval_count"], marker="o", markersize=3,
                    linewidth=1.5, color=color, label=label)
        ax_bot.plot(x, df["pct_genome_covered"], marker="o", markersize=3,
                    linewidth=1.5, color=color, label=label)

    # --- X axis: symlog so d=0 is visible but 100–100k is log-spaced ---
    linthresh = 100  # linear below 100 bp, log above
    ax_bot.set_xscale("symlog", linthresh=linthresh, linscale=0.3)
    ax_bot.xaxis.set_major_formatter(ticker.FuncFormatter(bp_formatter))
    ax_bot.set_xlabel("Left-side padding (bp)", fontsize=11)

    # Pick sensible major ticks manually
    major_ticks = [0, 100, 500, 1_000, 5_000, 10_000, 50_000, 100_000]
    ax_bot.set_xticks(major_ticks)
    ax_bot.tick_params(axis="x", which="minor", bottom=False)

    # --- Top panel: interval count (log Y) ---
    ax_top.set_yscale("log")
    ax_top.set_ylabel("Interval count", fontsize=11)
    ax_top.yaxis.set_major_formatter(
        ticker.FuncFormatter(lambda x, _: f"{x:,.0f}" if x >= 1 else "")
    )
    ax_top.grid(True, which="major", linestyle="--", linewidth=0.5, alpha=0.6)
    ax_top.grid(True, which="minor", linestyle=":", linewidth=0.3, alpha=0.4)
    ax_top.legend(loc="upper right", fontsize=9, framealpha=0.9)
    ax_top.tick_params(axis="x", which="both", bottom=False, labelbottom=False)

    # --- Bottom panel: % coverage ---
    ax_bot.set_ylabel("% genome covered", fontsize=11)
    ax_bot.grid(True, which="major", linestyle="--", linewidth=0.5, alpha=0.6)
    ax_bot.grid(True, which="minor", linestyle=":", linewidth=0.3, alpha=0.4)

    # Add a vertical dashed line at the linthresh boundary so reader knows
    # where the scale transitions
    for ax in (ax_top, ax_bot):
        ax.axvline(linthresh, color="gray", linewidth=0.8,
                   linestyle=":", alpha=0.7)

    fig.text(0.01, 0.5,
             f"← linear | log → (threshold: {linthresh} bp)",
             va="center", ha="left", fontsize=7, color="gray", rotation=90)

    plt.savefig(args.output, dpi=args.dpi, bbox_inches="tight")
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()
