#!/bin/python

import argparse
import os

import hail as hl


def parse_args():
    """Parse command line arguments.
    """
    parser = argparse.ArgumentParser()

    parser.add_argument('--vcf', help='Prepped vcf', required=True)
    parser.add_argument('--hail-table', help='Hail table output path', required=True)
    parser.add_argument('--reference', help='Reference genome', default="GRCh38")

    args = parser.parse_args()

    return args


def main():
    args = parse_args()
    hl.import_vcf(args.vcf, force_bgz=True, reference_genome=args.reference).write(args.hail_table, overwrite=True)


if __name__ == "__main__":
    main()