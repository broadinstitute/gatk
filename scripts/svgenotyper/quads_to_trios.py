#!/bin/python

import argparse
import sys


def create_family_cohort(ped_file_path, samples):
    families = {}
    sexes = {}
    phenotypes = {}
    with open(ped_file_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            tok = line.strip().split('\t')
            fid = tok[0]
            sample_id = tok[1]
            if sample_id not in samples:
                continue
            p1 = tok[2]
            p2 = tok[3]
            sex = tok[4]
            pheno = tok[5]
            phenotypes[sample_id] = pheno
            sexes[sample_id] = sex
            if p1 != '0' and p2 != '0':
                if fid not in families:
                    families[fid] = []
                families[fid].append((sample_id, p1, p2))

        for fid, famlist in families.items():
            for i, fam in enumerate(famlist):
                p1 = fam[0]
                p2 = fam[1]
                child = fam[2]
                if p1 in sexes and p2 in sexes and child in sexes:
                    fid_i = "{}_{}".format(fid, i)
                    sys.stdout.write("\t".join([fid_i, p1, '0', '0', sexes[p1], phenotypes[p1]]) + '\n')
                    sys.stdout.write("\t".join([fid_i, p2, '0', '0', sexes[p2], phenotypes[p2]]) + '\n')
                    sys.stdout.write("\t".join([fid_i, child, p1, p2, sexes[child], phenotypes[child]]) + '\n')


def parse_args():
    """Parse command line arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--samples', help='Samples list')
    parser.add_argument('--ped', help='Family structure file in PED format containing quads')
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    with open(args.samples, 'r') as f:
        samples = set([x.strip() for x in f.readlines()])
    create_family_cohort(args.ped, samples)


if __name__ == "__main__":
    main()
