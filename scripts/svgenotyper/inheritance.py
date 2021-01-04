import argparse
import os
import sys

from pysam import VariantFile
from matplotlib import pyplot as plt

INFOS_LIST = ['PPE', 'PSR1', 'PSR2', 'EPE', 'ESR1', 'ESR2', 'PHI_PE', 'PHI_SR1', 'PHI_SR2']

def is_depth_only(record):
    return record.info.ALGORITHMS == 'depth'

def max_carrier_cnlp(record, cutoff: int):
    if 'CNLP' not in record.samples[0]:
        return False
    svtype = record.info['SVTYPE']
    min_cnlps = [min(record.samples[called_sample]['CNLP']) for called_sample in [s for s in record.samples if is_carrier(record.samples[s], svtype)]]
    if len(min_cnlps) == 0:
        return False
    return min(min_cnlps) <= cutoff

def median(x: list):
    if len(x) == 0:
        raise ValueError('Empty list')
    midpoint = len(x) / 2
    return sorted(x)[midpoint]


def is_carrier(record_sample, svtype):
    if svtype == 'DUP':
        return record_sample['CN'] > record_sample['NCN']
    return sum(record_sample['GT']) > 0


class InheritanceSuite(object):
    def __init__(self,
                 contigs: list,
                 args):
        self.out_dir = args.out_dir
        self.out_name = args.out_name
        self.contigs = contigs
        self.vcf = VariantFile(args.vcf)
        self.samples = self.vcf.header.samples
        self.cohort = self.create_family_cohort(args.ped, set(self.samples))

    def create_family_cohort(self, ped_file_path, samples):
        individuals = {}
        with open(ped_file_path) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                tok = line.strip().split('\t')
                fid = tok[0]
                sample_id = tok[1]
                if sample_id not in samples:
                    continue
                sex = tok[4]
                phenotype = tok[5]
                individuals[sample_id] = Individual(id=sample_id, fid=fid, sex=sex, phenotype=phenotype)

        families = {}
        with open(ped_file_path) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                tok = line.strip().split('\t')
                fid = tok[0]
                sample_id = tok[1]
                if sample_id not in individuals:
                    continue
                p1 = tok[2]
                p2 = tok[3]
                if p1 == '0' or p2 == '0':
                    continue
                if p1 not in individuals or p2 not in individuals or sample_id not in individuals:
                    continue
                if fid in families:
                    sibs = tuple(list(families[fid].get_children()).append(individuals[sample_id]))
                else:
                    sibs = (individuals[sample_id],)
                families[fid] = Family(fid=fid, parent1=individuals[p1], parent2=individuals[p2], siblings=sibs)
        cohort = FamilyCohort()
        for fam in families.values():
            cohort.add_family(fam)
        return cohort

    def calculate_mendelian_violations(self):
        valid_proband_calls = []
        invalid_proband_calls = []
        valid_rate_by_vid = {}
        valid_rate_by_allele_count = {}
        num_zero_de_novo = {}
        svtype_counts = {}
        with open(os.path.join(self.out_dir, self.out_name + ".de_novo.bed"), 'w') as bed:
            for record in self.vcf:
                ac = sum([sum(record.samples[s]['GT']) for s in self.samples])
                valid, invalid = self.get_record_mendelian_stats(record)
                valid_proband_calls.extend(valid)
                invalid_proband_calls.extend(invalid)
                n_valid = len(valid)
                n_invalid = len(invalid)
                if n_valid + n_invalid == 0:
                    continue
                valid_rate = n_valid / float(n_valid + n_invalid)
                valid_rate_by_vid[record.id] = valid_rate
                n_de_novo = self.get_record_de_novos(record)
                if n_de_novo > 0:
                    bed.write("\t".join([record.chrom, str(record.pos), str(record.stop), record.info['SVTYPE']]) + '\n')
                svtype = record.info['SVTYPE']
                if svtype not in valid_rate_by_allele_count:
                    valid_rate_by_allele_count[svtype] = ([], [], [])
                    num_zero_de_novo[svtype] = 0
                    svtype_counts[svtype] = 0
                valid_rate_by_allele_count[svtype][0].append(ac)
                valid_rate_by_allele_count[svtype][1].append(valid_rate)
                valid_rate_by_allele_count[svtype][2].append(n_de_novo)
                if n_de_novo == 0:
                    num_zero_de_novo[svtype] += 1
                svtype_counts[svtype] += 1

        sys.stdout.write("Rate of sites lacking de novo calls:\n")
        for key in svtype_counts:
            rate = num_zero_de_novo[key] / float(svtype_counts[key])
            sys.stdout.write("{}: {} / {} ({} %)\n".format(key, num_zero_de_novo[key], svtype_counts[key], rate * 100))

        out_path = os.path.join(self.out_dir, self.out_name + ".valid_rate_by_ac.png")

        rows = 1
        cols = len(valid_rate_by_allele_count.keys())
        figure_width = 3 * len(valid_rate_by_allele_count.keys())
        figure_height = 2

        pf, ax = plt.subplots(rows, cols, figsize=(figure_width, figure_height))

        for i, key in enumerate(valid_rate_by_allele_count.keys()):
            ax[i].plot(valid_rate_by_allele_count[key][0], valid_rate_by_allele_count[key][1], '.k', alpha=0.1, markersize=1)
            ax[i].set_xlabel("AC ({})".format(key))
            if i == 0:
                ax[i].set_ylabel('valid_rate')

        plt.tight_layout()
        plt.savefig(out_path)
        plt.close(pf)

        out_path = os.path.join(self.out_dir, self.out_name + ".de_novo_by_ac.png")

        rows = 1
        cols = len(valid_rate_by_allele_count.keys())
        figure_width = 3 * len(valid_rate_by_allele_count.keys())
        figure_height = 2

        pf, ax = plt.subplots(rows, cols, figsize=(figure_width, figure_height))

        for i, key in enumerate(valid_rate_by_allele_count.keys()):
            ax[i].plot(valid_rate_by_allele_count[key][0], valid_rate_by_allele_count[key][2], '.k', alpha=0.1, markersize=1)
            ax[i].set_xlabel("AC ({})".format(key))
            if i == 0:
                ax[i].set_ylabel('de_novo_count')

        plt.tight_layout()
        plt.savefig(out_path)
        plt.close(pf)

        out_path = os.path.join(self.out_dir, self.out_name + ".de_novo_by_svtype_hist.png")

        rows = 1
        cols = len(valid_rate_by_allele_count.keys())
        figure_width = 3 * len(valid_rate_by_allele_count.keys())
        figure_height = 2

        pf, ax = plt.subplots(rows, cols, figsize=(figure_width, figure_height))

        for i, key in enumerate(valid_rate_by_allele_count.keys()):
            ax[i].hist(valid_rate_by_allele_count[key][2], bins=50)
            ax[i].set_xlabel("de novo count ({})".format(key))
            if i == 0:
                ax[i].set_ylabel('freq')

        plt.tight_layout()
        plt.savefig(out_path)
        plt.close(pf)

        out_path = os.path.join(self.out_dir, self.out_name + ".valid_rate_by_svtype_hist.png")

        rows = 1
        cols = len(valid_rate_by_allele_count.keys())
        figure_width = 3 * len(valid_rate_by_allele_count.keys())
        figure_height = 2

        pf, ax = plt.subplots(rows, cols, figsize=(figure_width, figure_height))

        for i, key in enumerate(valid_rate_by_allele_count.keys()):
            ax[i].hist(valid_rate_by_allele_count[key][1], bins=50)
            ax[i].set_xlabel("valid_rate ({})".format(key))
            if i == 0:
                ax[i].set_ylabel('freq')

        plt.tight_layout()
        plt.savefig(out_path)
        plt.close(pf)

    def get_record_de_novos(self, record):
        n_de_novo = 0
        for fid in self.cohort.get_family_ids():
            fam = self.cohort.get_family_by_family_id(fid)
            p_gt_sum = sum(record.samples[fam.parent1.id]['GT']) + sum(record.samples[fam.parent2.id]['GT'])
            for s in fam.get_siblings():
                s_sum = sum(record.samples[s.id]['GT'])
                if s_sum > p_gt_sum:
                    n_de_novo += 1
        return n_de_novo

    def get_record_mendelian_stats(self, record):
        valid = []
        invalid = []
        for fid in self.cohort.get_family_ids():
            fam = self.cohort.get_family_by_family_id(fid)
            p_gt_max = max(record.samples[fam.parent1.id]['GT']) + max(record.samples[fam.parent2.id]['GT'])
            p_gt_min = min(record.samples[fam.parent1.id]['GT']) + min(record.samples[fam.parent2.id]['GT'])
            for s in fam.get_siblings():
                s_sum = sum(record.samples[s.id]['GT'])
                if s_sum > p_gt_max or s_sum < p_gt_min:
                    invalid.append(s.id)
                else:
                    valid.append(s.id)
        return valid, invalid

class FamilyCohort(object):
    def __init__(self):
        self.families = {}
        self.individuals_to_family_ids = {}

    def add_family(self, family):
        self.families[family.id] = family
        self.individuals_to_family_ids[family.parent1] = family.id
        self.individuals_to_family_ids[family.parent2] = family.id
        for s in family.get_siblings():
            self.individuals_to_family_ids[s] = family.id

    def get_family_by_family_id(self, fid):
        return self.families[fid]

    def get_family_ids(self):
        return self.families.keys()

    def get_family_by_individual_id(self, id):
        return self.families[self.individuals_to_family_ids[id]]

    def get_children(self, parent):
        return self.individuals_to_family_ids[parent].get_siblings()

    def get_siblings(self, sib):
        return tuple(s for s in self.individuals_to_family_ids[sib].get_siblings() if s != sib)


class Individual(object):
    def __init__(self,
                 id: str,
                 fid: str,
                 sex: str,
                 phenotype: str):
        self.id = id
        self.fid = fid
        self.sex = sex
        self.phenotype = phenotype


class Family(object):
    def __init__(self,
                 fid: str,
                 parent1: Individual,
                 parent2: Individual,
                 siblings: tuple):
        self.id = fid
        self.parent1 = parent1
        self.parent2 = parent2
        self.siblings = siblings

    def get_parents(self):
        return (self.father, self.mother)

    def get_siblings(self):
        return tuple(self.siblings)


def parse_args():
    """Parse command line arguments.
    """
    parser = argparse.ArgumentParser()

    parser.add_argument('--vcf', help='Input vcf', required=True)
    parser.add_argument('--out-name', help='Output base name', required=True)
    parser.add_argument('--contigs', help='List of contigs', required=True)
    parser.add_argument('--out-dir', help='Output directory', default='.')
    parser.add_argument('--ped', help='Family structure file in PED format')
    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)

    with open(args.contigs, 'r') as f:
        contigs = [x.strip() for x in f.readlines()]

    suite = InheritanceSuite(contigs=contigs, args=args)

    print("Calculating Mendelian violations...")
    suite.calculate_mendelian_violations()


if __name__ == "__main__":
    main()
