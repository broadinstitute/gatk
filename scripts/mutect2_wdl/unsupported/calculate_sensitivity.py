"""
The input to this script is a tab-separated table with header columns:
STATUS	BAM_DEPTH	AF_EXP	TYPE

STATUS takes on values TP (true positive), FP (false positive), FN (false negative), FFN (filtered false negative)
as emitted by the GATK tool Concordance (the input table come from running VariantsToTable on an output vcf of Concordance).
BAM_DEPTH is the depth of the simulated tumor (ie the Hapmap mixture) at each site
AF_EXP is the expected allele fraction based on the component samples' genotypes and the mixing fractions
TYPE is either SNP or INDEL

This script divides the data by SNP/INDEL and by bins of depth and expected allele fraction and computes the sensitivity
(true positives / (true positives + false negatives + filtered false negatives)) for each bin, along with error bars
on the sensitivity and the raw counts of true and called variants.  It produces tables containing this information and plots of
the sensitivity with error bars.



"""

from __future__ import division
import numpy as np
import pandas as pd
import argparse
from itertools import chain
from statsmodels.stats.proportion import proportion_confint
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

SNP_AF_BINS = [0.05, 0.1, 0.2, 0.4, 0.8, 1.0]
INDEL_AF_BINS = [0.1, 0.2, 1]
AF_JITTER = 0.05
EPSILON = 0.0001 #regularizer for 0/0

EXPECTED_AF_COLUMN = 'AF_EXP'
DEPTH_COLUMN = 'BAM_DEPTH'
STATUS_COLUMN = 'STATUS' # 1 if called, 0 otherwise

def flatten(x):
    return list(chain.from_iterable(x))

def generate_sensitivity_table(df, snp_or_indel, depth_bins, depth_bin_width):

    def assign_bin(x, bins, jitter):
        """ Assigns 'x' into a bin +/- jitter. """
        for _bin in bins:
            lower, upper = _bin - jitter, _bin + jitter
            if (x > lower and x <= upper):
                return _bin
        return np.NaN

    def assign_depth_bin(depth):
        return assign_bin(depth, depth_bins, depth_bin_width)

    def assign_snp_af_bin(af):
        return assign_bin(af, SNP_AF_BINS, AF_JITTER)

    def assign_indel_af_bin(af):
        if (af > 0 and af <= 0.1):
            return 0.1
        if (af > 0.1 and af <= 0.2):
            return 0.2
        if (af > 0.2 and af <= 1):
            return 1
        else:
            return np.NaN

    # add the 'bin' columns for each row
    df.loc[:, 'af_bin'] = df[EXPECTED_AF_COLUMN].apply(assign_snp_af_bin if snp_or_indel == "SNP" else assign_indel_af_bin)
    df.loc[:, 'depth_bin'] = df[DEPTH_COLUMN].apply(assign_depth_bin)
    df.loc[:, 'is_tp'] = df['STATUS'] == 'TP'
    df = df[~np.isnan(df['af_bin'])]

    # make a matrix of binned AF vs binned depth where entries are the number of true positives
    called = pd.pivot_table(df, index = 'depth_bin', columns = 'af_bin',
                            values = ['is_tp'], aggfunc = sum, fill_value = 0)

    # make a matrix of binned AF vs binned depth where entries are the number of true positives and false negatives
    total = pd.pivot_table(df, index = 'depth_bin', columns = 'af_bin',
                           values = [STATUS_COLUMN], aggfunc = len, fill_value = 0)

    sensitivity = (called.values + EPSILON) / (total.values + EPSILON)

    # wrap in pandas dataframe
    af_bins = np.sort(df['af_bin'].unique())
    depth_bins = np.sort(df['depth_bin'].unique())
    sensitivity_table = pd.DataFrame(sensitivity, columns = [str(x) for x in af_bins], index = total.index)

    error_bars = [proportion_confint(x, max(y,1), alpha = 0.05, method = 'wilson') for x,y in zip(called.values.flat, total.values.flat)]

    #reshape confidence intervals
    num_rows, num_cols = called.shape
    new_shape = (num_rows, num_cols * 2)
    error_bars = np.array(error_bars).reshape(new_shape)

    #load into dataframe
    error_bars = pd.DataFrame(error_bars, index = total.index,
                              columns = flatten([(str(x) + "_LCI", str(x) + "_UCI") for x in af_bins]))

    sensitivity_with_confidence = sensitivity_table.merge(error_bars, left_index = True, right_index = True)

    # merge the sensitivity with the error bars, the counts total true variants, and the counts of called true variants into a single data frame
    total = pd.DataFrame(total.values, columns = [str(x) + "_total_count" for x in af_bins], index = total.index)
    called = pd.DataFrame(called.values, columns = [str(x) + "_called_count" for x in af_bins], index = total.index)
    sensitivity_with_confidence_and_truth_count = sensitivity_with_confidence.merge(total, left_index = True, right_index = True)
    sensitivity_with_confidence_and_counts = sensitivity_with_confidence_and_truth_count.merge(called, left_index = True, right_index = True)

    sensitivity_with_confidence_and_counts.to_csv( snp_or_indel + "_sensitivity.tsv", sep = '\t')

    return sensitivity_with_confidence_and_counts, af_bins

def draw_sensitivity_graph(df, af_bins, snp_or_indel):
    x = df.index    # depth

    for _bin in af_bins:
        y = df['{0}'.format(_bin)]
        low, high = y - df['{0}_LCI'.format(_bin)], df['{0}_UCI'.format(_bin)] - y
        plt.errorbar(x,y, yerr = [low, high], label = "{0}".format(_bin))

    axes = plt.gca()
    axes.set_ylim([0.0, 1.0])

    plt.legend(loc = 4)
    plt.title( snp_or_indel + " Sensitivity")
    plt.xlabel('Depth')
    plt.ylabel('Sensitivity')
    plt.savefig( snp_or_indel + "_sensitivity.png")
    plt.close()

def run():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file")
    parser.add_argument("--depth_bins", nargs='+', type=int)
    parser.add_argument("--depth_bin_width", type=int)
    args = parser.parse_args()

    df = pd.read_table(args.input_file, sep = '\t', dtype = {'CHROM': str, 'TYPE': str, 'STATUS': str})

    snps_df = df.loc[df['TYPE'] == 'SNP']
    indels_df = df.loc[df['TYPE'] == 'INDEL']

    snp_df, snp_af_bins = generate_sensitivity_table(snps_df, "SNP", args.depth_bins, args.depth_bin_width)
    indel_df, indel_af_bins = generate_sensitivity_table(indels_df, "Indel", args.depth_bins, args.depth_bin_width)

    draw_sensitivity_graph(snp_df, snp_af_bins, "SNP")
    draw_sensitivity_graph(indel_df, indel_af_bins, "Indel")

if __name__ == '__main__':
    run()