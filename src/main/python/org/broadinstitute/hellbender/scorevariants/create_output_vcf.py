#!/usr/bin/python3

from pysam import VariantFile
import re
import argparse
import sys

CONTIG_INDEX = 0;
POS_INDEX = 1;
REF_INDEX = 2;
ALT_INDEX = 3;
KEY_INDEX = 4;

def create_output_vcf(vcf_in, scores_file, vcf_out, label):
    variant_file = VariantFile(vcf_in)
    variant_file.reset()

    variant_file.header.info.add(id=label, number=1, type='Float', description='Log odds of being a true variant versus \
                    being false under the trained Convolutional Neural Network')
    header = variant_file.header.copy()
    vcfWriter = VariantFile(vcf_out, 'w', header=header)

    with open(scores_file) as scoredVariants:
        sv = next(scoredVariants)
        for variant in variant_file:
            scoredVariant = sv.split('\t')
            if variant.contig == scoredVariant[CONTIG_INDEX] and \
               variant.pos == int(scoredVariant[POS_INDEX]) and \
               variant.ref == scoredVariant[REF_INDEX] and \
               ', '.join(variant.alts or []) == re.sub('[\[\]]', '', scoredVariant[ALT_INDEX]):

                    if len(scoredVariant) > KEY_INDEX:
                        variant.info.update({label: float(scoredVariant[KEY_INDEX])})

                    vcfWriter.write(variant)

                    sv = next(scoredVariants, None)
            else:
                sys.exit("Score file out of sync with original VCF. Score file has: " + sv + "\nBut VCF has: " + str(variant))
