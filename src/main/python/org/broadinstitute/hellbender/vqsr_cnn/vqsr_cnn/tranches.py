import os
import pysam
from collections import Counter


# Package Imports
from vqsr_cnn import defines
from vqsr_cnn import arguments

def run():
    args = arguments.parse_args()
    if 'write_snp_tranches' == args.mode:
        write_tranches(args)
    elif 'write_indel_tranches' == args.mode:
            write_tranches(args, variant_type='INDEL')
    else:
        raise ValueError('Unknown tranche mode:', args.mode)


def write_tranches(args, variant_type='SNP'):
    score_key = args.score_keys[0]
    tranches = sorted(args.tranches, reverse=True)

    vcf_truth = pysam.VariantFile(args.train_vcf, 'rb')
    vcf_reader = pysam.VariantFile(args.input_vcf, 'r')

    print('Writing sensitivity tranches:', tranches, ' with score_key:', score_key, 'for', variant_type, 'variants.')
    stats = Counter()
    scores = []

    print('Iterate over tranche truth VCF.')
    for variant in vcf_truth:
        if not variant.alts:
            stats['Variant without alts'] += 1
            continue
        for allele_idx, allele in enumerate(variant.alts):
            v_scored = allele_in_vcf(allele, variant, vcf_reader)
            if not v_scored:
                stats['Variant not in input VCF'] += 1
                continue
            scores.append(v_scored.info[score_key])
            stats['count'] += 1
            if stats['count']%2000==0:
                for k in stats.keys():
                    print(k, ' has:', stats[k])
        if stats['count'] >= args.samples:
            break

    scores = sorted(scores, reverse=True)
    thresholds = []

    for i,t in enumerate(tranches):
        thresholds.append(scores[int(clamp(t*len(scores), 0, len(scores)-1))])
        if len(thresholds) == 1:
            filter_id = score_key + 'Tranche' + variant_type + format_tranche(t) + 'to100.0'
            vcf_reader.header.filters.add(filter_id, None, None, 'Tranche filtering from '+score_key)
        else:
            filter_id = score_key + 'Tranche' + variant_type + format_tranche(t) + 'to' + format_tranche(tranches[i-1])
            vcf_reader.header.filters.add(filter_id, None, None, 'Tranche filtering from '+score_key)
        print('At tranch:', t, ' score key:', score_key, ' cutoff is:', thresholds[-1], 'filter id is:', filter_id)

    vcf_writer = pysam.VariantFile(args.output_vcf, 'w', header=vcf_reader.header)
    for variant in vcf_reader:
        if score_key in variant.info:
            if variant_type == 'SNP' and is_snp(variant) or variant_type == 'INDEL' and is_indel(variant):
                variant.filter.add(score_to_filter_status(variant.info[score_key], score_key, variant_type, tranches, thresholds))

        vcf_writer.write(variant)

    for k in stats.keys():
        print(k, ' has:', stats[k])


def allele_in_vcf(allele, variant, vcf_ram):
    ''' Check if variant's allele is in a VCF file.

    Arguments
        allele: the allele from the provided variant that we are checking
        variant: the variant whose allele we are looking for
        vcf_ram: the VCF we look in, must have an index (tbi, or idx)

    Returns
        variant if it is found otherwise None
    '''
    variants = vcf_ram.fetch(variant.contig, variant.pos-1, variant.pos)

    for v in variants:
        if v.contig == variant.contig and v.pos == variant.pos and allele in v.alts:
            return v

    return None


def score_to_filter_status(score, score_key, variant_type, tranches, thresholds):

    for i, threshold in enumerate(thresholds):
        if score < threshold:
            if i == 0:
                high_thresh = '100.0'
            else:
                high_thresh = format_tranche(tranches[i-1])
            low_thresh = format_tranche(tranches[i])

            return score_key + 'Tranche' + variant_type + low_thresh + 'to' + high_thresh

    return 'PASS'


def format_tranche(t):
    return '{0:.1f}'.format(t*100)

def is_snp(variant):
    return len(variant.ref) == 1 and any(map(lambda x: len(x) == 1, variant.alts))

def is_indel(variant):
    return len(variant.ref) > 1 or any(map(lambda x: len(x) > 1, variant.alts))


def clamp(n, minn, maxn):
    return max(min(maxn, n), minn)


# Back to the top!
if "__main__" == __name__:
    run()