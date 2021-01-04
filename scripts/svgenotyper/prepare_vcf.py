from pysam import VariantFile, VariantHeader
import argparse
import sys

def parse_args():
    """Parse command line arguments.
    """
    parser = argparse.ArgumentParser()

    parser.add_argument('--vcf', help='Genotyped vcf', required=True)

    args = parser.parse_args()

    return args


def main():
    args = parse_args()
    vcf = VariantFile(args.vcf)

    header_str = str(vcf.header) \
        .replace('##FORMAT=<ID=CNLP,Number=1,Type=Integer,Description="Phred-scaled copy number posterior over the event region">',
                 '##FORMAT=<ID=CNLP,Number=5,Type=Integer,Description="Phred-scaled copy number posterior over the event region">') \
        .replace('##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">',
                 '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">')
    sys.stdout.write(header_str)

    for record in vcf.fetch():
        if record.info['SVTYPE'] == 'DUP':
            for s in record.samples:
                n_non_ref_alleles = int(record.samples[s]['CN']) - int(record.samples[s]['NCN'])
                if len(record.samples[s]['GT']) == 1:
                    if n_non_ref_alleles == 0:
                        record.samples[s]['GT'] = (record.alleles[1])
                    else:
                        record.samples[s]['GT'] = (record.alleles[0])
                else:
                    if n_non_ref_alleles == 0:
                        record.samples[s]['GT'] = (record.alleles[1], record.alleles[1])
                    elif n_non_ref_alleles == 1:
                        record.samples[s]['GT'] = (record.alleles[1], record.alleles[0])
                    else:
                        # TODO : greater copy numbers not supported
                        record.samples[s]['GT'] = (record.alleles[0], record.alleles[0])
        sys.stdout.write(str(record))


if __name__ == "__main__":
    main()
