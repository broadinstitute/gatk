import pandas as pd
from pysam import VariantFile

from . import constants
from .constants import SVTypes


def load_data(vcf_path, mean_coverage_path, logging=True):
    mean_count_df = pd.read_csv(mean_coverage_path, sep='\t', header=None)

    vcf = VariantFile(vcf_path)
    samples_list = list(vcf.header.samples)

    vids_list = []
    gt_list = []
    data_list = []
    cnlp_list = []
    svlen_list = []
    svtype_list = []

    num_processed = 0
    for record in vcf.fetch():
        # TODO: generate genotypes from depth posteriors
        if record.info['ALGORITHMS'] == 'depth':
            continue
        svtype_i = SVTypes[record.info['SVTYPE']]
        vids_list.append(record.id)
        gt_list.append([sum(record.samples[sample]['GT']) for sample in samples_list])
        data_list.append([[record.samples[sample][attr] for attr in constants.GENOTYPE_FIELDS] for sample in samples_list])
        cnlp_list.append([record.samples[sample]['CNLP'] for sample in samples_list])
        svlen_list.append(record.info['SVLEN'])
        if svtype_i != SVTypes.BND:
            svtype_list.append(svtype_i)
        else:
            strands = record.info['STRANDS']
            if strands == '+-':
                svtype_list.append(SVTypes.DEL)
            elif strands == '-+':
                svtype_list.append(SVTypes.DUP)
            elif strands == '++' or strands == '--':
                svtype_list.append(SVTypes.INV)
            else:
                raise ValueError('Invalid STRANDS value: ' + strands)
        num_processed += 1
        if logging and num_processed % 1000 == 0:
            print("Read {:d} vcf records...".format(num_processed))

    return vids_list, gt_list, data_list, cnlp_list, svlen_list, svtype_list, samples_list, mean_count_df
