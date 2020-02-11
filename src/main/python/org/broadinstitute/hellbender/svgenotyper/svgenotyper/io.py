import pandas as pd
from pysam import VariantFile

from . import constants
from .constants import SVTypes
from .preprocess import create_tensors, compute_preprocessed_tensors


def load_data(vcf_path: str, mean_coverage_path: str, svtype: SVTypes, num_states: int, depth_dilution_factor: float,
              device: str = 'cpu', logging=True):
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
        if svtype_i != svtype:
            continue
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

    gt_t, pe_t, sr1_t, sr2_t, ncn_t, cnlp_t, svlen_t, svtype_t, mean_count_t, samples_np, vids_np = \
        create_tensors(vids_list, gt_list, data_list, cnlp_list, svlen_list, svtype_list, samples_list, mean_count_df, device)

    return compute_preprocessed_tensors(num_states, svtype, depth_dilution_factor, pe_t, sr1_t, sr2_t, mean_count_t,
                                        cnlp_t, ncn_t, svtype_t, device)
