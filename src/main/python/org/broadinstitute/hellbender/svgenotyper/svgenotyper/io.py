import logging
import pandas as pd
from pysam import VariantFile

from . import constants
from .constants import SVTypes
from .preprocess import create_tensors, compute_preprocessed_tensors


def read_vcf(vcf_path: str, svtype: SVTypes):
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
        data_list.append([[record.samples[sample][attr] for attr in constants.INPUT_GENOTYPE_FIELDS] for sample in samples_list])
        cnlp_list.append([record.samples[sample]['CNLP'] for sample in samples_list])
        svlen_list.append(record.info['SVLEN'])
        if svtype_i != SVTypes.BND:
            svtype_list.append(svtype_i.value)
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
        if num_processed % 1000 == 0:
            logging.info("Read {:d} vcf records...".format(num_processed))
    vcf.close()
    return vids_list, samples_list, gt_list, data_list, cnlp_list, svlen_list, svtype_list


def load_data(vcf_path: str, mean_coverage_path: str, svtype: SVTypes, num_states: int, depth_dilution_factor: float,
              device: str = 'cpu'):
    mean_count_df = pd.read_csv(mean_coverage_path, sep='\t', header=None, index_col=0)
    vids_list, samples_list, gt_list, data_list, cnlp_list, svlen_list, svtype_list = read_vcf(vcf_path=vcf_path, svtype=svtype)
    if len(vids_list) == 0:
        return None, None
    gt_t, pe_t, sr1_t, sr2_t, ncn_t, cnlp_t, svlen_t, svtype_t, mean_count_t, samples_np, vids_np = create_tensors(
        vids_list, gt_list, data_list, cnlp_list, svlen_list, svtype_list, samples_list, mean_count_df, device)
    data = compute_preprocessed_tensors(num_states, svtype, depth_dilution_factor, pe_t, sr1_t, sr2_t, mean_count_t,
                                        cnlp_t, ncn_t, device)
    return vids_np, data


def write_output(input_vcf_path: str, output_vcf_path: str, output_data: dict):
    vcf_in = VariantFile(input_vcf_path)
    samples_list = list(vcf_in.header.samples)
    n_samples = len(samples_list)

    header = vcf_in.header
    header.info.add("PPE", "1", "Float", "Paired-end read support probability")
    header.info.add("PSR1", "1", "Float", "First split read support probability")
    header.info.add("PSR2", "1", "Float", "Second split read support probability")
    header.info.add("PRD", "1", "Float", "Read depth support probability")
    header.info.add("EPE", "1", "Float", "Paired-end read background")
    header.info.add("ESR1", "1", "Float", "First split read background")
    header.info.add("ESR2", "1", "Float", "Second split read background")
    header.formats.add("PL", "1", "Float", "Genotype probability")
    header.formats.add("GQ", "1", "Float", "Genotype quality (log odds ratio)")
    
    vcf_out = VariantFile(output_vcf_path, 'w', header=vcf_in.header)
    for record in vcf_in.fetch():
        # TODO: generate genotypes from depth posteriors
        if record.info['ALGORITHMS'] == 'depth':
            continue
        vid = record.id
        if vid not in output_data:
            logging.warning("Omitting missing VID: {:s}".format(vid))
            continue
        vid_data = output_data[vid]
        record.info['PPE'] = vid_data['p_m_pe']
        record.info['PSR1'] = vid_data['p_m_sr1']
        record.info['PSR2'] = vid_data['p_m_sr2']
        record.info['PRD'] = vid_data['p_m_rd']
        record.info['EPE'] = vid_data['eps_pe']
        record.info['ESR1'] = vid_data['eps_sr1']
        record.info['ESR2'] = vid_data['eps_sr2']
        gt = vid_data['gt']
        gt_p = vid_data['gt_p']
        gt_lod = vid_data['gt_lod']
        if n_samples != gt.size:
            raise ValueError("Samples list size {:d} did not match gt size {:d}".format(n_samples, gt.size))
        if n_samples != gt_p.size:
            raise ValueError("Samples list size {:d} did not match gt_p size {:d}".format(n_samples, gt_p.size))
        if n_samples != gt.size:
            raise ValueError("Samples list size {:d} did not match gt_lod size {:d}".format(n_samples, gt_lod.size))
        for i in range(n_samples):
            sample = samples_list[i]
            if gt[i] == 0:
                alleles = (0, 0)
            elif gt[i] == 1:
                alleles = (0, 1)
            else:
                alleles = (1, 1)
            record.samples[sample]['GT'] = alleles
            record.samples[sample]['PL'] = gt_p[i]
            record.samples[sample]['GQ'] = gt_lod[i]
        vcf_out.write(record)
    vcf_in.close()
    vcf_out.close()
