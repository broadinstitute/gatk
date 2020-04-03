import logging
import os
import numpy as np
import pandas as pd
from pysam import VariantFile
import torch

from gatktool import tool

from . import constants
from .constants import SVTypes
from .preprocess import compute_preprocessed_tensors
from .model import SVGenotyperData


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
        # TODO: handle depth-only variants
        if record.info['ALGORITHMS'] == 'depth':
            continue
        if SVTypes[record.info['SVTYPE']] != SVTypes.BND:
            svtype_i = SVTypes[record.info['SVTYPE']]
        else:
            strands = record.info['STRANDS']
            if strands == '+-':
                svtype_i = SVTypes.DEL
            elif strands == '-+':
                svtype_i = SVTypes.DUP
            elif strands == '++' or strands == '--':
                svtype_i = SVTypes.INV
            else:
                raise ValueError('Invalid STRANDS value: ' + strands)
        if svtype_i != svtype:
            continue
        vids_list.append(record.id)
        gt_list.append([sum(record.samples[sample]['GT']) for sample in samples_list])
        data_list.append([[record.samples[sample][attr] for attr in constants.INPUT_GENOTYPE_FIELDS] for sample in samples_list])
        cnlp_list.append([record.samples[sample]['CNLP'] for sample in samples_list])
        svlen_list.append(record.info['SVLEN'])
        svtype_list.append(svtype_i.value)
        num_processed += 1
        if num_processed % 1000 == 0:
            logging.info("Read {:d} vcf records...".format(num_processed))
    vcf.close()
    return vids_list, samples_list, gt_list, data_list, cnlp_list, svlen_list, svtype_list


def load_batch(batch_size: int,
               device: str):
    vid_list = []
    pe_list = []
    sr1_list = []
    sr2_list = []
    ncn_list = []
    cnlp_list = []
    for _ in range(batch_size):
        fifo_line = tool.readDataFIFO()
        fifo_data = fifo_line.split('\t')
        vid_list.append(fifo_data[0])
        pe_list.append([int(x)] for x in fifo_data[1].split(';'))
        sr1_list.append([int(x)] for x in fifo_data[2].split(';'))
        sr2_list.append([int(x)] for x in fifo_data[3].split(';'))
        ncn_list.append([int(x)] for x in fifo_data[4].split(';'))
        cnlp_list.append([int(x)] for x in fifo_data[5].split(';'))
    vid_np = np.asarray(vid_list)
    pe_t = torch.tensor(pe_list, device=device)
    sr1_t = torch.tensor(sr1_list, device=device)
    sr2_t = torch.tensor(sr2_list, device=device)
    ncn_t = torch.tensor(ncn_list, device=device)
    cnlp_t = torch.tensor(cnlp_list, device=device)
    return vid_np, pe_t, sr1_t, sr2_t, ncn_t, cnlp_t


def load_data(batch_size: int,
              mean_coverage_path: str,
              samples_path: str,
              svtype: SVTypes,
              num_states: int,
              depth_dilution_factor: float,
              device: str = 'cpu'):
    mean_count_df = pd.read_csv(mean_coverage_path, sep='\t', header=None, index_col=0)
    mean_count_t = torch.from_numpy(mean_count_df.values).to(device=device).squeeze(-1)
    samples_np = np.loadtxt(samples_path)
    vids_np, pe_t, sr1_t, sr2_t, ncn_t, cnlp_t = load_batch(batch_size=batch_size, device=device)
    if vids_np.shape[0] == 0:
        return None
    pe_t, sr1_t, sr2_t, depth_t, rd_gt_prob_t = compute_preprocessed_tensors(num_states, svtype, depth_dilution_factor,
                                                                             pe_t, sr1_t, sr2_t, mean_count_t, cnlp_t,
                                                                             ncn_t, device)
    return SVGenotyperData(svtype, vids_np, samples_np, pe_t, sr1_t, sr2_t, depth_t, rd_gt_prob_t)


def validate_mean_count(mean_count_df, samples_np):


def write_vcf(input_vcf_path: str, output_vcf_path: str, output_data: dict, global_stats: dict):
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
    header.info.add("PHI_PE", "1", "Float", "Paired-end mean bias")
    header.info.add("PHI_SR1", "1", "Float", "First split mean bias")
    header.info.add("PHI_SR2", "1", "Float", "Second split mean bias")
    if "PL" not in header.formats:
        header.formats.add("PL", "1", "Float", "Genotype probability")
    if "GQ" not in header.formats:
        header.formats.add("GQ", "1", "Float", "Genotype quality")

    for svtype in global_stats:
        items = {key: "{:0.4f}".format(float(global_stats[svtype][key])) for key in global_stats[svtype]}.items()
        header.add_meta(key="SVGenotyper_" + svtype.name, items=items)
    
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
        record.info['PHI_PE'] = vid_data['phi_pe']
        record.info['PHI_SR1'] = vid_data['phi_sr1']
        record.info['PHI_SR2'] = vid_data['phi_sr2']
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


def save_tensors(data: SVGenotyperData, base_path: str):
    data_vars = vars(data)
    for var in data_vars:
        torch.save(data_vars[var], base_path + "." + var + ".pt")


def save_list(data: list, path: str):
    with open(path, 'w') as f:
        f.writelines([x + '\n' for x in data])


def load_tensors(directory: str, model_name: str, svtype: SVTypes, device: str = 'cpu'):
    base_path = os.path.join(directory, model_name + "." + svtype.name)
    pe_t = torch.load(base_path + ".pe_t.pt", map_location=device)
    sr1_t = torch.load(base_path + ".sr1_t.pt", map_location=device)
    sr2_t = torch.load(base_path + ".sr2_t.pt", map_location=device)
    depth_t = torch.load(base_path + ".depth_t.pt", map_location=device)
    rd_gt_prob_t = torch.load(base_path + ".rd_gt_prob_t.pt", map_location=device)
    return SVGenotyperData(pe_t=pe_t, sr1_t=sr1_t, sr2_t=sr2_t, depth_t=depth_t, rd_gt_prob_t=rd_gt_prob_t)


def load_list(path: str):
    with open(path, 'r') as f:
        lines = f.readlines()
        if lines is None or len(lines) == 0:
            raise ValueError("Tried to read empty list: {:s}".format(path))
        return [x.strip() for x in lines]
