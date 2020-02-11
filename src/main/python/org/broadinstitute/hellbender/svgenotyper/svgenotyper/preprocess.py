import numpy as np
import pandas as pd
import pyro
import torch

from . import constants
from .constants import SVTypes
from .model import SVGenotyperData


def create_tensors(vids_list: list,
                   gt_list: list,
                   data_list: list,
                   cnlp_list: list,
                   svlen_list: list,
                   svtype_list: list,
                   samples_list: list,
                   mean_count_df: pd.DataFrame,
                   device_name: str,
                   default_dtype: torch.dtype):
    device = torch.device(device_name)
    torch.set_default_dtype(default_dtype)

    data_np = np.array(data_list)
    gt_np = np.array(gt_list)
    cnlp_np = np.array(cnlp_list)

    gt_t = torch.from_numpy(gt_np).long()
    pe_t = torch.from_numpy(data_np[..., constants.GENOTYPE_FIELDS.index('PE')]).to(dtype=torch.get_default_dtype(), device=device)
    sr1_t = torch.from_numpy(data_np[..., constants.GENOTYPE_FIELDS.index('SSR')]).to(dtype=torch.get_default_dtype(), device=device)
    sr2_t = torch.from_numpy(data_np[..., constants.GENOTYPE_FIELDS.index('ESR')]).to(dtype=torch.get_default_dtype(), device=device)
    ncn_t = torch.from_numpy(data_np[..., constants.GENOTYPE_FIELDS.index('NCN')]).long().to(device=device)
    cnlp_t = torch.from_numpy(cnlp_np).to(dtype=torch.get_default_dtype(), device=device)

    svlen_t = torch.from_numpy(np.array(svlen_list)).to(dtype=torch.get_default_dtype(), device=device)
    svtype_t = torch.from_numpy(np.array(svtype_list)).long().to(device=device).squeeze(-1)
    mean_count_t = torch.from_numpy(mean_count_df.values).to(dtype=torch.get_default_dtype(), device=device).squeeze(-1)

    samples_np = np.array(samples_list, dtype='str')
    vids_np = np.array(vids_list, dtype='str')

    return gt_t, pe_t, sr1_t, sr2_t, ncn_t, cnlp_t, svlen_t, svtype_t, mean_count_t, samples_np, vids_np


def compute_preprocessed_tensors(k: int,
                                 svtype: SVTypes,
                                 depth_dilution_factor: float,
                                 pe_t: torch.Tensor,
                                 sr1_t: torch.Tensor,
                                 sr2_t: torch.Tensor,
                                 mean_count_t: torch.Tensor,
                                 cnlp_t: torch.Tensor,
                                 ncn_t: torch.Tensor,
                                 svtype_t: torch.Tensor,
                                 device: str) -> SVGenotyperData:
    # Per-base depth
    depth_t = 0.5 * mean_count_t.squeeze(-1).to(dtype=torch.get_default_dtype(), device=device)

    # Copy state posterior probabilities
    cn_prob_t = torch.exp(cnlp_t / (-10. * torch.log10(torch.tensor(np.exp(1.), device=device))))

    # Filter PE for insertions
    pe_t = pe_t.where(svtype_t.ne(SVTypes.INS).unsqueeze(-1), torch.zeros(pe_t.shape, device=device))

    # Clamp values for numerical stability
    pe_t = pe_t.clone().clamp(0, constants.MAX_PE_COUNT)
    sr1_t = sr1_t.clone().clamp(0, constants.MAX_SR_COUNT)
    sr2_t = sr2_t.clone().clamp(0, constants.MAX_SR_COUNT)

    # TODO : we only use CNLP values from the VCF for CNVs
    rd_gt_prob_t = torch.ones((cn_prob_t.shape[0], cn_prob_t.shape[1], k), device=device) / float(k)

    if svtype == SVTypes.DEL:
        n = cn_prob_t.shape[0]
        m = cn_prob_t.shape[1]
        cn_prob_states = cn_prob_t.shape[2]
        filled_zeros = torch.zeros(n, m, cn_prob_states, device=device)
        filled_index = torch.arange(cn_prob_states, device=device).repeat(n, m, 1)
        hom_ref_prob = torch.where(filled_index >= ncn_t.unsqueeze(-1), cn_prob_t, filled_zeros).sum(dim=-1).unsqueeze(-1)
        probs_list = [hom_ref_prob]
        for i in range(1, k - 1):
            probs_list.append(torch.where(filled_index == ncn_t.unsqueeze(-1) - i, cn_prob_t, filled_zeros).sum(dim=-1).unsqueeze(-1))
        probs_list.append(torch.where(filled_index <= ncn_t.unsqueeze(-1) - k + 1, cn_prob_t, filled_zeros).sum(dim=-1).unsqueeze(-1))
        rd_gt_prob_t = torch.cat(probs_list, dim=-1)
        del filled_zeros, filled_index, probs_list

    if svtype == SVTypes.DUP:
        n = cn_prob_t.shape[0]
        m = cn_prob_t.shape[1]
        cn_prob_states = cn_prob_t.shape[2]
        filled_zeros = torch.zeros(n, m, cn_prob_states, device=device)
        filled_index = torch.arange(cn_prob_states, device=device).repeat(n, m, 1)
        hom_ref_prob = torch.where(filled_index <= ncn_t.unsqueeze(-1), cn_prob_t, filled_zeros).sum(dim=-1).unsqueeze(-1)
        probs_list = [hom_ref_prob]
        for i in range(1, k - 1):
            probs_list.append(torch.where(filled_index == ncn_t.unsqueeze(-1) + i, cn_prob_t, filled_zeros).sum(dim=-1).unsqueeze(-1))
        probs_list.append(torch.where(filled_index >= ncn_t.unsqueeze(-1) + k - 1, cn_prob_t, filled_zeros).sum(dim=-1).unsqueeze(-1))
        rd_gt_prob_t = torch.cat(probs_list, dim=-1)
        del filled_zeros, filled_index, probs_list

    rd_gt_prob_t += depth_dilution_factor
    rd_gt_prob_t = rd_gt_prob_t / rd_gt_prob_t.sum(dim=-1).unsqueeze(-1)

    return SVGenotyperData(pe_t, sr1_t, sr2_t, depth_t, rd_gt_prob_t)
