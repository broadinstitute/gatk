import numpy as np
import torch

from . import constants
from .constants import SVTypes
from .model import SVGenotyperData


def compute_preprocessed_tensors(num_states: int,
                                 svtype: SVTypes,
                                 depth_dilution_factor: float,
                                 pe_t: torch.Tensor,
                                 sr1_t: torch.Tensor,
                                 sr2_t: torch.Tensor,
                                 mean_count_t: torch.Tensor,
                                 cnlp_t: torch.Tensor,
                                 ncn_t: torch.Tensor,
                                 device: str) -> SVGenotyperData:
    # Per-base depth
    depth_t = mean_count_t.squeeze(-1).to(dtype=torch.get_default_dtype(), device=device)

    # Copy state posterior probabilities
    cn_prob_t = torch.exp(cnlp_t / (-10. * torch.log10(torch.tensor(np.exp(1.), device=device))))

    # Clamp values for numerical stability
    pe_t = pe_t.clamp(0, constants.MAX_PE_COUNT)
    sr1_t = sr1_t.clamp(0, constants.MAX_SR_COUNT)
    sr2_t = sr2_t.clamp(0, constants.MAX_SR_COUNT)

    if svtype == SVTypes.DEL:
        n = cn_prob_t.shape[0]
        m = cn_prob_t.shape[1]
        cn_prob_states = cn_prob_t.shape[2]
        filled_zeros = torch.zeros(n, m, cn_prob_states, device=device)
        filled_index = torch.arange(cn_prob_states, device=device).repeat(n, m, 1)
        hom_ref_prob = torch.where(filled_index >= ncn_t.unsqueeze(-1), cn_prob_t, filled_zeros).sum(dim=-1).unsqueeze(-1)
        probs_list = [hom_ref_prob]
        for i in range(1, num_states - 1):
            probs_list.append(torch.where(filled_index == ncn_t.unsqueeze(-1) - i, cn_prob_t, filled_zeros).sum(dim=-1).unsqueeze(-1))
        probs_list.append(torch.where(filled_index <= ncn_t.unsqueeze(-1) - num_states + 1, cn_prob_t, filled_zeros).sum(dim=-1).unsqueeze(-1))
        rd_gt_prob_t = torch.cat(probs_list, dim=-1)
        del filled_zeros, filled_index, probs_list
    elif svtype == SVTypes.DUP:
        n = cn_prob_t.shape[0]
        m = cn_prob_t.shape[1]
        cn_prob_states = cn_prob_t.shape[2]
        filled_zeros = torch.zeros(n, m, cn_prob_states, device=device)
        filled_index = torch.arange(cn_prob_states, device=device).repeat(n, m, 1)
        hom_ref_prob = torch.where(filled_index <= ncn_t.unsqueeze(-1), cn_prob_t, filled_zeros).sum(dim=-1).unsqueeze(-1)
        probs_list = [hom_ref_prob]
        for i in range(1, num_states - 1):
            probs_list.append(torch.where(filled_index == ncn_t.unsqueeze(-1) + i, cn_prob_t, filled_zeros).sum(dim=-1).unsqueeze(-1))
        probs_list.append(torch.where(filled_index >= ncn_t.unsqueeze(-1) + num_states - 1, cn_prob_t, filled_zeros).sum(dim=-1).unsqueeze(-1))
        rd_gt_prob_t = torch.cat(probs_list, dim=-1)
        del filled_zeros, filled_index, probs_list
    else:
        # TODO: we only use CNLP values from the VCF for CNVs
        rd_gt_prob_t = torch.ones((cn_prob_t.shape[0], cn_prob_t.shape[1], num_states), device=device) / float(num_states)

    rd_gt_prob_t += depth_dilution_factor
    rd_gt_prob_t = rd_gt_prob_t / rd_gt_prob_t.sum(dim=-1).unsqueeze(-1)

    return pe_t, sr1_t, sr2_t, depth_t, rd_gt_prob_t
