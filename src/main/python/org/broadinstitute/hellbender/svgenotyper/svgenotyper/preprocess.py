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
                                 device: str,
                                 tensor_dtype: torch.dtype) -> SVGenotyperData:
    # Per-base depth
    depth_t = mean_count_t.squeeze(-1).to(device=device, dtype=tensor_dtype)

    # Copy state posterior probabilities
    cn_prob_t = torch.exp(cnlp_t / (-10. * torch.log10(torch.tensor(np.exp(1.))))).to(device=device, dtype=tensor_dtype)

    # Clamp values for numerical stability
    pe_t = pe_t.clamp(0, constants.MAX_PE_COUNT)
    sr1_t = sr1_t.clamp(0, constants.MAX_SR_COUNT)
    sr2_t = sr2_t.clamp(0, constants.MAX_SR_COUNT)

    # Map from "local copy states" used in CN model to "genotype state" indexing in the final tensor dimension
    # e.g. for DELs on diploid chromosomes, map copy state 0 -> event state 2 (hom alt), 1 -> 1 (het), 2 -> 0 (hom ref)
    if svtype == SVTypes.DEL:
        n = cn_prob_t.shape[0]  # Number of sites
        m = cn_prob_t.shape[1]  # Number of samples
        cn_prob_states = cn_prob_t.shape[2]  # Number of states from copy number model
        filled_zeros = torch.zeros(n, m, cn_prob_states, device=device, dtype=tensor_dtype)  # All-zeros tensor
        filled_index = torch.arange(cn_prob_states, device=device, dtype=tensor_dtype).repeat(n, m, 1)  # Last dimension is [0, 1, ..., num_states-1]

        # CN model hom-ref probability where state is at least the local ploidy
        hom_ref_prob = torch.where(filled_index >= ncn_t.unsqueeze(-1), cn_prob_t, filled_zeros).sum(dim=-1).unsqueeze(-1)
        probs_list = [hom_ref_prob]
        for i in range(1, num_states - 1):
            # CN model probability where state is local ploidy minus i
            probs_list.append(torch.where(filled_index == ncn_t.unsqueeze(-1) - i, cn_prob_t, filled_zeros).sum(dim=-1).unsqueeze(-1))
        # CN model probability where state is at most local ploidy minus number genotyping max state
        probs_list.append(torch.where(filled_index <= ncn_t.unsqueeze(-1) - num_states + 1, cn_prob_t, filled_zeros).sum(dim=-1).unsqueeze(-1))
        rd_gt_prob_t = torch.cat(probs_list, dim=-1)
        del filled_zeros, filled_index, probs_list
    elif svtype == SVTypes.DUP:
        n = cn_prob_t.shape[0]  # Number of sites
        m = cn_prob_t.shape[1]  # Number of samples
        cn_prob_states = cn_prob_t.shape[2]  # Number of states from copy number model
        filled_zeros = torch.zeros(n, m, cn_prob_states, device=device, dtype=tensor_dtype)  # All-zeros tensor
        filled_index = torch.arange(cn_prob_states, device=device, dtype=tensor_dtype).repeat(n, m, 1)  # Last dimension is [0, 1, ..., num_states-1]

        # CN model hom-ref probability where state is at most the local ploidy
        hom_ref_prob = torch.where(filled_index <= ncn_t.unsqueeze(-1), cn_prob_t, filled_zeros).sum(dim=-1).unsqueeze(-1)
        probs_list = [hom_ref_prob]
        for i in range(1, num_states - 1):
            # CN model probability where state is local ploidy minus i
            probs_list.append(torch.where(filled_index == ncn_t.unsqueeze(-1) + i, cn_prob_t, filled_zeros).sum(dim=-1).unsqueeze(-1))
        # CN model probability where state is at least local ploidy plus number genotyping max state
        probs_list.append(torch.where(filled_index >= ncn_t.unsqueeze(-1) + num_states - 1, cn_prob_t, filled_zeros).sum(dim=-1).unsqueeze(-1))
        rd_gt_prob_t = torch.cat(probs_list, dim=-1)
        del filled_zeros, filled_index, probs_list
    else:
        # we only use CNLP values from the VCF for CNVs
        rd_gt_prob_t = None

    # Dilute and normalize
    if rd_gt_prob_t is not None:
        rd_gt_prob_t += depth_dilution_factor
        rd_gt_prob_t = rd_gt_prob_t / rd_gt_prob_t.sum(dim=-1).unsqueeze(-1)

    return pe_t, sr1_t, sr2_t, depth_t, rd_gt_prob_t
