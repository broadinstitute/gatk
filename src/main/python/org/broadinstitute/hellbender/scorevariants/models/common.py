import torch
from scorevariants.encoders import variant_label_map

def get_snp_indel_split(batch):
    indel_idx = (batch['type'] == 1)
    snp_idx = (batch['type'] == 0)
    return snp_idx, indel_idx

def predictions_to_snp_scores(predictions, eps=1e-7):
    snp = predictions[:, variant_label_map['SNP']]
    not_snp = predictions[:, variant_label_map['NOT_SNP']]
    return torch.log(eps + snp / (not_snp + eps))

def predictions_to_indel_scores(predictions, eps=1e-7):
    indel = predictions[:, variant_label_map['INDEL']]
    not_indel = predictions[:, variant_label_map['NOT_INDEL']]
    return torch.log(eps + indel / (not_indel + eps))

def predictions_to_score(predictions, types):
    indel_scores = predictions_to_indel_scores(predictions)
    snp_scores = predictions_to_snp_scores(predictions)
    indel_idx = (types == 1)
    snp_idx = (types == 0)
    scores = snp_scores
    scores[indel_idx] = indel_scores[indel_idx]
    scores[~snp_idx & ~indel_idx] = torch.maximum(indel_scores, snp_scores)[~snp_idx & ~indel_idx]
    return scores

