def variant_is_snp(variant):
    """Return whether or not the variant is a SNP

    Args:
        pysam.VariantRecord containing the variant of interest

    """
    if len(variant.ref) > 1:
        return False
    if not variant.alts:
        return False
    for alt in variant.alts:
        if alt is None:
            return False
        if alt not in ["A", "C", "G", "T", "N", "*"]:
            return False
    return True
