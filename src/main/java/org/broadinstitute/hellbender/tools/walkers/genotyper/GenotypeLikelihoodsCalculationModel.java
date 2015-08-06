package org.broadinstitute.hellbender.tools.walkers.genotyper;

/**
 * The model representing how we calculate genotype likelihoods
 */
public enum GenotypeLikelihoodsCalculationModel {
    SNP,
    INDEL,
    GENERALPLOIDYSNP,
    GENERALPLOIDYINDEL,
    BOTH;
}
