package org.broadinstitute.hellbender.tools.walkers.genotyper;

/**
 * Enumeration of possible genotyping modes.
 *
 * <p>
 *     A genotyping mode represents the general approach taken when detecting and calculate variant genotypes.
 * </p>
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public enum GenotypingOutputMode {

    /**
     * The genotyper will choose the most likely alternate allele
     */
    DISCOVERY,

    /**
     * Only the alleles passed by the user should be considered.
     */
    GENOTYPE_GIVEN_ALLELES
}
