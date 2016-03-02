package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.broadinstitute.hellbender.cmdline.Advanced;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.Hidden;
import org.broadinstitute.hellbender.utils.pairhmm.PairHMM;

/**
 * Set of arguments related to {@link ReadLikelihoodCalculationEngine} implementations
 */
public final class LikelihoodEngineArgumentCollection {
    
    @Advanced
    @Argument(fullName="gcpHMM", shortName="gcpHMM", doc="Flat gap continuation penalty for use in the Pair HMM", optional = true)
    public int gcpHMM = 10;

    /**
     * The PairHMM implementation to use for genotype likelihood calculations. The various implementations balance a tradeoff of accuracy and runtime.
     */
    @Hidden
    @Argument(fullName = "pair_hmm_implementation", shortName = "pairHMM", doc = "The PairHMM implementation to use for genotype likelihood calculations", optional = true)
    public PairHMM.Implementation pairHMM = PairHMM.Implementation.FASTEST_AVAILABLE;

    /**
     * The phredScaledGlobalReadMismappingRate reflects the average global mismapping rate of all reads, regardless of their
     * mapping quality.  This term effects the probability that a read originated from the reference haplotype, regardless of
     * its edit distance from the reference, in that the read could have originated from the reference haplotype but
     * from another location in the genome.  Suppose a read has many mismatches from the reference, say like 5, but
     * has a very high mapping quality of 60.  Without this parameter, the read would contribute 5 * Q30 evidence
     * in favor of its 5 mismatch haplotype compared to reference, potentially enough to make a call off that single
     * read for all of these events.  With this parameter set to Q30, though, the maximum evidence against any haplotype
     * that this (and any) read could contribute is Q30.
     *
     * Set this term to any negative number to turn off the global mapping rate.
     */
    @Advanced
    @Argument(fullName="phredScaledGlobalReadMismappingRate", shortName="globalMAPQ", doc="The global assumed mismapping rate for reads", optional = true)
    public int phredScaledGlobalReadMismappingRate = 45;

}
