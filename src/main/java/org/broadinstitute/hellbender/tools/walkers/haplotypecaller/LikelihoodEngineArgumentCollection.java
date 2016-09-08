package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.broadinstitute.hellbender.cmdline.Advanced;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.Hidden;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.pairhmm.PairHMM;

/**
 * Set of arguments related to {@link ReadLikelihoodCalculationEngine} implementations
 */
public final class LikelihoodEngineArgumentCollection {

    @Hidden
    @Advanced
    @Argument(fullName = "likelihoodCalculationEngine", shortName = "likelihoodEngine",
            doc= "What likelihood calculation engine to use to calculate the relative likelihood of reads vs haplotypes", optional = true)
    public ReadLikelihoodCalculationEngine.Implementation likelihoodEngineImplementation = ReadLikelihoodCalculationEngine.Implementation.PairHMM;

    /**
     * Bases with a quality below this threshold will reduced to the minimum usable qualiy score (6).
     */
    @Argument(fullName = "base_quality_score_threshold", shortName = "bqst", doc = "Base qualities below this threshold will be reduced to the minimum (" + QualityUtils.MIN_USABLE_Q_SCORE + ")", optional = true)
    public byte BASE_QUALITY_SCORE_THRESHOLD = PairHMM.BASE_QUALITY_SCORE_THRESHOLD;

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
     * When calculating the likelihood of variants, we can try to correct for PCR errors that cause indel artifacts.
     * The correction is based on the reference context, and acts specifically around repetitive sequences that tend
     * to cause PCR errors). The variant likelihoods are penalized in increasing scale as the context around a
     * putative indel is more repetitive (e.g. long homopolymer). The correction can be disabling by specifying
     * '-pcrModel NONE'; in that case the default base insertion/deletion qualities will be used (or taken from the
     * read if generated through the BaseRecalibrator). <b>VERY IMPORTANT: when using PCR-free sequencing data we
     * definitely recommend setting this argument to NONE</b>.
     */
    @Advanced
    @Argument(fullName = "pcr_indel_model", shortName = "pcrModel", doc = "The PCR indel model to use", optional = true)
    public PairHMMLikelihoodCalculationEngine.PCRErrorModel pcrErrorModel = PairHMMLikelihoodCalculationEngine.PCRErrorModel.CONSERVATIVE;


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
