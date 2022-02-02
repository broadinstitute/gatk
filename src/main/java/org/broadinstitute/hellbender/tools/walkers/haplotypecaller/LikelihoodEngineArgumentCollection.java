package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.Hidden;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.utils.dragstr.DragstrParams;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.pairhmm.PairHMM;

import java.io.Serializable;

/**
 * Set of arguments related to {@link ReadLikelihoodCalculationEngine} implementations
 */
public final class LikelihoodEngineArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    public static final String BASE_QUALITY_SCORE_THRESHOLD_FULLNAME = "base-quality-score-threshold";
    public static final String DRAGSTR_PARAMS_PATH_FULLNAME = "dragstr-params-path";
    public static final String DRAGSTR_HET_HOM_RATIO_FULLNAME = "dragstr-het-hom-ratio";
    public static final String DONT_USE_DRAGSTR_PAIRHMM_FULLNAME = "dont-use-dragstr-pair-hmm-scores";

    /**
     * Bases with a quality below this threshold will reduced to the minimum usable qualiy score (6).
     */
    @Argument(fullName = BASE_QUALITY_SCORE_THRESHOLD_FULLNAME, doc = "Base qualities below this threshold will be reduced to the minimum (" + QualityUtils.MIN_USABLE_Q_SCORE + ")", optional = true)
    public byte BASE_QUALITY_SCORE_THRESHOLD = PairHMM.BASE_QUALITY_SCORE_THRESHOLD;

    @Argument(fullName = DRAGSTR_PARAMS_PATH_FULLNAME, doc = "location of the DRAGstr model parameters for STR error correction used in the Pair HMM. When provided, it overrides other PCR error correcting mechanisms", optional = true)
    public GATKPath dragstrParams = null;

    @Argument(fullName = DRAGSTR_HET_HOM_RATIO_FULLNAME, doc="het to hom prior ratio use with DRAGstr on", optional = true)
    public int dragstrHetHomRatio = 2;

    @Argument(fullName = DONT_USE_DRAGSTR_PAIRHMM_FULLNAME, doc="disable DRAGstr pair-hmm score even when dragstr-params-path was provided", optional = false)
    public boolean dontUseDragstrPairHMMScores = false;

    @Advanced
    @Argument(fullName="pair-hmm-gap-continuation-penalty", doc="Flat gap continuation penalty for use in the Pair HMM", optional = true)
    public int gcpHMM = 10;

    @Advanced
    @Argument(fullName="expected-mismatch-rate-for-read-disqualification", doc="Error rate used to set expectation for post HMM read disqualification based on mismatches", optional = true)
    public double expectedErrorRatePerBase = PairHMMLikelihoodCalculationEngine.DEFAULT_EXPECTED_ERROR_RATE_PER_BASE;

    /**
     * The PairHMM implementation to use for genotype likelihood calculations. The various implementations balance a tradeoff of accuracy and runtime.
     */
    @Advanced
    @Argument(fullName = "pair-hmm-implementation", shortName = "pairHMM", doc = "The PairHMM implementation to use for genotype likelihood calculations", optional = true)
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
    @Argument(fullName = "pcr-indel-model", doc = "The PCR indel model to use", optional = true)
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
    @Argument(fullName="phred-scaled-global-read-mismapping-rate", doc="The global assumed mismapping rate for reads", optional = true)
    public int phredScaledGlobalReadMismappingRate = 45;

    @Advanced
    @Argument(fullName = "disable-symmetric-hmm-normalizing", doc="Toggle to revive legacy behavior of asymmetrically normalizing the arguments to the reference haplotype", optional = true)
    public boolean disableSymmetricallyNormalizeAllelesToReference = false;

    @Advanced
    @Argument(fullName ="disable-cap-base-qualities-to-map-quality", doc= "If false this disables capping of base qualities in the HMM to the mapping quality of the read", optional = true)
    public boolean disableCapReadQualitiesToMapQ = false;

    /**
     * If enabled, rather than disqualifying all reads over a threshold of minimum hmm scores we will instead choose a less strict
     * and less aggressive cap for disqualification based on the read length and base qualities.
     */
    @Argument(fullName="enable-dynamic-read-disqualification-for-genotyping", doc="Will enable less strict read disqualification low base quality reads")
    public boolean enableDynamicReadDisqualification = false;

    /**
     * Argument used to adjust the agressiveness of dynamic read disqualification
     */
    @Advanced
    @Hidden
    @Argument(fullName="dynamic-read-disqualification-threshold", doc="Constant used to scale the dynamic read disqualificaiton")
    public double readDisqualificationThresholdConstant = PairHMMLikelihoodCalculationEngine.DEFAULT_DYNAMIC_DISQUALIFICATION_SCALE_FACTOR;

    /**
     * Argument for generating a file of all of the inputs and ouptus for the pair hmm
     */
    @Advanced
    @Hidden
    @Argument(fullName="pair-hmm-results-file", doc="Constant used to scale the dynamic read disqualificaiton", optional = true)
    public GATKPath pairHmmResultsFile = null;

    @ArgumentCollection
    public PairHMMNativeArgumentCollection pairHMMNativeArgs = new PairHMMNativeArgumentCollection();
}
