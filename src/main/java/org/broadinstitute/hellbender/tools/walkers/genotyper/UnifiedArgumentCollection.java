package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.Hidden;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.utils.pairhmm.PairHMM;

public final class UnifiedArgumentCollection extends StandardCallerArgumentCollection {
    private static final long serialVersionUID = 1L;

    @Argument(fullName = "genotype_likelihoods_model", shortName = "glm", doc = "Genotype likelihoods calculation model to employ -- SNP is the default option, while INDEL is also available for calling indels and BOTH is available for calling both together", optional = true)
    public GenotypeLikelihoodsCalculationModel GLmodel = GenotypeLikelihoodsCalculationModel.SNP;

    /**
     * The PCR error rate is independent of the sequencing error rate, which is necessary because we cannot necessarily
     * distinguish between PCR errors vs. sequencing errors.  The practical implication for this value is that it
     * effectively acts as a cap on the base qualities.
     */
    @Argument(fullName = "pcr_error_rate", shortName = "pcr_error", doc = "The PCR error rate to be used for computing fragment-based likelihoods", optional = true)
    public Double PCR_error = 1e-4;

    /**
     * Note that calculating the SLOD increases the runtime by an appreciable amount.
     */
    @Argument(fullName = "computeSLOD", shortName = "slod", doc = "If provided, we will calculate the SLOD (SB annotation)", optional = true)
    public boolean COMPUTE_SLOD = false;

    /**
     * The PairHMM implementation to use for -glm INDEL genotype likelihood calculations. The various implementations balance a tradeoff of accuracy and runtime.
     */
    @Argument(fullName = "pair_hmm_implementation", shortName = "pairHMM", doc = "The PairHMM implementation to use for -glm INDEL genotype likelihood calculations", optional = true)
    public PairHMM.Implementation pairHMM = PairHMM.Implementation.FASTEST_AVAILABLE;

    /**
     * The minimum confidence needed in a given base for it to be used in variant calling.  Note that the base quality of a base
     * is capped by the mapping quality so that bases on reads with low mapping quality may get filtered out depending on this value.
     * Note too that this argument is ignored in indel calling. In indel calling, low-quality ends of reads are clipped off (with fixed threshold of Q20).
     */
    @Argument(fullName = "min_base_quality_score", shortName = "mbq", doc = "Minimum base quality required to consider a base for calling", optional = true)
    public int MIN_BASE_QUALTY_SCORE = 17;

    /**
     * If the fraction of reads with deletions spanning a locus is greater than this value, the site will not be considered callable and will be skipped.
     * To disable the use of this parameter, set its value to >1.
     */
    @Argument(fullName = "max_deletion_fraction", shortName = "deletions", doc = "Maximum fraction of reads with deletions spanning this locus for it to be callable", optional = true)
    public Double MAX_DELETION_FRACTION = 0.05;

    // indel-related arguments
    /**
     * A candidate indel is genotyped (and potentially called) if there are this number of reads with a consensus indel at a site.
     * Decreasing this value will increase sensitivity but at the cost of larger calling time and a larger number of false positives.
     */
    @Argument(fullName = "min_indel_count_for_genotyping", shortName = "minIndelCnt", doc = "Minimum number of consensus indels required to trigger genotyping run", optional = true)
    public int MIN_INDEL_COUNT_FOR_GENOTYPING = 5;

    /**
     * Complementary argument to minIndelCnt.  Only samples with at least this fraction of indel-containing reads will contribute
     * to counting and overcoming the threshold minIndelCnt.  This parameter ensures that in deep data you don't end
     * up summing lots of super rare errors up to overcome the 5 read default threshold.  Should work equally well for
     * low-coverage and high-coverage samples, as low coverage samples with any indel containing reads should easily over
     * come this threshold.
     */
    @Argument(fullName = "min_indel_fraction_per_sample", shortName = "minIndelFrac", doc = "Minimum fraction of all reads at a locus that must contain an indel (of any allele) for that sample to contribute to the indel count for alleles", optional = true)
    public double MIN_INDEL_FRACTION_PER_SAMPLE = 0.25;

    @Advanced
    @Argument(fullName = "indelGapContinuationPenalty", shortName = "indelGCP", doc = "Indel gap continuation penalty, as Phred-scaled probability.  I.e., 30 => 10^-30/10", optional = true)
    public byte INDEL_GAP_CONTINUATION_PENALTY = 10;

    @Advanced
    @Argument(fullName = "indelGapOpenPenalty", shortName = "indelGOP", doc = "Indel gap open penalty, as Phred-scaled probability.  I.e., 30 => 10^-30/10", optional = true)
    public byte INDEL_GAP_OPEN_PENALTY = 45;

    @Hidden
    @Argument(fullName = "indelHaplotypeSize", shortName = "indelHSize", doc = "Indel haplotype size", optional = true)
    public int INDEL_HAPLOTYPE_SIZE = 80;

    @Hidden
    @Argument(fullName = "indelDebug", shortName = "indelDebug", doc = "Output indel debug info", optional = true)
    public boolean OUTPUT_DEBUG_INDEL_INFO = false;

    @Hidden
    @Argument(fullName = "ignoreSNPAlleles", shortName = "ignoreSNPAlleles", doc = "expt", optional = true)
    public boolean IGNORE_SNP_ALLELES = false;

    /*
        Generalized ploidy argument (debug only): squash all reads into a single pileup without considering sample info
     */
    @Hidden
    @Argument(fullName = "allReadsSP", shortName = "dl", doc = "expt", optional = true)
    public boolean TREAT_ALL_READS_AS_SINGLE_POOL = false;

    /*
       Generalized ploidy argument (debug only): When building site error models, ignore lane information and build only
       sample-level error model
     */
    @Hidden
    @Argument(fullName = "ignoreLaneInfo", shortName = "ignoreLane", doc = "Ignore lane when building error model, error model is then per-site", optional = true)
    public boolean IGNORE_LANE_INFO = false;

    /*
        Generalized ploidy argument: VCF file that contains truth calls for reference sample. If a reference sample is included through argument -refsample,
        then this argument is required.
     */
    @Hidden
    @Argument(fullName="reference_sample_calls", shortName = "referenceCalls", doc="VCF file with the truth callset for the reference sample", optional=true)
    public FeatureInput<VariantContext> referenceSampleRod;

    /*
        Reference sample name: if included, a site-specific error model will be built in order to improve calling quality. This requires ideally
        that a bar-coded reference sample be included with the polyploid/pooled data in a sequencing experimental design.
        If argument is absent, no per-site error model is included and calling is done with a generalization of traditional statistical calling.
     */
    @Hidden
    @Argument(shortName="refsample", fullName="reference_sample_name", doc="Reference sample name.", optional=true)
    public String referenceSampleName;

    /**
     * The following argument are for debug-only tweaks when running generalized ploidy with a reference sample
     */
    @Hidden
    @Argument(shortName="minqs", fullName="min_quality_score", doc="Min quality score to consider. Smaller numbers process faster. Default: Q1.", optional=true)
    public byte minQualityScore= 1;

    @Hidden
    @Argument(shortName="maxqs", fullName="max_quality_score", doc="Max quality score to consider. Smaller numbers process faster. Default: Q40.", optional=true)
    public byte maxQualityScore= 40;

    @Hidden
    @Argument(shortName="site_prior", fullName="site_quality_prior", doc="Phred-Scaled prior quality of the site. Default: Q20.", optional=true)
    public byte phredScaledPrior = 20;

    @Hidden
    @Argument(shortName = "min_call_power", fullName = "min_power_threshold_for_calling", doc="The minimum confidence in the error model to make a call. Number should be between 0 (no power requirement) and 1 (maximum power required).", optional = true)
    public double minPower = 0.95;
}
