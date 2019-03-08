package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.tools.walkers.mutect.M2ArgumentCollection;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

public class M2FiltersArgumentCollection {
    private static final long serialVersionUID = 9345L;

    /**
     * meta-filtering parameters -- how we set the threshold on artifact probabilities
     */
    public static final String THRESHOLD_STRATEGY_LONG_NAME = "threshold-strategy";
    public static final String F_SCORE_BETA_LONG_NAME = "f-score-beta";
    public static final String FALSE_DISCOVERY_RATE_LONG_NAME = "false-discovery-rate";
    public static final String INITIAL_THRESHOLD_LONG_NAME = "initial-threshold";

    private static final ThresholdCalculator.Strategy DEFAULT_THRESHOLD_STRATEGY = ThresholdCalculator.Strategy.OPTIMAL_F_SCORE;
    private static final double DEFAULT_F_SCORE_BETA = 1.0;
    private static final double DEFAULT_MAX_FALSE_DISCOVERY_RATE = 0.05;
    private static final double DEFAULT_INITIAL_POSTERIOR_THRESHOLD = 0.1;

    @Argument(fullName = THRESHOLD_STRATEGY_LONG_NAME, optional = true, doc = "The method for optimizing the posterior probability threshold")
    public ThresholdCalculator.Strategy thresholdStrategy = DEFAULT_THRESHOLD_STRATEGY;

    @Argument(fullName = F_SCORE_BETA_LONG_NAME, optional = true, doc = "F score beta, the relative weight of recall to precision, used if OPTIMAL_F_SCORE strategy is chosen")
    public double fScoreBeta = DEFAULT_F_SCORE_BETA;

    @Argument(fullName = FALSE_DISCOVERY_RATE_LONG_NAME, optional = true, doc = "Maximum false discovery rate allowed if FALSE_DISCOVERY_RATE threshold strategy is chosen")
    public double maxFalsePositiveRate = DEFAULT_MAX_FALSE_DISCOVERY_RATE;

    @Argument(fullName = INITIAL_THRESHOLD_LONG_NAME, optional = true, doc = "Initial artifact probability threshold used in first iteration")
    public double initialPosteriorThreshold = DEFAULT_INITIAL_POSTERIOR_THRESHOLD;

    /**
     * Mitochondria mode includes the filters {@link LogOddsOverDepthFilter} and {@link ChimericOriginalAlignmentFilter}
     * and excludes the filters {@link ClusteredEventsFilter}, {@link MultiallelicFilter}, {@link PolymeraseSlippageFilter},
     * {@link FilteredHaplotypeFilter}, and {@link GermlineFilter}
     */
    @Argument(fullName = M2ArgumentCollection.MITOCHONDRIA_MODE_LONG_NAME, optional = true, doc = "Set filters to mitochondrial defaults")
    public boolean mitochondria = false;


    /**
     * Hard filter thresholds
     */
    public static final String MAX_EVENTS_IN_REGION_LONG_NAME = "max-events-in-region";
    public static final String MAX_ALT_ALLELE_COUNT_LONG_NAME = "max-alt-allele-count";
    public static final String UNIQUE_ALT_READ_COUNT_LONG_NAME = "unique-alt-read-count";
    public static final String MIN_MEDIAN_MAPPING_QUALITY_LONG_NAME = "min-median-mapping-quality";
    public static final String MIN_MEDIAN_BASE_QUALITY_LONG_NAME = "min-median-base-quality";
    public static final String MAX_MEDIAN_FRAGMENT_LENGTH_DIFFERENCE_LONG_NAME = "max-median-fragment-length-difference";
    public static final String MIN_MEDIAN_READ_POSITION_LONG_NAME = "min-median-read-position";
    public static final String MAX_N_RATIO_LONG_NAME = "max-n-ratio";
    public static final String MIN_LOG_10_ODDS_DIVIDED_BY_DEPTH = "lod-divided-by-depth";
    public static final String MIN_READS_ON_EACH_STRAND_LONG_NAME = "min-reads-per-strand";
    public static final String MAX_NUMT_FRACTION_LONG_NAME = "max-numt-fraction";

    private static final int DEFAULT_MAX_EVENTS_IN_REGION = 2;
    private static final int DEFAULT_MAX_ALT_ALLELES = 1;
    private static final int DEFAULT_MIN_UNIQUE_ALT_READS = 0;
    private static final int DEFAULT_MIN_MEDIAN_MAPPING_QUALITY = 30;
    private static final int DEFAULT_MIN_MEDIAN_BASE_QUALITY = 20;
    private static final int DEFAULT_MAX_MEDIAN_FRAGMENT_LENGTH_DIFFERENCE = 10000;
    private static final int DEFAULT_MIN_MEDIAN_READ_POSITION = 1;
    private static final double DEFAULT_MAX_N_RATIO = Double.POSITIVE_INFINITY;
    private static final double DEFAULT_MIN_LOG_10_ODDS_DIVIDED_BY_DEPTH = 0.0035;
    private static final int DEFAULT_MIN_READS_ON_EACH_STRAND = 0;
    private static final double DEFAULT_MAX_NUMT_FRACTION = 0.85;

    @Argument(fullName = MAX_EVENTS_IN_REGION_LONG_NAME, optional = true, doc = "Maximum events in a single assembly region.  Filter all variants if exceeded.")
    public int maxEventsInRegion = DEFAULT_MAX_EVENTS_IN_REGION;

    @Argument(fullName = MAX_ALT_ALLELE_COUNT_LONG_NAME, optional = true, doc = "Maximum alt alleles per site.")
    public int numAltAllelesThreshold = DEFAULT_MAX_ALT_ALLELES;

    @Argument(fullName = UNIQUE_ALT_READ_COUNT_LONG_NAME, shortName = "unique", optional = true, doc = "Minimum unique (i.e. deduplicated) reads supporting the alternate allele")
    public int uniqueAltReadCount = DEFAULT_MIN_UNIQUE_ALT_READS;

    @Argument(fullName = MIN_MEDIAN_MAPPING_QUALITY_LONG_NAME, optional = true, doc="Minimum median mapping quality of alt reads")
    public int minMedianMappingQuality = DEFAULT_MIN_MEDIAN_MAPPING_QUALITY;

    @Argument(fullName = MIN_MEDIAN_BASE_QUALITY_LONG_NAME, optional = true, doc="Minimum median base quality of alt reads")
    public int minMedianBaseQuality = DEFAULT_MIN_MEDIAN_BASE_QUALITY;

    @Argument(fullName = MAX_MEDIAN_FRAGMENT_LENGTH_DIFFERENCE_LONG_NAME, optional = true, doc="Maximum difference between median alt and ref fragment lengths")
    public int maxMedianFragmentLengthDifference = DEFAULT_MAX_MEDIAN_FRAGMENT_LENGTH_DIFFERENCE;

    @Argument(fullName = MIN_MEDIAN_READ_POSITION_LONG_NAME, optional = true, doc = "Minimum median distance of variants from the end of reads")
    public int minMedianReadPosition = DEFAULT_MIN_MEDIAN_READ_POSITION;

    @Argument(fullName = MAX_N_RATIO_LONG_NAME, optional = true, doc = "Maximum fraction of non-ref bases in the pileup that are N (unknown)")
    public double nRatio = DEFAULT_MAX_N_RATIO;

    @Argument(fullName = MIN_LOG_10_ODDS_DIVIDED_BY_DEPTH, doc="Minimum log10 odds divided by depth", optional = true)
    public double minLog10OddsDividedByDepth = DEFAULT_MIN_LOG_10_ODDS_DIVIDED_BY_DEPTH;

    @Argument(fullName = MIN_READS_ON_EACH_STRAND_LONG_NAME, optional = true, doc = "Minimum alt reads required on both forward and reverse strands")
    public int minReadsOnEachStrand = DEFAULT_MIN_READS_ON_EACH_STRAND;

    @Argument(fullName = MAX_NUMT_FRACTION_LONG_NAME, doc="Maximum fraction of alt reads that originally aligned outside the mitochondria.  These are due to NuMTs.", optional = true)
    public double maxNuMTFraction = DEFAULT_MAX_NUMT_FRACTION;


    /**
     * Input files and values to use if inputs are missing
     */
    public static final String CONTAMINATION_TABLE_LONG_NAME = "contamination-table";
    public static final String CONTAMINATION_ESTIMATE_LONG_NAME = "contamination-estimate";
    public static final String TUMOR_SEGMENTATION_LONG_NAME = "tumor-segmentation";

    private static final double DEFAULT_CONTAMINATION = 0.0;

    @Argument(fullName = CONTAMINATION_TABLE_LONG_NAME, optional = true, doc = "Tables containing contamination information.")
    public List<File> contaminationTables = new ArrayList<>();

    @Argument(fullName = CONTAMINATION_ESTIMATE_LONG_NAME, optional = true, doc = "Estimate of contamination.")
    public double contaminationEstimate = DEFAULT_CONTAMINATION;

    @Argument(fullName = TUMOR_SEGMENTATION_LONG_NAME, doc="Tables containing tumor segments' minor allele fractions for germline hets emitted by CalculateContamination", optional = true)
    public List<File> tumorSegmentationTables = new ArrayList<>();


    /**
     * Parameters
     */
    public static final String LOG_10_PRIOR_OF_SNV_LONG_NAME = "log10-somatic-snv-prior";
    public static final String LOG_10_PRIOR_OF_INDEL_LONG_NAME = "log10-somatic-indel-prior";
    public static final String INITIAL_LOG_10_PRIOR_OF_ARTIFACT_VERSUS_VARIANT_LONG_NAME = "initial-artifact-prior";
    public static final String NORMAL_P_VALUE_THRESHOLD_LONG_NAME = "normal-p-value-threshold";
    public static final String MIN_POLYMERASE_SLIPPAGE_LENGTH = "min-slippage-length";
    public static final String PCR_SLIPPAGE_RATE_LONG_NAME = "pcr-slippage-rate";
    public static final String MAX_DISTANCE_TO_FILTERED_CALL_ON_SAME_HAPLOTYPE_LONG_NAME = "distance-on-haplotype";
    public static final String LONG_INDEL_LENGTH_LONG_NAME = "long-indel-length";

    private static final double DEFAULT_LOG_10_PRIOR_OF_SNV = -6.0;
    private static final double DEFAULT_LOG_10_PRIOR_OF_INDEL = -7.0;
    private static final double DEFAULT_INITIAL_LOG_10_PRIOR_OF_VARIANT_VERSUS_ARTIFACT = -1;
    private static final double DEFAULT_NORMAL_P_VALUE_THRESHOLD = 0.001;
    private static final int DEFAULT_MIN_SLIPPAGE_LENGTH = 8;
    private static final double DEFAULT_SLIPPAGE_RATE = 0.1;
    private static final int DEFAULT_MAX_INTRA_HAPLOTYPE_DISTANCE = 100;
    private static final int DEFAULT_LONG_INDEL_SIZE = 5;

    @Argument(fullName= LOG_10_PRIOR_OF_SNV_LONG_NAME, doc="Log10 prior probability that a site has a somatic SNV", optional = true)
    public double log10SNVPrior = DEFAULT_LOG_10_PRIOR_OF_SNV;

    @Argument(fullName= LOG_10_PRIOR_OF_INDEL_LONG_NAME, doc="Log10 prior probability that a site has a somatic indel", optional = true)
    public double log10IndelPrior = DEFAULT_LOG_10_PRIOR_OF_INDEL;

    @Argument(fullName= INITIAL_LOG_10_PRIOR_OF_ARTIFACT_VERSUS_VARIANT_LONG_NAME, doc="Initial guess for log10 prior probability that a called site is not a technical artifact", optional = true)
    public double initialLog10PriorOfVariantVersusArtifact = DEFAULT_INITIAL_LOG_10_PRIOR_OF_VARIANT_VERSUS_ARTIFACT;

    @Argument(fullName = NORMAL_P_VALUE_THRESHOLD_LONG_NAME, optional = true, doc = "P value threshold for normal artifact filter")
    public static final double normalPileupPValueThreshold = DEFAULT_NORMAL_P_VALUE_THRESHOLD;

    @Argument(fullName = MIN_POLYMERASE_SLIPPAGE_LENGTH, optional = true, doc = "Minimum number of reference bases in an STR to suspect polymerase slippage")
    public int minSlippageLength = DEFAULT_MIN_SLIPPAGE_LENGTH;

    @Argument(fullName = PCR_SLIPPAGE_RATE_LONG_NAME, optional = true, doc = "The frequency of polymerase slippage in contexts where it is suspected")
    public double slippageRate = DEFAULT_SLIPPAGE_RATE;


    @Argument(fullName = MAX_DISTANCE_TO_FILTERED_CALL_ON_SAME_HAPLOTYPE_LONG_NAME, optional = true, doc = "On second filtering pass, variants with same PGT and PID tags as a filtered variant within this distance are filtered.")
    public int maxDistanceToFilteredCallOnSameHaplotype = DEFAULT_MAX_INTRA_HAPLOTYPE_DISTANCE;

    @Argument(fullName = LONG_INDEL_LENGTH_LONG_NAME, optional = true, doc = "Indels of this length or greater are treated specially by the mapping quality filter.")
    public int longIndelLength = DEFAULT_LONG_INDEL_SIZE;
}
