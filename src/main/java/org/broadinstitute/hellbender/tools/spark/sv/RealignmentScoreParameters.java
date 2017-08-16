package org.broadinstitute.hellbender.tools.spark.sv;

import org.broadinstitute.barclay.argparser.Argument;

import java.io.Serializable;

/**
 * Collects arguments that control how we score read re-alignments.
 */
@SuppressWarnings("WeakerAccess")
public final class RealignmentScoreParameters implements Serializable {

    private static final long serialVersionUID = -1L;

    public static final double DEFAULT_MATCH_COST = 0;
    public static final double DEFAULT_MISMATCH_COST = 40;
    public static final double DEFAULT_GAP_OPEN_COST = 60;
    public static final double DEFAULT_GAP_EXTEND_COST = 10;
    public static final int DEFAULT_MIN_MAPQ = 16;
    public static final double DEFAULT_DOMINANT_ALIGNED_BASE_COUNT_RATIO = 1.1;
    public static final double DEFAULT_UNMAPPED_FRAGMENT_PENALTY = 60;
    public static final double DEFAULT_IMPROPER_PAIR_PENALTY = 20;
    public static final double DEFAULT_MAXIMUM_LIKELIHOOD_DIFFERENCE_PER_TEMPLATE = 120;
    public static final double DEFAULT_INTER_HAPLOTYPE_PENALTY = 60;

    public static final String GAP_OPEN_COST_PARAM_FULL_NAME = "realignment-gap-open-penalty";
    public static final String GAP_OPEN_COST_PARAM_SHORT_NAME = GAP_OPEN_COST_PARAM_FULL_NAME;
    public static final String GAP_EXTEND_COST_PARAM_FULL_NAME = "realignment-gap-extend-penalty";
    public static final String GAP_EXTEND_COST_PARAM_SHORT_NAME = GAP_EXTEND_COST_PARAM_FULL_NAME;
    public static final String MATCH_COST_PARAM_FULL_NAME = "realignment-match-penalty";
    public static final String MATCH_COST_PARAM_SHORT_NAME = MATCH_COST_PARAM_FULL_NAME;
    public static final String MISMATCH_COST_PARAM_FULL_NAME = "realignment-mismatch-penalty";
    public static final String MISMATCH_COST_PARAM_SHORT_NAME = MISMATCH_COST_PARAM_FULL_NAME;
    public static final String MINIMUM_MAP_QUALITY_PARAM_FULL_NAME = "realignment-minimum-mq-for-base-call-counting";
    public static final String MINIMUM_MAP_QUALITY_PARAM_SHORT_NAME = MINIMUM_MAP_QUALITY_PARAM_FULL_NAME ;
    public static final String DOMINANT_ALIGNED_BASE_COUNT_RATIO_PARAM_FULL_NAME = "dominant-strand-aligned-base-count-ratio";
    public static final String DOMINANT_ALIGNED_BASE_COUNT_RATIO_PARAM_SHORT_NAME = DOMINANT_ALIGNED_BASE_COUNT_RATIO_PARAM_FULL_NAME;
    public static final String MAXIMUM_LIKELIHOOD_DIFFERENCE_PER_TEMPLATE_FULL_NAME = "maximum-likelihood-difference-per-template";
    public static final String MAXIMUM_LIKELIHOOD_DIFFERENCE_PER_TEMPLATE_SHORT_NAME = MAXIMUM_LIKELIHOOD_DIFFERENCE_PER_TEMPLATE_FULL_NAME;
    public static final String IMPROPER_PAIR_PENALTY_FULL_NAME = "improper-pair-penalty";
    public static final String IMPROPER_PAIR_PENALTY_SHORT_NAME = IMPROPER_PAIR_PENALTY_FULL_NAME;
    public static final String UNMAPPED_FRAGMENT_PENALTY_FULL_NAME = "unmapped-fragment-penalty";
    public static final String UNMAPPED_FRAGMENT_PENALITY_SHORT_NAME = UNMAPPED_FRAGMENT_PENALTY_FULL_NAME;
    public static final String INTER_HAPLOTYPE_PENALTY_FULL_NAME = "inter-haplotype-penalty";
    public static final String INTER_HAPLOTYPE_PENALTY_SHORT_NAME = INTER_HAPLOTYPE_PENALTY_FULL_NAME;

   @Argument(fullName = GAP_OPEN_COST_PARAM_FULL_NAME, shortName = GAP_OPEN_COST_PARAM_SHORT_NAME,
            doc = "Phred-scaled cost for a gap (indel) opening in a read realigned against a contig or haplotype",
            minValue = 0,
            optional =  true)
   public double gapOpenPenalty = DEFAULT_GAP_OPEN_COST;

   @Argument(fullName = GAP_EXTEND_COST_PARAM_FULL_NAME, shortName = GAP_EXTEND_COST_PARAM_SHORT_NAME,
            doc = "Phred-scaled cost for a gap (indel) extension in a read realigned against a contig or haplotype",
            minValue = 0,
            optional = true)
   public double gapExtendPenalty = DEFAULT_GAP_EXTEND_COST;

   @Argument(fullName = MATCH_COST_PARAM_FULL_NAME, shortName = MATCH_COST_PARAM_SHORT_NAME,
            doc = "Phred-scaled cost for base match in a read realigned against a contig or haplotype",
            minValue = 0,
            optional =  true)
   public double matchPenalty = DEFAULT_MATCH_COST;

   @Argument(fullName = MISMATCH_COST_PARAM_FULL_NAME, shortName = MISMATCH_COST_PARAM_SHORT_NAME,
            doc = "Phred-scaled cost for base mismatch in a read realigned against a contig or haplotype",
            minValue = 0,
            optional =  true)
   public double mismatchPenalty = DEFAULT_MISMATCH_COST;

   @Argument(fullName = MINIMUM_MAP_QUALITY_PARAM_FULL_NAME,
           shortName = MINIMUM_MAP_QUALITY_PARAM_SHORT_NAME,
            doc = "Minimum mapping quality of a local read realignment interval to be considered when comparing the " +
                    "number of base calls between strands to figure out which strand the sequence aligns to",
            minValue = 0,
            optional = true)
   public int minimumMappingQuality = DEFAULT_MIN_MAPQ;

   @Argument(fullName = DOMINANT_ALIGNED_BASE_COUNT_RATIO_PARAM_FULL_NAME,
            shortName = DOMINANT_ALIGNED_BASE_COUNT_RATIO_PARAM_SHORT_NAME,
            doc = "Minimum ratio of the number of bases aligned to the actual strand vs the other strand. A sequence is supposed to " +
                    "align to the strand whose number of base calls aligned is larger than the other strand's by at least this factor, " +
                    "regardless of actual editing distance",
            minValue = 1.0,
            optional = true)
   public double dominantStrandAlignedBaseCountRatio = DEFAULT_DOMINANT_ALIGNED_BASE_COUNT_RATIO;

   @Argument(fullName = UNMAPPED_FRAGMENT_PENALTY_FULL_NAME,
             shortName = UNMAPPED_FRAGMENT_PENALITY_SHORT_NAME,
             doc = "Phred scaled score penalty for an unmapped read/fragment",
             minValue = 0.0,
             optional = true)
   public double unmappedFragmentPenalty = DEFAULT_UNMAPPED_FRAGMENT_PENALTY;

   @Argument(fullName = IMPROPER_PAIR_PENALTY_FULL_NAME,
             shortName = IMPROPER_PAIR_PENALTY_SHORT_NAME,
             doc = "Phred scaled score penalty for a read pair to have an improper orientation",
             minValue = 0.0,
             optional = true)
   public double improperPairPenalty = DEFAULT_IMPROPER_PAIR_PENALTY;

   @Argument(fullName = MAXIMUM_LIKELIHOOD_DIFFERENCE_PER_TEMPLATE_FULL_NAME,
             shortName = MAXIMUM_LIKELIHOOD_DIFFERENCE_PER_TEMPLATE_SHORT_NAME,
             minValue = 0.0,
             optional = true)
   public double maximumLikelihoodDifferencePerTemplate = DEFAULT_MAXIMUM_LIKELIHOOD_DIFFERENCE_PER_TEMPLATE;

   @Argument(fullName = INTER_HAPLOTYPE_PENALTY_FULL_NAME,
             shortName = INTER_HAPLOTYPE_PENALTY_SHORT_NAME,
             doc = "score penalty imposed for a mapping of a haplotype vs another",
             minValue = 0.0,
             optional = true)
   public double interHaplotypePenalty = DEFAULT_INTER_HAPLOTYPE_PENALTY;

    public RealignmentScoreParameters() {}
}
