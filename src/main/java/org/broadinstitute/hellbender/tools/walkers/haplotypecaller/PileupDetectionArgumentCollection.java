package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;

/**
 * Set of arguments for configuring the pileup detection code
 */
public final class PileupDetectionArgumentCollection {

    public static final String PILEUP_DETECTION_LONG_NAME = "pileup-detection";
    public static final String PILEUP_DETECTION_ENABLE_INDELS = "pileup-detection-enable-indel-detection";
    public static final String PILEUP_DETECTION_SNP_THRESHOLD = "pileup-detection-enable-snp-alt-threshold";
    public static final String PILEUP_DETECTION_ABSOLUTE_ALT_DEPTH = "pileup-detection-absolute-alt-depth";
    public static final String PILEUP_DETECTION_INDEL_THRESHOLD = "pileup-detection-enable-indel-alt-threshold";
    public static final String PILEUP_DETECTION_FILTER_COVERAGE_LONG_NAME = "pileup-detection-filter-coverage-threshold";

    //TODO we currently don't see the same threshold from active region determination...
    public static final String PILEUP_DETECTION_ACETIVE_REGION_LOD_THRESHOLD_LONG_NAME = "pileup-detection-active-region-lod-threshold";


    // Arguments related to DRAGEN heuristics related to "read badness" intended to filter out false positives from the pileup deteciotn code
    public static final String PILEUP_DETECTION_BAD_READ_RATIO_LONG_NAME = "pileup-detection-bad-read-tolerance";
    public static final String PILEUP_DETECTION_PROPER_PAIR_READ_BADNESS_LONG_NAME = "pileup-detection-proper-pair-read-badness";
    public static final String PILEUP_DETECTION_EDIT_DISTANCE_BADNESS_LONG_NAME = "pileup-detection-edit-distance-read-badness-threshold";
    public static final String PILEUP_DETECTION_CHIMERIC_READ_BADNESS_LONG_NAME = "pileup-detection-chimeric-read-badness";
    //TODO these need to be implemented with some input from Illumina
    public static final String PILEUP_DETECTION_TLEN_MEAN_LONG_NAME = "pileup-detection-template-mean-badness-threshold";
    public static final String PILEUP_DETECTION_TLEN_STD_LONG_NAME = "pileup-detection-template-std-badness-threshold";

    // Argumetns related to filtering the list at the assembly graph level
    public static final String PILEUP_DETECTION_INDEL_SNP_BLOCKING_RANGE = "pileup-detection-snp-adjacent-to-assembled-indel-range";

    /**
     * Enables pileup-based haplotype creation and variant detection
     *
     * NOTE: --pileup-detection is a beta feature. Use this mode at your own risk.
     */
    @Argument(fullName= PILEUP_DETECTION_LONG_NAME, doc = "If enabled, the variant caller will create pileup-based haplotypes in addition to the assembly-based haplotype generation.", optional = true)
    public boolean usePileupDetection = false;

    /**
     * Enables detection of indels from the pileups in. (EXPERIMENTAL FEATURE)
     */
    @Advanced
    @Argument(fullName= PILEUP_DETECTION_ENABLE_INDELS, doc = "If enabled, pileup detection code will attempt to detect indels missing from assembly. (Requires `--pileup-detection` argument)", optional = true)
    public boolean detectIndels = false;


    /**
     * Percentage of reads required to supprot the alt for a variant to be considered
     */
    @Advanced
    @Argument(fullName= PILEUP_DETECTION_SNP_THRESHOLD, doc = "Percentage of alt supporting reads in order to consider alt SNP", optional = true)
    public double snpThreshold = 0.1;
    @Advanced
    @Argument(fullName= PILEUP_DETECTION_INDEL_THRESHOLD, doc = "Percentage of alt supporting reads in order to consider alt indel", optional = true)
    public double indelThreshold = 0.5;


    @Advanced
    @Argument(fullName= PILEUP_DETECTION_ABSOLUTE_ALT_DEPTH, doc = "Absolute number of alt reads necessary to be included in pileup events", optional = true)
    public double pileupAbsoluteDepth = 0;
    @Advanced
    @Argument(fullName= PILEUP_DETECTION_INDEL_SNP_BLOCKING_RANGE, doc = "Filters out pileup snps within this many bases of an assembled indel", optional = true)
    public int snpAdajacentToAssemblyIndel = 5;




    /**
     * Arguments related to the "bad read filtering" where alleles that are supported primarily by reads that fail at least one of a number of heuristics will be filtered out
     */
    @Advanced
    @Argument(fullName= PILEUP_DETECTION_BAD_READ_RATIO_LONG_NAME, doc = "Thresold of Alt reads rejected by bad reads heuristics to allow the variant", optional = true)
    public double badReadThreshold = 0.0;
    @Advanced
    @Argument(fullName= PILEUP_DETECTION_PROPER_PAIR_READ_BADNESS_LONG_NAME, doc = "Reject Alt reads not in proper-pairs", optional = true)
    public boolean badReadProperPair = true;
    @Advanced
    @Argument(fullName= PILEUP_DETECTION_EDIT_DISTANCE_BADNESS_LONG_NAME, doc = "Reject alt reads with greater than this percent edit distance from the reference", optional = true)
    public double badReadEditDistance = 0.08;
    @Advanced
    @Argument(fullName= PILEUP_DETECTION_CHIMERIC_READ_BADNESS_LONG_NAME, doc = "Reject reads that are chimeric or supplementary", optional = true)
    public boolean badReadSecondary = true;
    @Advanced
    @Argument(fullName= PILEUP_DETECTION_TLEN_MEAN_LONG_NAME, doc = "Mean template length to consider for read badness. Requires '--"+PILEUP_DETECTION_TLEN_STD_LONG_NAME+"' to also be set.", optional = true)
    public double templateLengthMean = 0.0;
    @Advanced
    @Argument(fullName= PILEUP_DETECTION_TLEN_STD_LONG_NAME, doc = "Standard deviation template to consider for read badness. Requires '--"+PILEUP_DETECTION_TLEN_MEAN_LONG_NAME+"' to also be set.", optional = true)
    public double templateLenghtStd = 0.0;
}
