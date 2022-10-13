package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.Hidden;
import org.checkerframework.checker.units.qual.A;

/**
 * Set of arguments for configuring the pileup detection code
 */
public final class PileupDetectionArgumentCollection {

    public static final String PILEUP_DETECTION_LONG_NAME = "pileup-detection";
    public static final String PILEUP_DETECTION_ENABLE_INDELS = "pileup-detection-enable-indel-pileup-calling";
    public static final String PILEUP_DETECTION_SNP_THRESHOLD = "pileup-detection-snp-alt-threshold";
    public static final String PILEUP_DETECTION_ABSOLUTE_ALT_DEPTH = "pileup-detection-absolute-alt-depth";
    public static final String PILEUP_DETECTION_INDEL_THRESHOLD = "pileup-detection-indel-alt-threshold";
    public static final String PILEUP_DETECTION_FILTER_COVERAGE_LONG_NAME = "pileup-detection-filter-coverage-threshold";
    public static final String PILEUP_DETECTION_SNP_BASEQUALITY_THRESHOLD = "pileup-detection-snp-basequlaity-filter";

    public static final String PILEUP_DETECTION_ACTIVE_REGION_LOD_THRESHOLD_LONG_NAME = "pileup-detection-active-region-phred-threshold";

    // Arguments related to DRAGEN heuristics related to "read badness" intended to filter out false positives from the pileup detection code
    public static final String PILEUP_DETECTION_BAD_READ_RATIO_LONG_NAME = "pileup-detection-bad-read-tolerance";
    public static final String PILEUP_DETECTION_PROPER_PAIR_READ_BADNESS_LONG_NAME = "pileup-detection-proper-pair-read-badness";
    public static final String PILEUP_DETECTION_EDIT_DISTANCE_BADNESS_LONG_NAME = "pileup-detection-edit-distance-read-badness-threshold";
    public static final String PILEUP_DETECTION_CHIMERIC_READ_BADNESS_LONG_NAME = "pileup-detection-chimeric-read-badness";
    //TODO these need to be implemented with some input from Illumina
    public static final String PILEUP_DETECTION_TLEN_MEAN_LONG_NAME = "pileup-detection-template-mean-badness-threshold";
    public static final String PILEUP_DETECTION_TLEN_STD_LONG_NAME = "pileup-detection-template-std-badness-threshold";

    // Argumetns related to filtering the list at the assembly graph level
    public static final String PILEUP_DETECTION_INDEL_SNP_BLOCKING_RANGE = "pileup-detection-snp-adjacent-to-assembled-indel-range";

    // Arguments related to post-assembly filtering heuristics.
    public static final String PILEUP_DETECTION_FILTER_ASSEMBLY_HAPS_THRESHOLD = "pileup-detection-filter-assembly-alt-bad-read-tolerance";
    public static final String PILEUP_DETECTION_EDIT_DISTANCE_BADNESS_FOR_ASSEMBLY_LONG_NAME = "pileup-detection-edit-distance-read-badness-for-assembly-filtering-threshold";


    public static final String GENERATE_PARTIALLY_DETERMINED_HAPLOTYPES_LONG_NAME = "use-pdhmm";
    public static final String DETERMINE_PD_HAPS = "make-determined-haps-from-pd-code";
    public static final String DEBUG_PILEUPCALLING_ARG = "print-pileupcalling-status";
    public static final String FALLBACK_TO_ALT_HAP_CONSTRUCITON_IF_ABORTED = "fallback-gga-if-pdhmm-fails";
    public static final String PDHMM_READ_OVERLAP_OPTIMIZATION = "use-pdhmm-overlap-optimization";

    /**
     * Enables pileup-based haplotype creation and variant detection
     *
     * NOTE: --pileup-detection is a beta feature. Use this mode at your own risk.
     */
    @Advanced
    @Argument(fullName= PILEUP_DETECTION_LONG_NAME, doc = "If enabled, the variant caller will create pileup-based haplotypes in addition to the assembly-based haplotype generation.", optional = true)
    public boolean usePileupDetection = false;
    /**
     * TODO
     */
    @Advanced
    @Argument(fullName= GENERATE_PARTIALLY_DETERMINED_HAPLOTYPES_LONG_NAME, doc = "", optional = true)
    public boolean generatePDHaplotypes = false;
    /**
     * A set of debug/testing arguments related to PDHMM and various off-label configurations for running pdhmm.
     */
    @Hidden
    @Argument(fullName= DETERMINE_PD_HAPS, doc = "", optional = true)
    public boolean determinePDHaps = false;
    @Hidden
    @Argument(fullName= DEBUG_PILEUPCALLING_ARG, doc = "", optional = true)
    public boolean debugPileupStdout = false;
    @Hidden
    @Argument(fullName= FALLBACK_TO_ALT_HAP_CONSTRUCITON_IF_ABORTED, doc = "", optional = true)
    public boolean useGGAFallback = true;
    @Hidden
    // NOTE: this optimization ASSUMES that we are not realigning reads as part of PDHMM (which is how DRAGEN works)
    @Argument(fullName= PDHMM_READ_OVERLAP_OPTIMIZATION, doc = "", optional = true)
    public boolean pdhmmOptimization = false;


    /**
     * Enables detection of indels from the pileups in. (EXPERIMENTAL FEATURE)
     */
    @Hidden
    @Argument(fullName= PILEUP_DETECTION_ENABLE_INDELS, doc = "Pileup Detection: If enabled, pileup detection code will attempt to detect indels missing from assembly. (Requires '--"+PILEUP_DETECTION_LONG_NAME+"' argument)", optional = true)
    public boolean detectIndels = false;
    @Advanced
    @Hidden
    @Argument(fullName= PILEUP_DETECTION_ACTIVE_REGION_LOD_THRESHOLD_LONG_NAME, doc = "Pileup Detection: This argument sets the minimum fast genotyper phred score necessary in active region determination to find a pileup variant. (Requires '--"+PILEUP_DETECTION_LONG_NAME+"'` argument)", optional = true)
    public double activeRegionPhredThreshold = 2.0;
    @Advanced
    @Hidden
    @Argument(fullName= "num-artificial-haplotypes-to-add-per-allele", doc = "Pileup Detection: This argument limits the maximum number of novel haplotypes to be added to the assembly haplotypes per pileup allele added. (Requires '--"+PILEUP_DETECTION_LONG_NAME+"'` argument)", optional = true)
    public int numHaplotypesToIterate = 5;
    @Advanced
    @Hidden
    @Argument(fullName= "artifical-haplotype-filtering-kmer-size", doc = "Pileup Detection: Controls what size to kmerize reads to in order to select best supported artificial haplotypes. (Requires '--"+PILEUP_DETECTION_LONG_NAME+"' argument)", optional = true)
    public int filteringKmerSize = 10;

    /**
     * Percentage of reads required to support the alt for a variant to be considered
     */
    @Hidden
    @Argument(fullName= PILEUP_DETECTION_SNP_THRESHOLD, doc = "Pileup Detection: Percentage of alt supporting reads in order to consider alt SNP. (Requires '--"+PILEUP_DETECTION_LONG_NAME+"' argument)", optional = true)
    public double snpThreshold = 0.1;
    @Hidden
    @Argument(fullName= PILEUP_DETECTION_INDEL_THRESHOLD, doc = "Pileup Detection: Percentage of alt supporting reads in order to consider alt indel. (Requires '--"+PILEUP_DETECTION_LONG_NAME+"' argument)", optional = true)
    public double indelThreshold = 0.1;

    @Hidden
    @Argument(fullName= PILEUP_DETECTION_ABSOLUTE_ALT_DEPTH, doc = "Pileup Detection: Absolute number of alt reads necessary to be included in pileup events. (Requires '--"+PILEUP_DETECTION_LONG_NAME+"' argument)", optional = true)
    public double pileupAbsoluteDepth = 0;
    @Hidden
    @Argument(fullName= PILEUP_DETECTION_INDEL_SNP_BLOCKING_RANGE, doc = "Pileup Detection: Filters out pileup snps within this many bases of an assembled indel. (Requires '--"+PILEUP_DETECTION_LONG_NAME+"' argument)", optional = true)
    public int snpAdajacentToAssemblyIndel = 5;
    @Hidden
    @Argument(fullName= PILEUP_DETECTION_SNP_BASEQUALITY_THRESHOLD, doc = "Pileup Detection: Filters out reads from pileup SNPs with base quality lower than this threshold. (Requires '--"+PILEUP_DETECTION_LONG_NAME+"' argument)", optional = true)
    public int qualityForSnpsInPileupDeteciton = 12; //WHY is this one different than the regular active region determination limit? Ask DRAGEN engineers.

    /**
     * Arguments related to the "bad read filtering" where alleles that are supported primarily by reads that fail at least one of a number of heuristics will be filtered out
     */
    @Hidden
    @Argument(fullName= PILEUP_DETECTION_BAD_READ_RATIO_LONG_NAME, doc = "Pileup Detection: Threshold of Alt reads rejected by bad reads heuristics to allow the variant. (Requires '--"+PILEUP_DETECTION_LONG_NAME+"' argument)", optional = true)
    public double badReadThreshold = 0.0;
    @Hidden
    @Argument(fullName= PILEUP_DETECTION_PROPER_PAIR_READ_BADNESS_LONG_NAME, doc = "Pileup Detection: Reject alt reads not in proper-pairs. (Requires '--"+PILEUP_DETECTION_LONG_NAME+"' argument)", optional = true)
    public boolean badReadProperPair = true;
    @Hidden
    @Argument(fullName= PILEUP_DETECTION_EDIT_DISTANCE_BADNESS_LONG_NAME, doc = "Pileup Detection: Reject alt reads with greater than this fraction of mismatching bases from the reference (proxied using the NM tag). (Requires '--"+PILEUP_DETECTION_LONG_NAME+"' argument)", optional = true)
    public double badReadEditDistance = 0.08;
    @Hidden
    @Argument(fullName= PILEUP_DETECTION_CHIMERIC_READ_BADNESS_LONG_NAME, doc = "Pileup Detection: Reject reads that are chimeric or supplementary. (Requires '--"+PILEUP_DETECTION_LONG_NAME+"' argument)", optional = true)
    public boolean badReadSecondaryOrSupplementary = true;
    @Hidden
    @Argument(fullName= PILEUP_DETECTION_TLEN_MEAN_LONG_NAME, doc = "Pileup Detection: Mean template length (T LEN) to consider for read badness. Requires '--"+PILEUP_DETECTION_TLEN_STD_LONG_NAME+"' to also be set.", optional = true)
    public double templateLengthMean = 0.0;
    @Hidden
    @Argument(fullName= PILEUP_DETECTION_TLEN_STD_LONG_NAME, doc = "Pileup Detection: Standard deviation template length (T LEN) to consider for read badness. Requires '--"+PILEUP_DETECTION_TLEN_MEAN_LONG_NAME+"' to also be set.", optional = true)
    public double templateLengthStd = 0.0;

    /**
     * Enables pileup-based haplotype creation and variant detection
     *
     * NOTE: --pileup-detection is a beta feature. Use this mode at your own risk.
     */
    @Advanced
    @Hidden
    @Argument(fullName= PILEUP_DETECTION_FILTER_ASSEMBLY_HAPS_THRESHOLD, doc = "If enabled (set to non-zero), will apply the \"badness\" filter to compatable assembled haplotypes.", optional = true)
    public double assemblyBadReadThreshold = 0.4;
    @Hidden
    @Argument(fullName= PILEUP_DETECTION_EDIT_DISTANCE_BADNESS_FOR_ASSEMBLY_LONG_NAME, doc = "Pileup Detection: Reject alt reads with greater than this fraction of mismatching bases from the reference (proxied using the NM tag). (Requires '--"+PILEUP_DETECTION_LONG_NAME+"' argument)", optional = true)
    public double assemblyBadReadEditDistance = 0.12;

}
