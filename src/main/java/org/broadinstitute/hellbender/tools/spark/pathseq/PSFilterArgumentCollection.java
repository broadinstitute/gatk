package org.broadinstitute.hellbender.tools.spark.pathseq;

import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.GATKPlugin.GATKReadFilterPluginDescriptor;
import org.broadinstitute.hellbender.engine.filters.AmbiguousBaseReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadLengthReadFilter;

import java.io.Serializable;
import java.util.List;

public final class PSFilterArgumentCollection implements Serializable {

    private static final long serialVersionUID = 1L;
    public final int bwaThreads = 1;

    public static final String KMER_FILE_PATH_LONG_NAME = "kmer-file";
    public static final String KMER_FILE_PATH_SHORT_NAME = "K";
    public static final String FILTER_METRICS_FILE_LONG_NAME = "filter-metrics";
    public static final String FILTER_METRICS_FILE_SHORT_NAME = "FM";
    public static final String FILTER_BWA_IMAGE_LONG_NAME = "filter-bwa-image";
    public static final String FILTER_BWA_IMAGE_SHORT_NAME = "FI";

    public static final String IS_HOST_ALIGNED_LONG_NAME = "is-host-aligned";
    public static final String IS_HOST_ALIGNED_SHORT_NAME = IS_HOST_ALIGNED_LONG_NAME;
    public static final String SKIP_FILTERS_LONG_NAME = "skip-quality-filters";
    public static final String SKIP_FILTERS_SHORT_NAME = SKIP_FILTERS_LONG_NAME;
    public static final String SKIP_PRE_BWA_REPARTITION_LONG_NAME = "skip-pre-bwa-repartition";
    public static final String SKIP_PRE_BWA_REPARTITION_SHORT_NAME = SKIP_PRE_BWA_REPARTITION_LONG_NAME;
    public static final String MIN_CLIPPED_READ_LENGTH_LONG_NAME = "min-clipped-read-length";
    public static final String MIN_CLIPPED_READ_LENGTH_SHORT_NAME = MIN_CLIPPED_READ_LENGTH_LONG_NAME;
    public static final String MAX_MASKED_BASES_LONG_NAME = "max-masked-bases";
    public static final String MAX_MASKED_BASES_SHORT_NAME = MAX_MASKED_BASES_LONG_NAME;
    public static final String MIN_BASE_QUALITY_LONG_NAME = "min-base-quality";
    public static final String MIN_BASE_QUALITY_SHORT_NAME = MIN_BASE_QUALITY_LONG_NAME;
    public static final String QUALITY_SCORE_THRESHOLD_LONG_NAME = "quality-threshold";
    public static final String QUALITY_SCORE_THRESHOLD_SHORT_NAME = QUALITY_SCORE_THRESHOLD_LONG_NAME;
    public static final String DUST_MASK_QUALITY_LONG_NAME = "dust-mask-quality";
    public static final String DUST_MASK_QUALITY_SHORT_NAME = DUST_MASK_QUALITY_LONG_NAME;
    public static final String DUST_WINDOW_SIZE_LONG_NAME = "dust-window";
    public static final String DUST_WINDOW_SIZE_SHORT_NAME = DUST_WINDOW_SIZE_LONG_NAME;
    public static final String DUST_T_SCORE_LONG_NAME = "dust-t";
    public static final String DUST_T_SCORE_SHORT_NAME = DUST_T_SCORE_LONG_NAME;
    public static final String HOST_MIN_IDENTITY_LONG_NAME = "host-min-identity";
    public static final String HOST_MIN_IDENTITY_SHORT_NAME = HOST_MIN_IDENTITY_LONG_NAME;
    public static final String HOST_KMER_COUNT_THRESHOLD_LONG_NAME = "host-kmer-thresh";
    public static final String HOST_KMER_COUNT_THRESHOLD_SHORT_NAME = HOST_KMER_COUNT_THRESHOLD_LONG_NAME;
    public static final String FILTER_BWA_SEED_LENGTH_LONG_NAME = "filter-bwa-seed-length";
    public static final String FILTER_BWA_SEED_LENGTH_SHORT_NAME = FILTER_BWA_SEED_LENGTH_LONG_NAME;
    public static final String FILTER_READS_PER_PARTITION_LONG_NAME = "filter-reads-per-partition";
    public static final String FILTER_READS_PER_PARTITION_SHORT_NAME = FILTER_READS_PER_PARTITION_LONG_NAME;
    public static final String FILTER_DUPLICATES_LONG_NAME = "filter-duplicates";
    public static final String FILTER_DUPLICATES_SHORT_NAME = FILTER_DUPLICATES_LONG_NAME;
    public static final String MAX_ADAPTER_MISMATCHES_LONG_NAME = "max-adapter-mismatches";
    public static final String MAX_ADAPTER_MISMATCHES_SHORT_NAME = MAX_ADAPTER_MISMATCHES_LONG_NAME;
    public static final String MIN_ADAPTER_LENGTH_LONG_NAME = "min-adapter-length";
    public static final String MIN_ADAPTER_LENGTH_SHORT_NAME = MIN_ADAPTER_LENGTH_LONG_NAME;

    /**
     * PathSeq will rapidly filter the reads if they are aligned to a host reference, thus reducing run time.
     */
    @Argument(doc = "Set if the input BAM is aligned to the host",
            fullName = IS_HOST_ALIGNED_LONG_NAME,
            shortName = IS_HOST_ALIGNED_SHORT_NAME,
            optional = true)
    public boolean alignedInput = false;

    @Argument(doc = "Path to host k-mer file generated with PathSeqBuildKmers. K-mer filtering is skipped if this is not specified.",
            fullName = KMER_FILE_PATH_LONG_NAME,
            shortName = KMER_FILE_PATH_SHORT_NAME,
            optional = true)
    public String kmerFilePath = null;

    @Argument(doc = "Skip low-quality and low-complexity read filtering",
            fullName = SKIP_FILTERS_LONG_NAME,
            shortName = SKIP_FILTERS_SHORT_NAME,
            optional = true)
    public boolean skipFilters = false;

    /**
     * <p>
     *     Advanced optimization option that should be used only in the case of inputs with a high proportion of microbial
     *     reads that are not host-aligned/coordinate-sorted.
     * </p>
     * <p>
         * In the filter tool, the input reads are initially divided up into smaller partitions (default size is usually
         * the size of one HDFS block, or ~64MB) that Spark works on in parallel. In samples with a low proportion of microbial
         * reads (e.g. < 1%), the steps leading up to the host BWA alignment will whittle these partitions down to a small
         * fraction of their original size. At that point, the distribution of reads across the partitions may be unbalanced.
     * </p>
     * <p>
     *     For example, say the input is 256MB and Spark splits this into 4 even partitions. It is possible that, after
     *     running through the quality filters and host kmer search, there are 5% remaining in partition #1, 8% in partition #2,
     *     2% in partition #3, and 20% in partition #4. Thus there is an imbalance of work across the partitions. To
     *     correct this, a "reparitioning" is invoked that distributes the reads evenly. Note this is especially important
     *     for host-aligned, coordinate-sorted inputs, in which unmapped reads would be concentrated in the last partitions.
     * </p>
     * <p>
     *     If, however, the proportion of microbial reads is higher, say 30%, then the partitions are generally more
     *     balanced (except for in the aforementioned coordinate-sorted case). In this case, the time spent doing
     *     the repartitioning is usually greater than the time saved by rebalancing, and this option should be invoked.
     * </p>
     */
    @Advanced
    @Argument(doc = "Skip pre-BWA repartition. Set to true for inputs with a high proportion of microbial reads that " +
            "are not host coordinate-sorted.",
            fullName = SKIP_PRE_BWA_REPARTITION_LONG_NAME,
            shortName = SKIP_PRE_BWA_REPARTITION_SHORT_NAME,
            optional = true)
    public boolean skipPreBwaRepartition = false;

    /**
     * Reads are trimmed based on base call quality and low-complexity content. Decreasing the value will enhance pathogen
     * detection (higher sensitivity) but also result in undesired false positives and ambiguous microbe alignments
     * (lower specificity).
     */
    @Argument(doc = "Minimum length of reads after quality trimming",
            fullName = MIN_CLIPPED_READ_LENGTH_LONG_NAME,
            shortName = MIN_CLIPPED_READ_LENGTH_SHORT_NAME,
            minValue = 1,
            minRecommendedValue = 31,
            optional = true)
    public int minReadLength = 31;

    /**
     * This is the threshold for filtering reads based on the number of 'N' values present in the sequence. Note that
     * the low-complexity DUST filter and quality filter mask using 'N' bases. Therefore, this parameter is the threshold
     * for the sum of:
     * <ul>
     *     <li>The number of N's in the original input read</li>
     *     <li>The number of low-quality base calls</li>
     *     <li>The number of low-complexity bases</li>
     * </ul>
     */
    @Argument(doc = "Max allowable number of masked bases per read",
            fullName = MAX_MASKED_BASES_LONG_NAME,
            shortName = MAX_MASKED_BASES_SHORT_NAME,
            minValue = 0,
            optional = true)
    public int maxAmbiguousBases = 2;

    @Argument(doc = "Bases below this call quality will be masked with 'N'",
            fullName = MIN_BASE_QUALITY_LONG_NAME,
            shortName = MIN_BASE_QUALITY_SHORT_NAME,
            minValue = 1,
            optional = true)
    public int qualPhredThresh = 15;

    /**
     * Controls the stingency of base call quality-based read trimming. Higher values result in more trimming.
     */
    @Argument(doc = "Quality score trimmer threshold",
            fullName = QUALITY_SCORE_THRESHOLD_LONG_NAME,
            shortName = QUALITY_SCORE_THRESHOLD_SHORT_NAME,
            minValue = 1,
            optional = true)
    public int readTrimThresh = 15;

    @Argument(doc = "Base quality to assign low-complexity bases",
            fullName = DUST_MASK_QUALITY_LONG_NAME,
            shortName = DUST_MASK_QUALITY_SHORT_NAME,
            optional = true)
    public int dustMask = 2;

    @Argument(doc = "DUST algorithm window size",
            fullName = DUST_WINDOW_SIZE_LONG_NAME,
            shortName = DUST_WINDOW_SIZE_SHORT_NAME,
            optional = true)
    public int dustW = 64;

    /**
     * Controls the stringency of low-complexity filtering.
     */
    @Argument(doc = "DUST algorithm score threshold",
            fullName = DUST_T_SCORE_LONG_NAME,
            shortName = DUST_T_SCORE_SHORT_NAME,
            optional = true)
    public double dustT = 20.0;

    /**
     * Controls the stringency of read filtering based on alignment to the host reference. The identity score is defined
     * as the number of matching bases less the number of deletions in the alignment.
     */
    @Argument(doc = "Host alignment identity score threshold, in bp",
            fullName = HOST_MIN_IDENTITY_LONG_NAME,
            shortName = HOST_MIN_IDENTITY_SHORT_NAME,
            minValue = 1,
            optional = true)
    public int minIdentity = 30;

    /**
     * Controls the stringency of read filtering based on host k-mer matching. Reads with at least this many matching
     * k-mers in the host reference will be filtered.
     */
    @Argument(doc = "Host kmer count threshold.",
            fullName = HOST_KMER_COUNT_THRESHOLD_LONG_NAME,
            shortName = HOST_KMER_COUNT_THRESHOLD_SHORT_NAME,
            minValue = 1,
            optional = true)
    public int hostKmerThresh = 1;

    /**
     * This file should be generated using BwaMemIndexImageCreator.
     */
    @Argument(doc = "The BWA image file of the host reference. This must be distributed to local disk on each node.",
            fullName = FILTER_BWA_IMAGE_LONG_NAME,
            shortName = FILTER_BWA_IMAGE_SHORT_NAME,
            optional = true)
    public String indexImageFile = null;

    /**
     * Controls the sensitivity of BWA alignment to the host reference. Shorter seed lengths will enhance detection of
     * host reads during the subtraction phase but will also increase run time.
     */
    @Argument(doc = "Minimum seed length for the host BWA alignment.",
            fullName = FILTER_BWA_SEED_LENGTH_LONG_NAME,
            shortName = FILTER_BWA_SEED_LENGTH_SHORT_NAME,
            minValue = 1,
            minRecommendedValue = 11,
            optional = true)
    public int minSeedLength = 19;

    /**
     * This is a parameter for fine-tuning memory performance. Lower values may result in less memory usage but possibly
     * at the expense of greater computation time.
     */
    @Advanced
    @Argument(doc = "Estimated reads per partition after quality, kmer, and BWA filtering",
            fullName = FILTER_READS_PER_PARTITION_LONG_NAME,
            shortName = FILTER_READS_PER_PARTITION_SHORT_NAME,
            minValue = 1,
            optional = true)
    public int filterReadsPerPartition = 200000;

    /**
     * If true, then for any two reads with identical sequences (or identical to the other's reverse complement), one
     * will be filtered.
     */
    @Argument(doc = "Filter duplicate reads",
            fullName = FILTER_DUPLICATES_LONG_NAME,
            shortName = FILTER_DUPLICATES_SHORT_NAME,
            optional = true)
    public boolean filterDuplicates = true;

    /**
     * Adapter trimming will require a match of at least this length to a known adapter.
     */
    @Argument(doc = "Minimum length of adapter sequence to trim",
            fullName = MIN_ADAPTER_LENGTH_LONG_NAME,
            shortName = MIN_ADAPTER_LENGTH_SHORT_NAME,
            minValue = 1,
            optional = true)
    public int minAdapterLength = 12;

    @Argument(doc = "Maximum number of mismatches for adapter trimming",
            fullName = MAX_ADAPTER_MISMATCHES_LONG_NAME,
            shortName = MAX_ADAPTER_MISMATCHES_SHORT_NAME,
            minValue = 0,
            optional = true)
    public int maxAdapterMismatches = 1;

    /**
     * If specified, records the number of reads remaining after each of the following steps:
     * <ul>
     *     <li>Pre-aligned host read filtering</li>
     *     <li>Low-quality and low-complexity sequence filtering</li>
     *     <li>Host read subtraction</li>
     *     <li>Read deduplication</li>
     * </ul>
     * <p>It also provides the following:</p>
     * <ul>
     *     <li>Number of low-quality and low-complexity reads removed</li>
     *     <li>Number of host reads removed</li>
     *     <li>Number of duplicate reads removed</li>
     *     <li>Final number of reads</li>
     *     <li>Final number of paired reads</li>
     *     <li>Final number of unpaired reads</li>
     * </ul>
     *<p>Note that using this option may substantially increase runtime.</p>
     */
    @Argument(doc = "Log counts of filtered reads to this file",
            fullName = FILTER_METRICS_FILE_LONG_NAME,
            shortName = FILTER_METRICS_FILE_SHORT_NAME,
            optional = true)
    public String filterMetricsFileUri = null;

    public void doReadFilterArgumentWarnings(final GATKReadFilterPluginDescriptor pluginDescriptor, final Logger logger) {
        final List<ReadFilter> readFilters = pluginDescriptor.getResolvedInstances();
        for (final ReadFilter filter : readFilters) {
            if (filter.getClass().isAssignableFrom(AmbiguousBaseReadFilter.class)) {
                logger.warn("Detected the use of AmbiguousBaseReadFilter, which is applied before the PathSeq " +
                        "base masking steps. Did you mean to use --max-masked-bases, which is applied after masking?");
            } else if (filter.getClass().isAssignableFrom(ReadLengthReadFilter.class)) {
                logger.warn("Detected the use of ReadLengthReadFilter, which is applied before the PathSeq " +
                        "clipping steps. Did you mean to use --min-clipped-read-length, which is applied after clipping?");
            }
        }
    }
}
