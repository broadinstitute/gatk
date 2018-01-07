package org.broadinstitute.hellbender.tools.spark.pathseq;

import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.GATKPlugin.GATKReadFilterPluginDescriptor;
import org.broadinstitute.hellbender.engine.filters.AmbiguousBaseReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadLengthReadFilter;

import java.io.Serializable;
import java.util.List;

public final class PSFilterArgumentCollection implements Serializable {

    private static final long serialVersionUID = 1L;

    @Argument(doc = "Set if the input BAM is aligned to the host",
            fullName = "isHostAligned",
            optional = true)
    public boolean alignedInput = false;
    @Argument(doc = "Path to host kmer library generated with PathSeqBuildKmers. Skipped if not specified.",
            fullName = "kmerLibraryPath",
            optional = true)
    public String kmerLibPath = null;
    @Argument(doc = "Skip read quality and low-complexity filtering",
            fullName = "skipFilters",
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
    @Argument(doc = "Skip pre-BWA repartition. Set to true for inputs with a high proportion of microbial reads that " +
            "are not host coordinate-sorted.",
            fullName = "skipPreBwaRepartition",
            optional = true)
    public boolean skipPreBwaRepartition = false;
    @Argument(doc = "Keep only clipped reads with length at least equal to the specified value",
            fullName = "minClippedReadLength",
            minValue = 0,
            minRecommendedValue = 31,
            optional = true)
    public int minReadLength = 31;
    @Argument(doc = "Max allowable number of masked bases per read",
            fullName = "maxMaskedBases",
            minValue = 0,
            optional = true)
    public int maxAmbiguousBases = 2;
    @Argument(doc = "Bases below this read quality will be transformed into 'N'",
            fullName = "minBaseQuality",
            minValue = 1,
            optional = true)
    public int qualPhredThresh = 15;
    @Argument(doc = "Quality score trimmer threshold",
            fullName = "readTrimmerThreshold",
            minValue = 1,
            optional = true)
    public int readTrimThresh = 15;
    @Argument(doc = "Base quality to assign DUST masked bases",
            fullName = "dustPhred",
            optional = true)
    public int dustMask = 2;
    @Argument(doc = "DUST window size",
            fullName = "dustWindowSize",
            optional = true)
    public int dustW = 64;
    @Argument(doc = "DUST score threshold",
            fullName = "dustTScore",
            optional = true)
    public double dustT = 20.0;
    @Argument(doc = "Host alignment identity score threshold, in bp",
            fullName = "filterMinIdentity",
            minValue = 1,
            optional = true)
    public int minIdentity = 30;
    @Argument(doc = "Host kmer count threshold.",
            fullName = "hostKmerThreshold",
            minValue = 1,
            optional = true)
    public int hostKmerThresh = 1;
    @Argument(doc = "The BWA image file of the host reference. This must be distributed to local disk on each node.",
            fullName = "filterBwaImage",
            optional = true)
    public String indexImageFile = null;
    public int bwaThreads = 1;
    @Argument(doc = "Minimum seed length for the host BWA alignment.",
            fullName = "filterMinSeedLength",
            minValue = 1,
            minRecommendedValue = 11,
            optional = true)
    public int minSeedLength = 19;
    @Argument(doc = "Estimated reads per partition after quality, kmer, and BWA filtering",
            fullName = "filterReadsPerPartition",
            minValue = 1,
            optional = true)
    public int readsPerPartition = 200000;
    @Argument(doc = "Filter duplicate reads",
            fullName = "filterDuplicates",
            optional = true)
    public boolean filterDuplicates = true;

    public void doReadFilterArgumentWarnings(final GATKReadFilterPluginDescriptor pluginDescriptor, final Logger logger) {
        final List<ReadFilter> readFilters = pluginDescriptor.getResolvedInstances();
        for (final ReadFilter filter : readFilters) {
            if (filter.getClass().isAssignableFrom(AmbiguousBaseReadFilter.class)) {
                logger.warn("Detected the use of AmbiguousBaseReadFilter, which is applied before the PathSeq " +
                        "base masking steps. Did you mean to use --maxMaskedBases, which is applied after masking?");
            } else if (filter.getClass().isAssignableFrom(ReadLengthReadFilter.class)) {
                logger.warn("Detected the use of ReadLengthReadFilter, which is applied before the PathSeq " +
                        "clipping steps. Did you mean to use --minClippedReadLength, which is applied after clipping?");
            }
        }
    }
    @Argument(doc = "Maximum number of mismatches for adapter trimming",
            fullName = "maxAdapterMismatches",
            minValue = 0,
            optional = true)
    public int maxAdapterMismatches = 1;
    @Argument(doc = "Minimum length of adapter sequence to trim",
            fullName = "minAdapterLength",
            minValue = 1,
            optional = true)
    public int minAdapterLength = 12;
}
