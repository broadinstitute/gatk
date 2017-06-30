package org.broadinstitute.hellbender.tools.spark.pathseq;

import org.broadinstitute.barclay.argparser.Argument;

import java.io.Serializable;

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
    @Argument(doc = "Keep only clipped reads with length at least equal to the specified value",
            fullName = "minClippedReadLength",
            minValue = 0,
            minRecommendedValue = 31,
            optional = true)
    public int minReadLength = 31;
    @Argument(doc = "Max allowable number of ambiguous bases",
            fullName = "filterMaxAmbiguousBases",
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
    @Argument(doc = "Host alignment coverage score threshold, in bp",
            fullName = "filterMinCoverage",
            minValue = 1,
            optional = true)
    public int minCoverage = 31;
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
    @Argument(doc = "Write remaining reads and time elapsed after each step to file. Substantially reduces performance.",
            fullName = "metricsFile",
            optional = true)
    public String metricsFileUri = null;
}
