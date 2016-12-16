package org.broadinstitute.hellbender.tools.spark.pathseq;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import htsjdk.samtools.SAMFileHeader;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.filters.AmbiguousBaseReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadLengthReadFilter;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.tools.spark.utils.ReadFilterSparkifier;
import org.broadinstitute.hellbender.tools.spark.utils.ReadTransformerSparkifier;
import org.broadinstitute.hellbender.transformers.BaseQualityClipReadTransformer;
import org.broadinstitute.hellbender.transformers.BaseQualityReadTransformer;
import org.broadinstitute.hellbender.transformers.DUSTReadTransformer;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.IOException;
import java.io.OutputStream;
import java.nio.charset.Charset;
import java.text.DecimalFormat;
import java.text.NumberFormat;

/**
 * This Spark tool is the first step in the PathSeq pipeline. The output is the set of high quality non-host reads.
 *
 * Filtering steps:
 * 1) Remove secondary and supplementary reads
 * 2) Mask repetitive sequences with 'N' and base quality --dustPhred using symmetric DUST
 * 3) Hard clip read ends using base qualities
 * 4) Remove reads shorter than --minClippedReadLength
 * 5) Mask bases whose Phred score is less than --minBaseQuality with 'N'
 * 6) Remove reads whose fraction of bases that are 'N' is greater than --maxAmbiguousBaseFraction
 * 7) If specified, Remove reads containing one or more kmers --kmerLibraryPath
 * 8) If specified, remove reads that align to the host BWA image --bwamemIndexImage with at least --minCoverage and --minIdentity
 * 9) If --filterDuplicates is set, remove exact duplicate reads (not using Mark Duplicates because it requires aligned reads)
 *
 * Notes:
 *
 * - Steps 2 - 6 can be skipped by setting --skipFilters.
 * - The tool assumes the BAM file is unaligned by default. If the BAM is aligned, use --isHostAligned to filter.
 * - Output will be two BAM files, OUTPUT_PATH.paired.bam and OUTPUT_PATH.unpaired.bam containing paired and unpaired reads.
 * - If the resulting set of reads is empty, the file will not be written.
 * - If --unpaired is set, pairedness flags will not be corrected after filtering, and all reads will be written to
 *     OUTPUT_PATH.unpaired.bam. This avoids two shuffles.
 *
 */
@CommandLineProgramProperties(summary = "Read preprocessing and host organism filtering on reads from a BAM file",
        oneLineSummary = "PathSeqFilter on Spark",
        programGroup = SparkProgramGroup.class)
public final class PathSeqFilterSpark extends GATKSparkTool {

    private static final long serialVersionUID = 1L;

    @Argument(doc = "Base uri for the output file(s). By default, two files are written appending '.paired.bam' " +
            "and '.unpaired.bam'",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    public String OUTPUT_PATH;
    @Argument(doc = "Set if the input BAM is aligned to the host",
            fullName = "isHostAligned",
            optional = true)
    public boolean ALIGNED_INPUT = false;
    @Argument(doc = "Path to kmer library generated with PathSeqKmerSpark. Skipped if not specified.",
            fullName = "kmerLibraryPath",
            optional = true)
    public String KMER_LIB_PATH = null;
    @Argument(doc = "Skip read quality and low-complexity filtering",
            fullName = "skipFilters",
            optional = true)
    public boolean SKIP_FILTERS = false;
    @Argument(doc = "Keep only clipped reads with length at least equal to the specified value",
            fullName = "minClippedReadLength",
            minValue = 0,
            minRecommendedValue = 31,
            optional = true)
    public int MIN_READ_LENGTH = 31;
    @Argument(doc = "Max allowable fraction of ambiguous bases",
            fullName = "maxAmbiguousBaseFraction",
            minValue = 0.0,
            maxValue = 1.0,
            optional = true)
    public double FRAC_N_THRESHOLD = 0.03;
    @Argument(doc = "Bases below this read quality will be transformed into 'N'",
            fullName = "minBaseQuality",
            minValue = 1,
            optional = true)
    public int QUAL_PHRED_THRESH = 15;
    @Argument(doc = "Size of kmers in the kmer library",
            fullName = "kmerSize",
            minValue = 1,
            maxValue = 31,
            minRecommendedValue = 25,
            optional = true)
    public int KMER_SIZE = 31;
    @Argument(doc = "Quality score trimmer threshold",
            fullName = "readTrimmerThreshold",
            minValue = 1,
            optional = true)
    public int READ_TRIM_THRESH = 15;
    @Argument(doc = "Base quality to assign DUST masked bases",
            fullName = "dustPhred",
            optional = true)
    public int DUST_MASK = 2;
    @Argument(doc = "DUST window size",
            fullName = "dustWindowSize",
            optional = true)
    public int DUST_W = 64;
    @Argument(doc = "DUST score threshold",
            fullName = "dustTScore",
            optional = true)
    public double DUST_T = 20.0;
    @Argument(doc = "Minimum seed length for the host BWA alignment.",
            fullName = "minSeedLength",
            minValue = 1,
            minRecommendedValue = 11,
            optional = true)
    public int MIN_SEED_LENGTH = 19;
    @Argument(doc = "Host alignment coverage threshold, in bp",
            fullName = "minCoverage",
            minValue = 1,
            optional = true)
    public int MIN_COVERAGE = 31;
    @Argument(doc = "Host alignment identity threshold, in bp",
            fullName = "minIdentity",
            minValue = 1,
            optional = true)
    public int MIN_IDENTITY = 30;
    @Argument(doc = "Comma-separated list of kmer base indices to mask (0-based), must have an even number of entries " +
            "and be in ascending order",
            fullName = "kmerMask",
            optional = true)
    public String KMER_MASK = "";
    @Argument(doc = "Host kmer count threshold.",
            fullName = "hostKmerThreshold",
            minValue = 1,
            optional = true)
    public int HOST_KMER_THRESH = 1;
    @Argument(doc = "The bwa mem index image file of the host reference, distributed to each executor. Must also " +
            "specify --referencePath",
            fullName = "bwamemIndexImage",
            optional = true)
    public String INDEX_IMAGE_FILE = null;
    @Argument(doc = "Number of bwa threads",
            fullName = "bwaThreads",
            minValue = 1,
            optional = true)
    public int BWA_THREADS = 1;
    @Argument(doc = "Filter duplicate reads",
            fullName = "filterDuplicates",
            optional = true)
    public boolean FILTER_DUPLICATES = false;
    @Argument(doc = "Ignore pairing information for deduplication and do not correct pairedness flags after filtering " +
            "(faster). All reads will be written to the unpaired read BAM.",
            fullName = "unpaired",
            optional = true)
    public boolean UNPAIRED = false;
    @Argument(doc = "Write remaining reads and time elapsed after each step to file. Substantially reduces performance.",
            fullName = "metricsFile",
            optional = true)
    public String METRICS_FILE_URI = null;

    protected long initialTimeMillis = 0;
    protected OutputStream metricsOutputStream = null;

    @Override
    public boolean requiresReads() {
        return true;
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {

        final byte[] kmerMask = PSUtils.parseMask(KMER_MASK, KMER_SIZE);
        final PipelineOptions options = getAuthenticatedGCSOptions();
        final SAMFileHeader header = getHeaderForReads();

        JavaRDD<GATKRead> reads;
        initializeMetrics(options);

        reads = getReads();
        recordMetrics("input", reads, options);

        reads = PSUtils.primaryReads(reads);
        recordMetrics("primary_reads", reads, options);

        if (ALIGNED_INPUT) {
            reads = reads.filter(new ReadFilterSparkifier(new HostAlignmentReadFilter(MIN_COVERAGE, MIN_IDENTITY)));
            recordMetrics("prealigned_filter", reads, options);
        }

        if (!SKIP_FILTERS) {

            //Apply DUST masking
            reads = reads.map(new ReadTransformerSparkifier(new DUSTReadTransformer(DUST_MASK, DUST_W, DUST_T)));

            //Apply base quality hard clipping
            reads = reads.map(new ReadTransformerSparkifier(new BaseQualityClipReadTransformer(READ_TRIM_THRESH)));

            //Filter reads with less than MIN_READ_LENGTH bases
            reads = reads.filter(new ReadFilterSparkifier(new ReadLengthReadFilter(MIN_READ_LENGTH, Integer.MAX_VALUE)));
            recordMetrics("quality_1", reads, options);

            //Change low-quality bases to 'N'
            reads = reads.map(new ReadTransformerSparkifier(new BaseQualityReadTransformer(QUAL_PHRED_THRESH)));

            //Filter reads with too many 'N's
            reads = reads.filter(new ReadFilterSparkifier(new AmbiguousBaseReadFilter(FRAC_N_THRESHOLD)));
            recordMetrics("quality_2", reads, options);
        }

        //Kmer filtering
        if (KMER_LIB_PATH != null) {
            reads = PSFilterUtils.doKmerFiltering(reads, options, KMER_LIB_PATH, kmerMask, KMER_SIZE, HOST_KMER_THRESH);
            recordMetrics("kmer", reads, options);
        }

        //Bwa host alignment filtering
        if (INDEX_IMAGE_FILE != null) {
            reads = PSFilterUtils.doBwaFilter(ctx, reads, header, INDEX_IMAGE_FILE, MIN_SEED_LENGTH, BWA_THREADS, MIN_COVERAGE, MIN_IDENTITY);
            recordMetrics("bwa", reads, options);
        }

        //Filter duplicates
        if (FILTER_DUPLICATES) {
            if (!UNPAIRED) {
                reads = PSFilterUtils.doSetPairFlags(reads);
            }
            reads = PSFilterUtils.filterDuplicateSequences(reads);
            recordMetrics("duplicates", reads, options);
        }

        if (!UNPAIRED) {
            //Unset paired read flags for reads that are not paired
            reads = PSFilterUtils.doSetPairFlags(reads);
            recordMetrics("set_flags", reads, options);
            final JavaRDD<GATKRead> pairedReads = PSUtils.pairedReads(reads);
            final JavaRDD<GATKRead> unpairedReads = PSUtils.unpairedReads(reads);
            if (!pairedReads.isEmpty()) {
                writeReads(ctx, OUTPUT_PATH + ".paired.bam", pairedReads);
            } else {
                logger.info("No paired reads to write - BAM will not be written.");
            }
            if (!unpairedReads.isEmpty()) {
                writeReads(ctx, OUTPUT_PATH + ".unpaired.bam", unpairedReads);
            } else {
                logger.info("No unpaired reads to write - BAM will not be written.");
            }
        } else {
            if (!reads.isEmpty()) {
                writeReads(ctx, OUTPUT_PATH + ".unpaired.bam", reads);
            } else {
                logger.info("No reads to write - BAM will not be written.");
            }
        }

        recordMetrics("output", reads, options);
        closeMetrics();
    }

    private void recordMetrics(final String name, final JavaRDD<GATKRead> reads, PipelineOptions options) {
        try {
            if (METRICS_FILE_URI != null) {
                if (metricsOutputStream == null) initializeMetrics(options);
                final long currentTimeMillis = System.currentTimeMillis();
                final double elapsedTimeSeconds = (currentTimeMillis - initialTimeMillis) / 1000.0;
                NumberFormat formatter = new DecimalFormat("#0.00");
                final String outputString = name + "\t" + reads.count() + "\t" + formatter.format(elapsedTimeSeconds) + "\n";
                if (metricsOutputStream != null) {
                    metricsOutputStream.write(outputString.getBytes(Charset.defaultCharset()));
                    metricsOutputStream.flush();
                } else {
                    logger.warn("Metrics output stream uninitialized. Metric named " + name + " will not be written.");
                }
            }
        } catch (final IOException e) {
            logger.warn("Could not write metrics to URI " + METRICS_FILE_URI, e);
        }
    }

    private void initializeMetrics(final PipelineOptions options) {
        try {
            if (METRICS_FILE_URI != null) {
                initialTimeMillis = System.currentTimeMillis();
                metricsOutputStream = BucketUtils.createFile(METRICS_FILE_URI, options);
                final String outputString = "step\treads\tseconds\n";
                metricsOutputStream.write(outputString.getBytes(Charset.defaultCharset()));
            }
        } catch (final IOException e) {
            logger.warn("Could not write metrics to URI " + METRICS_FILE_URI, e);
        }
    }

    private void closeMetrics() {
        if (metricsOutputStream != null) {
            try {
                metricsOutputStream.close();
            } catch (final IOException e) {
                logger.warn("Could not close metrics file stream", e);
            }
        }
    }

}
