package org.broadinstitute.hellbender.tools.spark.pathseq;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import htsjdk.samtools.SAMFileHeader;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSink;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadsWriteFormat;

import java.io.IOException;

@CommandLineProgramProperties(summary = "Aligns reads to a pathogen reference using Bwa Mem",
        oneLineSummary = "Aligns reads to pathogen reference",
        programGroup = SparkProgramGroup.class)
public final class PathSeqPathogenAlignSpark extends GATKSparkTool {

    private static final long serialVersionUID = 1L;
    private static final int BWA_MAX_ALT_HITS = 5000;
    private static final int BWA_SCORE_THRESHOLD = 30;

    @Argument(doc = "Output path URI", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    public String OUTPUT;

    @Argument(doc = "URI to paired reads bam",
            fullName = "pairedInput",
            optional = true)
    public String INPUT_PAIRED = null;

    @Argument(doc = "URI to unpaired reads bam",
            fullName = "unpairedInput",
            optional = true)
    public String INPUT_UNPAIRED = null;

    @Argument(doc = "URI of the bwa index image file",
            fullName = "bwamemIndexImage")
    public String BWA_IMAGE;

    @Argument(doc = "Reference URI corresponding to the image file",
            fullName = "referencePath")
    public String REFERENCE_PATH;

    @Argument(doc = "Minimum seed length for the BWA alignment.",
            fullName = "minSeedLength",
            minValue = 1,
            minRecommendedValue = 11,
            optional = true)
    public int SEED_LENGTH = 19;

    @Argument(doc = "Number of bwa threads",
            fullName = "bwaThreads",
            minValue = 1,
            optional = true)
    public int BWA_THREADS = 1;

    @Argument(doc = "Add new alignments to the reads in the input BAM using the SA tag, rather than overwriting them",
            fullName = "preserveInputAlignments",
            optional = true)
    public boolean PRESERVE_ALIGNMENTS = false;

    protected static JavaRDD<GATKRead> doBwaAlignment(final JavaRDD<GATKRead> reads,
                                                      final Broadcast<SAMFileHeader> headerBroadcast,
                                                      final String indexFileName,
                                                      final int minSeedLength,
                                                      final int numThreads,
                                                      final boolean preserveAlignments,
                                                      final boolean isPaired) {
        return reads.mapPartitions(itr -> (new PSBwaAligner(indexFileName,
                headerBroadcast.value(), minSeedLength, numThreads, BWA_MAX_ALT_HITS, BWA_SCORE_THRESHOLD,
                isPaired, preserveAlignments)).apply(itr));
    }

    private boolean doAlignmentAndWriteOutput(final String path, final boolean isPaired, final JavaSparkContext ctx,
                                           final PipelineOptions options, final ReadsSparkSource readsSource) {

        JavaRDD<GATKRead> reads = null;
        Broadcast<SAMFileHeader> headerBroadcast = null;

        if (path != null) {
            if (BucketUtils.fileExists(path, options)) {
                final SAMFileHeader header = readsSource.getHeader(path, null);
                reads = readsSource.getParallelReads(path, null, null, bamPartitionSplitSize);
                PSUtils.addReferenceSequencesToHeader(header, REFERENCE_PATH, getReferenceWindowFunction(), options);
                headerBroadcast = ctx.broadcast(header);
            } else {
                logger.warn("Could not find file " + path + ". Skipping...");
                return false;
            }
        } else {
            return false;
        }
        final JavaRDD<GATKRead> alignedReads = doBwaAlignment(reads, headerBroadcast, BWA_IMAGE, SEED_LENGTH,
                BWA_THREADS, PRESERVE_ALIGNMENTS, isPaired);

        try {
            final String outputSuffix;
            if (isPaired) {
                outputSuffix = ".paired.bam";
            } else {
                outputSuffix = ".unpaired.bam";
            }
            ReadsSparkSink.writeReads(ctx, OUTPUT + outputSuffix, REFERENCE_PATH, alignedReads, headerBroadcast.value(),
                    shardedOutput ? ReadsWriteFormat.SHARDED : ReadsWriteFormat.SINGLE,
                    PSUtils.pathseqGetRecommendedNumReducers(path, numReducers, options, getTargetPartitionSize()));
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(OUTPUT, "writing failed", e);
        }
        return true;
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {

        final ReadsSparkSource readsSource = new ReadsSparkSource(ctx, readArguments.getReadValidationStringency());
        final PipelineOptions options = getAuthenticatedGCSOptions();
        boolean bPairedSuccess = doAlignmentAndWriteOutput(INPUT_PAIRED, true, ctx, options, readsSource);
        boolean bUnpairedSuccess = doAlignmentAndWriteOutput(INPUT_UNPAIRED, false, ctx, options, readsSource);
        if (!bPairedSuccess && !bUnpairedSuccess) {
            throw new UserException.BadInput("No reads were loaded. Ensure --pairedInput and/or --unpairedInput are set and valid.");
        }
    }

}
