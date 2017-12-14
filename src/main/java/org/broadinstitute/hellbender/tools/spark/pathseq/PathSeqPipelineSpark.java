package org.broadinstitute.hellbender.tools.spark.pathseq;

import htsjdk.samtools.SAMFileHeader;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.apache.spark.storage.StorageLevel;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.GATKPlugin.GATKReadFilterPluginDescriptor;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.PathSeqProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSink;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.pathseq.loggers.*;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadsWriteFormat;
import scala.Tuple2;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

@DocumentedFeature
@CommandLineProgramProperties(summary = "Combined tool that performs all steps: read filtering, pathogen alignment, and abundance " +
        "scoring. Note that reducing reads per partition increases parallelism of pathogen alignment on Spark.",
        oneLineSummary = "Combined tool that performs all steps: read filtering, pathogen alignment, and abundance scoring",
        programGroup = PathSeqProgramGroup.class)
@BetaFeature
public class PathSeqPipelineSpark extends GATKSparkTool {

    private static final long serialVersionUID = 1L;

    @ArgumentCollection
    public PSFilterArgumentCollection filterArgs = new PSFilterArgumentCollection();

    @ArgumentCollection
    public PSBwaArgumentCollection bwaArgs = new PSBwaArgumentCollection();

    @ArgumentCollection
    public PSScoreArgumentCollection scoreArgs = new PSScoreArgumentCollection();

    @Argument(doc = "URI to the output BAM",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            optional = true)
    public String outputPath = null;

    @Argument(doc = "Log counts of filtered reads to this file",
            fullName = "filterMetricsFile",
            optional = true)
    public String filterMetricsFileUri = null;

    @Argument(doc = "Log counts of mapped and unmapped reads to this file",
            fullName = "scoreMetricsFile",
            optional = true)
    public String scoreMetricsFileUri = null;

    @Argument(doc = "Number of reads per partition to use for alignment and scoring.",
            fullName = "readsPerPartition",
            optional = true,
            minValue = 100)
    public int readsPerPartition = 5000;

    /**
     * Because numReducers is based on the input size, it causes too many partitions to be produced when the output size is much smaller.
     */
    @Argument(doc = "Number of reads per partition for output. Use this to control the number of sharded BAMs (not --numReducers).",
            fullName = "readsPerPartitionOutput",
            optional = true,
            minValue = 100,
            minRecommendedValue = 100000)
    public int readsPerPartitionOutput = 1000000;

    /**
     * Reduces number of partitions of paired reads, keeping pairs together.
     */
    private static JavaRDD<GATKRead> repartitionPairedReads(final JavaRDD<GATKRead> pairedReads, final int alignmentPartitions, final long numReads) {
        final int readsPerPartition = 1 + (int) (numReads / alignmentPartitions);
        return pairedReads.mapPartitions(iter -> PathSeqPipelineSpark.pairPartitionReads(iter, readsPerPartition))
                .repartition(alignmentPartitions)
                .flatMap(List::iterator);
    }

    /**
     * Maps partition of paired reads to a partition of Lists containing each pair. Assumes pairs are adjacent.
     */
    private static Iterator<List<GATKRead>> pairPartitionReads(final Iterator<GATKRead> iter, final int readsPerPartition) {
        final ArrayList<List<GATKRead>> readPairs = new ArrayList<>(readsPerPartition / 2);
        while (iter.hasNext()) {
            final List<GATKRead> list = new ArrayList<>(2);
            list.add(iter.next());
            if (!iter.hasNext()) throw new GATKException("Odd number of read pairs in paired reads partition");
            list.add(iter.next());
            if (!list.get(0).getName().equals(list.get(1).getName())) throw new GATKException("Pair did not have the same name in a paired reads partition");
            readPairs.add(list);
        }
        return readPairs.iterator();
    }

    @Override
    public boolean requiresReads() {
        return true;
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {

        filterArgs.doReadFilterArgumentWarnings(getCommandLineParser().getPluginDescriptor(GATKReadFilterPluginDescriptor.class), logger);
        SAMFileHeader header = PSUtils.checkAndClearHeaderSequences(getHeaderForReads(), filterArgs, logger);

        //Do not allow use of numReducers
        if (numReducers > 0) {
            throw new UserException.BadInput("Use --readsPerPartitionOutput instead of --numReducers.");
        }

        //Filter
        final Tuple2<JavaRDD<GATKRead>, JavaRDD<GATKRead>> filterResult;
        final PSFilter filter = new PSFilter(ctx, filterArgs, header);
        try (final PSFilterLogger filterLogger = filterMetricsFileUri != null ? new PSFilterFileLogger(getMetricsFile(), filterMetricsFileUri) : new PSFilterEmptyLogger()) {
            final JavaRDD<GATKRead> inputReads = getReads();
            filterResult = filter.doFilter(inputReads, filterLogger);
        }
        JavaRDD<GATKRead> pairedReads = filterResult._1;
        JavaRDD<GATKRead> unpairedReads = filterResult._2;

        //Counting forces an action on the RDDs to guarantee we're done with the Bwa image and kmer filter
        final long numPairedReads = pairedReads.count();
        final long numUnpairedReads = unpairedReads.count();
        final long numTotalReads = numPairedReads + numUnpairedReads;

        //Closes Bwa image, kmer filter, and metrics file if used
        //Note the host Bwa image before must be unloaded before trying to load the pathogen image
        filter.close();

        //Rebalance partitions using the counts
        final int numPairedPartitions = 1 + (int) (numPairedReads / readsPerPartition);
        final int numUnpairedPartitions = 1 + (int) (numUnpairedReads / readsPerPartition);
        pairedReads = repartitionPairedReads(pairedReads, numPairedPartitions, numPairedReads);
        unpairedReads = unpairedReads.repartition(numUnpairedPartitions);

        //Bwa pathogen alignment
        final PSBwaAlignerSpark aligner = new PSBwaAlignerSpark(ctx, bwaArgs);
        PSBwaUtils.addReferenceSequencesToHeader(header, bwaArgs.referencePath, getReferenceWindowFunction());
        final Broadcast<SAMFileHeader> headerBroadcast = ctx.broadcast(header);
        JavaRDD<GATKRead> alignedPairedReads = aligner.doBwaAlignment(pairedReads, true, headerBroadcast);
        JavaRDD<GATKRead> alignedUnpairedReads = aligner.doBwaAlignment(unpairedReads, false, headerBroadcast);

        //Cache this expensive result. Note serialization significantly reduces memory consumption.
        alignedPairedReads.persist(StorageLevel.MEMORY_AND_DISK_SER());
        alignedUnpairedReads.persist(StorageLevel.MEMORY_AND_DISK_SER());

        //Score pathogens
        final PSScorer scorer = new PSScorer(scoreArgs);
        final JavaRDD<GATKRead> readsFinal = scorer.scoreReads(ctx, alignedPairedReads, alignedUnpairedReads, header);

        //Clean up header
        header = PSBwaUtils.removeUnmappedHeaderSequences(header, readsFinal, logger);

        //Log read counts
        if (scoreMetricsFileUri != null) {
            try (final PSScoreLogger scoreLogger = new PSScoreFileLogger(getMetricsFile(), scoreMetricsFileUri)) {
                scoreLogger.logReadCounts(readsFinal);
            }
        }

        //Write reads to BAM, if specified
        if (outputPath != null) {
            try {
                //Reduce number of partitions since we previously went to ~5K reads per partition, which
                // is far too small for sharded output.
                final int numPartitions = Math.max(1, (int) (numTotalReads / readsPerPartitionOutput));
                final JavaRDD<GATKRead> readsFinalRepartitioned = readsFinal.coalesce(numPartitions, false);
                ReadsSparkSink.writeReads(ctx, outputPath, null, readsFinalRepartitioned, header,
                        shardedOutput ? ReadsWriteFormat.SHARDED : ReadsWriteFormat.SINGLE, numPartitions);
            } catch (final IOException e) {
                throw new UserException.CouldNotCreateOutputFile(outputPath, "writing failed", e);
            }
        }
        aligner.close();
    }

}
