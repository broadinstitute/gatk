package org.broadinstitute.hellbender.tools.spark.pathseq;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import htsjdk.samtools.SAMFileHeader;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.apache.spark.storage.StorageLevel;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.GATKPlugin.GATKReadFilterPluginDescriptor;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.PathSeqProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSink;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadsWriteFormat;
import scala.Tuple2;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

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

    @Argument(doc = "Number of reads per partition to use for alignment and scoring.",
            fullName = "readsPerPartition",
            optional = true,
            minValue = 100)
    public int readsPerPartition = 5000;

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

        //Filter
        final PSFilter filter = new PSFilter(ctx, filterArgs, getReads(), header);
        final Tuple2<JavaRDD<GATKRead>, JavaRDD<GATKRead>> result = filter.doFilter();
        JavaRDD<GATKRead> pairedReads = result._1;
        JavaRDD<GATKRead> unpairedReads = result._2;

        //Counting forces an action on the RDDs to guarantee we're done with the Bwa image and kmer filter
        final long numPairedReads = pairedReads.count();
        final long numUnpairedReads = unpairedReads.count();
        final long numTotalReads = numPairedReads + numUnpairedReads;
        logger.info("Number of paired reads after filtering: " + numPairedReads);
        logger.info("Number of unpaired reads after filtering: " + numUnpairedReads);
        logger.info("Number of total reads after filtering: " + numTotalReads);

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
        final PipelineOptions options = getAuthenticatedGCSOptions();
        PSBwaUtils.addReferenceSequencesToHeader(header, bwaArgs.referencePath, getReferenceWindowFunction(), options);
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

        //Write reads to BAM, if specified
        if (outputPath != null) {
            try {
                final int numReducers = Math.max(1, (int) (numTotalReads / 1000000));
                ReadsSparkSink.writeReads(ctx, outputPath, null, readsFinal, header,
                        shardedOutput ? ReadsWriteFormat.SHARDED : ReadsWriteFormat.SINGLE, numReducers);
            } catch (final IOException e) {
                throw new UserException.CouldNotCreateOutputFile(outputPath, "writing failed", e);
            }
        }
        aligner.close();
    }

}
