package org.broadinstitute.hellbender.tools.spark.longread;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.LongReadAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.UncheckedIOException;
import java.util.Comparator;
import java.util.Iterator;
import java.util.stream.StreamSupport;

/**
 * Shard PacBio sub-reads uBAM by ZMW (for running CCS algo.).
 *
 * Now that the output BAMs will NOT be the same as the number of ZMWs, as that will be too many files.
 * We cluster several ZMWs together in each output BAM.
 *
 * <h3>Usage example</h3>
 * <pre>
 *   gatk ShardPacBioSubReadsUBamByZMWClusterSpark \
 *     -I input_reads.bam \
 *     --read-index input_reads.bam.sbi \
 *     -O output_prefix_
 * </pre>
 */
@DocumentedFeature
@BetaFeature
@CommandLineProgramProperties(
        oneLineSummary = "Shard PacBio sub-reads uBAM by ZMW",
        summary =
                "Shard PacBio sub-reads uBAM by ZMW.",
        programGroup = LongReadAnalysisProgramGroup.class)
public class ShardPacBioSubReadsUBamByZMWClusterSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;


    private static final String ZMW_ATRRIBUTE_KEY = "zm";


    @Argument(doc = "the output prefix", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    protected String output;

    @Advanced
    @Argument(doc = "count of ZMWs in one shard output bam", fullName = "shard-size", optional = true)
    protected Integer shardSize = 20_000;

    @Advanced
    @Hidden
    @Argument(doc = "turn on debug mode (much longer runtime, and IO/memory burden", fullName = "debug", optional = true)
    protected Boolean debugMode = false;

    @Override
    public boolean requiresReads() {
        return true;
    }

    @Override
    protected void runTool( final JavaSparkContext ctx ) {

        // these redundant local variables are here just to avoid the whole tool having to be serialized...
        final int local_shard_size = shardSize;
        final String output_prefix = output;

        final SAMFileHeader headerForReads = getHeaderForReads();
        headerForReads.setSortOrder(SAMFileHeader.SortOrder.queryname); // each output will be sorted by read name

        if (! debugMode) {
            getUnfilteredReads()
                    .groupBy(read -> read.getAttributeAsInteger(ZMW_ATRRIBUTE_KEY) / local_shard_size) // shard index
                    .foreach( shard -> {
                        final Iterator<SAMRecord> readsInThisCluster =
                                StreamSupport.stream(shard._2.spliterator(), false)
                                        .map(read -> read.convertToSAMRecord(headerForReads))
                                        .sorted(Comparator.comparing(SAMRecord::getReadName))
                                        .iterator();

                        final Integer shardIdx = shard._1;
                        try ( SAMFileWriter writer = new SAMFileWriterFactory().setCreateIndex(false)
                                .makeBAMWriter(headerForReads, true, IOUtils.getPath(output_prefix + shardIdx + ".bam"))) {
                            readsInThisCluster.forEachRemaining(writer::addAlignment);
                        } catch ( final UncheckedIOException ie) {
                            throw new GATKException("Can't write BAM file for shard: " + shardIdx, ie);
                        }
                    });
        } else {
            final JavaPairRDD<Integer, Iterable<GATKRead>> cache = getUnfilteredReads()
                    .groupBy(read -> read.getAttributeAsInteger(ZMW_ATRRIBUTE_KEY) / local_shard_size).cache(); // shard index
            logger.info(String.format("Expect %d output BAMs.", cache.count()));

            cache.foreach( shard -> {
                final Iterator<SAMRecord> readsInThisCluster =
                        StreamSupport.stream(shard._2.spliterator(), false)
                                .map(read -> read.convertToSAMRecord(headerForReads))
                                .sorted(Comparator.comparing(SAMRecord::getReadName))
                                .iterator();

                final Integer shardIdx = shard._1;
                try ( SAMFileWriter writer = new SAMFileWriterFactory().setCreateIndex(false)
                        .makeBAMWriter(headerForReads, true, IOUtils.getPath(output_prefix + shardIdx + ".bam"))) {
                    readsInThisCluster.forEachRemaining(writer::addAlignment);
                } catch ( final UncheckedIOException ie) {
                    throw new GATKException("Can't write BAM file for shard: " + shardIdx, ie);
                }
            });
        }
    }
}
