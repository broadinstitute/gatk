package org.broadinstitute.hellbender.tools.spark.pathseq;

import htsjdk.samtools.SAMFileHeader;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
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
import scala.Tuple2;

import java.io.IOException;

/**
 * This Spark tool classifies reads taxonomically, assigning each taxon abundance scores. This is the last
 *   step of the PathSeq pipeline for abundance quantification. The output is a tab-delimited table with the
 *   following columns:
 * <p>
 * 1) taxonomic id
 * 2) taxonomic classification
 * 3) taxonomic type (order, family, genus, etc.)
 * 4) name
 * 5) raw abundance score:
 * <p>
 * score of species S = sum_{read r in reads aligned to S} 1 / ({# of species that r aligned to} x {length of species reference in Mbp})
 * <p>
 * score of higher taxa = sum of scores of child nodes
 * <p>
 * 6) abundance score normalized as a percentage of total pathogen-mapped reads
 * 7) total number of reads aligned to this taxon
 * 8) number of reads assigned uniquely to this taxon (i.e. mapped only to the node and/or its children)
 * 9) total taxon reference sequence length
 * <p>
 * For unpaired reads, a particular alignment is considered a hit if it aligned with minimal coverage score (read length
 *   less padding) given as a fraction of the full read length, and identity score (matches less deletions) given as a
 *   fraction of the read with padded bases removed.
 * For paired-end reads, an alignment to a particular taxon is considered a hit only if both reads aligned to that
 *   taxon with minimal coverage and identity.
 */

@CommandLineProgramProperties(summary = "Performs taxonomic classification of BAM file reads that have been aligned to a reference of pathogens",
        oneLineSummary = "PathSeq read classification tool",
        programGroup = SparkProgramGroup.class)
public class PathSeqScoreSpark extends GATKSparkTool {

    private static final long serialVersionUID = 1L;

    @Argument(doc = "URI to paired reads BAM",
            fullName = "pairedInput",
            optional = true)
    public String pairedInput = null;

    @Argument(doc = "URI to unpaired reads BAM",
            fullName = "unpairedInput",
            optional = true)
    public String unpairedInput = null;

    @Argument(doc = "URI to the output BAM",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            optional = true)
    public String outputPath = null;

    @ArgumentCollection
    public PSScoreArgumentCollection scoreArgs = new PSScoreArgumentCollection();

    private int recommendedNumReducers = 0;

    private Tuple2<JavaRDD<GATKRead>, SAMFileHeader> readInputWithHeader(final String path,
                                                                         final ReadsSparkSource readsSource) {
        if (path != null) {
            if (BucketUtils.fileExists(path)) {
                recommendedNumReducers += PSUtils.pathseqGetRecommendedNumReducers(path, numReducers, getTargetPartitionSize());
                final SAMFileHeader header = readsSource.getHeader(path, null);
                JavaRDD<GATKRead> reads = readsSource.getParallelReads(path, null, null, bamPartitionSplitSize);
                reads = PSUtils.primaryReads(reads);
                return new Tuple2<>(reads, header);
            } else {
                logger.warn("Could not find file " + path + ". Skipping...");
            }
        }
        return new Tuple2<>(null, null);
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {

        if (!readArguments.getReadFiles().isEmpty()) {
            throw new UserException.BadInput("Please use --pairedInput or --unpairedInput instead of --input");
        }

        final ReadsSparkSource readsSource = new ReadsSparkSource(ctx, readArguments.getReadValidationStringency());

        //Load reads
        final Tuple2<JavaRDD<GATKRead>, SAMFileHeader> pairedData = readInputWithHeader(pairedInput, readsSource);
        final Tuple2<JavaRDD<GATKRead>, SAMFileHeader> unpairedData = readInputWithHeader(unpairedInput, readsSource);

        final JavaRDD<GATKRead> pairedReads = pairedData._1;
        final SAMFileHeader pairedHeader = pairedData._2;
        final JavaRDD<GATKRead> unpairedReads = unpairedData._1;
        final SAMFileHeader unpairedHeader = unpairedData._2;

        if (pairedHeader != null && !pairedHeader.getSortOrder().equals(SAMFileHeader.SortOrder.queryname)) {
            throw new UserException.BadInput("Paired input BAM must be sorted by queryname");
        }

        //Join header sequences and read groups
        final SAMFileHeader header = PSScoreUtils.joinBamHeaders(pairedHeader, unpairedHeader);

        //Main tool routine
        final JavaRDD<GATKRead> readsFinal = PSScoreUtils.scoreReads(ctx, pairedReads, unpairedReads, header, scoreArgs);

        //Write reads to BAM, if specified
        //Note writeReads() is not used because we determine recommendedNumReducers differently with 2 input BAMs
        if (outputPath != null) {
            try {
                ReadsSparkSink.writeReads(ctx, outputPath, null, readsFinal, header,
                        shardedOutput ? ReadsWriteFormat.SHARDED : ReadsWriteFormat.SINGLE, recommendedNumReducers);
            } catch (final IOException e) {
                throw new UserException.CouldNotCreateOutputFile(outputPath, "writing failed", e);
            }
        }
    }


}
