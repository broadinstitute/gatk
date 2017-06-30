package org.broadinstitute.hellbender.tools.spark.pathseq;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.PathSeqProgramGroup;
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
 * step of the PathSeq pipeline for abundance quantification. The output is a tab-delimited table with the
 * following columns:
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
 * less padding) given as a fraction of the full read length, and identity score (matches less deletions) given as a
 * fraction of the read with padded bases removed.
 * For paired-end reads, an alignment to a particular taxon is considered a hit only if both reads aligned to that
 * taxon with minimal coverage and identity.
 */

@CommandLineProgramProperties(summary = "Performs taxonomic classification of reads that have been aligned to a pathogen " +
        "reference and generates a table of abundance scores. Briefly, an alignment is considered a 'hit' if it meets " +
        "minimum query coverage and identity criteria. For paired reads, hits that appear in only one mate and not the " +
        "other are discarded. Reads with a single hit add a score of 1 to the corresponding species or strain. For reads with " +
        "N hits, a score of 1/N is distributed to each. Scores are then normalized by genome length in megabases and added " +
        "to all taxonomic ancestors." +
        "\n\n" +
        "Users must supply the paired and/or unpaired read BAM(s) (do not use --input) and the taxonomic database " +
        "created using PathSeqBuildReferenceTaxonomy." +
        "\n\n" +
        "The output is a score table giving the abundance score, normalized score, number of hits, number of " +
        "unambiguous hits, and genome length of each taxon. Optionally, a BAM is also produced that attaches a " +
        PSScorer.HITS_TAG + " tag to each read listing its hit taxonomic IDs.",
        oneLineSummary = "Step 3: Classifies pathogen-aligned reads and generates abundance scores",
        programGroup = PathSeqProgramGroup.class)
@BetaFeature
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

    /**
     * Combines SAM file header sequences and read groups. If both headers are not null, additional header entries
     * from the unpaired header will not be copied over.
     */
    static SAMFileHeader joinBamHeaders(final SAMFileHeader pairedHeader, final SAMFileHeader unpairedHeader) {
        SAMFileHeader header;
        if (pairedHeader != null) {
            header = pairedHeader;
            if (unpairedHeader != null && !header.equals(unpairedHeader)) {
                //Add sequences
                for (final SAMSequenceRecord rec : unpairedHeader.getSequenceDictionary().getSequences()) {
                    if (header.getSequenceDictionary().getSequence(rec.getSequenceName()) == null) {
                        header.addSequence(rec);
                    }
                }
                //Add read groups
                for (final SAMReadGroupRecord rec : unpairedHeader.getReadGroups()) {
                    if (header.getReadGroup(rec.getReadGroupId()) == null) {
                        header.addReadGroup(rec);
                    }
                }
            }
        } else if (unpairedHeader != null) {
            header = unpairedHeader;
        } else {
            throw new UserException.BadInput("No headers were loaded");
        }
        return header;
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
        final SAMFileHeader header = joinBamHeaders(pairedHeader, unpairedHeader);

        //Main tool routine
        final PSScorer scorer = new PSScorer(scoreArgs);
        final JavaRDD<GATKRead> readsFinal = scorer.scoreReads(ctx, pairedReads, unpairedReads, header);

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
