package org.broadinstitute.hellbender.tools.spark.pathseq;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.MetagenomicsProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSink;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.pathseq.loggers.PSScoreFileLogger;
import org.broadinstitute.hellbender.tools.spark.pathseq.loggers.PSScoreLogger;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadsWriteFormat;
import scala.Tuple2;

import java.io.IOException;

/**
 * Classify reads and estimate abundances of each taxon in the reference. This is the third and final step of the PathSeq pipeline.
 *
 * <p>See PathSeqPipelineSpark for an overview of the PathSeq pipeline.</p>
 *
 * <p>This tool performs taxonomic classification of reads that have been aligned to a microbe reference. Using
 * alignments that are sufficiently identical to the reference, it scores each taxon in the reference based on aligned
 * reads. These scores can be used to detect and quantify microbe abundance.</p>
 *
 * <h3>Methods</h3>
 *
 * <p>Alignments with sufficient identity score (e.g. 90% of read length) are used to estimate read counts
 * and the relative abundance of microorganisms present in the sample at each level of the taxonomic tree (e.g. strain,
 * species, genus, family, etc.). If a read maps to more than one organism, only the best alignment and any others with
 * identity score within a margin of error of the best (e.g. 2%) are retained. For paired-end reads, alignments to organisms present in one read but not the other are
 * discarded. Reads with a single valid alignment add a score of 1 to the corresponding species or strain. For reads with N hits, a score of 1/N is distributed
 * to each organism. Scores are totaled for each taxon by summing the scores across all reads and the scores of any descendent taxa.</p>
 *
 * <h3>Input</h3>
 *
 * <ul>
 *     <li>Queryname-sorted BAM file containing only paired reads aligned to the microbe reference</li>
 *     <li>BAM file containing only unpaired reads aligned to the microbe reference</li>
 *     <li>*Taxonomy file generated using PathSeqBuildReferenceTaxonomy</li>
 * </ul>
 *
 * <p>*A standard microbe taxonomy file is available in the <a href="https://software.broadinstitute.org/gatk/download/bundle">GATK Resource Bundle</a>.</p>
 *
 * <h3>Output</h3>
 *
 * <p>Tab-delimited scores table with the following columns:</p>
 *
 * <ul>
 *      <li>NCBI taxonomic ID</li>
 *      <li>phylogenetic classification</li>
 *      <li>phylogenetic rank (order, family, genus, etc.)</li>
 *      <li>taxon name</li>
 *      <li>abundance score (described above)</li>
 *      <li>abundance score normalized as a percentage of total microbe-mapped reads</li>
 *      <li>total number of reads aligned to this taxon</li>
 *      <li>number of reads assigned unambiguously to the taxon (i.e. mapped only to the node and/or its children)</li>
 *      <li>total taxon reference sequence length</li>
 * </ul>
 *
 * <p>The tool may also optionally produce:</p>
 *
 * <ul>
 *     <li>BAM file of all reads annotated with the NCBI taxonomy IDs of mapped organisms</li>
 *     <li>Metrics file with the number of mapped and unmapped reads</li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <p>This tool can be run without explicitly specifying Spark options. That is to say, the given example command
 * without Spark options will run locally. See
 * <a href ="https://software.broadinstitute.org/gatk/documentation/article?id=10060">Tutorial#10060</a> for an example
 * of how to set up and run a Spark tool on a cloud Spark cluster.</p>
 *
 * <pre>
 * gatk PathSeqScoreSpark  \
 *   --paired-input input_reads_paired.bam \
 *   --unpaired-input input_reads_unpaired.bam \
 *   --taxonomy-file taxonomy.db \
 *   --scores-output scores.txt \
 *   --output output_reads.bam \
 *   --min-score-identity 0.90 \
 *   --identity-margin 0.02
 * </pre>
 *
 * @author Mark Walker &lt;markw@broadinstitute.org&gt;
 */

@CommandLineProgramProperties(summary = "Classify reads and estimate abundances of each taxon in the reference. This is the third and final step of the PathSeq pipeline.",
        oneLineSummary = "Step 3: Classifies pathogen-aligned reads and generates abundance scores",
        programGroup = MetagenomicsProgramGroup.class)
@DocumentedFeature
public class PathSeqScoreSpark extends GATKSparkTool {

    private static final long serialVersionUID = 1L;

    public static final String PAIRED_INPUT_LONG_NAME = "paired-input";
    public static final String PAIRED_INPUT_SHORT_NAME = "PI";
    public static final String UNPAIRED_INPUT_LONG_NAME = "unpaired-input";
    public static final String UNPAIRED_INPUT_SHORT_NAME = "UI";

    @Argument(doc = "Input queryname-sorted BAM containing only paired reads",
            fullName = PAIRED_INPUT_LONG_NAME,
            shortName = PAIRED_INPUT_SHORT_NAME,
            optional = true)
    public String pairedInput = null;

    @Argument(doc = "Input BAM containing only unpaired reads",
            fullName = UNPAIRED_INPUT_LONG_NAME,
            shortName = UNPAIRED_INPUT_SHORT_NAME,
            optional = true)
    public String unpairedInput = null;

    /**
     * Records have a "YP" tag that lists the NCBI taxonomy IDs of any mapped organisms. This tag is omitted if the
     * read is unmapped.
     */
    @Argument(doc = "Output BAM",
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
            throw new UserException.BadInput("Please use --paired-input or --unpaired-input instead of --input");
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
        if (scoreArgs.scoreMetricsFileUri != null) {
            try (final PSScoreLogger scoreLogger = new PSScoreFileLogger(getMetricsFile(), scoreArgs.scoreMetricsFileUri)) {
                scoreLogger.logReadCounts(readsFinal);
            }
        }

        //Write reads to BAM, if specified
        //Note writeReads() is not used because we determine recommendedNumReducers differently with 2 input BAMs
        if (outputPath != null) {
            try {
                ReadsSparkSink.writeReads(ctx, outputPath, null, readsFinal, header,
                        shardedOutput ? ReadsWriteFormat.SHARDED : ReadsWriteFormat.SINGLE, recommendedNumReducers, shardedPartsDir);
            } catch (final IOException e) {
                throw new UserException.CouldNotCreateOutputFile(outputPath, "writing failed", e);
            }
        }
    }


}
