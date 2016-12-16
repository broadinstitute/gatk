package org.broadinstitute.hellbender.tools.spark.pathseq;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceRecord;
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
import scala.Tuple2;

import java.io.IOException;
import java.util.Map;

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
 * 6) abundance score normalized as a percentage of total pathogen mapped reads
 * 7) total number of reads aligned
 * 8) number of reads alignmented assigned uniquely to this taxon (i.e. mapped only to the node and/or its children)
 * 9) reference sequence length
 * <p>
 * For unpaired reads, a particular alignment is considered a hit if it aligned with minimal coverage (read length
 * less padding) given as a fraction of the full read length and identity (matches less deletions) given as a
 * fraction of the read with padded bases removed.
 * For paired-end reads, a particular alignment is considered a hit only if both reads aligned to a reference with
 * the same taxonomic id with minimal coverage and identity.
 */

@CommandLineProgramProperties(summary = "Performs taxonomic classification of BAM file reads that have been aligned to a reference of pathogens",
        oneLineSummary = "PathSeq read classification tool",
        programGroup = SparkProgramGroup.class)
public class PathSeqClassifyReadsSpark extends GATKSparkTool {

    private static final long serialVersionUID = 1L;

    @Argument(doc = "URI to paired reads bam",
            fullName = "pairedInput",
            optional = true)
    public String INPUT_PAIRED = null;

    @Argument(doc = "URI to unpaired reads bam",
            fullName = "unpairedInput",
            optional = true)
    public String INPUT_UNPAIRED = null;

    @Argument(doc = "URI for the taxonomic scores output",
            fullName = "scoresOutputPath",
            optional = false)
    public String SCORES_PATH;

    @Argument(doc = "URI to the reference taxonomy database build using PathSeqBuildReferenceTaxonomy",
            fullName = "taxonomicDatabasePath",
            optional = false)
    public String TAXDB_PATH;

    @Argument(doc = "URI to the output BAM",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = true)
    public String OUTPUT_PATH = null;

    @Argument(doc = "Alignment coverage threshold, as a fraction of the read length (between 0 and 1)",
            fullName = "coverageThreshold",
            minValue = 0.0,
            maxValue = 1.0,
            optional = true)
    public double MIN_COVERAGE = 0.90;

    @Argument(doc = "Alignment identity threshold, as a fraction of the read length (between 0 and 1)",
            fullName = "identityThreshold",
            minValue = 0.0,
            maxValue = 1.0,
            optional = true)
    public double MIN_IDENTITY = 0.90;

    @Argument(doc = "Write accessions found in the reads header but not the taxonomy database to this file",
            fullName = "warningsFile",
            optional = true)
    public String HEADER_WARNING_FILE = null;

    int recommendedNumReducers = 0;

    private Tuple2<JavaRDD<GATKRead>, SAMFileHeader> readInputWithHeader(final String path, final PipelineOptions options, final ReadsSparkSource readsSource) {
        if (path != null) {
            if (BucketUtils.fileExists(path, options)) {
                recommendedNumReducers += PSUtils.pathseqGetRecommendedNumReducers(path, numReducers, options, getTargetPartitionSize());
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

        final PipelineOptions options = getAuthenticatedGCSOptions();
        final ReadsSparkSource readsSource = new ReadsSparkSource(ctx, readArguments.getReadValidationStringency());

        //Load reads
        final Tuple2<JavaRDD<GATKRead>, SAMFileHeader> pairedData = readInputWithHeader(INPUT_PAIRED, options, readsSource);
        final Tuple2<JavaRDD<GATKRead>, SAMFileHeader> unpairedData = readInputWithHeader(INPUT_UNPAIRED, options, readsSource);

        final JavaRDD<GATKRead> pairedReads = pairedData._1;
        final SAMFileHeader pairedHeader = pairedData._2;
        final JavaRDD<GATKRead> unpairedReads = unpairedData._1;
        final SAMFileHeader unpairedHeader = unpairedData._2;

        //Group reads into pairs
        JavaRDD<Iterable<GATKRead>> groupedReads;
        if (pairedReads != null) {
            groupedReads = PSUtils.groupReadPairs(pairedReads);
            if (unpairedReads != null) {
                groupedReads = groupedReads.union(unpairedReads.map(PSClassifyReadsUtils::wrapInIterable));
            }
        } else if (unpairedReads != null) {
            groupedReads = unpairedReads.map(PSClassifyReadsUtils::wrapInIterable);
        } else {
            throw new UserException.BadInput("No reads were loaded. Ensure --pairedInput and/or --unpairedInput are set and valid.");
        }

        //Join headers
        SAMFileHeader header;
        if (pairedHeader != null) {
            header = pairedHeader;
            if (unpairedHeader != null && !header.equals(unpairedHeader)) {
                for (final SAMSequenceRecord rec : unpairedHeader.getSequenceDictionary().getSequences()) {
                    if (header.getSequenceDictionary().getSequence(rec.getSequenceName()) == null) {
                        header.addSequence(rec);
                    }
                }
            }
        } else if (unpairedHeader != null) {
            header = unpairedHeader;
        } else {
            throw new UserException.BadInput("No headers were loaded");
        }

        //Load taxonomy database, created by running PathSeqBuildReferenceTaxonomy with this reference
        logger.info("Loading taxonomy database...");
        final PSTaxonomyDatabase taxDB = PSUtils.readTaxonomyDatabase(TAXDB_PATH, options);
        final Broadcast<Map<String, String>> contigToTaxBroadcast = ctx.broadcast(taxDB.contigToTaxIDMap);

        //Check header against database
        if (HEADER_WARNING_FILE != null) {
            logger.info("Checking header against database...");
            PSClassifyReadsUtils.doHeaderWarnings(HEADER_WARNING_FILE, getHeaderForReads(), taxDB, options, logger);
        }

        //Determine which alignments are hits and return their tax IDs
        final JavaRDD<Tuple2<Iterable<GATKRead>, PSHitInfo>> readHits = PSClassifyReadsUtils.mapGroupedReadsToTax(groupedReads, MIN_COVERAGE, MIN_IDENTITY, contigToTaxBroadcast);

        //Collect hit info objects
        logger.info("Collecting hits...");
        final Iterable<PSHitInfo> hitInfo = PSClassifyReadsUtils.collectTupleSecond(readHits);

        //Flatten the original reads, now with their hit tags set
        final JavaRDD<GATKRead> readsFinal = PSClassifyReadsUtils.flattenTupleFirstIterable(readHits);

        //Compute taxonomic scores
        logger.info("Computing scores...");
        final Map<String, PSScoreInfo> taxScores = PSClassifyReadsUtils.computeTaxScores(hitInfo, taxDB.tree);

        //Write scores to file
        logger.info("Writing output...");
        PSClassifyReadsUtils.writeScoresFile(taxScores, taxDB.tree, SCORES_PATH, options);

        //Write reads to BAM, if specified
        if (OUTPUT_PATH != null) {
            try {
                ReadsSparkSink.writeReads(ctx, OUTPUT_PATH, null, readsFinal, header,
                        shardedOutput ? ReadsWriteFormat.SHARDED : ReadsWriteFormat.SINGLE, recommendedNumReducers);
            } catch (final IOException e) {
                throw new UserException.CouldNotCreateOutputFile(OUTPUT_PATH, "writing failed", e);
            }
        }
    }


}
