package org.broadinstitute.hellbender.tools.spark.pathseq;

import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;

import java.io.Serializable;

public final class PSScoreArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    public static final String SCORES_OUTPUT_LONG_NAME = "scores-output";
    public static final String SCORES_OUTPUT_SHORT_NAME = "SO";
    public static final String TAXONOMIC_DATABASE_LONG_NAME = "taxonomy-file";
    public static final String TAXONOMIC_DATABASE_SHORT_NAME = "T";
    public static final String SCORE_METRICS_FILE_LONG_NAME = "score-metrics";
    public static final String SCORE_METRICS_FILE_SHORT_NAME = "SM";
    public static final String SCORE_WARNINGS_FILE_LONG_NAME = "score-warnings";
    public static final String SCORE_WARNINGS_FILE_SHORT_NAME = "SW";
    public static final String MIN_SCORE_IDENTITY_LONG_NAME = "min-score-identity";
    public static final String MIN_SCORE_IDENTITY_SHORT_NAME = MIN_SCORE_IDENTITY_LONG_NAME;
    public static final String IDENTITY_MARGIN_LONG_NAME = "identity-margin";
    public static final String IDENTITY_MARGIN_SHORT_NAME = IDENTITY_MARGIN_LONG_NAME;
    public static final String DIVIDE_BY_GENOME_LENGTH_LONG_NAME = "divide-by-genome-length";
    public static final String DIVIDE_BY_GENOME_LENGTH_SHORT_NAME = DIVIDE_BY_GENOME_LENGTH_LONG_NAME;
    public static final String NOT_NORMALIZED_BY_KINGDOM_LONG_NAME = "not-normalized-by-kingdom";
    public static final String NOT_NORMALIZED_BY_KINGDOM_SHORT_NAME = NOT_NORMALIZED_BY_KINGDOM_LONG_NAME;
    public static final String SCORE_READS_PER_PARTITION_LONG_NAME = "score-reads-per-partition-estimate";
    public static final String SCORE_READS_PER_PARTITION_SHORT_NAME = SCORE_READS_PER_PARTITION_LONG_NAME;

    @Argument(doc = "URI for the taxonomic scores output",
            fullName = SCORES_OUTPUT_LONG_NAME,
            shortName = SCORES_OUTPUT_SHORT_NAME)
    public String scoresPath;

    @Argument(doc = "URI to the microbe reference taxonomy database built using PathSeqBuildReferenceTaxonomy",
            fullName = TAXONOMIC_DATABASE_LONG_NAME,
            shortName = TAXONOMIC_DATABASE_SHORT_NAME)
    public String taxonomyDatabasePath;

    /**
     * This parameter controls the stringency of the microbe alignment. The identity score threshold is defined as the
     * number of matching bases minus number of deletions. Alignments below this threshold score will be ignored.
     */
    @Argument(doc = "Alignment identity score threshold, as a fraction of the read length (between 0 and 1).",
            fullName = MIN_SCORE_IDENTITY_LONG_NAME,
            shortName = MIN_SCORE_IDENTITY_SHORT_NAME,
            minValue = 0.0,
            maxValue = 1.0,
            optional = true)
    public double minIdentity = 0.90;

    /**
     * For reads having multiple alignments, the best hit is always counted as long as it is above the identity score
     * threshold. Any additional hits will be counted when its identity score is within this percentage of the best hit.
     *
     * <p>For example, consider a read that aligns to two different sequences, one with identity score 0.90 and the other with
     * 0.85. If the minimum identity score is 0.7, the best hit (with score 0.90) is counted. In addition, if the identity margin is 10%,
     * then any additional alignments at or above 0.90 * (1 - 0.10) = 0.81 would also be counted. Therefore in this example the second
     * alignment with score 0.85 would be counted.</p>
     */
    @Argument(doc = "Identity margin, as a fraction of the best hit (between 0 and 1). ",
            fullName = IDENTITY_MARGIN_LONG_NAME,
            shortName = IDENTITY_MARGIN_SHORT_NAME,
            minValue = 0.0,
            maxValue = 1.0,
            optional = true)
    public double identityMargin = 0.02;

    /**
     * If true, the score contributed by each read is divided by the mapped organism's genome length in the reference.
     */
    @Argument(doc = "Divide abundance scores by each taxon's reference genome length (in millions)",
            fullName = DIVIDE_BY_GENOME_LENGTH_LONG_NAME,
            shortName = DIVIDE_BY_GENOME_LENGTH_SHORT_NAME,
            optional = true)
    public boolean divideByGenomeLength = false;

    /**
     * Comparmentalizes the normalized abundance scores by kingdom.
     */
    @Argument(doc = "If true, normalized abundance scores will be reported as a percentage within their kingdom.",
            fullName = NOT_NORMALIZED_BY_KINGDOM_LONG_NAME,
            shortName = NOT_NORMALIZED_BY_KINGDOM_SHORT_NAME,
            optional = true)
    public boolean notNormalizedByKingdom = false;

    @Argument(doc = "Write accessions found in the reads header but not the taxonomy database to this file",
            fullName = SCORE_WARNINGS_FILE_LONG_NAME,
            shortName = SCORE_WARNINGS_FILE_SHORT_NAME,
            optional = true)
    public String headerWarningFile = null;

    /**
     * This parameter is for fine-tuning memory performance. Lower values may result in less memory usage but possibly
     * at the expense of greater computation time.
     */
    @Advanced
    @Argument(doc = "Estimated reads per Spark partition for scoring",
            fullName = SCORE_READS_PER_PARTITION_LONG_NAME,
            shortName = SCORE_READS_PER_PARTITION_SHORT_NAME,
            minValue = 1,
            optional = true)
    public int readsPerPartitionEstimate = 200000;

    /**
     * If specified, records the following metrics:
     * <ul>
     *     <li>Number of reads mapped to the microbial reference</li>
     *     <li>Number of unmapped reads</li>
     * </ul>
     *<p>Note that using this option may increase runtime.</p>
     */
    @Argument(doc = "Log counts of mapped and unmapped reads to this file",
            fullName = SCORE_METRICS_FILE_LONG_NAME,
            shortName = SCORE_METRICS_FILE_SHORT_NAME,
            optional = true)
    public String scoreMetricsFileUri = null;

}