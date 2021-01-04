package org.broadinstitute.hellbender.tools.spark.sv;

import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;

import java.io.Serializable;

public class DiscoverVariantsFromReadDepthArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    public static final int DEFAULT_MIN_EVENT_SIZE = 50;
    public static final int DEFAULT_INSERT_SIZE = 350;
    public static final int DEFAULT_PLOIDY = 2;
    public static final int DEFAULT_MIN_GRAPH_EVIDENCE_READ_PAIRS = 2;
    public static final double DEFAULT_GRAPH_PARTITION_RECIPROCAL_OVERLAP = 0.01;
    public static final int DEFAULT_MAX_GRAPH_BRANCHES = 1000;
    public static final double DEFAULT_MIN_EVENT_PROBABILITY = 0.95;
    public static final double DEFAULT_MAX_GRAPH_PATH_LENGTH_FACTOR = 10.0;
    public static final int DEFAULT_MAX_GRAPH_EDGE_VISITS = 3;
    public static final double DEFAULT_MIN_HAPLOTYPE_OUTPUT_PROBABILITY = 0.1;

    public static final String MIN_EVENT_SIZE_LONG_NAME = "min-size";
    public static final String INSERT_SIZE_LONG_NAME = "insert-size";
    public static final String PLOIDY_LONG_NAME = "ploidy";
    public static final String MIN_EVIDENCE_READ_PAIRS_LONG_NAME = "min-evidence-read-pairs";
    public static final String GRAPH_PARTITION_RECIPROCAL_LONG_NAME = "partition-ro";
    public static final String MAX_GRAPH_BRANCHES_LONG_NAME = "max-branches";
    public static final String MIN_EVENT_PROBABILITY_LONG_NAME = "min-event-prob";
    public static final String MAX_GRAPH_PATH_LENGTH_FACTOR_LONG_NAME = "max-path-length-factor";
    public static final String MAX_GRAPH_EDGE_VISITS_LONG_NAME = "max-edge-visits";
    public static final String MIN_HAPLOTYPE_OUTPUT_PROBABILITY_LONG_NAME = "min-haplotype-prob";

    @Argument(doc = "Min event size (bp)",
            fullName = MIN_EVENT_SIZE_LONG_NAME, optional = true)
    public int minEventSize = DEFAULT_MIN_EVENT_SIZE;

    @Argument(doc = "Library insert size (bp)",
            fullName = INSERT_SIZE_LONG_NAME, optional = true)
    public int insertSize = DEFAULT_INSERT_SIZE;

    @Argument(doc = "Autosomal ploidy",
            fullName = PLOIDY_LONG_NAME, optional = true)
    public int ploidy = DEFAULT_PLOIDY;

    @Argument(doc = "Min haplotype probability for output",
            fullName = MIN_HAPLOTYPE_OUTPUT_PROBABILITY_LONG_NAME, optional = true)
    public double minHaplotypeProb = DEFAULT_MIN_HAPLOTYPE_OUTPUT_PROBABILITY;

    @Argument(doc = "Min called event probability",
            fullName = MIN_EVENT_PROBABILITY_LONG_NAME, optional = true)
    public double minEventProb = DEFAULT_MIN_EVENT_PROBABILITY;

    @Advanced
    @Argument(doc = "Min number pairs for read pair evidence clusters",
            fullName = MIN_EVIDENCE_READ_PAIRS_LONG_NAME, optional = true)
    public int minLinkReadPairs = DEFAULT_MIN_GRAPH_EVIDENCE_READ_PAIRS;

    @Advanced
    @Argument(doc = "Min reciprocal overlap for fully overlapping subgraphs to be merged",
            fullName = GRAPH_PARTITION_RECIPROCAL_LONG_NAME, optional = true)
    public double paritionReciprocalOverlap = DEFAULT_GRAPH_PARTITION_RECIPROCAL_OVERLAP;

    @Advanced
    @Argument(doc = "Max number of queued path branches to hold in memory",
            fullName = MAX_GRAPH_BRANCHES_LONG_NAME, optional = true)
    public int maxBranches = DEFAULT_MAX_GRAPH_BRANCHES;

    @Advanced
    @Argument(doc = "Max path length as a multiple of the number of reference edges",
            fullName = MAX_GRAPH_PATH_LENGTH_FACTOR_LONG_NAME, optional = true)
    public double maxPathLengthFactor = DEFAULT_MAX_GRAPH_PATH_LENGTH_FACTOR;

    @Advanced
    @Argument(doc = "Max number of visits allowed per edge",
            fullName = MAX_GRAPH_EDGE_VISITS_LONG_NAME, optional = true)
    public int maxEdgeVisits = DEFAULT_MAX_GRAPH_EDGE_VISITS;

}
