package org.broadinstitute.hellbender.tools.spark.pipelines.metrics;

import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.metrics.MetricAccumulationLevel;
import org.broadinstitute.hellbender.metrics.MetricsArgumentCollection;

import java.util.EnumSet;
import java.util.Set;

// TODO: filter reads based on only isReverseStrand/mateIsReverseStrand (strand bias)
// TODO: filter reads based on {MATE_ON_SAME_CONTIG, MATE_DIFFERENT_STRAND, GOOD_CIGAR, NON_ZERO_REFERENCE_LENGTH_ALIGNMENT}
// TODO: filter reads based on length value (if too large), and/or minimum_pct like in Picard.
// TODO: case EITHER for enum EndToUser. For truncated genomic regions of putative SV breakpoints, not all reads have
//       both ends land in the region, so third case is possible: will use either end when only one end is available in
//       the region specified, and only first end if both are available.
// TODO: user argument validation (eg. maxMADTolerance)

/**
 * MetricsArgumentCollection argument collection for InsertSize metrics. All members should be
 * instantiable as command line arguments.
 */
public class InsertSizeMetricsArgumentCollection extends MetricsArgumentCollection {
    @Argument(doc = "A local path to PDF file where histogram plot will be saved in.",
            shortName = "HIST",
            fullName = "HistogramPlotPDF",
            optional = false)
    public String histogramPlotFile = null;

    @Argument(doc = "Generate mean, sd and plots by trimming the data down to MEDIAN + maxMADTolerance*MEDIAN_ABSOLUTE_DEVIATION. " +
            "This is done because insert size data typically includes enough anomalous values from chimeras and other " +
            "artifacts to make the mean and sd grossly misleading regarding the real distribution.",
            shortName = "TOL",
            fullName = "HistogramPlotDeviationsTolerance",
            optional = true)
    public double maxMADTolerance = 10.0;

    // read filtering criteria
    @Argument(doc = "If set to true, filter pairs of reads that are not properly--as judged by aligner--oriented.",
            shortName = "PP",
            fullName = "filterNonProperlyPairedReads",
            optional = true)
    public boolean filterNonProperlyPairedReads = false;

    @Argument(doc = "If set to true, include duplicated reads as well.",
            shortName = "Dup",
            fullName = "useDuplicateReads",
            optional = true)
    public boolean useDuplicateReads = false;

    @Argument(doc = "If set to true, include secondary alignments.",
            shortName = "S",
            fullName = "useSecondaryAlignments",
            optional = true)
    public boolean useSecondaryAlignments = false;

    @Argument(doc = "If set to true, include supplementary alignments.",
            shortName = "SS",
            fullName = "useSupplementaryAlignments",
            optional = true)
    public boolean useSupplementaryAlignments = false;

    @Argument(doc = "If set non-zero value, only include reads passing certain mapping quality threshold. " +
            "If set to zero, reads with zero mapping quality will be included in calculating metrics.",
            shortName = "MAPQ",
            fullName = "MAPQThreshold",
            optional = true)
    public int MQPassingThreshold = 0;

    @Argument(doc="The level(s) at which to accumulate metrics. Possible values are {ALL_READS, SAMPLE, LIBRARY, READ GROUP}.",
            shortName="LEVEL",
            fullName = "MetricsAccumulationLevel",
            optional = false)
    public Set<MetricAccumulationLevel> metricAccumulationLevel = EnumSet.of(MetricAccumulationLevel.ALL_READS);

    /**
     * Which end of a read pair to use for collecting insert size metrics.
     */
    public enum EndToUse {
        FIRST(1), SECOND(2);
        private final int value;
        EndToUse(int value){
            this.value = value;
        }
        public int getValue(){
            return value;
        }
    }

    @Argument(doc = "Which end of pairs to use for collecting information. " +
            "Possible values:{FIRST, SECOND}.",
            shortName = "E",
            fullName = "whichEndOfPairToUse",
            optional = true)
    public EndToUse useEnd = EndToUse.FIRST;
}
