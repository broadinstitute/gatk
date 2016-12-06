package org.broadinstitute.hellbender.metrics;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.MetricAccumulationLevelArgumentCollection;

import java.io.Serializable;

// TODO: filter reads based on only isReverseStrand/mateIsReverseStrand (strand bias)
// TODO: filter reads based on {MATE_ON_SAME_CONTIG, MATE_DIFFERENT_STRAND, GOOD_CIGAR, NON_ZERO_REFERENCE_LENGTH_ALIGNMENT}
// TODO: filter reads based on length value (if too large), and/or minimum_pct like in Picard.
// TODO: case EITHER for enum EndToUser. For truncated genomic regions of putative SV breakpoints, not all reads have
//       both ends land in the region, so third case is possible: will use either end when only one end is available in
//       the region specified, and only first end if both are available.
// TODO: user argument validation (eg. maxMADTolerance)

/**
 * ArgumentCollection for InsertSizeMetrics collectors.
 */
public class InsertSizeMetricsArgumentCollection extends MetricsArgumentCollection implements Serializable {

    private static final long serialVersionUID = 1L;

    @Argument(fullName = "histogramPlotFile",
            shortName="H",
            doc="File to write insert size Histogram chart to.")
    public String histogramPlotFile;

    @Argument(doc="Generate mean, sd and plots by trimming the data down to MEDIAN + maxMADTolerance*MEDIAN_ABSOLUTE_DEVIATION. " +
            "This is done because insert size data typically includes enough anomalous values from chimeras and other " +
            "artifacts to make the mean and sd grossly misleading regarding the real distribution.",
            shortName = "TOL",
            fullName = "HistogramPlotDeviationsTolerance",
            optional = true)
    public double maxMADTolerance = 10.0;

    @Argument(shortName="W", doc="Explicitly sets the Histogram width, overriding automatic truncation of Histogram tail. " +
            "Also, when calculating mean and standard deviation, only bins <= HISTOGRAM_WIDTH will be included.", optional=true)
    public Integer histogramWidth = null;

    @Argument(shortName="M", doc="When generating the Histogram, discard any data categories (out of FR, TANDEM, RF) that have fewer than this " +
            "percentage of overall reads. (Range: 0 to 1).")
    public float minimumPct = 0.05f;

    @Argument(doc = "Should an output plot be created")
    public boolean producePlot = false;

    @ArgumentCollection
    public MetricAccumulationLevelArgumentCollection metricAccumulationLevel = new MetricAccumulationLevelArgumentCollection();

}
