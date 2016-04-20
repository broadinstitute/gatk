package org.broadinstitute.hellbender.tools.spark.pipelines.metrics;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.metrics.Header;
import htsjdk.samtools.metrics.MetricsFile;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.engine.AuthHolder;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.metrics.MetricsUtils;
import org.broadinstitute.hellbender.metrics.QualityYieldMetrics;
import org.broadinstitute.hellbender.metrics.QualityYieldMetricsArgumentCollection;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.Serializable;
import java.util.List;

/**
 * QualityYieldMetricsCollector for Spark.
 */
public class QualityYieldMetricsCollectorSpark
        implements MetricsCollectorSpark<QualityYieldMetricsArgumentCollection>, Serializable {

    private static final long serialVersionUID = 1L;

    private QualityYieldMetricsArgumentCollection args = null;

    private MetricsFile<QualityYieldMetrics, Integer> metricsFile = null;

    /**
     * Initialize the collector with input arguments;
     */
    @Override
    public void initialize(final QualityYieldMetricsArgumentCollection inputArgs, final List<Header> defaultHeaders) {
        metricsFile = new MetricsFile<QualityYieldMetrics, Integer>();
        defaultHeaders.stream().forEach(h -> metricsFile.addHeader(h));
        this.args = inputArgs;
    }

    /**
     * Do the actual metrics collection on the provided RDD.
     * @param filteredReads The reads to be analyzed for this collector.
     * @param samHeader The SAMFileHeader associated with the reads in the input RDD.
     * @param authHolder authHolder
     * @param inputBaseName base name of the input file
     */
    public void collectMetrics(
            final JavaRDD<GATKRead> filteredReads,
            final SAMFileHeader samHeader,
            final String inputBaseName,
            final AuthHolder authHolder
    )
    {
        final QualityYieldMetrics metrics =
                filteredReads.aggregate(new QualityYieldMetrics().setUseOriginalQualities(args.useOriginalQualities),
                        (hgp, read) -> hgp.addRead(read),
                        (hgp1, hgp2) -> hgp1.merge(hgp2))
                        .finish();

        metricsFile.addMetric(metrics);

        // this could be moved to finishCollection, but it needs the autholder/basename
        MetricsUtils.saveMetrics(metricsFile, args.output.getAbsolutePath(), authHolder);
    }

    public void finishCollection() {}

}
