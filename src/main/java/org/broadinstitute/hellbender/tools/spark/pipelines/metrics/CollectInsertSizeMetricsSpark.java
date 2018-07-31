package org.broadinstitute.hellbender.tools.spark.pipelines.metrics;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.metrics.Header;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.metrics.InsertSizeMetricsArgumentCollection;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.List;

/**
 * Collects insert size distribution information in alignment data. The tool leverages the Spark framework
 * for faster operation.
 *
 * <h3>Usage example</h3>
 * <pre>
 * gatk CollectInsertSizeMetricsSpark \
 *   -I gs://cloud-bucket/input.bam \
 *   -H gs://cloud-bucket/insert_size_histogram.pdf \
 *   -O gs://cloud-bucket/insert_size_metrics.txt \
 *   -- \
 *   --spark-runner GCS \
 *   --cluster my-dataproc-cluster
 * </pre>
 * <p>
 * See <a href=http://broadinstitute.github.io/picard/picard-metric-definitions.html#InsertSizeMetrics>
 *     http://broadinstitute.github.io/picard/picard-metric-definitions.html#InsertSizeMetrics</a>
 * for an explanation of individual metrics. See <a href ="https://software.broadinstitute.org/gatk/documentation/article?id=10060">
 *     Tutorial#10060</a> for an example of how to set up and run a Spark tool on a cloud Spark cluster.
 * </p>
 */
@CommandLineProgramProperties(
        summary = "Collects insert size distribution information in SAM/BAM/CRAM file(s). The tool leverages " +
                "the Spark framework for faster operation.",
        oneLineSummary = "Collects insert size distribution information on alignment data",
        programGroup = DiagnosticsAndQCProgramGroup.class)
@DocumentedFeature
@BetaFeature
public final class CollectInsertSizeMetricsSpark
        extends MetricsCollectorSparkTool<InsertSizeMetricsArgumentCollection> {

    private static final long serialVersionUID = 1L;

    @ArgumentCollection
    private InsertSizeMetricsArgumentCollection insertSizeArgs = new InsertSizeMetricsArgumentCollection();

    private InsertSizeMetricsCollectorSpark insertSizeCollector = new InsertSizeMetricsCollectorSpark();

    @Override
    public InsertSizeMetricsArgumentCollection getInputArguments() {
        return insertSizeArgs;
    }

    @Override
    protected SortOrder getExpectedSortOrder() { return insertSizeCollector.getExpectedSortOrder(); }

    @Override
    protected void initialize(
            final InsertSizeMetricsArgumentCollection inputArgs,
            final SAMFileHeader samHeader,
            final List<Header> defaultHeaders)
    {
        insertSizeCollector.initialize(inputArgs, samHeader, defaultHeaders);
    }

    /**
     * Return the read filters required for this collector
     */
    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return insertSizeCollector.getDefaultReadFilters();
    }

    @Override
    protected void collectMetrics(
            final JavaRDD<GATKRead> filteredReads,
            final SAMFileHeader samHeader)
    {
        insertSizeCollector.collectMetrics(filteredReads, samHeader);
    }

    @Override
    protected void saveMetrics(final String inputName) {
        insertSizeCollector.saveMetrics(inputName);
    }
}
