package org.broadinstitute.hellbender.tools.spark.pipelines.metrics;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.metrics.Header;

import htsjdk.samtools.metrics.MetricsFile;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.engine.AuthHolder;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.metrics.*;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.Serializable;
import java.util.*;

/**
 * Worker class to collect insert size metrics, add metrics to file, and provides
 * accessors to stats of groups of different level.
 */
public class InsertSizeMetricsCollectorSpark
        implements MetricsCollectorSpark<InsertSizeMetricsArgumentCollection>,
        Serializable
{
    private static final long serialVersionUID = 1L;

    private InsertSizeMetricsArgumentCollection inputArgs = null;
    private InsertSizeMetricsCollector collector = null;

    private InsertSizeMetricsCollector resultMetrics = null;
    private MetricsFile<InsertSizeMetrics, Integer> metricsFile = null;

    /**
     * Initialize the collector with input arguments;
     */
    @Override
    public void initialize(
            final InsertSizeMetricsArgumentCollection inputArgs,
            final SAMFileHeader samHeader,
            final List<Header> defaultHeaders)
    {
        this.inputArgs = inputArgs;
        collector = new InsertSizeMetricsCollector(inputArgs, samHeader);
        metricsFile = new MetricsFile<>();
        if (defaultHeaders != null) {
            defaultHeaders.stream().forEach(h -> metricsFile.addHeader(h));
        }
    }

    /**
     * Return a read filter to be used for this collector.
     *
     * @param samFileHeader
     * @return ReadFilter
     */
    @Override
    public ReadFilter getReadFilter(final SAMFileHeader samFileHeader) {
        return collector.getReadFilter(samFileHeader);
    }

    @Override
    public void collectMetrics(
            final JavaRDD<GATKRead> filteredReads,
            final SAMFileHeader samHeader)
    {
        if (filteredReads.isEmpty()) {
            throw new GATKException("No valid reads found in input file.");
        }

        resultMetrics = filteredReads
                .mapPartitions(
                    it -> {
                        it.forEachRemaining(r -> collector.acceptRecord(r.convertToSAMRecord(samHeader), null));
                        return Arrays.asList(collector);
                    })
                .reduce(collector::combine);
    }

    public void saveMetrics(final String inputName, final AuthHolder authHolder) {
        resultMetrics.finish(metricsFile, inputName, authHolder);
    }

}