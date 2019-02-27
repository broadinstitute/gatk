package org.broadinstitute.hellbender.tools.examples.metrics.multi;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.metrics.*;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.io.Serializable;
import java.util.Collections;
import java.util.List;

/**
 * Example multi-level metrics collector for illustrating how to collect metrics on
 * specified accumulation levels.
 */
public final class ExampleMultiMetricsCollector
        extends MultiLevelReducibleCollector<
                    ExampleMultiMetrics,
                    Integer,
                    PerUnitExampleMultiMetricsCollectorArgs,
                    PerUnitExampleMultiMetricsCollector>
        implements Serializable
{
    private static final long serialVersionUID = 1L;
    private static final Logger log = LogManager.getLogger(ExampleMultiMetricsCollector.class);

    ExampleMultiMetricsArgumentCollection inputArgs = null;

    public ExampleMultiMetricsCollector(){}

    /**
     * @param inputArgs ExampleMultiMetricsArgumentCollection populated with argument values. May not be null.
     * @param samHeader samHeader for the input to be processed. May not be null.
     */
    public void initialize(
            final ExampleMultiMetricsArgumentCollection inputArgs,
            final SAMFileHeader samHeader)
    {
        Utils.nonNull(inputArgs);
        Utils.nonNull(samHeader);

        this.inputArgs = inputArgs;

        // call setup to create the level distributors
        setup(inputArgs.metricAccumulationLevel.accumulationLevels, samHeader.getReadGroups());
    }

    /**
     * Return the read filter for example metrics collector.
     * @return List of read filters to be used to filter records
     */
    public List<ReadFilter> getDefaultReadFilters() {
        return Collections.singletonList(new ReadFilterLibrary.AllowAllReadsReadFilter());
    }

    /**
     * Process the record and create a collector args object to hold the results of any
     * calculation for re-use across multiple levels.
     */
    @Override
    protected PerUnitExampleMultiMetricsCollectorArgs makeArg(
            final SAMRecord samRecord,
            final ReferenceSequence refSeq) {
        // Example collector doesn't use any cached values from the record so its a no-op.
        return new PerUnitExampleMultiMetricsCollectorArgs();
    }

    /**
     * Make an PerUnitExampleMultiMetricsCollector with the given arguments
     */
    @Override
    protected PerUnitExampleMultiMetricsCollector makeChildCollector(
            final String sample,
            final String library,
            final String readGroup) {
        return new PerUnitExampleMultiMetricsCollector(
                sample,
                library,
                readGroup);
    }

    /**
     * Combine two ExampleMultiMetricsCollector objects and return a single ExampleMultiMetricsCollector
     * object representing the combined results. NOTE: this implementation is destructive in that it
     * merges the source into the target and returns the target as the combined object.
     * @param target target destination of combined metrics. May not be null.
     * @param source source of metrics to be combined into target. May not be null.
     * @return single object representing the combined source and target objects
     */
    protected static ExampleMultiMetricsCollector combine(
            final ExampleMultiMetricsCollector target,
            final ExampleMultiMetricsCollector source) {
        Utils.nonNull(target);
        Utils.nonNull(source);
        target.combine(source);
        return target;
    }

    /**
     * Combine two PerUnitExampleMultiMetricsCollector objects and return a single PerUnitExampleMultiMetricsCollector
     * object representing the combined results.
     * @param collector1 source PerUnitExampleMultiMetricsCollector. May not be null.
     * @param collector2 target PerUnitExampleMultiMetricsCollector. May not be null.
     * @return single PerUnitExampleMultiMetricsCollector object representing the combined source and target objects
     */
    @Override
    public PerUnitExampleMultiMetricsCollector combineUnit(
            final PerUnitExampleMultiMetricsCollector collector1,
            final PerUnitExampleMultiMetricsCollector collector2) {
        Utils.nonNull(collector1);
        Utils.nonNull(collector2);
        return collector1.combine(collector2);
    }

    /**
     * Finish the metrics collection and save any results to the metrics file.
     * @param metricsFile a metricsFile where the collected metrics should be stored. May not be null.
     */
    public void saveMetrics(
            final MetricsFile<ExampleMultiMetrics, Integer> metricsFile)
    {
        Utils.nonNull(metricsFile);

        finish();

        addAllLevelsToFile(metricsFile);
        MetricsUtils.saveMetrics(metricsFile, inputArgs.output);
    }

}
