package org.broadinstitute.hellbender.metrics;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamPairUtil;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.IOUtil;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.AuthHolder;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.utils.R.RScriptExecutor;
import org.broadinstitute.hellbender.utils.R.RScriptExecutorException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.File;
import java.io.Serializable;

/**
 * Collects InsertSizeMetrics on the specified accumulationLevels
 */
public final class InsertSizeMetricsCollector
        extends MultiLevelReducibleCollector<InsertSizeMetrics, Integer, InsertSizeMetricsCollectorArgs, PerUnitInsertSizeMetricsCollector>
        implements Serializable
{
    private static final long serialVersionUID = 1L;
    private static final Logger log = LogManager.getLogger(InsertSizeMetricsCollector.class);

    InsertSizeMetricsArgumentCollection inputArgs = null;

    /**
     * @param inputArgs InsertSizeMetricsArgumentCollection populated with argument values. May not be null.
     * @param samHeader samHeader for the input to be processed. May not be null.
     */
    public InsertSizeMetricsCollector(
            final InsertSizeMetricsArgumentCollection inputArgs,
            final SAMFileHeader samHeader)
    {
        Utils.nonNull(inputArgs);
        Utils.nonNull(samHeader);

        this.inputArgs = inputArgs;
        setup(inputArgs.metricAccumulationLevel, samHeader.getReadGroups());
    }

    /**
     * Return the read filter for InsertSizeMetrics collector.
     * @param samHeader SAMFileHeader for the input file. May not be null
     * @return ReadFilter to be used to filter records
     */
    public ReadFilter getReadFilter(SAMFileHeader samHeader) {
        Utils.nonNull(samHeader);
        return new WellformedReadFilter(samHeader).and(
                new InsertSizeMetricsReadFilter(
                    inputArgs.useEnd,
                    inputArgs.filterNonProperlyPairedReads,
                    !inputArgs.useDuplicateReads,
                    !inputArgs.useSecondaryAlignments,
                    !inputArgs.useSupplementaryAlignments,
                    inputArgs.MQPassingThreshold
            )
        );
    }

    // We will pass insertSize and PairOrientation with the DefaultPerRecordCollectorArgs passed to
    // the record collectors. This method is called once per samRecord
    @Override
    protected InsertSizeMetricsCollectorArgs makeArg(SAMRecord samRecord, ReferenceSequence refSeq) {
        // inferred insert size is negative if the mate maps to lower position than the read, so use abs
        final int insertSize = Math.abs(samRecord.getInferredInsertSize());
        final SamPairUtil.PairOrientation orientation = SamPairUtil.getPairOrientation(samRecord);

        return new InsertSizeMetricsCollectorArgs(insertSize, orientation);
    }

    // Make an PerUnitInsertSizeMetricsCollector with the given arguments
    @Override
    protected PerUnitInsertSizeMetricsCollector makeChildCollector(
            final String sample,
            final String library,
            final String readGroup) {
        return new PerUnitInsertSizeMetricsCollector(
                sample,
                library,
                readGroup,
                inputArgs.minimumPct,
                inputArgs.maxMADTolerance,
                inputArgs.histogramWidth);
    }

    /**
     * Combine two InsertSizeMetricsCollector objects and return a single InsertSizeMetricsCollector
     * object representing the combined results. NOTE: this implementation is destructive in that it
     * merges the source into the target and returns the target as the combined object.
     * @param target target destination of combined metrics. May not be null.
     * @param source source of metrics to be combined into target. May not be null.
     * @return single object representing the combined source and target objects
     */
    public InsertSizeMetricsCollector combine(InsertSizeMetricsCollector target, InsertSizeMetricsCollector source) {
        Utils.nonNull(target);
        Utils.nonNull(source);
        target.combine(source);
        return target;
    }

    /**
     * Combine two PerUnitInsertSizeMetricsCollector objects and return a single PerUnitInsertSizeMetricsCollector
     * object representing the combined results.
     * @param collector1 source PerUnitInsertSizeMetricsCollector. May not be null.
     * @param collector2 target PerUnitInsertSizeMetricsCollector. May not be null.
     * @return single PerUnitInsertSizeMetricsCollector object representing the combined source and target objects
     */
    public PerUnitInsertSizeMetricsCollector combineUnit(
            PerUnitInsertSizeMetricsCollector collector1,
            PerUnitInsertSizeMetricsCollector collector2) {
        Utils.nonNull(collector1);
        Utils.nonNull(collector2);
        return collector1.combine(collector2);
    }

    /**
     * Finish the metrics collection by saving any results to a metrics file.
     * @param metricsFile a metricsFile where the collected metrics should be stored. May not be null.
     * @param inputName the name of the input, for optional inclusion in the metrics file. May not be null.
     * @param authHolder Authentication info for this context.
     */
    public void finish(
            final MetricsFile<InsertSizeMetrics, Integer> metricsFile,
            final String inputName,
            final AuthHolder authHolder)
    {
        Utils.nonNull(metricsFile);

        finish();
        addAllLevelsToFile(metricsFile);

        if (metricsFile.getNumHistograms() == 0) { // can happen if user sets MINIMUM_PCT = 0.5, etc.
            log.warn("All data categories were discarded because they contained < " + inputArgs.minimumPct +
                     " of the total aligned paired data.");
            final PerUnitInsertSizeMetricsCollector allReadsCollector = getAllReadsCollector();
            log.warn("Total mapped pairs in all categories: " + (allReadsCollector == null ?
                    allReadsCollector :
                    allReadsCollector.getTotalInserts()));
        }
        else {
            if (null != authHolder) {
                MetricsUtils.saveMetrics(metricsFile, inputArgs.output, authHolder);
            }
            else {
                metricsFile.write(new File(inputArgs.output));
            }
            if (inputArgs.producePlot) {
                writeHistogramPDF(inputName);
            }
        }
    }

    /**
     * Calls R script to plot histogram(s) in PDF.
     */
    private void writeHistogramPDF(final String inputName) throws RScriptExecutorException {
        // path to Picard R script for producing histograms in PDF files.
        final String R_SCRIPT = "insertSizeHistogram.R";

        File histFile = new File(inputArgs.histogramPlotFile);
        IOUtil.assertFileIsWritable(histFile);

        final RScriptExecutor executor = new RScriptExecutor();
        executor.addScript(new Resource(R_SCRIPT, InsertSizeMetricsCollector.class));
        executor.addArgs(
                inputArgs.output,           // text-based metrics file
                histFile.getAbsolutePath(), // PDF graphics file
                inputName                   // input bam file
        );
        if (inputArgs.histogramWidth != null) {
            executor.addArgs(String.valueOf(inputArgs.histogramWidth));
        }
        executor.exec();
    }

    /**
     * Custom read filter, paramaterized by InsertSizeMetricsArgumentCollection provided
     */
    private static final class InsertSizeMetricsReadFilter implements ReadFilter, Serializable {
        private static final long serialVersionUID = 1L;

        private final ReadFilter combinedReadFilter;

        public InsertSizeMetricsReadFilter(final InsertSizeMetricsArgumentCollection.EndToUse whichEnd,
                                           final boolean filterNonProperlyPairedReads,
                                           final boolean filterDuplicatedReads,
                                           final boolean filterSecondaryAlignments,
                                           final boolean filterSupplementaryAlignments,
                                           final int MQThreshold){

            final InsertSizeMetricsArgumentCollection.EndToUse endVal = whichEnd;

            ReadFilter tempFilter = ReadFilterLibrary.MAPPED;
            tempFilter = tempFilter.and(GATKRead::isPaired);
            tempFilter = tempFilter.and(read -> 0!=read.getFragmentLength());
            tempFilter = tempFilter.and(read -> endVal == (read.isFirstOfPair() ?
                    InsertSizeMetricsArgumentCollection.EndToUse.FIRST :
                    InsertSizeMetricsArgumentCollection.EndToUse.SECOND));

            if(filterNonProperlyPairedReads)  { tempFilter = tempFilter.and(GATKRead::isProperlyPaired); }
            if(filterDuplicatedReads)         { tempFilter = tempFilter.and(read -> !read.isDuplicate()); }
            if(filterSecondaryAlignments)     { tempFilter = tempFilter.and(read -> !read.isSecondaryAlignment()); }
            if(filterSupplementaryAlignments) { tempFilter = tempFilter.and(read -> !read.isSupplementaryAlignment()); }

            if(0!=MQThreshold)  { tempFilter = tempFilter.and(read -> read.getMappingQuality() >= MQThreshold);}

            combinedReadFilter = tempFilter;
        }

        @Override
        public boolean test(final GATKRead read){
            return combinedReadFilter.test(read);
        }
    }

}
