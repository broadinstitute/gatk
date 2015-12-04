package org.broadinstitute.hellbender.tools.picard.analysis.directed;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.IntervalList;
import org.broadinstitute.hellbender.metrics.MetricAccumulationLevel;

import java.io.File;
import java.util.*;

/**
 * Calculates HS metrics for a given SAM/BAM file. Requires the input of a list of
 * target intervals and a list of bait intervals. Can be invoked either on an entire
 * iterator of SAMRecords or be passed SAMRecords one at a time.
 *
 * @author Jonathan Burke
 */
public final class TargetedPcrMetricsCollector extends TargetMetricsCollector<TargetedPcrMetrics> {
    //maybe instead just inject this into the TargetedMetricCollector ->

    public TargetedPcrMetricsCollector(final Set<MetricAccumulationLevel> accumulationLevels, final List<SAMReadGroupRecord> samRgRecords, final ReferenceSequenceFile refFile, final File perTargetCoverage, final IntervalList targetIntervals, final IntervalList probeIntervals, final String probeSetName) {
        super(accumulationLevels, samRgRecords, refFile, perTargetCoverage, targetIntervals, probeIntervals, probeSetName);
    }

    @Override
    public TargetedPcrMetrics convertMetric(TargetMetrics targetMetrics) {
        final TargetedPcrMetrics pcrMetrics = new TargetedPcrMetrics();
        TargetMetricsCollector.reflectiveCopy(targetMetrics, pcrMetrics,
                new String[]{"PROBE_SET",           "PROBE_TERRITORY",    "ON_PROBE_BASES",    "NEAR_PROBE_BASES",    "OFF_PROBE_BASES",    "PCT_SELECTED_BASES",  "PCT_OFF_PROBE",    "ON_PROBE_VS_SELECTED",    "MEAN_PROBE_COVERAGE"},
                new String[]{"CUSTOM_AMPLICON_SET", "AMPLICON_TERRITORY", "ON_AMPLICON_BASES", "NEAR_AMPLICON_BASES", "OFF_AMPLICON_BASES", "PCT_AMPLIFIED_BASES", "PCT_OFF_AMPLICON", "ON_AMPLICON_VS_SELECTED", "MEAN_AMPLICON_COVERAGE"}
        );

        return pcrMetrics;
    }
}
