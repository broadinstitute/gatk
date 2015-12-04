package org.broadinstitute.hellbender.tools.picard.analysis.directed;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.IntervalList;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.QCProgramGroup;
import org.broadinstitute.hellbender.metrics.MetricAccumulationLevel;

import java.io.File;
import java.util.*;

/**
 * Collect metric information for target pcr metrics runs.  See CollectTargetedMetrics and TargetPcrMetricsCollector for
 * more information
 */
@CommandLineProgramProperties(
        summary = "Calculates a set of metrics to Illumina Truseq Custom Amplicon sequencing from an aligned SAM/BAM" +
                "file. If a reference sequence is provided, AT/GC dropout metrics will " +
                "be calculated, and the PER_TARGET_COVERAGE option can be used to output GC and " +
                "mean coverage information for every target.",
        oneLineSummary = "Produces Targeted PCR-related metrics given the provided SAM/BAM",
        programGroup = QCProgramGroup.class
)
public final class CollectTargetedPcrMetrics extends CollectTargetedMetrics<TargetedPcrMetrics, TargetedPcrMetricsCollector> {

    @Argument(shortName = "AI", doc = "An interval list file that contains the locations of the baits used.")
    public File AMPLICON_INTERVALS;

    @Argument(shortName = "N", doc = "Custom amplicon set name. If not provided it is inferred from the filename of the AMPLICON_INTERVALS intervals.", optional = true)
    public String CUSTOM_AMPLICON_SET_NAME;

    /**
     * @return AMPLICON_INTERVALS
     */
    @Override
    protected IntervalList getProbeIntervals() {
        return IntervalList.fromFile(AMPLICON_INTERVALS);
    }

    /**
     * @return CUSTOM_AMPLICON_SET_NAME
     */
    @Override
    protected String getProbeSetName() {
        return CUSTOM_AMPLICON_SET_NAME != null ? CUSTOM_AMPLICON_SET_NAME : CollectTargetedMetrics.renderProbeNameFromFile(AMPLICON_INTERVALS);
    }

    @Override
    protected TargetedPcrMetricsCollector makeCollector(final Set<MetricAccumulationLevel> accumulationLevels,
                                                        final List<SAMReadGroupRecord> samRgRecords,
                                                        final ReferenceSequenceFile refFile,
                                                        final File perTargetCoverage,
                                                        final IntervalList targetIntervals,
                                                        final IntervalList probeIntervals,
                                                        final String probeSetName) {
        return new TargetedPcrMetricsCollector(accumulationLevels, samRgRecords, refFile, perTargetCoverage, targetIntervals, probeIntervals, probeSetName);
    }
}
