package org.broadinstitute.hellbender.tools.picard.analysis.directed;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.PicardCommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.metrics.MetricAccumulationLevel;
import org.broadinstitute.hellbender.metrics.MultiLevelMetrics;
import org.broadinstitute.hellbender.utils.runtime.ProgressLogger;

import java.io.File;
import java.util.EnumSet;
import java.util.List;
import java.util.Set;

/**
 * Both CollectTargetedPCRMetrics and CalculateHybridSelection metrics share virtually identical program structures except
 * for the name of their targeting mechanisms (e.g. bait set or amplicon set).  The shared behavior of these programs
 * is encapsulated in CollectTargetedMetrics which is then subclassed by CalculateHsMetrics and CollectTargetedPcrMetrics.
 * <p/>
 * This program verifies the input parameters to TargetMetricsCollector and converts all files to
 * the format desired by TargetMetricsCollector.  Then it instantiates a TargetMetricsCollector and
 * collects metric information for all reads in the INPUT sam file.
 */
public abstract class CollectTargetedMetrics<METRIC extends MultiLevelMetrics, COLLECTOR extends TargetMetricsCollector<METRIC>> extends PicardCommandLineProgram {

    protected abstract IntervalList getProbeIntervals();

    protected abstract String getProbeSetName();

    /**
     * A factory method for the TargetMetricsCollector to use this time.  Examples of TargetMetricsCollector:
     * (TargetedPcrMetricsCollector, HsMetricsCalculator)
     *
     * @return A TargetMetricsCollector to which we will pass SAMRecords
     */
    protected abstract COLLECTOR makeCollector(final Set<MetricAccumulationLevel> accumulationLevels,
                                               final List<SAMReadGroupRecord> samRgRecords,
                                               final ReferenceSequenceFile refFile,
                                               final File perTargetCoverage,
                                               final IntervalList targetIntervals,
                                               final IntervalList probeIntervals,
                                               final String probeSetName);

    @Argument(shortName = "TI", doc = "An interval list file that contains the locations of the targets.")
    public List<File> TARGET_INTERVALS;

    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME, shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            doc = "An aligned SAM/BAM file.")
    public File INPUT;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "The output file to write the metrics to.")
    public File OUTPUT;

    @Argument(shortName = "LEVEL", doc = "The level(s) at which to accumulate metrics.")
    public Set<MetricAccumulationLevel> METRIC_ACCUMULATION_LEVEL = EnumSet.of(MetricAccumulationLevel.ALL_READS);

    @Argument(optional = true, doc = "An optional file to output per target coverage information to.")
    public File PER_TARGET_COVERAGE;

    /**
     * Asserts that files are readable and writable and then fires off an
     * HsMetricsCalculator instance to do the real work.
     */
    protected Object doWork() {
        for (final File targetInterval : TARGET_INTERVALS) IOUtil.assertFileIsReadable(targetInterval);
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        if (PER_TARGET_COVERAGE != null) IOUtil.assertFileIsWritable(PER_TARGET_COVERAGE);

        final SamReader reader = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT);
        final IntervalList targetIntervals = IntervalList.fromFiles(TARGET_INTERVALS);

        // Validate that the targets and baits have the same references as the reads file
        SequenceUtil.assertSequenceDictionariesEqual(
                reader.getFileHeader().getSequenceDictionary(),
                targetIntervals.getHeader().getSequenceDictionary());
        SequenceUtil.assertSequenceDictionariesEqual(
                reader.getFileHeader().getSequenceDictionary(),
                getProbeIntervals().getHeader().getSequenceDictionary()
        );

        ReferenceSequenceFile ref = null;
        if (REFERENCE_SEQUENCE != null) {
            IOUtil.assertFileIsReadable(REFERENCE_SEQUENCE);
            ref = ReferenceSequenceFileFactory.getReferenceSequenceFile(REFERENCE_SEQUENCE);
            SequenceUtil.assertSequenceDictionariesEqual(
                    reader.getFileHeader().getSequenceDictionary(), ref.getSequenceDictionary(),
                    INPUT, REFERENCE_SEQUENCE
            );
        }

        final COLLECTOR collector = makeCollector(
                METRIC_ACCUMULATION_LEVEL,
                reader.getFileHeader().getReadGroups(),
                ref,
                PER_TARGET_COVERAGE,
                targetIntervals,
                getProbeIntervals(),
                getProbeSetName()
        );

        final ProgressLogger progress = new ProgressLogger(logger);
        for (final SAMRecord record : reader) {
            collector.acceptRecord(record, null);
            progress.record(record);
        }

        // Write the output file
        final MetricsFile<METRIC, Integer> metrics = getMetricsFile();
        collector.finish();

        collector.addAllLevelsToFile(metrics);

        metrics.write(OUTPUT);

        CloserUtil.close(reader);
        return null;
    }

    /** Renders a probe name from the provided file, returning {@link File#getName()} with all extensions stripped. */
    static String renderProbeNameFromFile(final File probeIntervalFile) {
        final String name = probeIntervalFile.getName();
        final int firstPeriodIndex = name.indexOf('.');
        if (firstPeriodIndex == -1) {
            return name;
        } else {
            return name.substring(0, firstPeriodIndex);
        }
    }

    protected String[] customCommandLineValidation() {
        if (PER_TARGET_COVERAGE != null && (METRIC_ACCUMULATION_LEVEL.size() != 1 ||
                METRIC_ACCUMULATION_LEVEL.iterator().next() != MetricAccumulationLevel.ALL_READS)) {
            return new String[]{"PER_TARGET_COVERAGE can be specified only when METRIC_ACCUMULATION_LEVEL is set " +
                    "to ALL_READS."};
        }

        if (PER_TARGET_COVERAGE != null && REFERENCE_SEQUENCE == null) {
            return new String[]{"Must supply REFERENCE_SEQUENCE when supplying PER_TARGET_COVERAGE"};
        }

        return super.customCommandLineValidation();
    }
}
