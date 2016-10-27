package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.coveragestats.CalculateCoverageStats;
import org.broadinstitute.hellbender.tools.exome.coveragestats.SampleCoverageStats;
import org.broadinstitute.hellbender.tools.exome.coveragestats.SampleCoverageStatsReader;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;

/**
 * Tool to filter samples based on coverage statistics.
 * <p>
 *     This tool uses the output from {@link CalculateCoverageStats} in order
 *     to discard sample with extreme coverage mean or variance.
 * </p>
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Filter samples based on their target coverage statistics",
        oneLineSummary = "Filter samples",
        programGroup = CopyNumberProgramGroup.class
)
public class FilterSamples extends CommandLineProgram {

    public static final String REJECTION_SAMPLE_NAME_COLUMN = "SAMPLE";
    public static final String REJECTION_FILTER_NAME_COLUMN = "FILTER";
    public static final String REJECTION_REASON_COLUMN = "REASON";
    public static final String OUTPUT_NAME_COLUMN = "NAME";

    public static final String REJECTION_FILE_FULL_NAME = "rejectOutputFile";
    public static final String REJECTION_FILE_SHORT_NAME = "rejected";
    public static final String MINIMUM_MEAN_COVERAGE_FULL_NAME = "minimumCoverageMean";
    public static final String MINIMUM_MEAN_COVERAGE_SHORT_NAME = "minCovMean";
    public static final String MAXIMUM_MEAN_COVERAGE_FULL_NAME = "maximumCoverageMean";
    public static final String MAXIMUM_MEAN_COVERAGE_SHORT_NAME = "maxCovMean";
    public static final String MINIMUM_COVERAGE_VARIANCE_FULL_NAME = "minimumCoverageVariance";
    public static final String MINIMUM_COVERAGE_VARIANCE_SHORT_NAME = "minCovVar";
    public static final String MAXIMUM_COVERAGE_VARIANCE_FULL_NAME = "maximumCoverageVariance";
    public static final String MAXIMUM_COVERAGE_VARIANCE_SHORT_NAME = "maxCovVar";

    public static final double DEFAULT_MINIMUM_MEAN_COVERAGE = Double.NEGATIVE_INFINITY;
    public static final double DEFAULT_MAXIMUM_MEAN_COVERAGE = Double.POSITIVE_INFINITY;
    public static final double DEFAULT_MINIMUM_COVERAGE_VARIANCE = 0.0;
    public static final double DEFAULT_MAXIMUM_COVERAGE_VARIANCE = Double.POSITIVE_INFINITY;

    @Argument(
        doc = "Input sample statistics.",
        fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
        shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
        optional = false
    )
    protected File inputFile;

    @Argument(
        doc = "Output rejection file.",
        fullName = REJECTION_FILE_FULL_NAME,
        shortName = REJECTION_FILE_SHORT_NAME,
        optional = true
    )
    protected File rejectionFile;

    @Argument(
        doc = "Output file of samples to keep.",
        fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
        shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
        optional = true
    )
    protected File outputFile;

    @Argument(
        doc = "Minimum mean coverage; samples with a lower mean are discarded.",
        fullName = MINIMUM_MEAN_COVERAGE_FULL_NAME,
        shortName = MINIMUM_MEAN_COVERAGE_SHORT_NAME,
        optional = true
    )
    protected double minimumMeanCoverage = DEFAULT_MINIMUM_MEAN_COVERAGE;

    @Argument(
            doc = "Maximum mean coverage; samples with a lower mean are discarded.",
            fullName = MAXIMUM_MEAN_COVERAGE_FULL_NAME,
            shortName = MAXIMUM_MEAN_COVERAGE_SHORT_NAME,
            optional = true
    )
    protected double maximumMeanCoverage = DEFAULT_MAXIMUM_MEAN_COVERAGE;

    @Argument(
            doc = "Minimum mean coverage; samples with a lower mean are discarded.",
            fullName = MINIMUM_COVERAGE_VARIANCE_FULL_NAME,
            shortName = MINIMUM_COVERAGE_VARIANCE_SHORT_NAME,
            optional = true
    )
    protected double minimumCoverageVariance = DEFAULT_MINIMUM_COVERAGE_VARIANCE;

    @Argument(
            doc = "Maximum mean coverage; samples with a lower mean are discarded.",
            fullName = MAXIMUM_COVERAGE_VARIANCE_FULL_NAME,
            shortName = MAXIMUM_COVERAGE_VARIANCE_SHORT_NAME,
            optional = true
    )
    protected double maximumCoverageVariance = DEFAULT_MAXIMUM_COVERAGE_VARIANCE;

    public enum SampleFilter {
        ExtremeCoverageMean,
        ExtremeCoverageVariance
    }

    /**
     * Compose the list of filter predicates to apply to input samples
     * based on the values given to user arguments.
     * @return never {@code null} but perhaps empty.
     */
    private List<SampleFilterPredicate> composeFilterPredicateList() {
        final List<SampleFilterPredicate> result = new ArrayList<>(2);
        if (maximumMeanCoverage != Double.POSITIVE_INFINITY || minimumMeanCoverage != Double.NEGATIVE_INFINITY) {
            result.add(new SampleMeanCoverageFilterPredicate());
        }
        if (maximumCoverageVariance != Double.POSITIVE_INFINITY || minimumCoverageVariance > 0.0) {
            result.add(new SampleCoverageVarianceFilterPredicate());
        }
        return result;
    }

    @Override
    protected Object doWork() {
        final List<SampleFilterPredicate> filters = composeFilterPredicateList();
        try(final SampleCoverageStatsReader reader = new SampleCoverageStatsReader(inputFile);
            final SampleRejectionRecordWriter rejectionWriter = new SampleRejectionRecordWriter(rejectionFile);
            final SampleListWriter outputWriter = new SampleListWriter(outputFile)) {
            for (final SampleCoverageStats stats : reader) {
                final List<SampleFilterPredicate> failingFilters = filters.stream()
                        .filter(f -> !f.test(stats)).collect(Collectors.toList());
                if (failingFilters.isEmpty()) {
                    outputWriter.writeSampleName(stats.sample);
                } else {
                    failingFilters.stream().forEach(f -> rejectionWriter.writeRecord(
                            new SampleRejectionRecord(stats.sample, f.getFilter(), f.rejectReason(stats))));
                }
            }
            return "SUCCESS";
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile(inputFile, ex);
        }
    }

    /**
     * Common interface for all sample filter predicates.
     */
    interface SampleFilterPredicate extends Predicate<SampleCoverageStats> {
        String rejectReason(final SampleCoverageStats stats);
        SampleFilter getFilter();
    }

    final class SampleMeanCoverageFilterPredicate implements SampleFilterPredicate {
        @Override
        public boolean test(final SampleCoverageStats sampleCoverageStats) {
            return sampleCoverageStats.mean >= minimumMeanCoverage && sampleCoverageStats.mean <= maximumMeanCoverage;
        }

        @Override
        public String rejectReason(final SampleCoverageStats stats) {
            return String.format("the mean coverage (%g) is outside the non-extreme range [%g, %g]", stats.mean, minimumMeanCoverage, maximumMeanCoverage);
        }

        @Override
        public SampleFilter getFilter() {
            return SampleFilter.ExtremeCoverageMean;
        }
    }

    final class SampleCoverageVarianceFilterPredicate implements  SampleFilterPredicate {
        @Override
        public String rejectReason(final SampleCoverageStats stats) {
            return String.format("the coverage variance (%g) is outside the non-extreme range [%g, %g]", stats.variance, minimumCoverageVariance, maximumCoverageVariance);
        }

        @Override
        public SampleFilter getFilter() {
            return SampleFilter.ExtremeCoverageVariance;
        }

        @Override
        public boolean test(final SampleCoverageStats sampleCoverageStats) {
            return sampleCoverageStats.variance >= minimumCoverageVariance
                    && sampleCoverageStats.variance <= maximumCoverageVariance;
        }
    }

    static final class SampleRejectionRecord {
        public final String sample;
        public final SampleFilter filter;
        public final String reason;

        public SampleRejectionRecord(final String sample, final SampleFilter filter, final String reason) {
            this.sample = sample;
            this.filter = filter;
            this.reason = reason;
        }
    }

    /**
     * Rejection file writer.
     */
    static final class SampleRejectionRecordWriter implements AutoCloseable {

        private final TableWriter<SampleRejectionRecord> fileWriter;
        private final File file;

        /**
         * Creates a new writer.  A {@code null} input file name indicates that no rejection record is to be written.
         */
        public SampleRejectionRecordWriter(final File file) {
            this.file = file;
            if (file == null) {
                fileWriter = null;
            } else {
                try {
                    fileWriter = new TableWriter<SampleRejectionRecord>(file, new TableColumnCollection(
                            REJECTION_SAMPLE_NAME_COLUMN, REJECTION_FILTER_NAME_COLUMN, REJECTION_REASON_COLUMN)) {

                        @Override
                        protected void composeLine(SampleRejectionRecord record, DataLine dataLine) {
                            dataLine.append(record.sample)
                                    .append(record.filter.toString())
                                    .append(record.reason);
                        }
                    };
                } catch (final IOException ex) {
                    throw new UserException.CouldNotCreateOutputFile(file, ex);
                }
            }
        }

        public void writeRecord(final SampleRejectionRecord record) {
            if (file != null) {
                try {
                    fileWriter.writeRecord(record);
                } catch (final IOException ex) {
                    throw new UserException.CouldNotCreateOutputFile(file, ex);
                }
            }
        }

        @Override
        public void close() {
            if (fileWriter != null) {
                try {
                    fileWriter.close();
                } catch (final IOException ex) {
                    throw new UserException.CouldNotCreateOutputFile(file, ex);
                }
            }
        }
    }

    /**
     * Main output writer.
     */
    static final class SampleListWriter implements AutoCloseable {

        private TableWriter<String> writer;
        private File file;

        SampleListWriter(final File outputFile) {
           file = outputFile;
            try {
                writer = new TableWriter<String>(outputFile, new TableColumnCollection(OUTPUT_NAME_COLUMN)) {

                    @Override
                    protected void composeLine(final String record, final DataLine dataLine) {
                        dataLine.append(record);
                    }
                };
            } catch (final IOException ex) {
                throw new UserException.CouldNotCreateOutputFile(file, ex);
            }
        }

        void writeSampleName(final String name) {
            try {
                writer.writeRecord(name);
            } catch (final IOException ex) {
                throw new UserException.CouldNotCreateOutputFile(file, ex);
            }
        }

        @Override
        public void close() {
            try {
                writer.close();
            } catch (final IOException ex) {
                throw new UserException.CouldNotCreateOutputFile(file, ex);
            }
        }
    }
}

