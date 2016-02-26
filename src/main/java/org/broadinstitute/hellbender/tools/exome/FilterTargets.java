package org.broadinstitute.hellbender.tools.exome;

import org.apache.commons.io.output.NullWriter;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Predicate;

/**
 * Tool to determine which targets to include in analysis and which to exclude
 * based on several properties, annotations and coverage statistics.
 * <p>
 *     This tool accepts several inputs with the original list of targets
 *     to consider and annotations used to do the filtering.
 * </p>
 * <p>
 *     The main output ({@link #outputFile} argument), is a target table file with the targets that make the cut.
 * </p>
 * <p>
 *     Additionally the user can specify an additional output file (using {@link #rejectedOutputFile}) that explains
 *     the reason for targets to be excluded from the output.
 * </p>
 * <p>
 *     By default no filtering is done; the user can activate filters by selecting excluding thresholds through
 *     user argument values.
 * </p>
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "filters targets based on coverage and annotations",
        oneLineSummary = "filters targets on coverage and annotations",
        programGroup = CopyNumberProgramGroup.class
)
public class FilterTargets extends CommandLineProgram {

    public static String REJECT_FILE_TARGET_COLUMN_NAME = "TARGET";
    public static String REJECT_FILE_FILTER_COLUMN_NAME = "FILTER";
    public static String REJECT_FILE_REASON_COLUMN_NAME = "REASON";

    /**
     * Enum with possible filter targets.
     */
    public enum TargetFilter {
        ExtremeTargetSize, ExtremeGCContent, ExtremeRepeatContent,
        ExtremeCoverageMean, ExtremeCoverageVariance, ExtremeCoverageInterquartileRange;
    }

    public static final int DEFAULT_MINIMUM_TARGET_SIZE = 0;
    public static final int DEFAULT_MAXIMUM_TARGET_SIZE = Integer.MAX_VALUE;
    public static final double DEFAULT_MAXIMUM_GC_CONTENT = 1.0;
    public static final double DEFAULT_MINIMUM_GC_CONTENT = 0.0;
    public static final double DEFAULT_MAXIMUM_REPEAT_CONTENT = 1.0;
    public static final double DEFAULT_MAXIMUM_COVERAGE_MEAN = Double.POSITIVE_INFINITY;
    public static final double DEFAULT_MINIMUM_COVERAGE_MEAN = Double.NEGATIVE_INFINITY;
    public static final double DEFAULT_MAXIMUM_COVERAGE_VARIANCE = Double.POSITIVE_INFINITY;
    public static final double DEFAULT_MINIMUM_COVERAGE_VARIANCE = Double.NEGATIVE_INFINITY;
    public static final double DEFAULT_MAXIMUM_COVERAGE_INTERQUARTILE_RANGE = Double.POSITIVE_INFINITY;

    public static final String MINIMUM_TARGET_SIZE_FULL_NAME = "minimumSize";
    public static final String MINIMUM_TARGET_SIZE_SHORT_NAME = "minSize";
    public static final String MAXIMUM_TARGET_SIZE_FULL_NAME = "maximumSize";
    public static final String MAXIMUM_TARGET_SIZE_SHORT_NAME = "maxSize";
    public static final String MINIMUM_GC_CONTENT_FULL_NAME = "minimumGCContent";
    public static final String MINIMUM_GC_CONTENT_SHORT_NAME = "minGC";
    public static final String MAXIMUM_GC_CONTENT_FULL_NAME = "maximumGCContent";
    public static final String MAXIMUM_GC_CONTENT_SHORT_NAME = "maxGC";
    public static final String MAXIMUM_REPEAT_CONTENT_FULL_NAME = "maximumRepeatContent";
    public static final String MAXIMUM_REPEAT_CONTENT_SHORT_NAME = "maxRepeat";
    public static final String MINIMUM_COVERAGE_MEAN_FULL_NAME = "minimumCoverageMean";
    public static final String MINIMUM_COVERAGE_MEAN_SHORT_NAME = "minCovMean";
    public static final String MAXIMUM_COVERAGE_MEAN_FULL_NAME = "maximumCoverageMean";
    public static final String MAXIMUM_COVERAGE_MEAN_SHORT_NAME = "maxCovMean";
    public static final String MINIMUM_COVERAGE_VARIANCE_FULL_NAME = "minimumCoverageVariance";
    public static final String MINIMUM_COVERAGE_VARIANCE_SHORT_NAME = "minCovVar";
    public static final String MAXIMUM_COVERAGE_VARIANCE_FULL_NAME = "maximumCoverageVariance";
    public static final String MAXIMUM_COVERAGE_VARIANCE_SHORT_NAME = "maxCovVar";
    public static final String MAXIMUM_COVERAGE_INTERQUARTILE_RANGE_FULL_NAME = "maximumCoverageIQR";
    public static final String MAXIMUM_COVERAGE_INTERQUARTILE_RANGE_SHORT_NAME = "maxCovIQR";

    public static final String REJECT_OUTPUT_FILE_FULL_NAME = "rejectedOutputFile";
    public static final String REJECT_OUTPUT_FILE_SHORT_NAME = "rejected";

    @Argument(
            doc = "Minimum target size; targets spanning fewer base-pairs are filtered out",
            fullName = MINIMUM_TARGET_SIZE_FULL_NAME,
            shortName = MINIMUM_TARGET_SIZE_SHORT_NAME,
            optional = true
    )
    protected int minimumTargetSize = DEFAULT_MINIMUM_TARGET_SIZE;

    @Argument(
            doc = "Maximum target size; targets spanning more base-pairs are filtered out",
            fullName = MAXIMUM_TARGET_SIZE_FULL_NAME,
            shortName = MAXIMUM_TARGET_SIZE_SHORT_NAME,
            optional = true
    )
    protected int maximumTargetSize = DEFAULT_MAXIMUM_TARGET_SIZE;

    @Argument(
            doc = "Minimum GC content; targets with less GC content will be excluded",
            fullName = MINIMUM_GC_CONTENT_FULL_NAME,
            shortName = MINIMUM_GC_CONTENT_SHORT_NAME,
            optional = true
    )
    protected double minimumGCContent = DEFAULT_MINIMUM_GC_CONTENT;

    @Argument(
            doc = "Maximum GC content; targets with more GC content will be excluded",
            fullName = MAXIMUM_GC_CONTENT_FULL_NAME,
            shortName = MAXIMUM_GC_CONTENT_SHORT_NAME,
            optional = true
    )
    protected double maximumGCContent = DEFAULT_MAXIMUM_GC_CONTENT;

    @Argument(
            doc = "Maximum repeat content; targets with a larger fraction of repeated content will be excluded",
            fullName = MAXIMUM_REPEAT_CONTENT_FULL_NAME,
            shortName = MAXIMUM_REPEAT_CONTENT_SHORT_NAME,
            optional = true
    )
    protected double maximumRepeatContent = DEFAULT_MAXIMUM_REPEAT_CONTENT;

    @Argument(
            doc = "Minimum target coverage mean; targets with a lower mean will be discarded",
            fullName = MINIMUM_COVERAGE_MEAN_FULL_NAME,
            shortName = MINIMUM_COVERAGE_MEAN_SHORT_NAME,
            optional = true
    )
    protected double minimumCoverageMean = DEFAULT_MINIMUM_COVERAGE_MEAN;

    @Argument(
            doc = "Maximum target coverage mean; targets with a larger mean with be discarded",
            fullName = MAXIMUM_COVERAGE_MEAN_FULL_NAME,
            shortName = MAXIMUM_COVERAGE_MEAN_SHORT_NAME,
            optional = true
    )
    protected double maximumCoverageMean = DEFAULT_MAXIMUM_COVERAGE_MEAN;

    @Argument(
            doc = "Minimum target coverage variance; targets with a lower variance will be discarded",
            fullName = MINIMUM_COVERAGE_VARIANCE_FULL_NAME,
            shortName = MINIMUM_COVERAGE_VARIANCE_SHORT_NAME,
            optional = true
    )
    protected double minimumCoverageVariance = DEFAULT_MINIMUM_COVERAGE_VARIANCE;

    @Argument(
            doc = "Maximum target coverage variance; targets with a larger variance with be discarded",
            fullName = MAXIMUM_COVERAGE_VARIANCE_FULL_NAME,
            shortName = MAXIMUM_COVERAGE_VARIANCE_SHORT_NAME,
            optional = true
    )
    protected double maximumCoverageVariance = DEFAULT_MAXIMUM_COVERAGE_VARIANCE;

    @Argument(
            doc = "Maximum target coverage interquartile range; targets with a larger interquartile range with be discarded",
            fullName = MAXIMUM_COVERAGE_INTERQUARTILE_RANGE_FULL_NAME,
            shortName = MAXIMUM_COVERAGE_INTERQUARTILE_RANGE_SHORT_NAME,
            optional = true
    )
    protected double maximumCoverageInterquartileRange = DEFAULT_MAXIMUM_COVERAGE_INTERQUARTILE_RANGE;


    @ArgumentCollection
    protected TargetArgumentCollection targetArguments
            = new TargetArgumentCollection();

    @Argument(
            doc = "Output left in targets",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            optional = false)
    protected File outputFile;

    @Argument(
            doc = "Output for left out targets",
            fullName = REJECT_OUTPUT_FILE_FULL_NAME,
            shortName = REJECT_OUTPUT_FILE_SHORT_NAME,
            optional = true)
    protected File rejectedOutputFile;

    @Override
    protected Object doWork() {

        final List<TargetFilterPredicate> filters = composeFilterPredicateList();

        final TargetCollection<Target> targets = targetArguments.readTargetCollection(false);

        try (final TargetTableWriter outputWriter = new TargetTableWriter(outputFile);
             final TargetRejectWriter rejectWriter = new TargetRejectWriter(rejectedOutputFile)) {
            for (final Target target : targets.targets()) {
                if (filters.stream().anyMatch(filter -> !filter.test(target))) {
                   filters.stream()
                           .filter(filter -> !filter.test(target))
                           .forEach((filter) -> rejectWriter.writeReason(target, filter));
                } else {
                   outputWriter.writeRecord(target);
                }
            }
            rejectWriter.logStats(logger);
            return "SUCCESS";
        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, ex);
        }
    }

    /**
     * Composes the list of target-filter predicates to be applied given
     * the user argument values.
     * @return never {@code null} but empty if no filter applies.
     */
    private List<TargetFilterPredicate> composeFilterPredicateList() {
        final List<TargetFilterPredicate> result = new ArrayList<>(3);
        if (minimumTargetSize > 0 || maximumTargetSize < Integer.MAX_VALUE) {
            result.add(new ExtremeTargetSize());
        }
        if (minimumGCContent > 0.0 || maximumGCContent < 1.0) {
            result.add(new ExtremeTargetGCContent());
        }
        if (maximumRepeatContent < 1.0) {
            result.add(new ExtremeTargetRepeatContent());
        }
        if (maximumCoverageMean < Double.POSITIVE_INFINITY ||
            minimumCoverageMean > Double.NEGATIVE_INFINITY) {
            result.add(new ExtremeCoverageMean());
        }
        if (maximumCoverageVariance < Double.POSITIVE_INFINITY ||
                minimumCoverageVariance > Double.NEGATIVE_INFINITY) {
            result.add(new ExtremeCoverageVariance());
        }
        if (maximumCoverageInterquartileRange < Double.POSITIVE_INFINITY) {
            result.add(new ExtremeCoverageInterquartileRange());
        }
        return result;
    }

    /**
     * Common interface for target filters.
     */
    private interface TargetFilterPredicate extends Predicate<Target> {

        /**
         * A name that identifies the filter.
         * @return never {@code null}
         */
        TargetFilter getFilter();

        /**
         * Compose a string explaining the reason to filter out the given
         * target.
         * <p>
         *     The result is undefined if {@link #test} does not return
         *     {@code false} of the input target.
         * </p>
         *
         * @return never null.
         */
        String reasonToFilter(final Target target);
    }


    /**
     * Common code for filter predicates that depend on a single
     * double typed target annotation.
     * <p>
     * It takes care of converting the annotation value into a double
     * throwing any pertinent {@link UserException} if the input values are
     * not valid doubles ore are absent.
     * </p>
     */
    private abstract class DoubleAnnotationBasedFilterPredicate implements TargetFilterPredicate {

        /**
         * The annotation the filter depends on.
         * @return never {@code null}.
         */
        abstract TargetAnnotation inputAnnotation();

        /**
         * Test the predicate given the double value for the dependent annotation.
         *
         * @param value the double value for the dependent annotation.
         * @return {@code true} iff the value passes the filter.
         */
        abstract boolean test(final double value);

        /**
         * User friendly explanation test as to why the annotation value does not pass the filter.
         * @param value the tested annotation value.
         * @return never {@code null}.
         */
        abstract String reasonToFilter(final double value);

        @Override
        public boolean test(final Target target) {
            if (!target.getAnnotations().hasAnnotation(inputAnnotation())) {
                throw new UserException.BadInput(String.format("requested to filter targets based on %s but this annotation is not available for target %s", inputAnnotation(), target.getName()));
            } else {
                final String value = target.getAnnotations().get(inputAnnotation());
                try {
                    final double doubleValue = Double.parseDouble(value);
                    return test(doubleValue);
                } catch (final NumberFormatException ex) {
                    throw new UserException.BadInput(String.format("%s annotation for target (%s) is not a valid double: (%s)", inputAnnotation(), target.getName(), value));
                }
            }
        }

        @Override
        public String reasonToFilter(final Target target) {
            final double doubleValue = target.getAnnotations().getDouble(inputAnnotation());
            return reasonToFilter(doubleValue);
        }
    }

    /**
     * Filter for extreme target coverage mean across samples.
     */
    private class ExtremeCoverageMean extends DoubleAnnotationBasedFilterPredicate {

        @Override
        TargetAnnotation inputAnnotation() {
            return TargetAnnotation.MEAN_COVERAGE;
        }

        @Override
        boolean test(final double value) {
            return value <= maximumCoverageMean && value >= minimumCoverageMean;
        }

        @Override
        String reasonToFilter(final double value) {
            return String.format("the coverage mean %g is outside the valid range [%g, %g]", value, minimumCoverageMean, maximumCoverageMean);
        }

        @Override
        public TargetFilter getFilter() {
            return TargetFilter.ExtremeCoverageMean;
        }
    }

    /**
     * Extreme coverage variance across sample target filter predicate.
     */
    private class ExtremeCoverageVariance extends DoubleAnnotationBasedFilterPredicate {

        @Override
        TargetAnnotation inputAnnotation() {
            return TargetAnnotation.COVERAGE_VARIANCE;
        }

        @Override
        boolean test(final double value) {
            return value <= maximumCoverageVariance && value >= minimumCoverageVariance;
        }

        @Override
        String reasonToFilter(final double value) {
            return String.format("the coverage variance %g is outside the valid range [%g, %g]", value, minimumCoverageVariance, maximumCoverageVariance);
        }

        @Override
        public TargetFilter getFilter() {
            return TargetFilter.ExtremeCoverageVariance;
        }
    }

    /**
     * Extreme coverage variance across sample target filter predicate.
     */
    private class ExtremeCoverageInterquartileRange extends DoubleAnnotationBasedFilterPredicate {

        @Override
        TargetAnnotation inputAnnotation() {
            return TargetAnnotation.COVERAGE_INTERQUARTILE_RANGE;
        }

        @Override
        boolean test(final double value) {
            return value <= maximumCoverageInterquartileRange;
        }

        @Override
        String reasonToFilter(final double value) {
            return String.format("the coverage interquartile range %g is greater than the maximum %g", value, maximumCoverageInterquartileRange);
        }

        @Override
        public TargetFilter getFilter() {
            return TargetFilter.ExtremeCoverageInterquartileRange;
        }
    }

    /**
     * Extreme target repeat fraction filter predicate.
     */
    private class ExtremeTargetRepeatContent extends DoubleAnnotationBasedFilterPredicate {

        @Override
        public TargetFilter getFilter() {
            return TargetFilter.ExtremeRepeatContent;
        }

        @Override
        public String reasonToFilter(final double repeatContent) {
            return String.format("the repeat content fraction %g is larger than the maximum permitted (%g)", repeatContent, maximumRepeatContent);
        }

        @Override
        public boolean test(final double repeatContent) {
            return repeatContent <= maximumRepeatContent;
        }

        @Override
        public TargetAnnotation inputAnnotation() {
            return TargetAnnotation.REPEAT_FRACTION;
        }
    }

    /**
     * Extreme GC Content filter predicate.
     */
    private class ExtremeTargetGCContent extends DoubleAnnotationBasedFilterPredicate {

        @Override
        public TargetFilter getFilter() {
            return TargetFilter.ExtremeGCContent;
        }

        @Override
        public String reasonToFilter(final double gcContent) {
            return String.format("the GC content fraction %g is outside the range [%g, %g]", gcContent, minimumGCContent, maximumGCContent);
        }

        @Override
        public boolean test(final double gcContent) {
            return gcContent <= maximumGCContent && gcContent >= minimumGCContent;
        }

        @Override
        public TargetAnnotation inputAnnotation() {
            return TargetAnnotation.GC_CONTENT;
        }
    }

    /**
     * Filters targets that are too small as per
     * the {@link #minimumTargetSize} user argument.
     */
    private class ExtremeTargetSize implements TargetFilterPredicate {

        /**
         * Creates the filter.
         */
        public ExtremeTargetSize() {}

        @Override
        public TargetFilter getFilter() {
            return TargetFilter.ExtremeTargetSize;
        }

        @Override
        public String reasonToFilter(final Target target) {
            return String.format("the target size (%d) is out of range [%d, %d]",
                    target.getInterval().size(), minimumTargetSize, maximumTargetSize);
        }

        @Override
        public boolean test(final Target target) {
            final SimpleInterval interval = target.getInterval();
            return interval.size() >= minimumTargetSize && interval.size() <= maximumTargetSize;
        }
    }

    /**
     * Writes filter failure explanatory records into a tab-table file.
     * <p>
     * Keeps counters on how many targets have failed each filter.
     * </p>
     */
    private class TargetRejectWriter implements  AutoCloseable {

        /**
         * Reference to the underlying reason table file writer.
         */
        private final TableWriter<Pair<Target, TargetFilterPredicate>> outputWriter;

        /**
         * Holds counter for how many times we have rejected a target
         * based on each filter.
         * <p>Initialized to all zeros</p>
         */
        private final int[] filterCounters;

        /**
         * Output file.
         * <p>It can be null, in which case we don't output reasons</p>
         */
        private final File outputFile;

        /**
         * Creates a target filter rejection writer.
         * <p>
         * The output file can be {@code null} in which case nothing will
         * be output. Nevertheless this writer will keep collecting some
         * stats that might be dumped to a logger using {@link #logStats}.
         * </p>
         *
         * @param file the output file. {@code null} if nothing is to be
         *             outputted.
         * @throws IOException
         */
        public TargetRejectWriter(final File file) throws IOException {
            this.outputFile = file;
            this.outputWriter = composeWriter(file);
            this.filterCounters = new int[TargetFilter.values().length];
        }

        /**
         * Constructs the output table writer given the output file-name.
         * <p>The output file can be {@code null}, in which case
         * reason are dump into no file.</p>
         *
         * @param file the output-file.
         * @return never {@code null}
         */
        private TableWriter<Pair<Target, TargetFilterPredicate>> composeWriter(final File file) {
            final TableColumnCollection columns =
                    new TableColumnCollection(
                            REJECT_FILE_TARGET_COLUMN_NAME,
                            REJECT_FILE_FILTER_COLUMN_NAME,
                            REJECT_FILE_REASON_COLUMN_NAME);
            try {
                final Writer outputWriter = (file == null) ? new NullWriter() : new FileWriter(file);
                return new TableWriter<Pair<Target, TargetFilterPredicate>>(outputWriter, columns) {
                    @Override
                    protected void composeLine(final Pair<Target, TargetFilterPredicate> record, final DataLine dataLine) {
                        dataLine.append(record.getLeft().getName())
                                .append(record.getRight().getFilter().toString())
                                .append(record.getRight().reasonToFilter(record.getLeft()));
                    }
                };
            } catch (final IOException ex) {
                if (file == null) {
                    throw new GATKException("unexpected IO exception", ex);
                } else {
                    throw new UserException.CouldNotCreateOutputFile(file, ex);
                }
            }
        }

        /**
         * Dumps the rejection reason for a target due to a filter.
         *
         * <p>
         *     This method assumes that the target has already failed the filter.
         * </p>
         *
         * @param target the target.
         * @param filter the target filter.
         */
        public void writeReason(final Target target, final TargetFilterPredicate filter) {
            Utils.nonNull(target);
            Utils.nonNull(filter);
            try {
                filterCounters[filter.getFilter().ordinal()]++;
                outputWriter.writeRecord(new ImmutablePair<>(target, filter));
            } catch (final IOException ex) {
                throw new UserException.CouldNotCreateOutputFile(outputFile, ex);
            }
        }

        @Override
        public void close()  {
            try {
                outputWriter.close();
            } catch (final IOException ex) {
                throw new UserException.CouldNotCreateOutputFile(outputFile, ex);
            }
        }

        /**
         * Dump filtering stats to a logger as a set of info messages.
         *
         * @param logger the logger where to dump the stats.
         */
        public void logStats(final Logger logger) {
            for (final TargetFilter filter : TargetFilter.values()) {
                final int targetCount = this.filterCounters[filter.ordinal()];
                logger.info(String.format("%s : %d targets filtered out", filter, targetCount));
            }
        }
    }
}