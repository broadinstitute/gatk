package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ExomeAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Collections;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Calculate read-counts on an exome given its intervals.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        usage = "Count overlapping reads exon by exon",
        usageShort = "Count overlapping reads exon by exon",
        programGroup = ExomeAnalysisProgramGroup.class
)

public final class ExomeReadCounts extends ReadWalker {

    /**
     * Default cohort name used unless one is provided using argument {@link #cohortName}.
     */
    public static final String DEFAULT_COHORT_NAME = "<ALL>";

    /**
     * Full name for the {@link #groupBy} argument.
     */
    protected static final String GROUP_BY_FULL_NAME = "groupBy";

    /**
     * Short name for the {@link #groupBy} argument.
     */
    protected static final String GROUP_BY_SHORT_NAME = GROUP_BY_FULL_NAME;

    /**
     * Short name for the {@link #columnSummaryOutput} argument.
     */
    protected static final String COLUMN_SUMMARY_OUTPUT_SHORT_NAME = "CSO";

    /**
     * Full name for the {@link #columnSummaryOutput} argument.
     */
    protected static final String COLUMN_SUMMARY_OUTPUT_FULL_NAME = "columnSummaryOutput";

    /**
     * Short name for the {@link #rowSummaryOutput} argument.
     */
    protected static final String ROW_SUMMARY_OUTPUT_SHORT_NAME = "RSO";

    /**
     * Full name for the {@link #rowSummaryOutput} argument.
     */
    protected static final String ROW_SUMMARY_OUTPUT_FULL_NAME = "rowSummaryOutput";

    /**
     * Full name for the {@link #cohortName} argument.
     */
    protected static final String COHORT_FULL_NAME = "cohortName";

    /**
     * Short name for the {@link #cohortName} argument.
     */
    protected static final String COHORT_SHORT_NAME = "cohort";

    /**
     * Format string for the column total average per bp.
     */
    protected static final String AVERAGE_DOUBLE_FORMAT = "%.4f";

    /**
     * Format string for the pcov output.
     */
    private static final String PCOV_OUTPUT_DOUBLE_FORMAT = "%.4g";

    /**
     * Transform argument full name.
     */
    protected static final String TRANSFORM_FULL_NAME = "transform";

    /**
     * Transform argument short name.
     */
    protected static final String TRANSFORM_SHORT_NAME = TRANSFORM_FULL_NAME;

    @Argument(
            doc = "output tabular file with the counts",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = false)
    protected File output = null;

    @Argument(
            doc = "summary per column output tabular file",
            shortName = COLUMN_SUMMARY_OUTPUT_SHORT_NAME,
            fullName = COLUMN_SUMMARY_OUTPUT_FULL_NAME,
            optional = true)
    protected File columnSummaryOutput = null;

    @Argument(
            doc = "summary per row output tabular file",
            shortName = ROW_SUMMARY_OUTPUT_SHORT_NAME,
            fullName = ROW_SUMMARY_OUTPUT_FULL_NAME,
            optional = true)
    protected File rowSummaryOutput = null;

    @Argument(
            doc = "group counts by Cohort (all samples in one), Sample or ReadGroup",
            shortName = GROUP_BY_SHORT_NAME,
            fullName = GROUP_BY_FULL_NAME,
            optional = true)
    protected GroupBy groupBy = GroupBy.COHORT;

    @Argument(
            doc = "Cohort name used to name the read count column when the groupBy " +
                    "argument is set to COHORT (default)",
            shortName = COHORT_SHORT_NAME,
            fullName = COHORT_FULL_NAME,
            optional = true)
    protected String cohortName = DEFAULT_COHORT_NAME;


    @Argument(
            doc = "Transformation to perform to the individual count output values",
            shortName = TRANSFORM_SHORT_NAME,
            fullName = TRANSFORM_FULL_NAME,
            optional = true
    )
    protected Transform transform = Transform.RAW;

    /**
     * Writer to the main output file indicated by {@link #output}.
     */
    private PrintWriter outputWriter;

    /**
     * Writer to the per-column summary output file indicated by {@link #columnSummaryOutput}.
     */
    private PrintWriter columnSummaryOutputWriter;

    /**
     * Writer to the per-row summary output file indicated by {@link #rowSummaryOutput}.
     */
    private PrintWriter rowSummaryOutputWriter;

    /**
     * Sample and read-group correspondence database.
     */
    private SampleCollection sampleCollection;

    /**
     * Reference to the logger.
     */
    private static final Logger logger = LogManager.getLogger(ExomeReadCounts.class);

    /**
     * Exon database reference.
     */
    private ExonCollection<SimpleInterval> exonCollection;

    /**
     * Counts table.
     */
    private CountColumns countColumns;

    /**
     * Count matrix indexed by count column and then exon.
     */
    private int[][] counts;

    @Override
    public ReadFilter makeReadFilter() {
        return super.makeReadFilter()
                .and(ReadFilterLibrary.MAPPED);
    }

    @Override
    public void onTraversalStart() {

        // Initializing meta-data structures (about samples, exons and so forth):
        final SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();
        if (sequenceDictionary == null) {
            throw new UserException("missing reference sequence dictionary, as no input contains such a dictionary... please consider to provide fail-over dictionary using " +
                    "a reference (e.g. --reference my_reference.fasta)");
        }
        if (!intervalArgumentCollection.intervalsSpecified()) {
            throw new UserException("you must indicate the set of exon as input intervals");
        }
        sampleCollection = new SampleCollection(getHeaderForReads());

        logger.log(Level.INFO, "Reading exons locations from intervals...");

        exonCollection = new IntervalBackedExonCollection(
                intervalArgumentCollection.getIntervals(sequenceDictionary));

        // Initializing count and count column management member fields:
        countColumns = groupBy.countColumns(this);
        final int columnCount = countColumns.columnCount();
        counts = new int[columnCount][exonCollection.exonCount()];

        // Open output files and write headers:
        outputWriter = openOutputWriter(output, composeMatrixOutputHeader(getCommandLine(), groupBy, countColumns.columnNames()));
        if (columnSummaryOutput != null) {
            columnSummaryOutputWriter = openOutputWriter(columnSummaryOutput,
                    composeColumnSummaryHeader(getCommandLine(), groupBy, exonCollection.exonCount(), exonCollection.exomeSize()));
        }
        if (rowSummaryOutput != null) {
            rowSummaryOutputWriter = openOutputWriter(rowSummaryOutput,
                    composeRowOutputHeader(getCommandLine(), groupBy, countColumns.columnCount()));
        }

        // Next we start the traversal:
        logger.log(Level.INFO, "Collecting read counts ...");
    }

    /**
     * Opens the output file for writing with a print-writer.
     *
     * @param output     the output file.
     * @param headerText to be printed immediately after opening the writer.
     * @return never {@code null}.
     * @throws UserException.CouldNotCreateOutputFile if there was some problem creating or overwriting {@code output}.
     */
    private PrintWriter openOutputWriter(final File output, final String headerText) {
        try {
            final PrintWriter result = new PrintWriter(output);
            result.println(headerText);
            result.flush();
            return result;
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(output, e);
        }
    }

    @Override
    public void apply(final SAMRecord read, final ReferenceContext referenceContext, final FeatureContext featureContext) {

        final SimpleInterval readLocation = new SimpleInterval(read);

        final int columnIndex = countColumns.columnIndex(read);
        if (columnIndex >= 0) { // < 0 would means that the read is to be ignored.
            exonCollection.indexRange(readLocation).forEach(i -> counts[columnIndex][i]++);
        }
    }

    @Override
    public Object onTraversalDone() {
        super.onTraversalDone();
        logger.log(Level.INFO, "Collecting read counts done.");

        logger.log(Level.INFO, "Writing counts ...");
        final int exonCount = exonCollection.exonCount();
        final int columnCount = counts.length;
        final long[] columnTotals = calculateColumnTotals();

        for (int i = 0; i < exonCount; i++) {
            final int[] countBuffer = new int[columnCount];
            final SimpleInterval exonInterval = exonCollection.exon(i);
            for (int j = 0; j < columnCount; j++) {
                countBuffer[j] = counts[j][i];
            }
            writeOutputRows(countBuffer, columnTotals, exonInterval);
        }
        logger.log(Level.INFO, "Writing counts done.");

        writeColumnSummaryOutput();
        closeOutputs();
        return "SUCCESS";
    }

    /**
     * Calculates the column totals.
     *
     * @return never {@code null}.
     */
    private long[] calculateColumnTotals() {
        final long[] result = new long[counts.length];

        for (int i = 0; i < counts.length; i++) {
            result[i] = IntStream.of(counts[i]).sum();
        }
        return result;
    }

    /**
     * Writes the column summary output table.
     */
    private void writeColumnSummaryOutput() {
        if (columnSummaryOutputWriter == null) {
            return;
        }

        final long exomeSize = exonCollection.exomeSize();

        final List<String> columnNames = countColumns.columnNames();
        for (int i = 0; i < columnNames.size(); i++) {
            final long sum = IntStream.of(counts[i]).sum();
            columnSummaryOutputWriter.println(
                    String.join("\t",
                            columnNames.get(i),
                            String.valueOf(sum),
                            String.format(AVERAGE_DOUBLE_FORMAT, sum / (double) exomeSize)));
        }
        columnSummaryOutputWriter.close();
    }

    /**
     * Close all output writers.
     */
    private void closeOutputs() {
        outputWriter.close();
        if (rowSummaryOutputWriter != null) {
            rowSummaryOutputWriter.close();
        }
        if (columnSummaryOutputWriter != null) {
            columnSummaryOutputWriter.close();
        }
    }

    /**
     * Writes the row in the main matrix output file for an exon and, if requested,
     * the corresponding row in the row summary output file.
     *
     * @param countBuffer  the counts for the target exon.
     * @param exonInterval genomic location of the target exon.
     */
    private void writeOutputRows(final int[] countBuffer, final long[] columnTotals,
                                 final SimpleInterval exonInterval) {
        final String countString = IntStream.range(0, countBuffer.length).mapToObj(
                i -> transform.apply(countBuffer[i], columnTotals[i])).collect(Collectors.joining("\t"));
        final String coordinateString = String.join("\t",
                exonInterval.getContig(),
                Integer.toString(exonInterval.getStart()),
                Integer.toString(exonInterval.getEnd()));

        outputWriter.println(String.join("\t", coordinateString, countString));

        if (rowSummaryOutputWriter != null) {
            final long sum = MathUtils.sum(countBuffer);
            rowSummaryOutputWriter.println(String.join("\t",
                    coordinateString,
                    Long.toString(sum), String.format(AVERAGE_DOUBLE_FORMAT,
                            sum / ((float) countColumns.columnCount() * exonInterval.size()))));
        }
    }

    /**
     * Composes the column summary output header.
     * <p>
     * Returns a multiline header prepended with the comment scape '##' sequence for general information an
     * final uncommented single line with the column headers.
     * </p>
     *
     * @param commandLine the command-line.
     * @param exonCount   number of exons processed.
     * @param exomeSize   size of the exome processed in bps.
     * @return never {@code null}.
     */
    private static String composeColumnSummaryHeader(final String commandLine, final GroupBy groupBy,
                                                     final int exonCount, final long exomeSize) {
        return String.format(
                String.join("\n",
                        "##fileFormat  = tsv",
                        "##commandLine = %s",
                        "##title       = Summary counts per %s",
                        "##metaData    = {",
                        "##    exonCount = %d,",
                        "##    exomeSize = %d (bp)",
                        "##}",
                        String.join("\t", "GROUP", "SUM", "AVG.BP")),
                commandLine, groupBy.toString(), exonCount, exomeSize);
    }

    /**
     * Composes the row summary output header.
     *
     * @param commandLine the execution command line.
     * @param groupBy     the value of the group-by argument used.
     * @param columnCount number of count-column involved in the analysis.
     * @return never {@code null}.
     */
    private static String composeRowOutputHeader(final String commandLine, final GroupBy groupBy,
                                                 final int columnCount) {
        return String.format(
                String.join("\n",
                        "##fileFormat  = tsv",
                        "##commandLine = %s",
                        "##title       = Summary counts per exon",
                        "##metaData = {",
                        "##    groupBy     = %s,",
                        "##    columnCount = %d",
                        "##}",
                        String.join("\t", "CHROM", "START", "STOP", "SUM", "AVG.COL.BP")),
                commandLine, groupBy, columnCount);
    }

    /**
     * Composes the main output header.
     *
     * @param commandLine      the tool command line.
     * @param groupBy          the group-by argument used.
     * @param countColumnNames the column names.
     * @return never {@code null}.
     */
    private static String composeMatrixOutputHeader(final String commandLine, final GroupBy groupBy,
                                                    final List<String> countColumnNames) {
        final String countColumnHeaderString = String.join("\t", countColumnNames).replace("%", "%%");
        final String formatString = String.join("\n",
                "##fileFormat  = tsv",
                "##commandLine = %s",
                "##title       = Read counts per exon and %s",
                String.join("\t", "CHROM", "START", "STOP", countColumnHeaderString));
        return String.format(formatString, commandLine, groupBy.toString());
    }

    /////////////////////////////////
    // Count column holder classes //
    /////////////////////////////////

    /**
     * Count-column manager.
     * <p>
     * Instances implementing this class are
     * responsible to determine the number
     * of count-columns, the header names and
     * how reads distribute across columns.
     * </p>
     */
    private abstract static class CountColumns {

        /**
         * Reference to the counting tool.
         */
        protected final ExomeReadCounts tool;

        /**
         * Creates a new count-counts instance given the reference to the embedding tool instance.
         *
         * @param tool the read-count tool instance.
         */
        private CountColumns(final ExomeReadCounts tool) {
            this.tool = tool;
        }

        /**
         * Composes the list of count column names.
         * <p>
         * <p>
         * The result list would have exactly <code>{@link #columnCount()}</code> elements.
         * </p>
         *
         * @return never {@code null}.
         */
        protected abstract List<String> columnNames();

        /**
         * Returns the number of count columns headers.
         *
         * @return never {@code null}.
         */
        public final int columnCount() {
            return columnNames().size();
        }

        /**
         * Returns the count column index corresponding to a read.
         *
         * @param read the query read.
         * @return a value from 0 to <code>{@link #columnCount()} - 1</code>, or -1 if the value is to be discarded.
         */
        protected abstract int columnIndex(SAMRecord read);
    }

    /**
     * Counts class to be use when {@link #groupBy} is {@code GroupBy#COHORT}.
     */
    protected static final class CohortCountColumns extends CountColumns {

        /**
         * Creates a new count-counts instance given the reference to the embedding tool instance.
         */
        private CohortCountColumns(final ExomeReadCounts tool) {
            super(tool);
        }

        @Override
        protected List<String> columnNames() {
            return Collections.singletonList(tool.cohortName);
        }

        @Override
        protected int columnIndex(final SAMRecord read) {
            return 0;
        }
    }

    /**
     * Counts class to be used with {@link GroupBy#SAMPLE}.
     * <p>
     * Keeps independent counts for each sample.
     * </p>
     */
    private static final class SampleCountColumns extends CountColumns {

        /**
         * Creates a new per sample count-counts manager given the reference to the embedding tool instance.
         */
        private SampleCountColumns(final ExomeReadCounts tool) {
            super(tool);
        }

        @Override
        protected List<String> columnNames() {
            return tool.sampleCollection.sampleIds();
        }

        @Override
        protected int columnIndex(final SAMRecord read) {
            return tool.sampleCollection.sampleIndexByRead(read);
        }
    }

    /**
     * Count-column manager when grouping by read-groups.
     * <p>
     * Keeps independent counts per read-group.
     * </p>
     */
    private static final class ReadGroupCountColumns extends CountColumns {

        /**
         * Creates a new per sample count-counts manager given the reference to the embedding tool instance.
         */
        private ReadGroupCountColumns(final ExomeReadCounts tool) {
            super(tool);
        }

        @Override
        protected List<String> columnNames() {
            return tool.sampleCollection.readGroups();
        }

        @Override
        protected int columnIndex(final SAMRecord read) {
            return tool.sampleCollection.readGroupIndexByRead(read);
        }
    }

    /**
     * Possible read count grouping.
     */
    protected enum GroupBy {

        /**
         * Single read count for all input.
         */
        COHORT(CohortCountColumns::new),

        /**
         * Count per sample.
         */
        SAMPLE(SampleCountColumns::new),

        /**
         * Count per read-group.
         */
        READ_GROUP(ReadGroupCountColumns::new);

        /**
         * Corresponding count-columns manager factory.
         */
        private final Function<ExomeReadCounts, ? extends CountColumns> countColumnsFactory;

        /**
         * Creates a new group-by option.
         *
         * @param countColumnsFactory the count-columns manager factory to be used with
         *                            this group-by option. Assumed not to be {@code null}.
         */
        GroupBy(final Function<ExomeReadCounts, ? extends CountColumns> countColumnsFactory) {
            this.countColumnsFactory = countColumnsFactory;
        }

        /**
         * Returns a lowercase user friendly name for this group-by option.
         * <p>This is to be used for some user targeted messages (error, log, etc)</p>
         *
         * @return never {@code null}.
         */
        @Override
        public String toString() {
            return name().toLowerCase().replace("_", " ");
        }

        /**
         * Creates a column count manager for a given read-counts tool.
         *
         * @param tool the targeted tool.
         * @return never {@code null}.
         * @throws IllegalArgumentException if {@code tool} is {@code null}.
         * @throws IllegalStateException    if {@code tool} is not in a compatible state for this type of count-columns.
         */
        protected CountColumns countColumns(final ExomeReadCounts tool) {
            return countColumnsFactory.apply(tool);
        }
    }

    /**
     * Matrix value transformation function to apply to the main output matrix values.
     */
    protected enum Transform {

        /**
         * Default read-count raw value (non-)transformation.
         */
        RAW((count, columnTotal) -> Integer.toString(count)),

        /**
         * Proportional coverage transformation.
         * <p>Individual counts are transformed into the fraction of the total
         * count across the enclosing column.</p>
         */
        PCOV((count, columnTotal) ->
                String.format(PCOV_OUTPUT_DOUBLE_FORMAT, count / (double) columnTotal));

        /**
         * Functional interface for the count transformation.
         */
        @FunctionalInterface
        protected interface Operator {

            /**
             * Output matrix value transformer method.
             * <p>It takes exactly two double arguments: the individual count and the column total sum.</p>
             * <p>Implementation of this method can assume that individual input {@code count} is 0 or greater and
             * that is not greater than {@code columnTotal}.</p>
             *
             * @param count       the individual count for an exon and count group
             * @param columnTotal the total count for the enclosing count group.
             * @return never {@code null}.
             */
            String apply(final int count, final long columnTotal);
        }

        /**
         * Holds a reference to the transformation operator.
         */
        private final Operator operator;

        /**
         * Creates a {@link Transform} instance given the corresponding transformation operator.
         *
         * @param operator the value transformation operator.
         */
        Transform(final Operator operator) {
            this.operator = operator;
        }

        /**
         * Transforms and composes the string representation of an individual count.
         * <p>The output string must be fully formatted human friendly representation of the
         * transformed value.</p>
         *
         * @param count       the individual count value.
         * @param columnTotal the corresponding column total sum.
         * @return never {@code null}.
         * @throws IllegalArgumentException if {@code count} is less than 0 or greater than {@code columnTotal}.
         */
        protected String apply(final int count, final long columnTotal) {
            if (count < 0) {
                throw new IllegalArgumentException("the count count cannot less than 0");
            }
            if (count > columnTotal) {
                throw new IllegalArgumentException("the count count cannot be larger than the column total");
            }
            return operator.apply(count, columnTotal);
        }
    }
}
