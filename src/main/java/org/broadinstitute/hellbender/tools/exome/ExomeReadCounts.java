package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.bed.BEDFeature;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ExomeAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureManager;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Calculate read-counts on an exome given its intervals.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Count overlapping reads exon by exon",
        oneLineSummary = "Count overlapping reads exon by exon",
        programGroup = ExomeAnalysisProgramGroup.class
)
public final class ExomeReadCounts extends ReadWalker {

    /**
     * Default cohort name used unless one is provided using argument {@link #cohortName}.
     */
    public static final String DEFAULT_COHORT_NAME = "<ALL>";

    /**
     * Header name for the column that contains the exon contig name.
     */
    public static final String EXON_CONTIG_COLUMN_NAME = "CONTIG";

    /**
     * Header name for the column that contains the exon start position.
     */
    public static final String EXON_START_COLUMN_NAME = "START";

    /**
     * Header name for the column that contains the exon end position.
     */
    public static final String EXON_END_COLUMN_NAME = "END";

    /**
     * Header name for the column that contains the exon name.
     */
    public static final String EXON_NAME_COLUMN_NAME = "NAME";

    /**
     * Header name for the column that contains a count sum (across columns or rows).
     */
    public static final String SUM_COLUMN_NAME = "SUM";

    /**
     * Header name for the column that contains an average per base-pair per column.
     */
    public static final String AVG_BP_COLUMN_NAME = "AVG.BP";

    /**
     * Header name for the column that contains an average per base-pair per row.
     */
    public static final String AVG_COL_BP_COLUMN_NAME = "AVG.COL.BP";

    /**
     * Output value separator.
     */
    public static final String COLUMN_SEPARATOR = "\t";

    /**
     * Output line separator.
     */
    public static final String LINE_SEPARATOR = "\n";

    /**
     * Output no-value string.
     */
    public static final String NO_VALUE_STRING = ".";

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

    /**
     * Exome file argument full name.
     */
    protected static final String EXOME_FILE_FULL_NAME = "exome";

    /**
     * Exome file argument short name.
     */
    protected static final String EXOME_FILE_SHORT_NAME = EXOME_FILE_FULL_NAME;

    /**
     * Exon output info argument full name
     */
    protected static final String EXON_OUT_INFO_FULL_NAME = "exonInformationColumns";

    /**
     * Exon output info argument short name
     */
    protected static final String EXON_OUT_INFO_SHORT_NAME = "exonInfo";


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

    @Argument(
            doc = "File containing the exons for analysis",
            shortName = EXOME_FILE_SHORT_NAME,
            fullName = EXOME_FILE_FULL_NAME,
            optional = true
    )
    protected File exomeFile = null;

    @Argument(
            doc = "Exon identification information to be show in outputs",
            shortName = EXON_OUT_INFO_SHORT_NAME,
            fullName = EXON_OUT_INFO_FULL_NAME,
            optional = true
    )
    protected ExonOutInfo exonOutInfo = ExonOutInfo.COORDS;

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
    private ExonCollection<? extends Locatable> exonCollection;

    /**
     * Counts table.
     */
    private CountColumns countColumns;

    /**
     * Count matrix indexed by count column and then exon.
     */
    private int[][] counts;

    @Override
    public CountingReadFilter makeReadFilter() {
        return super.makeReadFilter()
                .and(new CountingReadFilter("Mapped", ReadFilterLibrary.MAPPED))
                .and(new CountingReadFilter("Not_Duplicate", ReadFilterLibrary.NOT_DUPLICATE))
                .and(new CountingReadFilter("Non_Zero_Reference_Length", ReadFilterLibrary.NON_ZERO_REFERENCE_LENGTH_ALIGNMENT));
    }


    @Override
    public void onTraversalStart() {

        sampleCollection = new SampleCollection(getHeaderForReads());

        logger.log(Level.INFO, "Reading exons locations from intervals...");

        exonCollection = resolveExonCollection();

        // Initializing count and count column management member fields:
        countColumns = groupBy.countColumns(this);
        final int columnCount = countColumns.columnCount();
        counts = new int[columnCount][exonCollection.exonCount()];

        // Open output files and write headers:
        outputWriter = openOutputWriter(output, composeMatrixOutputHeader(getCommandLine(), exonOutInfo, groupBy, countColumns.columnNames()));
        if (columnSummaryOutput != null) {
            columnSummaryOutputWriter = openOutputWriter(columnSummaryOutput,
                    composeColumnSummaryHeader(getCommandLine(), groupBy, exonCollection.exonCount(), exonCollection.exomeSize()));
        }
        if (rowSummaryOutput != null) {
            rowSummaryOutputWriter = openOutputWriter(rowSummaryOutput,
                    composeRowOutputHeader(getCommandLine(), exonOutInfo, groupBy, countColumns.columnCount()));
        }

        // Next we start the traversal:
        logger.log(Level.INFO, "Collecting read counts ...");
    }

    /**
     * Builds the exon collection given the values of user arguments.
     *
     * @throws UserException if there is some inconsistency in user arguments and inputs.
     * @throws GATKException if there was any problem parsing the content of the exome file.
     * @return never {@code null}.
     */
    private ExonCollection<? extends Locatable> resolveExonCollection() {
        final ExonCollection<? extends Locatable> result;
        if (exomeFile != null) {
            result = resolveExonCollectionFromExomeFile();
        } else if (hasIntervals()) {
            final SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();
            result = ExonCollections.fromSimpleIntervalList(intervalArgumentCollection.getIntervals(sequenceDictionary));
        } else {
            throw new UserException(String.format("You must indicate the set of exon as input intervals (e.g. -L target-intervals.list) or a exome feature file (e.g. -%s my-targets.bed) ",EXOME_FILE_SHORT_NAME));
        }
        if (exonOutInfo.requiresUniqueExonName()) {
            checkAllExonsHaveName(result);
        }
        return result;
    }

    /**
     * Checks whether all exons in the input exon collection have a designated name.
     * @param result the query exon collection.
     * @throws UserException if there are some exons with no name in {@code result}.
     */
    private <T extends Locatable> void checkAllExonsHaveName(ExonCollection<T> result) {
        if (result.exons().stream().anyMatch(e -> result.name(e) == null)) {
            throw new UserException(String.format("Exon output info requested '%s' requires that each exon has a designated unique name/id but there are some with no names: %s", exonOutInfo.name(),
                    result.exons().stream().filter(e -> result.name(e) == null).limit(10).map(e -> result.location(e).toString()).collect(Collectors.joining(", "))));
        }
    }

    /**
     * Constructs the exon collection from an exome-file passed by the user.
     *
     * @return never {@code null}.
     */
    private ExonCollection<? extends Locatable> resolveExonCollectionFromExomeFile() {
        Utils.regularReadableUserFile(exomeFile);
        final FeatureCodec<? extends Feature, ?> codec = FeatureManager.getCodecForFile(exomeFile);
        logger.log(Level.INFO,String.format("Reading exon intervals from exome file '%s' ...",exomeFile.getAbsolutePath()));
        final Class<? extends Feature> featureType = codec.getFeatureType();
        final ExonCollection<? extends Locatable> result;
        if (BEDFeature.class.isAssignableFrom(featureType)) {
            @SuppressWarnings("unchecked")
            final FeatureCodec<? extends BEDFeature, ?> bedFeatureCodec = (FeatureCodec<? extends BEDFeature, ?>) codec;
            result = ExonCollections.fromBEDFeatureFile(exomeFile, bedFeatureCodec);
        } else {
            throw new UserException.BadInput(String.format("currently only BED formatted exome file are supported. '%s' does not seem to be a BED file",exomeFile.getAbsolutePath()));
        }
        logger.log(Level.INFO,String.format("Found %d exons to analyze.",result.exonCount()));
        return result;
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
    public void apply(final GATKRead read, final ReferenceContext referenceContext, final FeatureContext featureContext) {

        final SimpleInterval readLocation = referenceContext.getInterval();

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
            for (int j = 0; j < columnCount; j++) {
                countBuffer[j] = counts[j][i];
            }
            writeOutputRows(countBuffer, columnTotals, i);
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
                    String.join(COLUMN_SEPARATOR,
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
     * @param index the index of exon within the exon collection.
     */
    private void writeOutputRows(final int[] countBuffer, final long[] columnTotals,
                                 final int index) {
        final String countString = IntStream.range(0, countBuffer.length).mapToObj(
                i -> transform.apply(countBuffer[i], columnTotals[i])).collect(Collectors.joining(COLUMN_SEPARATOR));
        final String exonInfoString = exonOutInfo.composeExonOutInfoString(index, exonCollection);

        outputWriter.println(String.join(COLUMN_SEPARATOR, exonInfoString, countString));

        if (rowSummaryOutputWriter != null) {
            final long sum = MathUtils.sum(countBuffer);
            final SimpleInterval location = exonCollection.location(index);
            final int exonSize = location.size();
            rowSummaryOutputWriter.println(String.join(COLUMN_SEPARATOR,
                    exonInfoString,
                    Long.toString(sum), String.format(AVERAGE_DOUBLE_FORMAT,
                            sum / ((float) countColumns.columnCount() * exonSize))));
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
                String.join(LINE_SEPARATOR,
                        "##fileFormat  = tsv",
                        "##commandLine = %s",
                        "##title       = Summary counts per %s",
                        "##metaData    = {",
                        "##    exonCount = %d,",
                        "##    exomeSize = %d (bp)",
                        "##}",
                        String.join(COLUMN_SEPARATOR, groupBy.name(), SUM_COLUMN_NAME, AVG_BP_COLUMN_NAME)),
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
    private static String composeRowOutputHeader(final String commandLine, final ExonOutInfo exonOutInfo, final GroupBy groupBy,
                                                 final int columnCount) {
        return String.format(
                String.join(LINE_SEPARATOR,
                        "##fileFormat  = tsv",
                        "##commandLine = %s",
                        "##title       = Summary counts per exon",
                        "##metaData = {",
                        "##    groupBy     = %s,",
                        "##    columnCount = %d",
                        "##}",
                        String.join(COLUMN_SEPARATOR,exonOutInfo.headerString(),SUM_COLUMN_NAME,AVG_COL_BP_COLUMN_NAME)),
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
    private static String composeMatrixOutputHeader(final String commandLine, final ExonOutInfo exonOutInfo, final GroupBy groupBy,
                                                    final List<String> countColumnNames) {
        final String countColumnHeaderString = String.join(COLUMN_SEPARATOR, countColumnNames).replace("%", "%%");
        final String formatString = String.join(LINE_SEPARATOR,
                "##fileFormat  = tsv",
                "##commandLine = %s",
                "##title       = Read counts per exon and %s",
                String.join(COLUMN_SEPARATOR, exonOutInfo.headerString(), countColumnHeaderString));
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
        protected abstract int columnIndex(GATKRead read);
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
        protected int columnIndex(final GATKRead read) {
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
        protected int columnIndex(final GATKRead read) {
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
        protected int columnIndex(final GATKRead read) {
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

    /**
     * Possible per-exon information included in the output.
     */
    protected enum ExonOutInfo {

        /**
         * Only exon coordinates (contig name, start and end positions).
         */
        COORDS(ExonOutInfo::coordinateComposer,
                EXON_CONTIG_COLUMN_NAME,EXON_START_COLUMN_NAME, EXON_END_COLUMN_NAME),

        /**
         * Only exon name/id.
         */
        NAME(ExonOutInfo::nameComposer,
                EXON_NAME_COLUMN_NAME),

        /**
         * Exon coordinates and name.
         */
        FULL(ExonOutInfo::coordinateAndNameComposer,
                EXON_CONTIG_COLUMN_NAME,EXON_START_COLUMN_NAME, EXON_END_COLUMN_NAME,EXON_NAME_COLUMN_NAME);

        /**
         * Exon information string composer for the genomic coordinate part of the exon.
         *
         * @param index      the target exon within the collection.
         * @param collection the containing exon collection.
         * @param <T>        the exon type.
         * @return never {@code null}.
         */
        private static <T> String coordinateComposer(final int index, final ExonCollection<T> collection) {
            final SimpleInterval location = collection.location(index);
            if (location == null) {
                return String.join(COLUMN_SEPARATOR, NO_VALUE_STRING, NO_VALUE_STRING, NO_VALUE_STRING);
            } else {
                return String.format(String.join(COLUMN_SEPARATOR, "%s", "%d", "%d"),
                        location.getContig(), location.getStart(), location.getEnd());
            }
        }

        /**
         * Exon information string composer for the name part of the exon information.
         *
         * @param index the target exon index.
         * @param collection the containing exon collection.
         * @param <T> the exon type.
         * @return never {@code null}.
         */
        private static <T> String nameComposer(final int index, final ExonCollection<T> collection) {
            final String name = collection.name(collection.exon(index));
            return name == null ? NO_VALUE_STRING : name;
        }

        /**
         * Exon information string composer for the combined coordinate name part of the exon information.
         *
         * @param index the target exon index in the input collection.
         * @param collection the containing exon collection.
         * @param <T> the exon type.
         * @return never {@code null}.
         */
        private static <T> String coordinateAndNameComposer(final int index, final ExonCollection<T> collection) {
            return String.join(COLUMN_SEPARATOR,coordinateComposer(index,collection),nameComposer(index,collection));
        }

        /**
         * Holds a reference to a unmodifiable list of column names for this exon output info.
         */
        private final List<String> headerNames;


        /**
         * Common interface for the lamba responsible to compose the output
         * exon information given the exon index in the enclosing exon collection.
         */
        @FunctionalInterface
        private interface Composer {
            /**
             * Generates the exon output inforamtion string given its index in the
             * enclosing exon collection.
             * @param index the subject exon index in the collection.
             * @param exonCollection the enclosing exon collection.
             * @param <T> the exon data type.
             * @throws IllegalArgumentException if {@code index} is not a valid index in {@code exonCollection} or
             *   {@code exonCollection} is {@code null}.
             * @return never {@code null}.
             *
             */
            <T> String apply(final int index, ExonCollection<T> exonCollection);
        }

        /**
         * Holds a reference to the composer function for this exon-out-info instance.
         */
        private final Composer composer;

        /**
         * Creates a new exon-out-info enum value given the composer function reference
         * a the list of output column names.
         * @param composer the composer lambda reference.
         * @param headerNames list of info column names in the order these are going to be lay out
         *                    in output files.
         * @throws IllegalArgumentException if {@code composer} or {@code headerNames} are {@code null},
         *   or {@code headerNames} contains a {@code null}.
         */
        ExonOutInfo(final Composer composer, final String... headerNames) {
            this.composer = Utils.nonNull(composer, "the info string composer cannot be null");
            this.headerNames = Collections.unmodifiableList(Arrays.asList(Utils.nonNull(headerNames, "the header name list provided cannot be null")));
            if (this.headerNames.stream().anyMatch(Objects::isNull)) {
                throw new IllegalArgumentException("the input header-name cannot contain nulls");
            }
        }

        /**
         * Number of columns spanned by the exon-info.
         * @return 1 or greater.
         */
        protected int columnCount() {
            return headerNames.size();
        }

        /**
         * Returns the header names
         * @return never {@code null}, an unmodifiable list with the header names.
         */
        protected List<String> headerNames() {
            return headerNames;
        }

        /**
         * Returns a string with the tab-separated header column names.
         * @return never {@code null}.
         */
        protected String headerString() {
            return String.join(COLUMN_SEPARATOR,headerNames());
        }


        /**
         * Checks whether this exon output info choice requires that all exon have
         * a unique name.
         */
        protected boolean requiresUniqueExonName() {
            return this == NAME;
        }

        /**
         * Composes the exon information output string.
         *
         * @param index of the exon in the collection.
         * @param collection the exon containing collection.
         * @throws IllegalArgumentException if either {@code exon} or {@code collection} is {@code null}.
         */
        protected <T> String composeExonOutInfoString(final int index, final ExonCollection<T> collection) {
            Utils.nonNull(collection,"the collection cannot be null");
            Utils.validIndex(index, collection.exonCount());
            return composer.apply(index,collection);
        }
    }
}
