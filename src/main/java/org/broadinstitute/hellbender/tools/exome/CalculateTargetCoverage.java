package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Calculate read-counts across targets given their reference coordinates.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 *
 * <h3>Examples</h3>
 *
 * <p>
 *     The command encompasses empirically determined parameters for TCGA project data.
 *     You may obtain better results with different parameters.
 * </p>
 *
 * <p>For whole exome sequencing (WES) data: </p>
 *
 * <pre>
 * java -Xmx4g -jar $gatk_jar CalculateTargetCoverage \
 *   --input bam.bam \
 *   --targets padded_targets.tsv \
 *   --disableReadFilter NotDuplicateReadFilter \
 *   --output base_filename.coverage.tsv
 * </pre>
 *
 * <p>
 *     The interval targets are exome target intervals padded, e.g. with 250 bases on either side.
 *     Target intervals do NOT overlap. Use the PadTargets tool to generate non-overlapping padded intervals from exome targets.
 *     Do NOT use BED format. See ConvertBedToTargetFile.
 * </p>
 *
 * <p>For whole genome sequencing (WGS) data, use SparkGenomeReadCounts instead.</p>
 *
 */
@CommandLineProgramProperties(
        summary = "Count overlapping reads target by target",
        oneLineSummary = "Count overlapping reads target by target",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
public final class CalculateTargetCoverage extends ReadWalker {

    public static final String DEFAULT_COHORT_NAME = "<ALL>";

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

    public static final String COLUMN_SEPARATOR = "\t";
    public static final String LINE_SEPARATOR = "\n";
    public static final String NO_VALUE_STRING = ".";
    protected static final String GROUP_BY_FULL_NAME = "groupBy";
    protected static final String GROUP_BY_SHORT_NAME = GROUP_BY_FULL_NAME;
    protected static final String COLUMN_SUMMARY_OUTPUT_SHORT_NAME = "CSO";
    protected static final String COLUMN_SUMMARY_OUTPUT_FULL_NAME = "columnSummaryOutput";
    protected static final String ROW_SUMMARY_OUTPUT_SHORT_NAME = "RSO";
    protected static final String ROW_SUMMARY_OUTPUT_FULL_NAME = "rowSummaryOutput";
    protected static final String COHORT_FULL_NAME = "cohortName";
    protected static final String COHORT_SHORT_NAME = "cohort";
    protected static final String AVERAGE_DOUBLE_FORMAT = "%.4f";
    protected static final String TRANSFORM_FULL_NAME = "transform";
    protected static final String TRANSFORM_SHORT_NAME = TRANSFORM_FULL_NAME;
    protected static final String TARGET_FILE_FULL_NAME = "targets";
    protected static final String TARGET_FILE_SHORT_NAME = "T";
    protected static final String TARGET_OUT_INFO_FULL_NAME = "targetInformationColumns";
    protected static final String TARGET_OUT_INFO_SHORT_NAME = "targetInfo";

    private static final String PCOV_OUTPUT_DOUBLE_FORMAT = "%.4g";

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
    protected GroupBy groupBy = GroupBy.SAMPLE;

    @Argument(
            doc = "Cohort name used to name the read count column when the groupBy " +
                    "argument is set to COHORT",
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
    protected Transform transform = Transform.PCOV;

    @Argument(
            doc = "TSV file listing 1-based genomic intervals with specific column headers. Do NOT use BED format.",
            shortName = TARGET_FILE_SHORT_NAME,
            fullName = TARGET_FILE_FULL_NAME,
            optional = true
    )
    protected File targetsFile = null;

    @Argument(
            doc = "Target identification information to be showed in outputs",
            shortName = TARGET_OUT_INFO_SHORT_NAME,
            fullName = TARGET_OUT_INFO_FULL_NAME,
            optional = true
    )
    protected TargetOutInfo targetOutInfo = TargetOutInfo.FULL;

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
    private static final Logger logger = LogManager.getLogger(CalculateTargetCoverage.class);

    /**
     * Target database reference.
     */
    private TargetCollection<Target> targetCollection;

    /**
     * Counts table.
     */
    private CountColumns countColumns;

    /**
     * Count matrix indexed by count column and then target.
     */
    private int[][] counts;

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        final List<ReadFilter> filters = new ArrayList<>(super.getDefaultReadFilters());
        filters.add(ReadFilterLibrary.MAPPED);
        filters.add(ReadFilterLibrary.NON_ZERO_REFERENCE_LENGTH_ALIGNMENT);
        filters.add(ReadFilterLibrary.NOT_DUPLICATE);

        return filters;
    }

    @Override
    public void onTraversalStart() {

        sampleCollection = new SampleCollection(getHeaderForReads());

        logger.log(Level.INFO, "Reading targets locations from intervals...");

        targetCollection = resolveTargetCollection();

        // Initializing count and count column management member fields:
        countColumns = groupBy.countColumns(this);
        final int columnCount = countColumns.columnCount();
        counts = new int[columnCount][targetCollection.targetCount()];

        // Open output files and write headers:
        outputWriter = openOutputWriter(output, composeMatrixOutputHeader(getCommandLine(), targetOutInfo, groupBy, countColumns.columnNames()));
        if (columnSummaryOutput != null) {
            columnSummaryOutputWriter = openOutputWriter(columnSummaryOutput,
                    composeColumnSummaryHeader(getCommandLine(), groupBy, targetCollection.targetCount(), targetCollection.totalSize()));
        }
        if (rowSummaryOutput != null) {
            rowSummaryOutputWriter = openOutputWriter(rowSummaryOutput,
                    composeRowOutputHeader(getCommandLine(), targetOutInfo, groupBy, countColumns.columnCount()));
        }

        // Next we start the traversal:
        logger.log(Level.INFO, "Collecting read counts ...");
    }

    /**
     * Builds the target collection given the values of user arguments.
     *
     * @throws UserException if there is some inconsistency in user arguments and inputs.
     * @throws GATKException if there was any problem parsing the content of the targets file.
     * @return never {@code null}.
     */
    private TargetCollection<Target> resolveTargetCollection() {
        final TargetCollection<Target> result;
        if (targetsFile != null) {
            result = resolveTargetsFromFile();
        } else if (hasIntervals()) {
            final SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();
            final List<SimpleInterval> intervals = intervalArgumentCollection.getIntervals(sequenceDictionary);

            //this constructor automatically generates target names from the intervals
            final List<Target> targets = intervals.stream().map(Target::new).collect(Collectors.toList());
            result = new HashedListTargetCollection<>(targets);
        } else {
            throw new UserException(String.format("You must indicate the set of target as input intervals (e.g. -L target-intervals.list) or a target feature file (e.g. -%s my-targets.tsv) ", TARGET_FILE_SHORT_NAME));
        }
        if (targetOutInfo.requiresUniqueTargetName()) {
            checkAllTargetsHaveName(result);
        }
        return result;
    }

    /**
     * Checks whether all targets in the input target collection have a designated name.
     * @param result the query target collection.
     * @throws UserException if there are some targets with no name in {@code result}.
     */
    private void checkAllTargetsHaveName(TargetCollection<Target> result) {
        if (result.targets().stream().anyMatch(t -> t.getName() == null || t.getName().equals(""))) {
            throw new UserException(String.format("Target output info requested '%s' requires that each target has a designated unique name/id but there are some with no names: %s", targetOutInfo.name(),
                    result.targets().stream().filter(t -> t.getName() == null || t.getName().equals("")).limit(10).map(e -> result.location(e).toString()).collect(Collectors.joining(", "))));
        }
    }

    /**
     * Constructs the target collection from an target-file passed by the user.
     *
     * @return never {@code null}.
     */
    private TargetCollection<Target> resolveTargetsFromFile() {
        IOUtils.canReadFile(targetsFile);
        logger.log(Level.INFO,String.format("Reading target intervals from targets file '%s' ...", targetsFile.getAbsolutePath()));
        final List<Target> targets = TargetTableReader.readTargetFile(targetsFile);
        return new HashedListTargetCollection<>(targets);
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
            targetCollection.indexRange(readLocation).forEach(i -> counts[columnIndex][i]++);
        }
    }

    @Override
    public Object onTraversalSuccess() {
        logger.log(Level.INFO, "Collecting read counts done.");
        logger.log(Level.INFO, "Writing counts ...");
        final long[] columnTotals = calculateColumnTotals();

        IntStream.range(0, targetCollection.targetCount()).forEach(target -> {
            final int[] countBuffer = IntStream.range(0, counts.length).map(column -> counts[column][target]).toArray();
            writeOutputRows(countBuffer, columnTotals, target);
        });
        logger.log(Level.INFO, "Writing counts done.");

        writeColumnSummaryOutput();
        return "SUCCESS";
    }

    @Override
    public void closeTool() {
        if (columnSummaryOutputWriter != null) {
            columnSummaryOutputWriter.close();
        }
        if (outputWriter != null){
            outputWriter.close();
        }
        if (rowSummaryOutputWriter != null) {
            rowSummaryOutputWriter.close();
        }
        if (columnSummaryOutputWriter != null) {
            columnSummaryOutputWriter.close();
        }
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

        final long totalSize = targetCollection.totalSize();

        final List<String> columnNames = countColumns.columnNames();
        for (int i = 0; i < columnNames.size(); i++) {
            final long sum = IntStream.of(counts[i]).sum();
            columnSummaryOutputWriter.println(
                    String.join(COLUMN_SEPARATOR,
                            columnNames.get(i),
                            String.valueOf(sum),
                            String.format(AVERAGE_DOUBLE_FORMAT, sum / (double) totalSize)));
        }
    }

    /**
     * Writes the row in the main matrix output file for a target and, if requested,
     * the corresponding row in the row summary output file.
     *
     * @param countBuffer  the counts for the target.
     * @param index the index of target within the target collection.
     */
    private void writeOutputRows(final int[] countBuffer, final long[] columnTotals,
                                 final int index) {
        final String countString = IntStream.range(0, countBuffer.length).mapToObj(
                i -> transform.apply(countBuffer[i], columnTotals[i])).collect(Collectors.joining(COLUMN_SEPARATOR));
        final String targetInfoString = targetOutInfo.composeTargetOutInfoString(index, targetCollection);

        outputWriter.println(String.join(COLUMN_SEPARATOR, targetInfoString, countString));

        if (rowSummaryOutputWriter != null) {
            final long sum = MathUtils.sum(countBuffer);
            final SimpleInterval location = targetCollection.location(index);
            final int targetSize = location.size();
            rowSummaryOutputWriter.println(String.join(COLUMN_SEPARATOR,
                    targetInfoString,
                    Long.toString(sum), String.format(AVERAGE_DOUBLE_FORMAT,
                            sum / ((float) countColumns.columnCount() * targetSize))));
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
     * @param targetCount   number of targets processed.
     * @param totalSize   combined total size of all targets processed in bps.
     * @return never {@code null}.
     */
    private static String composeColumnSummaryHeader(final String commandLine, final GroupBy groupBy,
                                                     final int targetCount, final long totalSize) {
        return String.format(
                String.join(LINE_SEPARATOR,
                        "##fileFormat  = tsv",
                        "##commandLine = %s",
                        "##title       = Summary counts per %s",
                        "##metaData    = {",
                        "##    targetCount = %d,",
                        "##    totalSize = %d (bp)",
                        "##}",
                        String.join(COLUMN_SEPARATOR, groupBy.name(), SUM_COLUMN_NAME, AVG_BP_COLUMN_NAME)),
                commandLine, groupBy.toString(), targetCount, totalSize);
    }

    /**
     * Composes the row summary output header.
     *
     * @param commandLine the execution command line.
     * @param groupBy     the value of the group-by argument used.
     * @param columnCount number of count-column involved in the analysis.
     * @return never {@code null}.
     */
    private static String composeRowOutputHeader(final String commandLine, final TargetOutInfo targetOutInfo, final GroupBy groupBy,
                                                 final int columnCount) {
        return String.format(
                String.join(LINE_SEPARATOR,
                        "##fileFormat  = tsv",
                        "##commandLine = %s",
                        "##title       = Summary counts per target",
                        "##metaData = {",
                        "##    groupBy     = %s,",
                        "##    columnCount = %d",
                        "##}",
                        String.join(COLUMN_SEPARATOR, targetOutInfo.headerString(),SUM_COLUMN_NAME,AVG_COL_BP_COLUMN_NAME)),
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
    private static String composeMatrixOutputHeader(final String commandLine, final TargetOutInfo targetOutInfo, final GroupBy groupBy,
                                                    final List<String> countColumnNames) {
        final String countColumnHeaderString = String.join(COLUMN_SEPARATOR, countColumnNames).replace("%", "%%");
        final String formatString = String.join(LINE_SEPARATOR,
                "##fileFormat  = tsv",
                "##commandLine = %s",
                "##title       = Read counts per target and %s",
                String.join(COLUMN_SEPARATOR, targetOutInfo.headerString(), countColumnHeaderString));
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
        protected final CalculateTargetCoverage tool;

        /**
         * Creates a new count-counts instance given the reference to the embedding tool instance.
         *
         * @param tool the read-count tool instance.
         */
        private CountColumns(final CalculateTargetCoverage tool) {
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
        private CohortCountColumns(final CalculateTargetCoverage tool) {
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
        private SampleCountColumns(final CalculateTargetCoverage tool) {
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
        private ReadGroupCountColumns(final CalculateTargetCoverage tool) {
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
        private final Function<CalculateTargetCoverage, ? extends CountColumns> countColumnsFactory;

        /**
         * Creates a new group-by option.
         *
         * @param countColumnsFactory the count-columns manager factory to be used with
         *                            this group-by option. Assumed not to be {@code null}.
         */
        GroupBy(final Function<CalculateTargetCoverage, ? extends CountColumns> countColumnsFactory) {
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
        protected CountColumns countColumns(final CalculateTargetCoverage tool) {
            return countColumnsFactory.apply(tool);
        }
    }

    /**
     * Matrix value transformation function to apply to the main output matrix values.
     */
    protected enum Transform {

        /**
         * Read-count raw value (non-)transformation.
         */
        RAW((count, columnTotal) -> Integer.toString(count)),

        /**
         * Default proportional coverage transformation.
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
             * @param count       the individual count for a target and count group
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
            ParamUtils.isPositiveOrZero(count, "the count cannot less than 0");
            Utils.validateArg(count <= columnTotal, "the count cannot be larger than the column total");
            return operator.apply(count, columnTotal);
        }
    }

    /**
     * Possible per-target information included in the output.
     */
    protected enum TargetOutInfo {

        /**
         * Only target coordinates (contig name, start and end positions).
         */
        COORDS(TargetOutInfo::coordinateComposer,
                TargetTableColumn.CONTIG.toString(), TargetTableColumn.START.toString(), TargetTableColumn.END.toString()),

        /**
         * Only target name/id.
         */
        NAME(TargetOutInfo::nameComposer, TargetTableColumn.NAME.toString()),

        /**
         * Target coordinates and name.
         */
        FULL(TargetOutInfo::coordinateAndNameComposer,
                TargetTableColumn.CONTIG.toString(), TargetTableColumn.START.toString(), TargetTableColumn.END.toString(), TargetTableColumn.NAME.toString());

        /**
         * Target information string composer for the genomic coordinate part of the target.
         *
         * @param index      the index of a target within the collection.
         * @param collection the containing target collection.
         * @return never {@code null}.
         */
        private static String coordinateComposer(final int index, final TargetCollection<Target> collection) {
            final SimpleInterval location = collection.location(index);
            if (location == null) {
                return String.join(COLUMN_SEPARATOR, NO_VALUE_STRING, NO_VALUE_STRING, NO_VALUE_STRING);
            } else {
                return String.format(String.join(COLUMN_SEPARATOR, "%s", "%d", "%d"),
                        location.getContig(), location.getStart(), location.getEnd());
            }
        }

        /**
         * Target information string composer for the name part of the target information.
         *
         * @param index the target index.
         * @param collection the containing target collection.
         * @return never {@code null}.
         */
        private static String nameComposer(final int index, final TargetCollection<Target> collection) {
            final String name = collection.target(index).getName();
            return name == null ? NO_VALUE_STRING : name;
        }

        /**
         * Target information string composer for the combined coordinate name part of the target information.
         *
         * @param index the target index in the input collection.
         * @param collection the containing target collection.
         * @return never {@code null}.
         */
        private static String coordinateAndNameComposer(final int index, final TargetCollection<Target> collection) {
            return String.join(COLUMN_SEPARATOR,coordinateComposer(index,collection),nameComposer(index,collection));
        }

        /**
         * Holds a reference to a unmodifiable list of column names for this target output info.
         */
        private final List<String> headerNames;


        /**
         * Common interface for the lamba responsible to compose the output
         * target information given the target index in the enclosing target collection.
         */
        @FunctionalInterface
        private interface Composer {
            /**
             * Generates the target output inforamtion string given its index in the
             * enclosing target collection.
             * @param index the subject target index in the collection.
             * @param targetCollection the enclosing target collection.
             * @throws IllegalArgumentException if {@code index} is not a valid index in {@code targetCollection} or
             *   {@code targetCollection} is {@code null}.
             * @return never {@code null}.
             *
             */
            String apply(final int index, TargetCollection<Target> targetCollection);
        }

        /**
         * Holds a reference to the composer function for this target-out-info instance.
         */
        private final Composer composer;

        /**
         * Creates a new target-out-info enum value given the composer function reference
         * a the list of output column names.
         * @param composer the composer lambda reference.
         * @param headerNames list of info column names in the order these are going to be lay out
         *                    in output files.
         * @throws IllegalArgumentException if {@code composer} or {@code headerNames} are {@code null},
         *   or {@code headerNames} contains a {@code null}.
         */
        TargetOutInfo(final Composer composer, final String... headerNames) {
            this.composer = Utils.nonNull(composer, "the info string composer cannot be null");
            this.headerNames = Collections.unmodifiableList(Arrays.asList(Utils.nonNull(headerNames, "the header name list provided cannot be null")));
            Utils.validateArg(this.headerNames.stream().noneMatch(Objects::isNull), "the input header-name cannot contain nulls");
        }

        /**
         * Number of columns spanned by the target-info.
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
         * Checks whether this target output info choice requires that all target have
         * a unique name.
         */
        protected boolean requiresUniqueTargetName() {
            return this == NAME;
        }

        /**
         * Composes the target information output string.
         *
         * @param index of the target in the collection.
         * @param collection the target containing collection.
         * @throws IllegalArgumentException if either {@code target} or {@code collection} is {@code null}.
         */
        protected String composeTargetOutInfoString(final int index, final TargetCollection<Target> collection) {
            Utils.nonNull(collection,"the collection cannot be null");
            Utils.validIndex(index, collection.targetCount());
            return composer.apply(index,collection);
        }
    }
}
