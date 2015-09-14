package org.broadinstitute.hellbender.utils.recalibration;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import org.apache.commons.lang3.tuple.MutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.recalibration.covariates.Covariate;
import org.broadinstitute.hellbender.utils.recalibration.covariates.StandardCovariateList;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.R.RScriptExecutor;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.collections.NestedIntegerArray;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.report.GATKReport;
import org.broadinstitute.hellbender.utils.report.GATKReportTable;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.*;
import java.util.*;

/**
 * This helper class holds the data HashMap as well as submaps that represent the marginal distributions collapsed over all needed dimensions.
 * It also has static methods that are used to perform the various solid recalibration modes that attempt to correct the reference bias.
 * This class holds the parsing methods that are shared between BaseRecalibrator and PrintReads.
 */
public final class RecalUtils {
    public static final String ARGUMENT_REPORT_TABLE_TITLE = "Arguments";
    public static final String QUANTIZED_REPORT_TABLE_TITLE = "Quantized";
    public static final String READGROUP_REPORT_TABLE_TITLE = "RecalTable0";
    public static final String QUALITY_SCORE_REPORT_TABLE_TITLE = "RecalTable1";
    public static final String ALL_COVARIATES_REPORT_TABLE_TITLE = "RecalTable2";

    public static final String ARGUMENT_COLUMN_NAME = "Argument";
    public static final String ARGUMENT_VALUE_COLUMN_NAME = "Value";
    public static final String QUANTIZED_VALUE_COLUMN_NAME = "QuantizedScore";
    public static final String QUANTIZED_COUNT_COLUMN_NAME = "Count";
    public static final String READGROUP_COLUMN_NAME = "ReadGroup";
    public static final String EVENT_TYPE_COLUMN_NAME = "EventType";
    public static final String EMPIRICAL_QUALITY_COLUMN_NAME = "EmpiricalQuality";
    public static final String ESTIMATED_Q_REPORTED_COLUMN_NAME = "EstimatedQReported";
    public static final String QUALITY_SCORE_COLUMN_NAME = "QualityScore";
    public static final String COVARIATE_VALUE_COLUMN_NAME = "CovariateValue";
    public static final String COVARIATE_NAME_COLUMN_NAME = "CovariateName";
    public static final String NUMBER_OBSERVATIONS_COLUMN_NAME = "Observations";
    public static final String NUMBER_ERRORS_COLUMN_NAME = "Errors";

    private static final String COLOR_SPACE_ATTRIBUTE_TAG = "CS"; // The tag that holds the color space for SOLID bams
    private static final String COLOR_SPACE_INCONSISTENCY_TAG = "ZC"; // A new tag made up for the recalibrator which will hold an array of ints which say if this base is inconsistent with its color
    private static boolean warnUserNullPlatform = false;

    private static final String SCRIPT_FILE = "BQSR.R";

    private static final Pair<String, String> covariateValue     = new MutablePair<>(RecalUtils.COVARIATE_VALUE_COLUMN_NAME, "%s");
    private static final Pair<String, String> covariateName      = new MutablePair<>(RecalUtils.COVARIATE_NAME_COLUMN_NAME, "%s");
    private static final Pair<String, String> eventType          = new MutablePair<>(RecalUtils.EVENT_TYPE_COLUMN_NAME, "%s");
    private static final Pair<String, String> empiricalQuality   = new MutablePair<>(RecalUtils.EMPIRICAL_QUALITY_COLUMN_NAME, "%.4f");
    private static final Pair<String, String> estimatedQReported = new MutablePair<>(RecalUtils.ESTIMATED_Q_REPORTED_COLUMN_NAME, "%.4f");
    private static final Pair<String, String> nObservations      = new MutablePair<>(RecalUtils.NUMBER_OBSERVATIONS_COLUMN_NAME, "%d");
    private static final Pair<String, String> nErrors            = new MutablePair<>(RecalUtils.NUMBER_ERRORS_COLUMN_NAME, "%.2f");

    /**
     * Component used to print out csv representation of the reports that can be use to perform analysis in
     * external tools. E.g. generate plots using R scripts.
     * <p/>
     * A header is always printed into the output stream (or file) when the printer is created. Then you only need
     * to call {@link #print(RecalibrationReport, String) print} for each report you want to include in the csv file.
     * Once finished, you close the printer calling {@link #close() close}
     *
     */
    private static class CsvPrinter {

        private final PrintStream ps;
        private final StandardCovariateList covariates;

        /**
         * Constructs a printer redirected to an output file.
         * @param out the output file.
         * @param covs covariates to print out.
         * @throws FileNotFoundException if the file could not be created anew.
         */
        protected CsvPrinter(final File out, final StandardCovariateList covs)
                throws FileNotFoundException {
            this(new FileOutputStream(out), covs);
        }

        /**
         * Constructs a printer redirected to an output stream
         * @param os the output.
         * @param covs  covariates to print out.
         */
        protected CsvPrinter(final OutputStream os, final StandardCovariateList covs) {
            covariates = covs;
            ps = new PrintStream(os);
            printHeader();
        }

        /**
         * Prints the header out.
         * <p/>
         * Should only be invoked at creation.
         */
        protected void printHeader() {
            RecalUtils.printHeader(ps);
        }

        /**
         * Prints out a report into the csv file.
         *
         *
         * @param report the report to print out.
         * @param mode  the report associated mode. (typically ORIGINAL, RECALIBRATED
         */
        public void print(final RecalibrationReport report, final String mode) {
            RecalUtils.writeCsv(ps, report.getRecalibrationTables(), mode, covariates, false);
        }

        /**
         * Close the csv printer.
         *
         * No further output will be allowed or take place after calling this method.
         */
        public void close() {
            ps.close();
        }

    }

    /**
     * Returns a csv output printer.
     *
     * @param out the output file. It will be overridden
     * @param covs list of covariates to print out.
     *
     * @throws FileNotFoundException if <code>out</code> could not be created anew.
     *
     * @return never <code>null</code>
     */
    protected static CsvPrinter csvPrinter(final File out, final StandardCovariateList covs)
        throws FileNotFoundException
    {
        if (covs == null) {
            throw new IllegalArgumentException("the input covariate array cannot be null");
        }
        return new CsvPrinter(out,covs);
    }

    /**
     * Prints out a collection of reports into a file in Csv format in a way
     * that can be used by R scripts (such as the plot generator script).
     * <p/>
     * The set of covariates is take as the minimum common set from all reports.
     *
     * @param out the output file. It will be overridden.
     * @param reports map where keys are the unique 'mode' (ORIGINAL, RECALIBRATED, ...)
     *                of each report and the corresponding value the report itself.
     * @throws FileNotFoundException if <code>out</code> could not be created anew.
     */
    public static void generateCsv(final File out, final Map<String, RecalibrationReport> reports) throws FileNotFoundException {
        if (reports.isEmpty()) {
            throw new GATKException("no reports");
        }
        final RecalibrationReport firstReport = reports.values().iterator().next();
        final StandardCovariateList covariates = firstReport.getCovariates();
        writeCsv(out, reports, covariates);
    }

    /**
     * Print out a collection of reports into a file in Csv format in a way
     * that can be used by R scripts (such as the plot generator script).
     *
     * @param out
     * @param reports map where keys are the unique 'mode' (ORIGINAL, RECALIBRATED, ...)
     *                of each report and the corresponding value the report itself.
     * @param covs the covariates to print out.
     * @throws FileNotFoundException if <code>out</code> could not be created anew.
     */
    private static void writeCsv(final File out, final Map<String, RecalibrationReport> reports, final StandardCovariateList covs)
        throws FileNotFoundException {
        final CsvPrinter p = csvPrinter(out, covs);
        for (final Map.Entry<String,RecalibrationReport> e : reports.entrySet()) {
            p.print(e.getValue(),e.getKey());
        }
        p.close();
    }

    public enum SOLID_RECAL_MODE {
        /**
         * Treat reference inserted bases as reference matching bases. Very unsafe!
         */
        DO_NOTHING,
        /**
         * Set reference inserted bases and the previous base (because of color space alignment details) to Q0. This is the default option.
         */
        SET_Q_ZERO,
        /**
         * In addition to setting the quality scores to zero, also set the base itself to 'N'. This is useful to visualize in IGV.
         */
        SET_Q_ZERO_BASE_N,
        /**
         * Look at the color quality scores and probabilistically decide to change the reference inserted base to be the base which is implied by the original color space instead of the reference.
         */
        REMOVE_REF_BIAS;
        
        public static SOLID_RECAL_MODE recalModeFromString(String recalMode) {
            if (recalMode.equals("DO_NOTHING"))
                return SOLID_RECAL_MODE.DO_NOTHING;
            if (recalMode.equals("SET_Q_ZERO"))
                return SOLID_RECAL_MODE.SET_Q_ZERO;
            if (recalMode.equals("SET_Q_ZERO_BASE_N"))
                return SOLID_RECAL_MODE.SET_Q_ZERO_BASE_N;
            if (recalMode.equals("REMOVE_REF_BIAS"))
                return SOLID_RECAL_MODE.REMOVE_REF_BIAS;

            throw new UserException.BadArgumentValue(recalMode, "is not a valid SOLID_RECAL_MODE value");
        }
    }

    public enum SOLID_NOCALL_STRATEGY {
        /**
         * When a no call is detected throw an exception to alert the user that recalibrating this SOLiD data is unsafe. This is the default option.
         */
        THROW_EXCEPTION,
        /**
         * Leave the read in the output bam completely untouched. This mode is only okay if the no calls are very rare.
         */
        LEAVE_READ_UNRECALIBRATED,
        /**
         * Mark these reads as failing vendor quality checks so they can be filtered out by downstream analyses.
         */
        PURGE_READ;

        public static SOLID_NOCALL_STRATEGY nocallStrategyFromString(String nocallStrategy) {
            if (nocallStrategy.equals("THROW_EXCEPTION"))
                return SOLID_NOCALL_STRATEGY.THROW_EXCEPTION;
            if (nocallStrategy.equals("LEAVE_READ_UNRECALIBRATED"))
                return SOLID_NOCALL_STRATEGY.LEAVE_READ_UNRECALIBRATED;
            if (nocallStrategy.equals("PURGE_READ"))
                return SOLID_NOCALL_STRATEGY.PURGE_READ;

            throw new UserException.BadArgumentValue(nocallStrategy, "is not a valid SOLID_NOCALL_STRATEGY value");
        }
    }

    private static List<GATKReportTable> generateReportTables(final RecalibrationTables recalibrationTables, final StandardCovariateList covariates, boolean sortByCols) {
        List<GATKReportTable> result = new LinkedList<>();
        int rowIndex = 0;

        GATKReportTable allCovsReportTable = null;

        for (NestedIntegerArray<RecalDatum> table : recalibrationTables){
            final ArrayList<Pair<String, String>> columnNames = new ArrayList<>(); // initialize the array to hold the column names
            columnNames.add(new MutablePair<>(covariates.getReadGroupCovariate().parseNameForReport(), "%s")); // save the required covariate name so we can reference it in the future
            if (!recalibrationTables.isReadGroupTable(table)) {
                columnNames.add(new MutablePair<>(covariates.getQualityScoreCovariate().parseNameForReport(), "%s")); // save the required covariate name so we can reference it in the future
                if (recalibrationTables.isAdditionalCovariateTable(table)) {
                    columnNames.add(covariateValue);
                    columnNames.add(covariateName);
                }
            }

            columnNames.add(eventType); // the order of these column names is important here
            columnNames.add(empiricalQuality);
            if (recalibrationTables.isReadGroupTable(table)) {
                columnNames.add(estimatedQReported); // only the read group table needs the estimated Q reported
            }
            columnNames.add(nObservations);
            columnNames.add(nErrors);

            final String reportTableName = getReportTableName(recalibrationTables, table);

            final GATKReportTable.TableSortingWay sort = sortByCols ? GATKReportTable.TableSortingWay.SORT_BY_COLUMN : GATKReportTable.TableSortingWay.DO_NOT_SORT;

            final GATKReportTable reportTable;
            final boolean addToList;

            //XXX this "if" implicitly uses the knowledge about the ordering of tables.
            if (!recalibrationTables.isAdditionalCovariateTable(table)) {
                reportTable = makeNewTableWithColumns(columnNames, reportTableName, sort);
                rowIndex = 0; // reset the row index since we're starting with a new table
                addToList = true;
            } else if (allCovsReportTable == null && recalibrationTables.isAdditionalCovariateTable(table)) {
                reportTable = makeNewTableWithColumns(columnNames, reportTableName, sort);
                rowIndex = 0; // reset the row index since we're starting with a new table
                allCovsReportTable =  reportTable;
                addToList = true;
            } else {
                reportTable = allCovsReportTable;
                addToList = false;
            }

            for (final NestedIntegerArray.Leaf<RecalDatum> row : table.getAllLeaves()) {
                final RecalDatum datum = row.value;
                final int[] keys = row.keys;

                int columnIndex = 0;
                int keyIndex = 0;
                reportTable.set(rowIndex, columnNames.get(columnIndex++).getLeft(), covariates.getReadGroupCovariate().formatKey(keys[keyIndex++]));
                if (!recalibrationTables.isReadGroupTable(table)) {
                    reportTable.set(rowIndex, columnNames.get(columnIndex++).getLeft(), covariates.getQualityScoreCovariate().formatKey(keys[keyIndex++]));
                    if (recalibrationTables.isAdditionalCovariateTable(table)){
                        final Covariate covariate = recalibrationTables.getCovariateForTable(table);
                        reportTable.set(rowIndex, columnNames.get(columnIndex++).getLeft(), covariate.formatKey(keys[keyIndex++]));
                        reportTable.set(rowIndex, columnNames.get(columnIndex++).getLeft(), covariate.parseNameForReport());
                    }
                }

                final EventType event = EventType.eventFrom(keys[keyIndex]);
                reportTable.set(rowIndex, columnNames.get(columnIndex++).getLeft(), event.toString());

                reportTable.set(rowIndex, columnNames.get(columnIndex++).getLeft(), datum.getEmpiricalQuality());
                if (recalibrationTables.isReadGroupTable(table)) {
                    reportTable.set(rowIndex, columnNames.get(columnIndex++).getLeft(), datum.getEstimatedQReported()); // we only add the estimated Q reported in the RG table
                }
                reportTable.set(rowIndex, columnNames.get(columnIndex++).getLeft(), datum.getNumObservations());
                reportTable.set(rowIndex, columnNames.get(columnIndex).getLeft(), datum.getNumMismatches());

                rowIndex++;
            }
            if (addToList){
                result.add(reportTable);      //XXX using a set would be slow because the equals method on GATKReportTable is expensive.
            }
        }

        return result;
    }

    private static String getReportTableName(RecalibrationTables recalibrationTables, NestedIntegerArray<RecalDatum> table) {
        final String reportTableName;
        if (recalibrationTables.isReadGroupTable(table)){
            return READGROUP_REPORT_TABLE_TITLE;
        } else if (recalibrationTables.isQualityScoreTable(table)){
            return QUALITY_SCORE_REPORT_TABLE_TITLE;
        } else {
            return ALL_COVARIATES_REPORT_TABLE_TITLE;
        }
    }

    private static GATKReportTable makeNewTableWithColumns(ArrayList<Pair<String, String>> columnNames, String reportTableName, GATKReportTable.TableSortingWay sort) {
        GATKReportTable rt = new GATKReportTable(reportTableName, "", columnNames.size(), sort);
        for (final Pair<String, String> columnName : columnNames) {
            rt.addColumn(columnName.getLeft(), columnName.getRight());
        }
        return rt;
    }

    /**
     * Outputs the GATK report to RAC.RECAL_TABLE.
     *
     * @param RAC The list of shared command line arguments
     * @param quantizationInfo Quantization info
     * @param recalibrationTables Recalibration tables
     * @param covariates The list of requested covariates
     * @param sortByCols True to use GATKReportTable.TableSortingWay.SORT_BY_COLUMN, false to use GATKReportTable.TableSortingWay.DO_NOT_SORT
     */
    public static void outputRecalibrationReport(final PrintStream recalTableStream, final RecalibrationArgumentCollection RAC, final QuantizationInfo quantizationInfo, final RecalibrationTables recalibrationTables, final StandardCovariateList covariates, boolean sortByCols) {
        final GATKReport report = createRecalibrationGATKReport(RAC.generateReportTable(covariates.covariateNames()), quantizationInfo.generateReportTable(sortByCols), generateReportTables(recalibrationTables, covariates, sortByCols));
        report.print(recalTableStream);
    }

    /**
     * Creates a consolidated GATK report, first generating report tables. Report can then be written to a stream via GATKReport.print(PrintStream).
     *
     * @param argumentTable Argument table
     * @param quantizationInfo Quantization info
     * @param recalibrationTables Recalibration tables
     * @param covariates The list of covariates
     * @param sortByCols True to use GATKReportTable.TableSortingWay.SORT_BY_COLUMN, false to use GATKReportTable.TableSortingWay.DO_NOT_SORT
     * @return GATK report
     */
    public static GATKReport createRecalibrationGATKReport(final GATKReportTable argumentTable, final QuantizationInfo quantizationInfo, final RecalibrationTables recalibrationTables, final StandardCovariateList covariates, final boolean sortByCols) {
        return createRecalibrationGATKReport(argumentTable, quantizationInfo.generateReportTable(sortByCols), generateReportTables(recalibrationTables, covariates, sortByCols));
    }

    /**
     * Creates a consolidated GATK report from the tables. Report can then be written to a stream via GATKReport.print(PrintStream).
     *
     * @param argumentTable Argument table
     * @param quantizationTable Quantization Table
     * @param recalTables Other recal tables
     * @return GATK report
     */
    private static GATKReport createRecalibrationGATKReport(final GATKReportTable argumentTable, final GATKReportTable quantizationTable, final List<GATKReportTable> recalTables) {
        final GATKReport report = new GATKReport();
        report.addTable(argumentTable);
        report.addTable(quantizationTable);
        report.addTables(recalTables);
        return report;
    }

    /**                                                s
     * Write recalibration plots into a file
     *
     * @param csvFile location of the intermediary file
     * @param maybeGzipedExampleReportFile where the report arguments are collected from.
     * @param output result plot file name.
     */
    public static void generatePlots(final File csvFile, final File maybeGzipedExampleReportFile, final File output) {
        final File exampleReportFile = IOUtils.gunzipToTempIfNeeded(maybeGzipedExampleReportFile);
        final RScriptExecutor executor = new RScriptExecutor();
        executor.addScript(new Resource(SCRIPT_FILE, RecalUtils.class));
        executor.addArgs(csvFile.getAbsolutePath());
        executor.addArgs(exampleReportFile.getAbsolutePath());
        executor.addArgs(output.getAbsolutePath());
        LogManager.getLogger(RecalUtils.class).debug("R command line: " + executor.getApproximateCommandLine());
        executor.exec();
    }

    private static void writeCsv(final PrintStream deltaTableFile, final RecalibrationTables recalibrationTables, final String recalibrationMode, final StandardCovariateList covariates, final boolean printHeader) {

        final NestedIntegerArray<RecalDatum> deltaTable = createDeltaTable(recalibrationTables, covariates.size());

        // add the quality score table to the delta table
        final NestedIntegerArray<RecalDatum> qualTable = recalibrationTables.getQualityScoreTable();
        for (final NestedIntegerArray.Leaf<RecalDatum> leaf : qualTable.getAllLeaves()) { // go through every element in the covariates table to create the delta table
            final int[] newCovs = new int[4];
            newCovs[0] = leaf.keys[0];
            newCovs[1] = covariates.size(); // replace the covariate name with an arbitrary (unused) index for QualityScore. This is a HACK.
            newCovs[2] = leaf.keys[1];
            newCovs[3] = leaf.keys[2];
            addToDeltaTable(deltaTable, newCovs, leaf.value); // add this covariate to the delta table
        }

        // add the optional covariates to the delta table
        for(final NestedIntegerArray<RecalDatum> covTable : recalibrationTables.getAdditionalTables()){
            for (final NestedIntegerArray.Leaf<RecalDatum> leaf : covTable.getAllLeaves()) {
                final int[] covs = new int[4];
                covs[0] = leaf.keys[0];
                covs[1] = covariates.indexByClass(recalibrationTables.getCovariateForTable(covTable).getClass()); // reset the quality score covariate to 0 from the keyset (so we aggregate all rows regardless of QS)
                covs[2] = leaf.keys[2];
                covs[3] = leaf.keys[3];
                addToDeltaTable(deltaTable, covs, leaf.value); // add this covariate to the delta table
            }
        }

        // output the csv file
        if (printHeader) {
            printHeader(deltaTableFile);
        }

        // print each data line
        for (final NestedIntegerArray.Leaf<RecalDatum> leaf : deltaTable.getAllLeaves()) {
            final List<Object> deltaKeys = generateValuesFromKeys(leaf.keys, covariates);
            final RecalDatum deltaDatum = leaf.value;
            deltaTableFile.print(Utils.join(",", deltaKeys));
            deltaTableFile.print(',' + deltaDatum.stringForCSV());
            deltaTableFile.println(',' + recalibrationMode);
        }
    }

    private static void printHeader(PrintStream out) {
        final List<String> header = new LinkedList<>();
        header.add("ReadGroup");
        header.add("CovariateValue");
        header.add("CovariateName");
        header.add("EventType");
        header.add("Observations");
        header.add("Errors");
        header.add("EmpiricalQuality");
        header.add("AverageReportedQuality");
        header.add("Accuracy");
        header.add("Recalibration");
        out.println(Utils.join(",", header));
    }

    /*
     * Return an initialized nested integer array with appropriate dimensions for use with the delta tables
     *
     * @param recalibrationTables     the recal tables
     * @param numCovariates           the total number of covariates being used
     * @return a non-null nested integer array
     */
    private static NestedIntegerArray<RecalDatum> createDeltaTable(final RecalibrationTables recalibrationTables, final int numCovariates) {
        final int[] dimensionsForDeltaTable = new int[4];

        // initialize the dimensions with those of the qual table to start with
        final NestedIntegerArray<RecalDatum> qualTable = recalibrationTables.getQualityScoreTable();
        final int[] dimensionsOfQualTable = qualTable.getDimensions();
        dimensionsForDeltaTable[0] = dimensionsOfQualTable[0];    // num read groups
        dimensionsForDeltaTable[1] = numCovariates + 1;           // num covariates
        dimensionsForDeltaTable[2] = dimensionsOfQualTable[1];
        dimensionsForDeltaTable[3] = dimensionsOfQualTable[2];

        // now, update the dimensions based on the optional covariate tables as needed
        for ( final NestedIntegerArray<RecalDatum> covTable : recalibrationTables.getAdditionalTables()) {
            final int[] dimensionsOfCovTable = covTable.getDimensions();
            dimensionsForDeltaTable[2] = Math.max(dimensionsForDeltaTable[2], dimensionsOfCovTable[2]);
            dimensionsForDeltaTable[3] = Math.max(dimensionsForDeltaTable[3], dimensionsOfCovTable[3]);
        }

        return new NestedIntegerArray<>(dimensionsForDeltaTable);
    }

    static List<Object> generateValuesFromKeys(final int[] keys, final StandardCovariateList covariates) {
        final List<Object> values = new ArrayList<>(4);
        values.add(covariates.getReadGroupCovariate().formatKey(keys[0]));

        final int covariateIndex = keys[1];
        final int covariateKey = keys[2];
        final Covariate covariate = (covariateIndex == covariates.size()) ? covariates.getQualityScoreCovariate() : covariates.get(covariateIndex);
        values.add(covariate.formatKey(covariateKey));
        values.add(covariate.parseNameForReport());
        values.add(EventType.eventFrom(keys[3]).prettyPrint());

        return values;
    }

    /**
     * Updates the current RecalDatum element in the delta table.
     *
     * If it doesn't have an element yet, it creates an RecalDatum element and adds it to the delta table.
     *
     * @param deltaTable the delta table
     * @param deltaKey the key to the table
     * @param recalDatum the recal datum to combine with the accuracyDatum element in the table
     */
    private static void addToDeltaTable(final NestedIntegerArray<RecalDatum> deltaTable, final int[] deltaKey, final RecalDatum recalDatum) {
        final RecalDatum deltaDatum = deltaTable.get(deltaKey); // check if we already have a RecalDatum for this key
        if (deltaDatum == null)
            // if we don't have a key yet, create a new one with the same values as the current datum
            deltaTable.put(new RecalDatum(recalDatum), deltaKey);
        else
            // if we do have a datum, combine it with this one
            deltaDatum.combine(recalDatum);
    }

    /**
     * Section of code shared between the two recalibration walkers which uses the command line arguments to adjust attributes of the read such as quals or platform string
     *
     * @param read The read to adjust
     * @param RAC  The list of shared command line arguments
     */
    public static void parsePlatformForRead(final GATKRead read, final SAMFileHeader header, final RecalibrationArgumentCollection RAC) {
        SAMReadGroupRecord readGroup = ReadUtils.getSAMReadGroupRecord(read, header);

        if (RAC.FORCE_PLATFORM != null && (readGroup.getPlatform() == null || !readGroup.getPlatform().equals(RAC.FORCE_PLATFORM))) {
            readGroup.setPlatform(RAC.FORCE_PLATFORM);
        }

        if (readGroup.getPlatform() == null) {
            if (RAC.DEFAULT_PLATFORM != null) {
                if (!warnUserNullPlatform) {
                    Utils.warnUser("The input .bam file contains reads with no platform information. " +
                            "Defaulting to platform = " + RAC.DEFAULT_PLATFORM + ". " +
                            "First observed at read with name = " + read.getName());
                    warnUserNullPlatform = true;
                }
                readGroup.setPlatform(RAC.DEFAULT_PLATFORM);
            }
            else {
                throw new UserException.MalformedRead(read, "The input .bam file contains reads with no platform information. First observed at read with name = " + read.getName());
            }
        }
    }

    private static boolean hasNoCallInColorSpace(final byte[] colorSpace) {
        final int length = colorSpace.length;
        for (int i = 1; i < length; i++) {  // skip the sentinal
            final byte color = colorSpace[i];
            if (color != (byte) '0' && color != (byte) '1' && color != (byte) '2' && color != (byte) '3') {
                return true; // There is a bad color in this SOLiD read
            }
        }

        return false; // There aren't any color no calls in this SOLiD read
    }

    /**
     * Given the base and the color calculate the next base in the sequence
     *
     * @param read     the read
     * @param prevBase The base
     * @param color    The color
     * @return The next base in the sequence
     */
    private static byte getNextBaseFromColor(GATKRead read, final byte prevBase, final byte color) {
        switch (color) {
            case '0':
                return prevBase;
            case '1':
                return performColorOne(prevBase);
            case '2':
                return performColorTwo(prevBase);
            case '3':
                return performColorThree(prevBase);
            default:
                throw new UserException.MalformedRead(read, "Unrecognized color space in SOLID read, color = " + (char) color +
                        " Unfortunately this bam file can not be recalibrated without full color space information because of potential reference bias.");
        }
    }

    /**
     * Parse through the color space of the read and add a new tag to the read that says which bases are
     * inconsistent with the color space. If there is a no call in the color space, this method returns false meaning
     * this read should be skipped
     *
     * @param strategy the strategy used for SOLID no calls
     * @param read     The GATKRead to parse
     * @param header   SAM header for the read
     * @return true if this read is consistent or false if this read should be skipped
     */
    public static boolean isColorSpaceConsistent(final SOLID_NOCALL_STRATEGY strategy, final GATKRead read, final SAMFileHeader header) {
        if (!ReadUtils.isSOLiDRead(read, header)) // If this is a SOLID read then we have to check if the color space is inconsistent. This is our only sign that SOLID has inserted the reference base
            return true;

        // Haven't calculated the inconsistency array yet for this read
        if (! read.hasAttribute(RecalUtils.COLOR_SPACE_INCONSISTENCY_TAG)) {
            if (read.hasAttribute(RecalUtils.COLOR_SPACE_ATTRIBUTE_TAG)) {
                byte[] colorSpace = read.getAttributeAsString(RecalUtils.COLOR_SPACE_ATTRIBUTE_TAG).getBytes();

                final boolean badColor = hasNoCallInColorSpace(colorSpace);
                if (badColor) {
                    if (strategy == SOLID_NOCALL_STRATEGY.LEAVE_READ_UNRECALIBRATED) {
                        return false; // can't recalibrate a SOLiD read with no calls in the color space, and the user wants to skip over them
                    }
                    else if (strategy == SOLID_NOCALL_STRATEGY.PURGE_READ) {
                        read.setFailsVendorQualityCheck(true);
                        return false;
                    }
                }

                byte[] readBases = read.getBases(); // Loop over the read and calculate first the inferred bases from the color and then check if it is consistent with the read
                if (read.isReverseStrand())
                    readBases = BaseUtils.simpleReverseComplement(read.getBases());

                final byte[] inconsistency = new byte[readBases.length];
                int i;
                byte prevBase = colorSpace[0]; // The sentinel
                for (i = 0; i < readBases.length; i++) {
                    final byte thisBase = getNextBaseFromColor(read, prevBase, colorSpace[i + 1]);
                    inconsistency[i] = (byte) (thisBase == readBases[i] ? 0 : 1);
                    prevBase = readBases[i];
                }
                read.setAttribute(RecalUtils.COLOR_SPACE_INCONSISTENCY_TAG, inconsistency);
            }
            else if (strategy == SOLID_NOCALL_STRATEGY.THROW_EXCEPTION) // if the strategy calls for an exception, throw it
                throw new UserException.MalformedRead(read, "Unable to find color space information in SOLiD read. First observed at read with name = " + read.getName() + " Unfortunately this .bam file can not be recalibrated without color space information because of potential reference bias.");

            else
                return false; // otherwise, just skip the read
        }

        return true;
    }

    /**
     * Check if this base is inconsistent with its color space. If it is then SOLID inserted the reference here and we should reduce the quality
     *
     * @param read   The read which contains the color space to check against
     * @param offset The offset in the read at which to check
     * @return Returns true if the base was inconsistent with the color space
     */
    public static boolean isColorSpaceConsistent(final GATKRead read, final int offset) {
        if (read.hasAttribute(RecalUtils.COLOR_SPACE_INCONSISTENCY_TAG)) {
            final byte[] inconsistency = read.getAttributeAsByteArray(RecalUtils.COLOR_SPACE_INCONSISTENCY_TAG);
            // NOTE: The inconsistency array is in the direction of the read, not aligned to the reference!
            if (read.isReverseStrand()) { // Negative direction
                return inconsistency[inconsistency.length - offset - 1] == (byte) 0;
            }
            else { // Forward direction
                return inconsistency[offset] == (byte) 0;
            }

            // This block of code is for if you want to check both the offset and the next base for color space inconsistency
            //if( read.isReverseStrand() ) { // Negative direction
            //    if( offset == 0 ) {
            //        return inconsistency[0] != 0;
            //    } else {
            //        return (inconsistency[inconsistency.length - offset - 1] != 0) || (inconsistency[inconsistency.length - offset] != 0);
            //    }
            //} else { // Forward direction
            //    if( offset == inconsistency.length - 1 ) {
            //        return inconsistency[inconsistency.length - 1] != 0;
            //    } else {
            //        return (inconsistency[offset] != 0) || (inconsistency[offset + 1] != 0);
            //    }
            //}

        }
        else { // No inconsistency array, so nothing is inconsistent
            return true;
        }
    }

    /**
     * Computes all requested covariates for every offset in the given read
     * by calling covariate.getValues(..).
     *
     * It populates an array of covariate values where result[i][j] is the covariate
     * value for the ith position in the read and the jth covariate in
     * reqeustedCovariates list.
     *
     * @param read                The read for which to compute covariate values.
     * @param header              SAM header for the read
     * @param covariates The list of requested covariates.
     * @return a matrix with all the covariates calculated for every base in the read
     */
    public static ReadCovariates computeCovariates(final GATKRead read, final SAMFileHeader header, final StandardCovariateList covariates) {
        final ReadCovariates readCovariates = new ReadCovariates(read.getLength(), covariates.size());
        computeCovariates(read, header, covariates, readCovariates);
        return readCovariates;
    }

    /**
     * Computes all requested covariates for every offset in the given read
     * by calling covariate.getValues(..).
     *
     * It populates an array of covariate values where result[i][j] is the covariate
     * value for the ith position in the read and the jth covariate in
     * covariates list.
     *
     * @param read                The read for which to compute covariate values.
     * @param header              SAM header for the read
     * @param covariates          The list of covariates.
     * @param resultsStorage      The object to store the covariate values
     */
    public static void computeCovariates(final GATKRead read, final SAMFileHeader header, final StandardCovariateList covariates, final ReadCovariates resultsStorage) {
        covariates.recordAllValuesInStorage(read, header, resultsStorage);
    }

    /**
     * Perform a certain transversion (A <-> C or G <-> T) on the base.
     *
     * @param base the base [AaCcGgTt]
     * @return the transversion of the base, or the input base if it's not one of the understood ones
     */
    private static byte performColorOne(byte base) {
        switch (base) {
            case 'A':
            case 'a':
                return 'C';
            case 'C':
            case 'c':
                return 'A';
            case 'G':
            case 'g':
                return 'T';
            case 'T':
            case 't':
                return 'G';
            default:
                return base;
        }
    }

    /**
     * Perform a transition (A <-> G or C <-> T) on the base.
     *
     * @param base the base [AaCcGgTt]
     * @return the transition of the base, or the input base if it's not one of the understood ones
     */
    private static byte performColorTwo(byte base) {
        switch (base) {
            case 'A':
            case 'a':
                return 'G';
            case 'C':
            case 'c':
                return 'T';
            case 'G':
            case 'g':
                return 'A';
            case 'T':
            case 't':
                return 'C';
            default:
                return base;
        }
    }

    /**
     * Return the complement (A <-> T or C <-> G) of a base.
     *
     * @param base the base [AaCcGgTt]
     * @return the complementary base, or the input base if it's not one of the understood ones
     */
    private static byte performColorThree(byte base) {
        switch (base) {
            case 'A':
            case 'a':
                return 'T';
            case 'C':
            case 'c':
                return 'G';
            case 'G':
            case 'g':
                return 'C';
            case 'T':
            case 't':
                return 'A';
            default:
                return base;
        }
    }

    /**
     * Combines the recalibration data for table1 and table2 into table1
     *
     * Note that table1 is the destination, so it is modified
     *
     * @param table1 the destination table to merge table2 into
     * @param table2 the source table to merge into table1
     */
    public static void combineTables(final NestedIntegerArray<RecalDatum> table1, final NestedIntegerArray<RecalDatum> table2) {
        if ( table1 == null ) throw new IllegalArgumentException("table1 cannot be null");
        if ( table2 == null ) throw new IllegalArgumentException("table2 cannot be null");
        if ( ! Arrays.equals(table1.getDimensions(), table2.getDimensions()))
            throw new IllegalArgumentException("Table1 " + Utils.join(",", table1.getDimensions()) + " not equal to " + Utils.join(",", table2.getDimensions()));

        for (final NestedIntegerArray.Leaf<RecalDatum> row : table2.getAllLeaves()) {
            final RecalDatum myDatum = table1.get(row.keys);

            if (myDatum == null)
                table1.put(row.value, row.keys);
            else
                myDatum.combine(row.value);
        }
    }

    /**
     * Increments the RecalDatum at the specified position in the specified table, or put a new item there
     * if there isn't already one.
     *
     * @param table the table that holds/will hold our item
     * @param qual qual for this event
     * @param isError error value for this event
     * @param keys location in table of our item
     */
    public static void incrementDatumOrPutIfNecessary( final NestedIntegerArray<RecalDatum> table,
                                                          final byte qual,
                                                          final double isError,
                                                          final int... keys ) {
        final RecalDatum existingDatum = table.get(keys);

        if ( existingDatum == null ) {
            // No existing item, put a new one
            table.put(createDatumObject(qual, isError), keys);
        } else {
            // Easy case: already an item here, so increment it
            existingDatum.increment(1L, isError);
        }
    }

    /**
     * creates a datum object with one observation and one or zero error
     *
     * @param reportedQual  the quality score reported by the instrument for this base
     * @param isError       whether or not the observation is an error
     * @return a new RecalDatum object with the observation and the error
     */
    private static RecalDatum createDatumObject(final byte reportedQual, final double isError) {
        return new RecalDatum(1, isError, reportedQual);
    }
}
