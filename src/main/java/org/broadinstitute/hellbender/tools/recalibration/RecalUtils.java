package org.broadinstitute.hellbender.tools.recalibration;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import org.apache.commons.lang3.tuple.MutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.recalibration.covariates.*;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.R.RScriptExecutor;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.collections.NestedIntegerArray;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.recalibration.EventType;
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
    public final static String ARGUMENT_REPORT_TABLE_TITLE = "Arguments";
    public final static String QUANTIZED_REPORT_TABLE_TITLE = "Quantized";
    public final static String READGROUP_REPORT_TABLE_TITLE = "RecalTable0";
    public final static String QUALITY_SCORE_REPORT_TABLE_TITLE = "RecalTable1";
    public final static String ALL_COVARIATES_REPORT_TABLE_TITLE = "RecalTable2";

    public final static String ARGUMENT_COLUMN_NAME = "Argument";
    public final static String ARGUMENT_VALUE_COLUMN_NAME = "Value";
    public final static String QUANTIZED_VALUE_COLUMN_NAME = "QuantizedScore";
    public static final String QUANTIZED_COUNT_COLUMN_NAME = "Count";
    public final static String READGROUP_COLUMN_NAME = "ReadGroup";
    public final static String EVENT_TYPE_COLUMN_NAME = "EventType";
    public final static String EMPIRICAL_QUALITY_COLUMN_NAME = "EmpiricalQuality";
    public final static String ESTIMATED_Q_REPORTED_COLUMN_NAME = "EstimatedQReported";
    public final static String QUALITY_SCORE_COLUMN_NAME = "QualityScore";
    public final static String COVARIATE_VALUE_COLUMN_NAME = "CovariateValue";
    public final static String COVARIATE_NAME_COLUMN_NAME = "CovariateName";
    public final static String NUMBER_OBSERVATIONS_COLUMN_NAME = "Observations";
    public final static String NUMBER_ERRORS_COLUMN_NAME = "Errors";

    private final static String COLOR_SPACE_ATTRIBUTE_TAG = "CS"; // The tag that holds the color space for SOLID bams
    private final static String COLOR_SPACE_INCONSISTENCY_TAG = "ZC"; // A new tag made up for the recalibrator which will hold an array of ints which say if this base is inconsistent with its color
    private static boolean warnUserNullPlatform = false;

    private static final String SCRIPT_FILE = "BQSR.R";

    private static final Pair<String, String> covariateValue     = new MutablePair<>(RecalUtils.COVARIATE_VALUE_COLUMN_NAME, "%s");
    private static final Pair<String, String> covariateName      = new MutablePair<>(RecalUtils.COVARIATE_NAME_COLUMN_NAME, "%s");
    private static final Pair<String, String> eventType          = new MutablePair<>(RecalUtils.EVENT_TYPE_COLUMN_NAME, "%s");
    private static final Pair<String, String> empiricalQuality   = new MutablePair<>(RecalUtils.EMPIRICAL_QUALITY_COLUMN_NAME, "%.4f");
    private static final Pair<String, String> estimatedQReported = new MutablePair<>(RecalUtils.ESTIMATED_Q_REPORTED_COLUMN_NAME, "%.4f");
    private static final Pair<String, String> nObservations      = new MutablePair<>(RecalUtils.NUMBER_OBSERVATIONS_COLUMN_NAME, "%d");
    private static final Pair<String, String> nErrors            = new MutablePair<>(RecalUtils.NUMBER_ERRORS_COLUMN_NAME, "%.2f");

    public static List<Class<? extends Covariate>> getStandardCovariateClasses(){
        return Arrays.asList(ContextCovariate.class, CycleCovariate.class);
    }

    public static List<Class<? extends Covariate>> getRequiredCovariateClasses(){
        return Arrays.asList(ReadGroupCovariate.class, QualityScoreCovariate.class);
    }

    public static List<Class<? extends Covariate>> getExperimentalCovariateClasses() {
        return Arrays.asList(RepeatLengthCovariate.class, RepeatUnitAndLengthCovariate.class, RepeatUnitCovariate.class);
    }

    public static List<Class<? extends Covariate>> getAllCovariateClasses(){
        List<Class<? extends Covariate>> ALL = new ArrayList<>();
        ALL.addAll(getStandardCovariateClasses());
        ALL.addAll(getRequiredCovariateClasses());
        ALL.addAll(getExperimentalCovariateClasses());
        return Collections.unmodifiableList(ALL);
    }
    /**
     * Generates two lists : required covariates and optional covariates based on the user's requests.
     *
     * Performs the following tasks in order:
     *  1. Adds all requierd covariates in order
     *  2. Check if the user asked to use the standard covariates and adds them all if that's the case
     *  3. Adds all covariates requested by the user that were not already added by the two previous steps
     *
     * @param argumentCollection the argument collection object for the recalibration walker
     * @return a pair of ordered lists : required covariates (first) and optional covariates (second)
     */
    public static Pair<ArrayList<Covariate>, ArrayList<Covariate>> initializeCovariates(RecalibrationArgumentCollection argumentCollection) {
        final List<Class<? extends Covariate>> covariateClasses = getAllCovariateClasses();
        final List<Class<? extends Covariate>> requiredClasses = getRequiredCovariateClasses();
        final List<Class<? extends Covariate>> standardClasses = getStandardCovariateClasses();

        final ArrayList<Covariate> requiredCovariates = addRequiredCovariatesToList(requiredClasses); // add the required covariates
        ArrayList<Covariate> optionalCovariates = new ArrayList<>();
        if (!argumentCollection.DO_NOT_USE_STANDARD_COVARIATES)
            optionalCovariates = addStandardCovariatesToList(standardClasses); // add the standard covariates if -standard was specified by the user

        // parse the -cov arguments that were provided, skipping over the ones already specified
        if (argumentCollection.COVARIATES != null) {
            for (String requestedCovariateString : argumentCollection.COVARIATES) {
                // help the transition from BQSR v1 to BQSR v2
                if ( requestedCovariateString.equals("DinucCovariate") )
                    throw new UserException.CommandLineException("DinucCovariate has been retired.  Please use its successor covariate " +
                            "ContextCovariate instead, which includes the 2 bp (dinuc) substitution model of the retired DinucCovariate " +
                            "as well as an indel context to model the indel error rates");

                boolean foundClass = false;
                for (Class<? extends Covariate> covClass : covariateClasses) {
                    if (requestedCovariateString.equalsIgnoreCase(covClass.getSimpleName())) { // -cov argument matches the class name for an implementing class
                        foundClass = true;
                        if (!requiredClasses.contains(covClass) &&
                                (argumentCollection.DO_NOT_USE_STANDARD_COVARIATES || !standardClasses.contains(covClass))) {
                            try {
                                final Covariate covariate = covClass.newInstance(); // now that we've found a matching class, try to instantiate it
                                optionalCovariates.add(covariate);
                            } catch (Exception e) {
                                throw new UserException.DynamicClassResolutionException(covClass, e);
                            }
                        }
                    }
                }

                if (!foundClass) {
                    throw new UserException.CommandLineException("The requested covariate type (" + requestedCovariateString + ") isn't a valid covariate option. Use --list to see possible covariates.");
                }
            }
        }
        return new MutablePair<>(requiredCovariates, optionalCovariates);
    }

    /**
     * Adds the required covariates to a covariate list
     *
     * Note: this method really only checks if the classes object has the expected number of required covariates, then add them by hand.
     *
     * @param classes list of classes to add to the covariate list
     * @return the covariate list
     */
    private static ArrayList<Covariate> addRequiredCovariatesToList(List<Class<? extends Covariate>> classes) {
        //FIXME this is bogus: this method is not checking the classes in the List it's passed at all
        //FIXME (just the number of classes). It then proceeds to add a new ReadGroupCovariate()
        //FIXME and new QualityScoreCovariate() unconditionally.
        //FIXME Should either be a parameterless method, or should actually
        //FIXME check and instantiate the classes in the List it's passed.
        ArrayList<Covariate> dest = new ArrayList<>(classes.size());
        if (classes.size() != 2)
            throw new GATKException("The number of required covariates has changed, this is a hard change in the code and needs to be inspected");

        dest.add(new ReadGroupCovariate()); // enforce the order with RG first and QS next.
        dest.add(new QualityScoreCovariate());
        return dest;
    }

    /**
     * Adds the standard covariates to a covariate list
     *
     * @param classes list of classes to add to the covariate list
     * @return the covariate list
     */
    private static ArrayList<Covariate> addStandardCovariatesToList(List<Class<? extends Covariate>> classes) {
        ArrayList<Covariate> dest = new ArrayList<>(classes.size());
        for (Class<?> covClass : classes) {
            try {
                final Covariate covariate = (Covariate) covClass.newInstance();
                dest.add(covariate);
            } catch (Exception e) {
                throw new UserException.DynamicClassResolutionException(covClass, e);
            }
        }
        return dest;
    }

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
        private final Covariate[] covariates;

        /**
         * Constructs a printer redirected to an output file.
         * @param out the output file.
         * @param c covariates to print out.
         * @throws java.io.FileNotFoundException if the file could not be created anew.
         */
        protected CsvPrinter(final File out, final Covariate ... c)
                throws FileNotFoundException {
            this(new FileOutputStream(out), c);
        }

        /**
         * Constructs a printer redirected to an output stream
         * @param os the output.
         * @param c  covariates to print out.
         */
        protected CsvPrinter(final OutputStream os, final Covariate ... c) {
            covariates = c == null ? new Covariate[0] : Arrays.copyOf(c, c.length);
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
            RecalUtils.writeCSV(ps,report.getRecalibrationTables(),mode,covariates,false);
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
     * @param c list of covariates to print out.
     *
     * @throws java.io.FileNotFoundException if <code>out</code> could not be created anew.
     *
     * @return never <code>null</code>
     */
    protected static CsvPrinter csvPrinter(final File out, final Covariate ... c)
        throws FileNotFoundException
    {
        if (c == null) {
            throw new IllegalArgumentException("the input covariate array cannot be null");
        }
        return new CsvPrinter(out,c);
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
     * @throws java.io.FileNotFoundException if <code>out</code> could not be created anew.
     */
    public static void generateCsv(final File out, final Map<String, RecalibrationReport> reports)
            throws FileNotFoundException {
        if (reports.size() == 0) {
            writeCsv(out, reports, new Covariate[0]);
        } else {
            final Iterator<RecalibrationReport> rit = reports.values().iterator();
            final RecalibrationReport first = rit.next();
            final Covariate[] firstCovariates = first.getRequestedCovariates();
            final Set<Covariate> covariates = new LinkedHashSet<>();
            Utils.addAll(covariates,firstCovariates);
            while (rit.hasNext() && covariates.size() > 0) {
                final Covariate[] nextCovariates = rit.next().getRequestedCovariates();
                final Set<String> nextCovariateNames = new LinkedHashSet<>(nextCovariates.length);
                for (final Covariate nc : nextCovariates) {
                    nextCovariateNames.add(nc.getClass().getSimpleName());
                }
                final Iterator<Covariate> cit = covariates.iterator();
                while (cit.hasNext()) {
                    if (!nextCovariateNames.contains(cit.next().getClass().getSimpleName())) {
                        cit.remove();
                    }
                }
            }
            writeCsv(out, reports, covariates.toArray(new Covariate[covariates.size()]));
        }
    }

    /**
     * Print out a collection of reports into a file in Csv format in a way
     * that can be used by R scripts (such as the plot generator script).
     *
     * @param out
     * @param reports map where keys are the unique 'mode' (ORIGINAL, RECALIBRATED, ...)
     *                of each report and the corresponding value the report itself.
     * @param c the covariates to print out.
     * @throws java.io.FileNotFoundException if <code>out</code> could not be created anew.
     */
    private static void writeCsv(final File out,
            final Map<String, RecalibrationReport> reports, final Covariate[] c)
        throws FileNotFoundException {
        final CsvPrinter p = csvPrinter(out,c);
        for (Map.Entry<String,RecalibrationReport> e : reports.entrySet()) {
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

    private static List<GATKReportTable> generateReportTables(final RecalibrationTables recalibrationTables, final Covariate[] requestedCovariates, boolean sortByCols) {
        List<GATKReportTable> result = new LinkedList<>();
        int reportTableIndex = 0;
        int rowIndex = 0;
        final Map<Covariate, String> covariateNameMap = new HashMap<>(requestedCovariates.length);
        for (final Covariate covariate : requestedCovariates)
            covariateNameMap.put(covariate, parseCovariateName(covariate));

        for (int tableIndex = 0; tableIndex < recalibrationTables.numTables(); tableIndex++) {

            final ArrayList<Pair<String, String>> columnNames = new ArrayList<>(); // initialize the array to hold the column names
            columnNames.add(new MutablePair<>(covariateNameMap.get(requestedCovariates[0]), "%s")); // save the required covariate name so we can reference it in the future
            if (tableIndex != RecalibrationTables.TableType.READ_GROUP_TABLE.ordinal()) {
                columnNames.add(new MutablePair<>(covariateNameMap.get(requestedCovariates[1]), "%s")); // save the required covariate name so we can reference it in the future
                if (tableIndex >= RecalibrationTables.TableType.OPTIONAL_COVARIATE_TABLES_START.ordinal()) {
                    columnNames.add(covariateValue);
                    columnNames.add(covariateName);
                }
            }

            columnNames.add(eventType); // the order of these column names is important here
            columnNames.add(empiricalQuality);
            if (tableIndex == RecalibrationTables.TableType.READ_GROUP_TABLE.ordinal())
                columnNames.add(estimatedQReported); // only the read group table needs the estimated Q reported
            columnNames.add(nObservations);
            columnNames.add(nErrors);

            final GATKReportTable reportTable;
            if (tableIndex <= RecalibrationTables.TableType.OPTIONAL_COVARIATE_TABLES_START.ordinal()) {
                if(sortByCols) {
                    reportTable = new GATKReportTable("RecalTable" + reportTableIndex++, "", columnNames.size(), GATKReportTable.TableSortingWay.SORT_BY_COLUMN);
                } else {
                    reportTable = new GATKReportTable("RecalTable" + reportTableIndex++, "", columnNames.size(), GATKReportTable.TableSortingWay.DO_NOT_SORT);
                }
                for (final Pair<String, String> columnName : columnNames)
                    reportTable.addColumn(columnName.getLeft(), columnName.getRight());
                rowIndex = 0; // reset the row index since we're starting with a new table
            } else {
                reportTable = result.get(RecalibrationTables.TableType.OPTIONAL_COVARIATE_TABLES_START.ordinal());
            }

            final NestedIntegerArray<RecalDatum> table = recalibrationTables.getTable(tableIndex);
            for (final NestedIntegerArray.Leaf<RecalDatum> row : table.getAllLeaves()) {
                final RecalDatum datum = row.value;
                final int[] keys = row.keys;

                int columnIndex = 0;
                int keyIndex = 0;
                reportTable.set(rowIndex, columnNames.get(columnIndex++).getLeft(), requestedCovariates[0].formatKey(keys[keyIndex++]));
                if (tableIndex != RecalibrationTables.TableType.READ_GROUP_TABLE.ordinal()) {
                    reportTable.set(rowIndex, columnNames.get(columnIndex++).getLeft(), requestedCovariates[1].formatKey(keys[keyIndex++]));
                    if (tableIndex >= RecalibrationTables.TableType.OPTIONAL_COVARIATE_TABLES_START.ordinal()) {
                        final Covariate covariate = requestedCovariates[tableIndex];

                        reportTable.set(rowIndex, columnNames.get(columnIndex++).getLeft(), covariate.formatKey(keys[keyIndex++]));
                        reportTable.set(rowIndex, columnNames.get(columnIndex++).getLeft(), covariateNameMap.get(covariate));
                    }
                }

                final EventType event = EventType.eventFrom(keys[keyIndex]);
                reportTable.set(rowIndex, columnNames.get(columnIndex++).getLeft(), event.toString());

                reportTable.set(rowIndex, columnNames.get(columnIndex++).getLeft(), datum.getEmpiricalQuality());
                if (tableIndex == RecalibrationTables.TableType.READ_GROUP_TABLE.ordinal())
                    reportTable.set(rowIndex, columnNames.get(columnIndex++).getLeft(), datum.getEstimatedQReported()); // we only add the estimated Q reported in the RG table
                reportTable.set(rowIndex, columnNames.get(columnIndex++).getLeft(), datum.getNumObservations());
                reportTable.set(rowIndex, columnNames.get(columnIndex).getLeft(), datum.getNumMismatches());

                rowIndex++;
            }
            result.add(reportTable);
        }

        return result;
    }

    private static String parseCovariateName(final Covariate covariate) {
        return covariate.getClass().getSimpleName().split("Covariate")[0];
    }

    /**
     * Return a human-readable string representing the used covariates
     *
     * @param requestedCovariates a vector of covariates
     * @return a non-null comma-separated string
     */
    public static String covariateNames(final Covariate[] requestedCovariates) {
        final List<String> names = new ArrayList<>(requestedCovariates.length);
        for ( final Covariate cov : requestedCovariates )
            names.add(cov.getClass().getSimpleName());
        return Utils.join(",", names);
    }

    /**
     * Outputs the GATK report to RAC.RECAL_TABLE.
     *
     * @param RAC The list of shared command line arguments
     * @param quantizationInfo Quantization info
     * @param recalibrationTables Recalibration tables
     * @param requestedCovariates The list of requested covariates
     * @param sortByCols True to use GATKReportTable.TableSortingWay.SORT_BY_COLUMN, false to use GATKReportTable.TableSortingWay.DO_NOT_SORT
     */
    public static void outputRecalibrationReport(final RecalibrationArgumentCollection RAC, final QuantizationInfo quantizationInfo, final RecalibrationTables recalibrationTables, final Covariate[] requestedCovariates, boolean sortByCols) {
        final GATKReport report = createRecalibrationGATKReport(RAC.generateReportTable(covariateNames(requestedCovariates)), quantizationInfo.generateReportTable(sortByCols), generateReportTables(recalibrationTables, requestedCovariates, sortByCols));
        report.print(RAC.RECAL_TABLE);
    }

    /**
     * Creates a consolidated GATK report, first generating report tables. Report can then be written to a stream via GATKReport.print(PrintStream).
     *
     * @param argumentTable Argument table
     * @param quantizationInfo Quantization info
     * @param recalibrationTables Recalibration tables
     * @param requestedCovariates The list of requested covariates
     * @param sortByCols True to use GATKReportTable.TableSortingWay.SORT_BY_COLUMN, false to use GATKReportTable.TableSortingWay.DO_NOT_SORT
     * @return GATK report
     */
    public static GATKReport createRecalibrationGATKReport(final GATKReportTable argumentTable, final QuantizationInfo quantizationInfo, final RecalibrationTables recalibrationTables, final Covariate[] requestedCovariates, final boolean sortByCols) {
        return createRecalibrationGATKReport(argumentTable, quantizationInfo.generateReportTable(sortByCols), generateReportTables(recalibrationTables, requestedCovariates, sortByCols));
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

    private static void writeCSV(final PrintStream deltaTableFile, final RecalibrationTables recalibrationTables, final String recalibrationMode, final Covariate[] requestedCovariates, final boolean printHeader) {

        final NestedIntegerArray<RecalDatum> deltaTable = createDeltaTable(recalibrationTables, requestedCovariates.length);

        // add the quality score table to the delta table
        final NestedIntegerArray<RecalDatum> qualTable = recalibrationTables.getQualityScoreTable();
        for (final NestedIntegerArray.Leaf<RecalDatum> leaf : qualTable.getAllLeaves()) { // go through every element in the covariates table to create the delta table
            final int[] newCovs = new int[4];
            newCovs[0] = leaf.keys[0];
            newCovs[1] = requestedCovariates.length; // replace the covariate name with an arbitrary (unused) index for QualityScore
            newCovs[2] = leaf.keys[1];
            newCovs[3] = leaf.keys[2];
            addToDeltaTable(deltaTable, newCovs, leaf.value); // add this covariate to the delta table
        }

        // add the optional covariates to the delta table
        for (int i = RecalibrationTables.TableType.OPTIONAL_COVARIATE_TABLES_START.ordinal(); i < requestedCovariates.length; i++) {
            final NestedIntegerArray<RecalDatum> covTable = recalibrationTables.getTable(i);
            for (final NestedIntegerArray.Leaf<RecalDatum> leaf : covTable.getAllLeaves()) {
                final int[] covs = new int[4];
                covs[0] = leaf.keys[0];
                covs[1] = i; // reset the quality score covariate to 0 from the keyset (so we aggregate all rows regardless of QS)
                covs[2] = leaf.keys[2];
                covs[3] = leaf.keys[3];
                addToDeltaTable(deltaTable, covs, leaf.value); // add this covariate to the delta table
            }
        }

        // output the csv file
        if (printHeader) {
            printHeader(deltaTableFile);
        }

        final Map<Covariate, String> covariateNameMap = new HashMap<>(requestedCovariates.length);
        for (final Covariate covariate : requestedCovariates)
            covariateNameMap.put(covariate, parseCovariateName(covariate));

        // print each data line
        for (final NestedIntegerArray.Leaf<RecalDatum> leaf : deltaTable.getAllLeaves()) {
            final List<Object> deltaKeys = generateValuesFromKeys(leaf.keys, requestedCovariates, covariateNameMap);
            final RecalDatum deltaDatum = leaf.value;
            deltaTableFile.print(Utils.join(",", deltaKeys));
            deltaTableFile.print("," + deltaDatum.stringForCSV());
            deltaTableFile.println("," + recalibrationMode);
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
        for ( int i = RecalibrationTables.TableType.OPTIONAL_COVARIATE_TABLES_START.ordinal(); i < numCovariates; i++ ) {
            final NestedIntegerArray<RecalDatum> covTable = recalibrationTables.getTable(i);
            final int[] dimensionsOfCovTable = covTable.getDimensions();
            dimensionsForDeltaTable[2] = Math.max(dimensionsForDeltaTable[2], dimensionsOfCovTable[2]);
            dimensionsForDeltaTable[3] = Math.max(dimensionsForDeltaTable[3], dimensionsOfCovTable[3]);
        }

        return new NestedIntegerArray<>(dimensionsForDeltaTable);
    }

    protected static List<Object> generateValuesFromKeys(final int[] keys, final Covariate[] covariates, final Map<Covariate, String> covariateNameMap) {
        final List<Object> values = new ArrayList<>(4);
        values.add(covariates[RecalibrationTables.TableType.READ_GROUP_TABLE.ordinal()].formatKey(keys[0]));

        final int covariateIndex = keys[1];
        final int covariateKey = keys[2];
        final Covariate covariate = covariateIndex == covariates.length ? covariates[RecalibrationTables.TableType.QUALITY_SCORE_TABLE.ordinal()] : covariates[covariateIndex];
        values.add(covariate.formatKey(covariateKey));
        values.add(covariateNameMap.get(covariate));
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
    public static void parsePlatformForRead(final SAMRecord read, final RecalibrationArgumentCollection RAC) {
        SAMReadGroupRecord readGroup = read.getReadGroup();

        if (RAC.FORCE_PLATFORM != null && (readGroup.getPlatform() == null || !readGroup.getPlatform().equals(RAC.FORCE_PLATFORM))) {
            readGroup.setPlatform(RAC.FORCE_PLATFORM);
        }

        if (readGroup.getPlatform() == null) {
            if (RAC.DEFAULT_PLATFORM != null) {
                if (!warnUserNullPlatform) {
                    Utils.warnUser("The input .bam file contains reads with no platform information. " +
                            "Defaulting to platform = " + RAC.DEFAULT_PLATFORM + ". " +
                            "First observed at read with name = " + read.getReadName());
                    warnUserNullPlatform = true;
                }
                readGroup.setPlatform(RAC.DEFAULT_PLATFORM);
            }
            else {
                throw new UserException.MalformedBAM(read, "The input .bam file contains reads with no platform information. First observed at read with name = " + read.getReadName());
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
    private static byte getNextBaseFromColor(SAMRecord read, final byte prevBase, final byte color) {
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
                throw new UserException.MalformedBAM(read, "Unrecognized color space in SOLID read, color = " + (char) color +
                        " Unfortunately this bam file can not be recalibrated without full color space information because of potential reference bias.");
        }
    }

    /**
     * Parse through the color space of the read and add a new tag to the SAMRecord that says which bases are
     * inconsistent with the color space. If there is a no call in the color space, this method returns false meaning
     * this read should be skipped
     *
     * @param strategy the strategy used for SOLID no calls
     * @param read     The SAMRecord to parse
     * @return true if this read is consistent or false if this read should be skipped
     */
    public static boolean isColorSpaceConsistent(final SOLID_NOCALL_STRATEGY strategy, final SAMRecord read) {
        if (!ReadUtils.isSOLiDRead(read)) // If this is a SOLID read then we have to check if the color space is inconsistent. This is our only sign that SOLID has inserted the reference base
            return true;

        // Haven't calculated the inconsistency array yet for this read
        if (read.getAttribute(RecalUtils.COLOR_SPACE_INCONSISTENCY_TAG) == null) {
            final Object attr = read.getAttribute(RecalUtils.COLOR_SPACE_ATTRIBUTE_TAG);
            if (attr != null) {
                byte[] colorSpace;
                if (attr instanceof String)
                    colorSpace = ((String) attr).getBytes();
                else
                    throw new UserException.MalformedBAM(read, String.format("Value encoded by %s in %s isn't a string!", RecalUtils.COLOR_SPACE_ATTRIBUTE_TAG, read.getReadName()));

                final boolean badColor = hasNoCallInColorSpace(colorSpace);
                if (badColor) {
                    if (strategy == SOLID_NOCALL_STRATEGY.LEAVE_READ_UNRECALIBRATED) {
                        return false; // can't recalibrate a SOLiD read with no calls in the color space, and the user wants to skip over them
                    }
                    else if (strategy == SOLID_NOCALL_STRATEGY.PURGE_READ) {
                        read.setReadFailsVendorQualityCheckFlag(true);
                        return false;
                    }
                }

                byte[] readBases = read.getReadBases(); // Loop over the read and calculate first the inferred bases from the color and then check if it is consistent with the read
                if (read.getReadNegativeStrandFlag())
                    readBases = BaseUtils.simpleReverseComplement(read.getReadBases());

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
                throw new UserException.MalformedBAM(read, "Unable to find color space information in SOLiD read. First observed at read with name = " + read.getReadName() + " Unfortunately this .bam file can not be recalibrated without color space information because of potential reference bias.");

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
    public static boolean isColorSpaceConsistent(final SAMRecord read, final int offset) {
        final Object attr = read.getAttribute(RecalUtils.COLOR_SPACE_INCONSISTENCY_TAG);
        if (attr != null) {
            final byte[] inconsistency = (byte[]) attr;
            // NOTE: The inconsistency array is in the direction of the read, not aligned to the reference!
            if (read.getReadNegativeStrandFlag()) { // Negative direction
                return inconsistency[inconsistency.length - offset - 1] == (byte) 0;
            }
            else { // Forward direction
                return inconsistency[offset] == (byte) 0;
            }

            // This block of code is for if you want to check both the offset and the next base for color space inconsistency
            //if( read.getReadNegativeStrandFlag() ) { // Negative direction
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
     * @param requestedCovariates The list of requested covariates.
     * @return a matrix with all the covariates calculated for every base in the read
     */
    public static ReadCovariates computeCovariates(final SAMRecord read, final Covariate[] requestedCovariates) {
        final ReadCovariates readCovariates = new ReadCovariates(read.getReadLength(), requestedCovariates.length);
        computeCovariates(read, requestedCovariates, readCovariates);
        return readCovariates;
    }

    /**
     * Computes all requested covariates for every offset in the given read
     * by calling covariate.getValues(..).
     *
     * It populates an array of covariate values where result[i][j] is the covariate
     * value for the ith position in the read and the jth covariate in
     * requestedCovariates list.
     *
     * @param read                The read for which to compute covariate values.
     * @param requestedCovariates The list of requested covariates.
     * @param resultsStorage      The object to store the covariate values
     */
    public static void computeCovariates(final SAMRecord read, final Covariate[] requestedCovariates, final ReadCovariates resultsStorage) {
        // Loop through the list of requested covariates and compute the values of each covariate for all positions in this read
        for (int i = 0; i < requestedCovariates.length; i++) {
            resultsStorage.setCovariateIndex(i);
            requestedCovariates[i].recordValues(read, resultsStorage);
        }
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
