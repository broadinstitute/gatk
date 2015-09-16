package org.broadinstitute.hellbender.tools.recalibration;


import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.recalibration.covariates.Covariate;
import org.broadinstitute.hellbender.tools.recalibration.covariates.StandardCovariateList;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.collections.NestedIntegerArray;
import org.broadinstitute.hellbender.utils.recalibration.EventType;
import org.broadinstitute.hellbender.utils.report.GATKReport;
import org.broadinstitute.hellbender.utils.report.GATKReportTable;

import java.io.File;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

/**
 * This class has all the static functionality for reading a recalibration report file into memory. 
 */
public final class RecalibrationReport {
    private QuantizationInfo quantizationInfo; // histogram containing the counts for qual quantization (calculated after recalibration is done)
    private final RecalibrationTables recalibrationTables; // quick access reference to the tables
    private final StandardCovariateList covariates; // list of all covariates to be used in this calculation

    private final GATKReportTable argumentTable; // keep the argument table untouched just for output purposes
    private final RecalibrationArgumentCollection RAC; // necessary for quantizing qualities with the same parameter

    public RecalibrationReport(final File recalFile) {
        this(new GATKReport(recalFile));
    }

    public RecalibrationReport(final InputStream recalibrationTableStream){
        this(new GATKReport(recalibrationTableStream));
    }

    public RecalibrationReport(final GATKReport report){
        this(report, getReadGroups(report));
    }

    public RecalibrationReport(final GATKReport report, final SortedSet<String> allReadGroups) {
        argumentTable = report.getTable(RecalUtils.ARGUMENT_REPORT_TABLE_TITLE);
        RAC = initializeArgumentCollectionTable(argumentTable);

        GATKReportTable quantizedTable = report.getTable(RecalUtils.QUANTIZED_REPORT_TABLE_TITLE);
        quantizationInfo = initializeQuantizationTable(quantizedTable);

        covariates = new StandardCovariateList(RAC, new ArrayList<>(allReadGroups));

        recalibrationTables = new RecalibrationTables(covariates, allReadGroups.size());

        initializeReadGroupCovariates(allReadGroups);

        parseReadGroupTable(report.getTable(RecalUtils.READGROUP_REPORT_TABLE_TITLE), recalibrationTables.getReadGroupTable());

        parseQualityScoreTable(report.getTable(RecalUtils.QUALITY_SCORE_REPORT_TABLE_TITLE), recalibrationTables.getQualityScoreTable());

        parseAllCovariatesTable(report.getTable(RecalUtils.ALL_COVARIATES_REPORT_TABLE_TITLE), recalibrationTables);

    }

    /**
     * Gets the unique read groups in the table
     *
     * @param report the GATKReport containing the table with RecalUtils.READGROUP_REPORT_TABLE_TITLE
     * @return the unique read groups
     */
    private static SortedSet<String> getReadGroups(final GATKReport report) {
        final GATKReportTable reportTable = report.getTable(RecalUtils.READGROUP_REPORT_TABLE_TITLE);
        final SortedSet<String> readGroups = new TreeSet<>();
        for ( int i = 0; i < reportTable.getNumRows(); i++ ) {
            readGroups.add(reportTable.get(i, RecalUtils.READGROUP_COLUMN_NAME).toString());
        }
        return readGroups;
    }

    /**
    * Combines two recalibration reports by adding all observations and errors
    *
    * Note: This method DOES NOT recalculate the empirical qualities and quantized qualities. You have to recalculate
    * them after combining. The reason for not calculating it is because this function is intended for combining a
    * series of recalibration reports, and it only makes sense to calculate the empirical qualities and quantized
    * qualities after all the recalibration reports have been combined. Having the user recalculate when appropriate,
    * makes this method faster
    *
    * Note2: The empirical quality reported, however, is recalculated given its simplicity.
    *
    * @param other the recalibration report to combine with this one
    */
    public void combine(final RecalibrationReport other) {
        for ( int tableIndex = 0; tableIndex < recalibrationTables.numTables(); tableIndex++ ) {
            final NestedIntegerArray<RecalDatum> myTable = recalibrationTables.getTable(tableIndex);
            final NestedIntegerArray<RecalDatum> otherTable = other.recalibrationTables.getTable(tableIndex);
            RecalUtils.combineTables(myTable, otherTable);
        }
    }

    public QuantizationInfo getQuantizationInfo() {
        return quantizationInfo;
    }

    public RecalibrationTables getRecalibrationTables() {
        return recalibrationTables;
    }

    public StandardCovariateList getCovariates() {
        return covariates;
    }

    /**
     * Initialize read group keys using the shared list of all the read groups.
     *
     * By using the same sorted set of read groups across all recalibration reports, even if
     * one report is missing a read group, all the reports use the same read group keys.
     *
     * @param allReadGroups The list of all possible read groups
     */
    private void initializeReadGroupCovariates(final SortedSet<String> allReadGroups) {
        for (final String readGroup: allReadGroups) {
            covariates.getReadGroupCovariate().keyFromValue(readGroup);
        }
    }

    /**
     * Compiles the list of keys for the Covariates table and uses the shared parsing utility to produce the actual table
     *
     * @param reportTable            the GATKReport table containing data for this table
     * @param recalibrationTables    the recalibration tables
     */
    private void parseAllCovariatesTable(final GATKReportTable reportTable, final RecalibrationTables recalibrationTables) {
        final int[] tempCOVarray = new int[4];

        for ( int i = 0; i < reportTable.getNumRows(); i++ ) {
            final Object rg = reportTable.get(i, RecalUtils.READGROUP_COLUMN_NAME);
            tempCOVarray[0] = covariates.getReadGroupCovariate().keyFromValue(rg);
            final Object qual = reportTable.get(i, RecalUtils.QUALITY_SCORE_COLUMN_NAME);
            tempCOVarray[1] = covariates.getQualityScoreCovariate().keyFromValue(qual);

            final String covName = (String)reportTable.get(i, RecalUtils.COVARIATE_NAME_COLUMN_NAME);
            final Object covValue = reportTable.get(i, RecalUtils.COVARIATE_VALUE_COLUMN_NAME);
            final Covariate cov = covariates.getCovariateByParsedName(covName);
            tempCOVarray[2] = cov.keyFromValue(covValue);

            final EventType event = EventType.eventFrom((String)reportTable.get(i, RecalUtils.EVENT_TYPE_COLUMN_NAME));
            tempCOVarray[3] = event.ordinal();

            recalibrationTables.getTableForCovariate(cov).put(getRecalDatum(reportTable, i, false), tempCOVarray);
        }
    }

    /**
     *
     * Compiles the list of keys for the QualityScore table and uses the shared parsing utility to produce the actual table
     * @param reportTable            the GATKReport table containing data for this table
     * @param qualTable               the map representing this table
     */
    private void parseQualityScoreTable(final GATKReportTable reportTable, final NestedIntegerArray<RecalDatum> qualTable) {
        final int[] tempQUALarray = new int[3];
        for ( int i = 0; i < reportTable.getNumRows(); i++ ) {
            final Object rg = reportTable.get(i, RecalUtils.READGROUP_COLUMN_NAME);
            tempQUALarray[0] = covariates.getReadGroupCovariate().keyFromValue(rg);
            final Object qual = reportTable.get(i, RecalUtils.QUALITY_SCORE_COLUMN_NAME);
            tempQUALarray[1] = covariates.getQualityScoreCovariate().keyFromValue(qual);
            final EventType event = EventType.eventFrom((String)reportTable.get(i, RecalUtils.EVENT_TYPE_COLUMN_NAME));
            tempQUALarray[2] = event.ordinal();

            qualTable.put(getRecalDatum(reportTable, i, false), tempQUALarray);
        }
    }

    /**
     * Compiles the list of keys for the ReadGroup table and uses the shared parsing utility to produce the actual table
     *
     * @param reportTable            the GATKReport table containing data for this table
     * @param rgTable                the map representing this table
     */
    private void parseReadGroupTable(final GATKReportTable reportTable, final NestedIntegerArray<RecalDatum> rgTable) {
        final int[] tempRGarray = new int[2];
        for ( int i = 0; i < reportTable.getNumRows(); i++ ) {
            final Object rg = reportTable.get(i, RecalUtils.READGROUP_COLUMN_NAME);
            tempRGarray[0] = covariates.getReadGroupCovariate().keyFromValue(rg);
            final EventType event = EventType.eventFrom((String)reportTable.get(i, RecalUtils.EVENT_TYPE_COLUMN_NAME));
            tempRGarray[1] = event.ordinal();

            rgTable.put(getRecalDatum(reportTable, i, true), tempRGarray);
        }
    }

    private double asDouble(final Object o) {
        if ( o instanceof Double)
            return (Double)o;
        else if ( o instanceof Integer)
            return (Integer)o;
        else if ( o instanceof Long)
            return (Long)o;
        else
            throw new GATKException("Object " + o + " is expected to be either a double, long or integer but it's not either: " + o.getClass());
    }

    private long asLong(final Object o) {
        if ( o instanceof Long)
            return (Long)o;
        else if ( o instanceof Integer)
            return ((Integer)o).longValue();
        else if ( o instanceof Double)
            return ((Double)o).longValue();
        else
            throw new GATKException("Object " + o + " is expected to be a long but it's not: " + o.getClass());
    }

    private RecalDatum getRecalDatum(final GATKReportTable reportTable, final int row, final boolean hasEstimatedQReportedColumn) {
        final long nObservations = asLong(reportTable.get(row, RecalUtils.NUMBER_OBSERVATIONS_COLUMN_NAME));
        final double nErrors = asDouble(reportTable.get(row, RecalUtils.NUMBER_ERRORS_COLUMN_NAME));

        // the estimatedQreported column only exists in the ReadGroup table
        final double estimatedQReported = hasEstimatedQReportedColumn ?
                (Double) reportTable.get(row, RecalUtils.ESTIMATED_Q_REPORTED_COLUMN_NAME) : // we get it if we are in the read group table
                Byte.parseByte((String) reportTable.get(row, RecalUtils.QUALITY_SCORE_COLUMN_NAME)); // or we use the reported quality if we are in any other table

        final RecalDatum datum = new RecalDatum(nObservations, nErrors, (byte)1);
        datum.setEstimatedQReported(estimatedQReported);
        //datum.setEmpiricalQuality(empiricalQuality); // don't set the value here because we will want to recompute with a different conditional Q score prior value
        return datum;
    }

    /**
     * Parses the quantization table from the GATK Report and turns it into a map of original => quantized quality scores
     *
     * @param table the GATKReportTable containing the quantization mappings
     * @return an ArrayList with the quantization mappings from 0 to MAX_SAM_QUAL_SCORE
     */
    private QuantizationInfo initializeQuantizationTable(GATKReportTable table) {
        final Byte[] quals  = new Byte[QualityUtils.MAX_SAM_QUAL_SCORE + 1];
        final Long[] counts = new Long[QualityUtils.MAX_SAM_QUAL_SCORE + 1];
        for ( int i = 0; i < table.getNumRows(); i++ ) {
            final byte originalQual = (byte)i;
            final Object quantizedObject = table.get(i, RecalUtils.QUANTIZED_VALUE_COLUMN_NAME);
            final Object countObject = table.get(i, RecalUtils.QUANTIZED_COUNT_COLUMN_NAME);
            final byte quantizedQual = Byte.parseByte(quantizedObject.toString());
            final long quantizedCount = Long.parseLong(countObject.toString());
            quals[originalQual] = quantizedQual;
            counts[originalQual] = quantizedCount;
        }
        return new QuantizationInfo(new ArrayList<>(Arrays.asList(quals)), new ArrayList<>(Arrays.asList(counts)));
    }

    /**
     * Parses the arguments table from the GATK Report and creates a RAC object with the proper initialization values
     *
     * @param table the GATKReportTable containing the arguments and its corresponding values
     * @return a RAC object properly initialized with all the objects in the table
     */
    private RecalibrationArgumentCollection initializeArgumentCollectionTable(GATKReportTable table) {
        final RecalibrationArgumentCollection RAC = new RecalibrationArgumentCollection();

        final List<String> standardCovariateClassNames = new StandardCovariateList(RAC, Collections.emptyList()).getStandardCovariateClassNames();

        for ( int i = 0; i < table.getNumRows(); i++ ) {
            final String argument = table.get(i, "Argument").toString();
            Object value = table.get(i, RecalUtils.ARGUMENT_VALUE_COLUMN_NAME);
            if (value.equals("null")) {
                value = null; // generic translation of null values that were printed out as strings | todo -- add this capability to the GATKReport
            }

            if (argument.equals("covariate") && value != null) {
                final List<String> covs = new ArrayList<>(Arrays.asList(value.toString().split(",")));
                if (!covs.equals(standardCovariateClassNames)) {
                    throw new UserException("Non-standard covariates are not supported. Only the following are supported " + standardCovariateClassNames + " but was " + covs);
                }
            } else if (argument.equals("no_standard_covs")) {
                final boolean no_standard_covs = Boolean.parseBoolean((String) value);
                if (no_standard_covs){
                    throw new UserException("Non-standard covariates are not supported. Only the following are supported " + standardCovariateClassNames + " but no_standard_covs was true");
                }
            } else if (argument.equals("solid_recal_mode"))
                RAC.SOLID_RECAL_MODE = RecalUtils.SOLID_RECAL_MODE.recalModeFromString((String) value);

            else if (argument.equals("solid_nocall_strategy"))
                RAC.SOLID_NOCALL_STRATEGY = RecalUtils.SOLID_NOCALL_STRATEGY.nocallStrategyFromString((String) value);

            else if (argument.equals("mismatches_context_size"))
                RAC.MISMATCHES_CONTEXT_SIZE = Integer.parseInt((String) value);

            else if (argument.equals("indels_context_size"))
                RAC.INDELS_CONTEXT_SIZE = Integer.parseInt((String) value);

            else if (argument.equals("mismatches_default_quality"))
                RAC.MISMATCHES_DEFAULT_QUALITY = Byte.parseByte((String) value);

            else if (argument.equals("insertions_default_quality"))
                RAC.INSERTIONS_DEFAULT_QUALITY = Byte.parseByte((String) value);

            else if (argument.equals("deletions_default_quality"))
                RAC.DELETIONS_DEFAULT_QUALITY = Byte.parseByte((String) value);

            else if (argument.equals("maximum_cycle_value"))
                RAC.MAXIMUM_CYCLE_VALUE = Integer.parseInt((String) value);

            else if (argument.equals("low_quality_tail"))
                RAC.LOW_QUAL_TAIL = Byte.parseByte((String) value);

            else if (argument.equals("default_platform"))
                RAC.DEFAULT_PLATFORM = (String) value;

            else if (argument.equals("force_platform"))
                RAC.FORCE_PLATFORM = (String) value;

            else if (argument.equals("quantizing_levels"))
                RAC.QUANTIZING_LEVELS = Integer.parseInt((String) value);

            else if (argument.equals("recalibration_report"))
                RAC.existingRecalibrationReport = (value == null) ? null : new File((String) value);

            else if (argument.equals("binary_tag_name"))
                RAC.BINARY_TAG_NAME = (value == null) ? null : (String) value;

            else if (argument.equals("sort_by_all_columns"))
                RAC.SORT_BY_ALL_COLUMNS = Boolean.parseBoolean((String) value);
        }

        return RAC;
    }

    /**
     * Creates the recalibration report.  Report can then be written to a stream via GATKReport.print(PrintStream).
     *
     * @return newly created recalibration report
     */
    public GATKReport createGATKReport() {
        return RecalUtils.createRecalibrationGATKReport(argumentTable, quantizationInfo, recalibrationTables, covariates, RAC.SORT_BY_ALL_COLUMNS);
    }

    public RecalibrationArgumentCollection getRAC() {
        return RAC;
    }

    /**
     * @return true if the report has no data
     */
    public boolean isEmpty() {
        return recalibrationTables.isEmpty();
    }
}
