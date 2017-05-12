package org.broadinstitute.hellbender.utils.report;

import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.recalibration.RecalUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.io.Reader;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.NavigableMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

/**
 * Container class for GATK report tables
 */
public final class GATKReport {
    public static final String RECAL_FILE = "input covariates table file for base quality score recalibration";
    public static final String GATKREPORT_HEADER_PREFIX = "#:GATKReport.";
    public static final GATKReportVersion LATEST_REPORT_VERSION = GATKReportVersion.V1_1;
    private static final String SEPARATOR = ":";
    private GATKReportVersion version = LATEST_REPORT_VERSION;

    private final NavigableMap<String, GATKReportTable> tables = new TreeMap<>();

    /**
     * Create a new, empty GATKReport.
     */
    public GATKReport() {
    }

    /**
     * Create a new GATKReport with the contents of a GATKReport on disk.
     *
     * @param filename the path to the file to load
     */
    public GATKReport(String filename) {
        this(BucketUtils.openFile(filename));
    }

    /**
     * Create a new GATKReport with the contents of a GATKReport on disk.
     *
     * @param file the file to load
     */
    public GATKReport(File file) {
        this(file.getPath());
    }

    public GATKReport(InputStream in){
        loadReport(new InputStreamReader(in));
    }

    /**
     * Create a new GATK report from GATK report tables
     * @param tables Any number of tables that you want to add to the report
     */
    public GATKReport(GATKReportTable... tables) {
        for( GATKReportTable table: tables)
            addTable(table);
    }

    /**
     * Gets the unique read groups in the table
     *
     * @return the unique read groups
     */
    public SortedSet<String> getReadGroups() {
        final GATKReportTable reportTable = getTable(RecalUtils.READGROUP_REPORT_TABLE_TITLE);
        final SortedSet<String> readGroups = new TreeSet<>();
        for ( int i = 0; i < reportTable.getNumRows(); i++ ) {
            readGroups.add(reportTable.get(i, RecalUtils.READGROUP_COLUMN_NAME).toString());
        }
        return readGroups;
    }

    /**
     * Load a GATKReport from a {@link Reader}
     *
     * @param in the reader to load from
     */
    private void loadReport(Reader in) {
        BufferedReader reader = new BufferedReader(in);
        String reportHeader;
        try {
            reportHeader = reader.readLine();
        } catch (IOException e) {
            throw new UserException("Could not read " + RECAL_FILE, e);
        }

        if ( reportHeader == null ) {
            throw new UserException(RECAL_FILE + " is empty.");
        }

        // Read the first line for the version and number of tables.
        version = GATKReportVersion.fromHeader(reportHeader);
        if (version.equals(GATKReportVersion.V0_1) ||
                version.equals(GATKReportVersion.V0_2))
            throw new UserException("The GATK no longer supports reading legacy GATK Reports. Please use v1.0 or newer.");

        int nTables = Integer.parseInt(reportHeader.split(":")[2]);

        // Read each table according ot the number of tables
        for (int i = 0; i < nTables; i++) {
            addTable(new GATKReportTable(reader, version));
        }
    }


    /**
     * Add a new, empty table to the report
     *
     * @param tableName        the name of the table
     * @param tableDescription the description of the table
     * @param numColumns       the number of columns in this table
     */
    public void addTable(final String tableName, final String tableDescription, final int numColumns) {
        addTable(tableName, tableDescription, numColumns, GATKReportTable.Sorting.DO_NOT_SORT);
    }

    /**
     * Add a new, empty table to the report
     *
     * @param tableName        the name of the table
     * @param tableDescription the description of the table
     * @param numColumns       the number of columns in this table
     * @param sortingWay       way to sort table
     */
    public void addTable(final String tableName, final String tableDescription, final int numColumns, final GATKReportTable.Sorting sortingWay) {
        GATKReportTable table = new GATKReportTable(tableName, tableDescription, numColumns, sortingWay);
        tables.put(tableName, table);
    }

    /**
     * Adds a table, empty or populated, to the report
     *
     * @param table the table to add
     */
    public void addTable(GATKReportTable table) {
        tables.put(table.getTableName(), table);
    }

    public void addTables(List<GATKReportTable> gatkReportTableV2s) {
        for ( GATKReportTable table : gatkReportTableV2s )
            addTable(table);
    }

    /**
     * Return a table with a given name
     *
     * @param tableName the name of the table
     * @return the table object
     */
    public GATKReportTable getTable(String tableName) {
        GATKReportTable table = tables.get(tableName);
        if (table == null)
            throw new GATKException("Table is not in GATKReport: " + tableName);
        return table;
    }

    /**
     * Print all tables contained within this container to a PrintStream
     *
     * @param out the PrintStream to which the tables should be written
     */
    public void print(PrintStream out) {
        out.println(GATKREPORT_HEADER_PREFIX + getVersion() + SEPARATOR + getTables().size());
        for (GATKReportTable table : tables.values()) {
            table.write(out);
        }
    }

    /**
     * Print all tables contained within this container to a PrintStream
     *
     * @param out the PrintStream to which the tables should be written
     */
    public void print(PrintStream out, GATKReportTable.Sorting sortingWay) {
        out.println(GATKREPORT_HEADER_PREFIX + getVersion() + SEPARATOR + getTables().size());
        for (GATKReportTable table : tables.values()) {
            table.write(out, sortingWay);
        }
    }



    public Collection<GATKReportTable> getTables() {
        return tables.values();
    }

    /**
     * This is the main function is charge of gathering the reports. It checks that the reports are compatible and then
     * calls the table gathering functions.
     *
     * @param input another GATKReport of the same format
     */
    public void concat(GATKReport input) {

        if ( !isSameFormat(input) ) {
            throw new GATKException("Failed to combine GATKReport, format doesn't match!");
        }

        for ( Map.Entry<String, GATKReportTable> table : tables.entrySet() ) {
            table.getValue().concat(input.getTable(table.getKey()));
        }
    }

    public GATKReportVersion getVersion() {
        return version;
    }

    /**
     * Returns whether or not the two reports have the same format, from columns, to tables, to reports, and everything
     * in between. This does not check if the data inside is the same. This is the check to see if the two reports are
     * gatherable or reduceable.
     *
     * @param report another GATK report
     * @return true if the the reports are gatherable
     */
    public boolean isSameFormat(GATKReport report) {
        if (!version.equals(report.version)) {
            return false;
        }
        if (!tables.keySet().equals(report.tables.keySet())) {
            return false;
        }
        for (String tableName : tables.keySet()) {
            if (!getTable(tableName).isSameFormat(report.getTable(tableName)))
                return false;
        }
        return true;
    }

    /**
     * Checks that the reports are exactly the same.
     *
     * @param report another GATK report
     * @return true if all field in the reports, tables, and columns are equal.
     */
    public boolean equals(GATKReport report) {
        if (!version.equals(report.version)) {
            return false;
        }
        if (!tables.keySet().equals(report.tables.keySet())) {
            return false;
        }
        for (String tableName : tables.keySet()) {
            if (!getTable(tableName).equals(report.getTable(tableName)))
                return false;
        }
        return true;
    }

    /**
     * The constructor for a simplified GATK Report. Simplified GATK report are designed for reports that do not need
     * the advanced functionality of a full GATK Report.
     * <p/>
     * A simple GATK Report consists of:
     * <p/>
     * - A single table
     * - No primary key ( it is hidden )
     * <p/>
     * Optional:
     * - Only untyped columns. As long as the data is an Object, it will be accepted.
     * - Default column values being empty strings.
     * <p/>
     * Limitations:
     * <p/>
     * - A simple GATK report cannot contain multiple tables.
     * - It cannot contain typed columns, which prevents arithmetic gathering.
     *
     * @param tableName The name of your simple GATK report table
     * @param columns   The names of the columns in your table
     * @return a simplified GATK report
     */
    public static GATKReport newSimpleReport(final String tableName, GATKReportTable.Sorting sorting, final String... columns) {
        return newSimpleReportWithDescription(tableName, "A simplified GATK table report", sorting, columns);
    }

    /**
     * @see #newSimpleReport(String, GATKReportTable.Sorting, String...) but with a customized description
     */
    public static GATKReport newSimpleReportWithDescription(final String tableName, final String desc, GATKReportTable.Sorting sorting, final String... columns) {
        GATKReportTable table = new GATKReportTable(tableName, desc, columns.length, sorting);

        for (String column : columns) {
            table.addColumn(column, "");
        }

        GATKReport output = new GATKReport();
        output.addTable(table);

        return output;
    }

    /**
     * This method provides an efficient way to populate a simplified GATK report. This method will only work on reports
     * that qualify as simplified GATK reports. See the newSimpleReport() constructor for more information.
     *
     * @param values     the row of data to be added to the table.
     *               Note: the number of arguments must match the columns in the table.
     */
    public void addRow(final Object... values) {
        // Must be a simple report
        if ( tables.size() != 1 )
            throw new GATKException("Cannot write a row to a complex GATK Report");

        GATKReportTable table = tables.firstEntry().getValue();
        if ( table.getNumColumns() != values.length )
            throw new GATKException("The number of arguments in writeRow (" + values.length + ") must match the number of columns in the table (" + table.getNumColumns() + ")" );

        final int rowIndex = table.getNumRows();
        for ( int i = 0; i < values.length; i++ )
            table.set(rowIndex, i, values[i]);
    }

}
