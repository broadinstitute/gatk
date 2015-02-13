/*
* Copyright (c) 2012 The Broad Institute
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.hellbender.utils.report;

import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.*;
import java.util.*;

/**
 * Container class for GATK report tables
 */
public class GATKReport {
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
        this(new File(filename));
    }

    /**
     * Create a new GATKReport with the contents of a GATKReport on disk.
     *
     * @param file the file to load
     */
    public GATKReport(File file) {
        loadReport(file);
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
     * Load a GATKReport file from disk
     *
     * @param file the file to load
     */
    private void loadReport(File file) {
        BufferedReader reader;
        String reportHeader;
        try {
            reader = new BufferedReader(IOUtils.makeReaderMaybeGzipped(file));
            reportHeader = reader.readLine();
        } catch (FileNotFoundException e) {
            throw new UserException.CouldNotReadInputFile(file, "it does not exist");
        } catch (IOException e) {
            throw new UserException.CouldNotReadInputFile(file, e);
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
        addTable(tableName, tableDescription, numColumns, GATKReportTable.TableSortingWay.DO_NOT_SORT);
    }

    /**
     * Add a new, empty table to the report
     *
     * @param tableName        the name of the table
     * @param tableDescription the description of the table
     * @param numColumns       the number of columns in this table
     * @param sortingWay       way to sort table
     */
    public void addTable(final String tableName, final String tableDescription, final int numColumns, final GATKReportTable.TableSortingWay sortingWay) {
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
     * Return true if table with a given name exists
     *
     * @param tableName the name of the table
     * @return true if the table exists, false otherwise
     */
    public boolean hasTable(String tableName) {
        return tables.containsKey(tableName);
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
        out.println(GATKREPORT_HEADER_PREFIX + getVersion().toString() + SEPARATOR + getTables().size());
        for (GATKReportTable table : tables.values())
            table.write(out);
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
    public static GATKReport newSimpleReport(final String tableName, final String... columns) {
        return newSimpleReportWithDescription(tableName, "A simplified GATK table report", columns);
    }

    /**
     * @see #newSimpleReport(String, String...) but with a customized description
     * @param tableName
     * @param desc
     * @param columns
     * @return
     */
    public static GATKReport newSimpleReportWithDescription(final String tableName, final String desc, final String... columns) {
        GATKReportTable table = new GATKReportTable(tableName, desc, columns.length);

        for (String column : columns) {
            table.addColumn(column, "");
        }

        GATKReport output = new GATKReport();
        output.addTable(table);

        return output;
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
    public static GATKReport newSimpleReport(final String tableName, final List<String> columns) {
        GATKReportTable table = new GATKReportTable(tableName, "A simplified GATK table report", columns.size());

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

    /**
     * This method provides an efficient way to populate a simplified GATK report. This method will only work on reports
     * that qualify as simplified GATK reports. See the newSimpleReport() constructor for more information.
     *
     * @param values     the row of data to be added to the table.
     *               Note: the number of arguments must match the columns in the table.
     */
    public void addRowList(final List<Object> values) {
        if ( tables.size() != 1 )
            throw new GATKException("Cannot write a row to a complex GATK Report");

        GATKReportTable table = tables.firstEntry().getValue();
        if ( table.getNumColumns() != values.size() )
            throw new GATKException("The number of arguments in writeRow() must match the number of columns in the table");

        final int rowIndex = table.getNumRows();
        int idx = 0;
        for ( Object value : values ) {
            table.set(rowIndex,idx,value);
            idx++;
        }
    }
}
