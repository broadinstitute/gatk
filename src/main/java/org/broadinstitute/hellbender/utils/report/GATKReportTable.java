package org.broadinstitute.hellbender.utils.report;

import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.text.TextFormattingUtils;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintStream;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public final class GATKReportTable {
    /**
     * REGEX that matches any table with an invalid name
     */
    public static final String INVALID_TABLE_NAME_REGEX = "[^a-zA-Z0-9_\\-\\.]";
    private static final String GATKTABLE_HEADER_PREFIX = "#:GATKTable";
    private static final String SEPARATOR = ":";
    private static final String ENDLINE = ":;";

    private static final Comparator<? super Object[]> ROW_COMPARATOR = (objectArr1, objectArr2) -> {
        final int EQUAL = 0;
        int result = EQUAL;
        for (int x = 0; x < objectArr1.length; x++) {
            if (!Objects.equals(objectArr1[x].getClass(), objectArr2[x].getClass())){
                result = objectArr1[x].toString().compareTo(objectArr2[x].toString()); //compare by toString if types different
            } else if (objectArr1[x] instanceof Byte) {
                result = ((Byte) objectArr1[x]).compareTo((Byte) objectArr2[x]);
            } else if (objectArr1[x] instanceof Short) {
                result = ((Short) objectArr1[x]).compareTo((Short) objectArr2[x]);
            } else if (objectArr1[x] instanceof Integer) {
                result = ((Integer) objectArr1[x]).compareTo((Integer) objectArr2[x]);
            } else if (objectArr1[x] instanceof Long) {
                result = ((Long) objectArr1[x]).compareTo((Long) objectArr2[x]);
            } else if (objectArr1[x] instanceof BigInteger) {
                result = ((BigInteger) objectArr1[x]).compareTo((BigInteger) objectArr2[x]);
            } else if (objectArr1[x] instanceof Float) {
                result = ((Float) objectArr1[x]).compareTo((Float) objectArr2[x]);
            } else if (objectArr1[x] instanceof Double) {
                result = ((Double) objectArr1[x]).compareTo((Double) objectArr2[x]);
            } else if (objectArr1[x] instanceof BigDecimal) {
                result = ((BigDecimal) objectArr1[x]).compareTo((BigDecimal) objectArr2[x]);
            } else { // default uses String comparison
                result = objectArr1[x].toString().compareTo(objectArr2[x].toString());
            }
            if (result != EQUAL) {
                return result;
            }
        }
        return result;
    };

    private final String tableName;
    private final String tableDescription;

    private final Sorting sortingWay;

    private final List<Object[]> underlyingData;
    private final List<GATKReportColumn> columnInfo;
    private final Map<Object, Integer> columnNameToIndex;
    private final Map<Object, Integer> rowIdToIndex;

    private static final String COULD_NOT_READ_HEADER = "Could not read the header of this file -- ";
    private static final String COULD_NOT_READ_COLUMN_NAMES = "Could not read the column names of this file -- ";
    private static final String COULD_NOT_READ_DATA_LINE = "Could not read a data line of this table -- ";
    private static final String COULD_NOT_READ_EMPTY_LINE = "Could not read the last empty line of this table -- ";
    private static final String OLD_GATK_TABLE_VERSION = "We no longer support older versions of the GATK Tables";

    private static final int INITITAL_ARRAY_SIZE = 10000;

    protected enum TableDataHeaderFields {
        COLS(2),
        ROWS(3),
        FORMAT_START(4);

        private final int index;
        TableDataHeaderFields(int index) { this.index = index; }
        public int index() { return index; }
    }

    public enum Sorting {
        SORT_BY_ROW,
        SORT_BY_COLUMN,
        DO_NOT_SORT
    }

    protected enum TableNameHeaderFields {
        NAME(2),
        DESCRIPTION(3);

        private final int index;
        TableNameHeaderFields(int index) { this.index = index; }
        public int index() { return index; }
    }

    /**
     * Construct a new GATK report table from the reader
     * Note that the row ID mappings are just the index -> index
     *
     * @param reader        the reader
     * @param version       the GATK report version
     */
    public GATKReportTable(BufferedReader reader, GATKReportVersion version) {

        switch ( version ) {
            case V1_1:
                // read in the header lines
                final String[] tableData, tableNameData;
                try {
                    tableData = reader.readLine().split(SEPARATOR);
                    tableNameData = reader.readLine().split(SEPARATOR);
                } catch (IOException e) {
                    throw new GATKException(COULD_NOT_READ_HEADER + e.getMessage());
                }

                // parse the header fields
                tableName = tableNameData[TableNameHeaderFields.NAME.index()];
                tableDescription = (tableNameData.length <= TableNameHeaderFields.DESCRIPTION.index()) ? "" : tableNameData[TableNameHeaderFields.DESCRIPTION.index()];                                           // table may have no description! (and that's okay)

                // when reading from a file, we do not re-sort the rows
                sortingWay = Sorting.DO_NOT_SORT;

                // initialize the data
                final int nColumns = Integer.parseInt(tableData[TableDataHeaderFields.COLS.index()]);
                final int nRows = Integer.parseInt(tableData[TableDataHeaderFields.ROWS.index()]);
                underlyingData = new ArrayList<>(nRows);
                columnInfo = new ArrayList<>(nColumns);
                columnNameToIndex = new LinkedHashMap<>(nColumns);

                // when reading from a file, the row ID mapping is just the index
                rowIdToIndex = new LinkedHashMap<>();
                for ( int i = 0; i < nRows; i++ )
                    rowIdToIndex.put(i, i);

                // read the column names
                final String columnLine;
                try {
                    columnLine = reader.readLine();
                } catch (IOException e) {
                    throw new GATKException(COULD_NOT_READ_COLUMN_NAMES);
                }

                final List<Integer> columnStarts = TextFormattingUtils.getWordStarts(columnLine);
                final String[] columnNames = TextFormattingUtils.splitFixedWidth(columnLine, columnStarts);

                // Put in columns using the format string from the header
                for ( int i = 0; i < nColumns; i++ ) {
                    final String format = tableData[TableDataHeaderFields.FORMAT_START.index() + i];
                    addColumn(columnNames[i], format);
                }

                // fill in the table
                try {
                    for ( int i = 0; i < nRows; i++ ) {
                        // read a data line
                        final String dataLine = reader.readLine();
                        final List<String> lineSplits = Arrays.asList(TextFormattingUtils.splitFixedWidth(dataLine, columnStarts));

                        underlyingData.add(new Object[nColumns]);
                        for ( int columnIndex = 0; columnIndex < nColumns; columnIndex++ ) {

                            final GATKReportDataType type = columnInfo.get(columnIndex).getDataType();
                            final String columnName = columnNames[columnIndex];
                            set(i, columnName, type.Parse(lineSplits.get(columnIndex)));

                        }
                    }
                } catch (IOException e) {
                    throw new GATKException(COULD_NOT_READ_DATA_LINE + e.getMessage());
                }

                try {
                    reader.readLine();
                } catch (IOException e) {
                    throw new GATKException(COULD_NOT_READ_EMPTY_LINE + e.getMessage());
                }
            break;

            default:
                throw new GATKException(OLD_GATK_TABLE_VERSION);
        }
    }

    /**
     * Construct a new GATK report table with the specified name and description
     *
     * Sorting will be done with {@link Sorting#SORT_BY_ROW}
     * @param tableName        the name of the table
     * @param tableDescription the description of the table
     * @param numColumns       the number of columns in this table
     *
     *
     */
    public GATKReportTable(final String tableName, final String tableDescription, final int numColumns) {
        this(tableName, tableDescription, numColumns, Sorting.SORT_BY_ROW);
    }

    /**
     * Construct a new GATK report table with the specified name and description and whether to sort rows by the row ID.
     *
     * @param tableName          the name of the table
     * @param tableDescription   the description of the table
     * @param numColumns         the number of columns in this table
     * @param sortingWay         in what way to sort rows (instead of the order in which they were added)
     */
    public GATKReportTable(final String tableName, final String tableDescription, final int numColumns, final Sorting sortingWay) {
        if ( !isValidName(tableName) ) {
            throw new GATKException("Attempted to set a GATKReportTable name of '" + tableName + "'.  GATKReportTable names must be purely alphanumeric - no spaces or special characters are allowed.");
        }

        if ( !isValidDescription(tableDescription) ) {
            throw new GATKException("Attempted to set a GATKReportTable description of '" + tableDescription + "'.  GATKReportTable descriptions must not contain newlines.");
        }

        this.tableName = tableName;
        this.tableDescription = tableDescription;
        this.sortingWay = sortingWay;

        underlyingData = new ArrayList<>(INITITAL_ARRAY_SIZE);
        columnInfo = new ArrayList<>(numColumns);
        columnNameToIndex = new LinkedHashMap<>(numColumns);
        rowIdToIndex = new LinkedHashMap<>();
    }

    /**
    * Verifies that a table or column name has only alphanumeric characters - no spaces or special characters allowed
    *
    * @param name the name of the table or column
    * @return true if the name is valid, false if otherwise
    */
    private boolean isValidName(String name) {
        Pattern p = Pattern.compile(INVALID_TABLE_NAME_REGEX);
        Matcher m = p.matcher(name);

        return !m.find();
    }

    /**
     * Verifies that a table or column name has only alphanumeric characters - no spaces or special characters allowed
     *
     * @param description the name of the table or column
     * @return true if the name is valid, false if otherwise
     */
    private boolean isValidDescription(String description) {
        Pattern p = Pattern.compile("\\r|\\n");
        Matcher m = p.matcher(description);

        return !m.find();
    }

    /**
     * Add a mapping from ID to the index of a new row added to the table.
     *
     * @param ID                    the unique ID
     * @param populateFirstColumn   should we automatically populate the first column with the row's ID?
     */
    public void addRowID(final String ID, final boolean populateFirstColumn) {
        addRowIDMapping(ID, underlyingData.size(), populateFirstColumn);
    }

    /**
     * Add a mapping from ID to row index.
     *
     * @param ID                    the unique ID
     * @param index                 the index associated with the ID
     */
    public void addRowIDMapping(final String ID, final int index) {
        addRowIDMapping(ID, index, false);
    }

    /**
     * Add a mapping from ID to row index.
     *
     * @param ID                    the unique ID
     * @param index                 the index associated with the ID
     * @param populateFirstColumn   should we automatically populate the first column with the row's ID?
     */
    public void addRowIDMapping(final Object ID, final int index, final boolean populateFirstColumn) {
        expandTo(index, false);
        rowIdToIndex.put(ID, index);

        if ( populateFirstColumn )
            set(index, 0, ID);
    }

    /**
     * Add a column to the report and the format string used to display the data.
     *
     * @param columnName   the name of the column
     * @param format       the format string used to display data
     */
    public void addColumn(String columnName, String format) {
        columnNameToIndex.put(columnName, columnInfo.size());
        columnInfo.add(new GATKReportColumn(columnName, format));
    }

    /**
     * Check if the requested cell is valid and expand the table if necessary
     *
     * @param rowIndex    the row index
     * @param colIndex    the column index
     */
    private void verifyEntry(final int rowIndex, final int colIndex) {
        if ( rowIndex < 0 || colIndex < 0 || colIndex >= getNumColumns() )
            throw new GATKException("attempted to access a cell that does not exist in table '" + tableName + "'");
    }

    /**
     * expand the underlying table if needed to include the given row index
     *
     * @param rowIndex        the row index
     * @param updateRowIdMap  should we update the row ID map?
     */
    private void expandTo(final int rowIndex, final boolean updateRowIdMap) {
        int currentSize = underlyingData.size();
        if ( rowIndex >= currentSize ) {
            final int numNewRows = rowIndex - currentSize + 1;
            for ( int i = 0; i < numNewRows; i++ ) {
                if ( updateRowIdMap )
                    rowIdToIndex.put(currentSize, currentSize);
                underlyingData.add(new Object[getNumColumns()]);
                currentSize++;
            }
        }
    }

    /**
     * Set the value for a given position in the table.
     * If the row ID doesn't exist, it will create a new row in the table with the given ID.
     *
     * @param rowID        the row ID
     * @param columnName   the name of the column
     * @param value        the value to set
     */
    public void set(final Object rowID, final String columnName, final Object value) {
        if ( !rowIdToIndex.containsKey(rowID) ) {
            rowIdToIndex.put(rowID, underlyingData.size());
            expandTo(underlyingData.size(), false);
        }
        set(rowIdToIndex.get(rowID), columnNameToIndex.get(columnName), value);
    }

    /**
     * Set the value for a given position in the table.
     * If the row index doesn't exist, it will create new rows in the table accordingly.
     *
     * @param rowIndex     the row index
     * @param colIndex     the column index
     * @param value        the value to set
     */
    public void set(final int rowIndex, final int colIndex, Object value) {
        expandTo(rowIndex, true);
        verifyEntry(rowIndex, colIndex);
        GATKReportColumn column = columnInfo.get(colIndex);

        // We do not accept internal null values
        if (value == null)
            value = "null";
        else
            value = fixType(value, column);

        if ( column.getDataType().equals(GATKReportDataType.fromObject(value)) || column.getDataType().equals(GATKReportDataType.Unknown) ) {
            underlyingData.get(rowIndex)[colIndex] = value;
            column.updateFormatting(value);
        } else {
            throw new GATKException(String.format("Tried to add an object of type: %s to a column of type: %s", GATKReportDataType.fromObject(value).name(), column.getDataType().name()));
        }
    }

    /**
    * Increment the value for a given position in the table.
    * Throws an exception if the value in the cell is not an integer.
    *
    * @param rowID        the row ID
    * @param columnName   the name of the column
    */
    public void increment(final Object rowID, final String columnName) {
        int prevValue;
        if ( !rowIdToIndex.containsKey(rowID) ) {
            rowIdToIndex.put(rowID, underlyingData.size());
            underlyingData.add(new Object[getNumColumns()]);
            prevValue = 0;
        } else {
            Object obj = get(rowID, columnName);
            if ( !(obj instanceof Integer) )
                throw new GATKException("Attempting to increment a value in a cell that is not an integer");
            prevValue = (Integer)obj;
        }

        set(rowIdToIndex.get(rowID), columnNameToIndex.get(columnName), prevValue + 1);
    }

    private Object fixType(final Object value, final GATKReportColumn column) {
        // Below is some code to convert a string into its appropriate type.

        // todo -- Types have to be more flexible. For example, %d should accept Integers, Shorts and Bytes.

        Object newValue = null;
        if ( value instanceof String && !column.getDataType().equals(GATKReportDataType.String) ) {
            // Integer case
            if ( column.getDataType().equals(GATKReportDataType.Integer) ) {
                try {
                    newValue = Long.parseLong((String) value);
                } catch (Exception e) {
                    /** do nothing */
                }
            }
            if ( column.getDataType().equals(GATKReportDataType.Decimal) ) {
                try {
                    newValue = Double.parseDouble((String) value);
                } catch (Exception e) {
                    /** do nothing */
                }
            }
            if ( column.getDataType().equals(GATKReportDataType.Character) && ((String) value).length() == 1 ) {
                newValue = ((String) value).charAt(0);
            }
        }

        return  (newValue != null) ? newValue : value;
    }

    /**
     * Get a value from the given position in the table
     *
     * @param rowID       the row ID
     * @param columnName  the name of the column
     * @return the value stored at the specified position in the table
     */
    public Object get(final Object rowID, final String columnName) {
        return get(rowIdToIndex.get(rowID), columnNameToIndex.get(columnName));
    }

    /**
     * Get a value from the given position in the table
     *
     * @param rowIndex       the row ID
     * @param columnName  the name of the column
     * @return the value stored at the specified position in the table
     */
    public Object get(final int rowIndex, final String columnName) {
        return get(rowIndex, columnNameToIndex.get(columnName));
    }

    /**
     * Get a value from the given position in the table
     *
     * @param rowIndex    the index of the row
     * @param columnIndex the index of the column
     * @return the value stored at the specified position in the table
     */
    public Object get(int rowIndex, int columnIndex) {
        verifyEntry(rowIndex, columnIndex);
        return underlyingData.get(rowIndex)[columnIndex];
    }

    /**
     * Write the table to the PrintStream, formatted nicely to be human-readable, AWK-able, and R-friendly.
     *
     * @param out the PrintStream to which the table should be written
     */
     void write(final PrintStream out) {
         write(out, sortingWay);
     }

    /**
     * Write the table to the PrintStream, formatted nicely to be human-readable, AWK-able, and R-friendly.
     *
     * @param out the PrintStream to which the table should be written
     * @param sortingWay how to sort the table
     */
     void write(final PrintStream out, Sorting sortingWay) {

         /*
          * Table header:
          * #:GATKTable:nColumns:nRows:(DataType for each column):;
          * #:GATKTable:TableName:Description :;
          * key   colA  colB
          * row1  xxxx  xxxxx
         */

         // write the table definition
         out.printf(GATKTABLE_HEADER_PREFIX + ":%d:%d", getNumColumns(), getNumRows());

         // write the formats for all the columns
         for ( final GATKReportColumn column : columnInfo )
             out.print(SEPARATOR + column.getFormat());
         out.println(ENDLINE);

         // write the table name & description
         out.printf(GATKTABLE_HEADER_PREFIX + ":%s:%s\n", tableName, tableDescription);

         // write the column names
         boolean needsPadding = false;
         for ( final GATKReportColumn column : columnInfo ) {
             if ( needsPadding )
                 out.printf("  ");
             needsPadding = true;

             out.printf(column.getColumnFormat().getNameFormat(), column.getColumnName());
         }
         out.println();

         // write the table body
         switch (sortingWay) {
             case SORT_BY_COLUMN:
                 Collections.sort(underlyingData, ROW_COMPARATOR);
                 for ( final Object[] row : underlyingData ) {
                     writeRow(out, row);
                 }
                 break;
             case SORT_BY_ROW:
                 // make sure that there are exactly the correct number of ID mappings
                 if ( rowIdToIndex.size() != underlyingData.size() )
                     throw new GATKException("There isn't a 1-to-1 mapping from row ID to index; this can happen when rows are not created consistently");

                 final TreeMap<Object, Integer> sortedMap;
                 try {
                     sortedMap = new TreeMap<>(rowIdToIndex);
                 } catch (ClassCastException e) {
                     throw new GATKException("Unable to sort the rows based on the row IDs because the ID Objects are of different types");
                 }
                 for ( final Map.Entry<Object, Integer> rowKey : sortedMap.entrySet() )
                     writeRow(out, underlyingData.get(rowKey.getValue()));
                 break;
             case DO_NOT_SORT:
                 for ( final Object[] row : underlyingData )
                     writeRow(out, row);
         }
         out.println();
     }

    private void writeRow(final PrintStream out, final Object[] row) {
        boolean needsPadding = false;
        for ( int i = 0; i < row.length; i++ ) {
            if ( needsPadding )
                out.printf("  ");
            needsPadding = true;

            final Object obj = row[i];
            final String value;

            final GATKReportColumn info = columnInfo.get(i);

            if ( obj == null )
                value = "null";
            else if ( info.getDataType().equals(GATKReportDataType.Unknown) && (obj instanceof Double || obj instanceof Float) )
                value = String.format("%.8f", obj);
            else
                value = String.format(info.getFormat(), obj);

            out.printf(info.getColumnFormat().getValueFormat(), value);
        }

        out.println();
    }

    public int getNumRows() {
        return underlyingData.size();
    }

    public int getNumColumns() {
        return columnInfo.size();
    }

    public List<GATKReportColumn> getColumnInfo() {
        return columnInfo;
    }

    public String getTableName() {
        return tableName;
    }

    /**
     * Concatenates the rows from the table to this one
     *
     * @param table another GATK table
     */
    public void concat(final GATKReportTable table) {
        if ( !isSameFormat(table) )
            throw new GATKException("Error trying to concatenate tables with different formats");

        // add the data
        underlyingData.addAll(table.underlyingData);

        // update the row index map
        final int currentNumRows = getNumRows();
        for ( Map.Entry<Object, Integer> entry : table.rowIdToIndex.entrySet() )
            rowIdToIndex.put(entry.getKey(), entry.getValue() + currentNumRows);
    }

    /**
     * Returns whether or not the two tables have the same format including columns and everything in between. This does
     * not check if the data inside is the same. This is the check to see if the two tables are gatherable or
     * reduceable
     *
     * @param table another GATK table
     * @return true if the the tables are gatherable
     */
    public boolean isSameFormat(final GATKReportTable table) {
        if ( !tableName.equals(table.tableName) ||
                !tableDescription.equals(table.tableDescription) ||
                columnInfo.size() != table.columnInfo.size() )
            return false;

        for ( int i = 0; i < columnInfo.size(); i++ ) {
            if ( !columnInfo.get(i).getFormat().equals(table.columnInfo.get(i).getFormat()) ||
                    !columnInfo.get(i).getColumnName().equals(table.columnInfo.get(i).getColumnName()) )
                return false;
        }

        return true;
    }

    /**
     * Checks that the tables are exactly the same.
     *
     * @param table another GATK report
     * @return true if all field in the reports, tables, and columns are equal.
     */
    public boolean equals(final GATKReportTable table) {
        if ( !isSameFormat(table) ||
                underlyingData.size() != table.underlyingData.size() )
            return false;

        final List<Object[]> myOrderedRows = getOrderedRows();
        final List<Object[]> otherOrderedRows = table.getOrderedRows();

        for ( int i = 0; i < underlyingData.size(); i++ ) {
            final Object[] myData = myOrderedRows.get(i);
            final Object[] otherData = otherOrderedRows.get(i);
            for ( int j = 0; j < myData.length; j++ ) {
                if ( !myData[j].toString().equals(otherData[j].toString()) )       // need to deal with different typing (e.g. Long vs. Integer)
                    return false;
            }
        }

        return true;
    }

    private List<Object[]> getOrderedRows() {

        switch (sortingWay) {
            case SORT_BY_COLUMN:
                Collections.sort(underlyingData, ROW_COMPARATOR);
                return underlyingData;
            case SORT_BY_ROW:
                final TreeMap<Object, Integer> sortedMap;
                try {
                    sortedMap = new TreeMap<>(rowIdToIndex);
                } catch (ClassCastException e) {
                    return underlyingData;
                }

                final List<Object[]> orderedData = new ArrayList<>(underlyingData.size());
                for ( final int rowKey : sortedMap.values() )
                    orderedData.add(underlyingData.get(rowKey));

                return orderedData;
            default:
                return underlyingData;
        }
    }
}
