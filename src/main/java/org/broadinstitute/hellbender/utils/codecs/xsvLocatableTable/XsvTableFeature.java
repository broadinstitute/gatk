package org.broadinstitute.hellbender.utils.codecs.xsvLocatableTable;

import htsjdk.tribble.Feature;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * A feature to represent a line in an arbitrarily delimited (XSV) file (i.e. TSV, CSV, etc.).
 * Created by jonn on 12/4/17.
 */
public class XsvTableFeature implements Feature {

    //==================================================================================================================
    // Public Static Members:

    //==================================================================================================================
    // Private Static Members:

    //==================================================================================================================
    // Private Members:

    /** The column index from which to read the contig for this {@link XsvTableFeature}. */
    private final int contigColumn;
    /** The column index from which to read the start position for this {@link XsvTableFeature}. */
    private final int startColumn;
    /** The column index from which to read the end position for this {@link XsvTableFeature}. */
    private final int endColumn;
    /** The indices for location information in descending order (to be used to remove these columns). */
    private final List<Integer> locationColumnRemoveIndiciesInOrder;
    /** Name of the source from which this feature was derived. */
    private final String dataSourceName;
    /** The contig for this feature. */
    private final String contig;
    /** The genomic start position for this feature. */
    private final int start;
    /** The genomic start position for this feature. */
    private final int end;
    /** The names of the columns in this {@link XsvTableFeature}. */
    private final List<String> columnNames;
    /** The values for each of the columns in this {@link XsvTableFeature}. */
    private final List<String> columnValues;

    //==================================================================================================================
    // Constructors:

    /**
     * Create an {@link XsvTableFeature}.
     * @param contigColumn The column index from which to read the contig for this {@link XsvTableFeature}.
     * @param startColumn The column index from which to read the start position for this {@link XsvTableFeature}.
     * @param endColumn The column index from which to read the end position for this {@link XsvTableFeature}.
     * @param columnNames The names of the columns in this {@link XsvTableFeature}.  Must have same number of entries as {@code columnValues}.
     * @param columnValues The values for each of the columns in this {@link XsvTableFeature}.  Must have same number of entries as {@code columnNames}.
     * @param dataSourceName The name of the source from which this {@link XsvTableFeature} was derived.
     */
    public XsvTableFeature(final int contigColumn, final int startColumn, final int endColumn,
                           final List<String> columnNames, final List<String> columnValues, final String dataSourceName) {

        this.contigColumn = contigColumn;
        this.startColumn = startColumn;
        this.endColumn = endColumn;
        this.dataSourceName = dataSourceName;

        contig = columnValues.get(contigColumn);
        try {
            start = Integer.valueOf(columnValues.get(startColumn));
            end = Integer.valueOf(columnValues.get(endColumn));
        }
        catch ( final NumberFormatException ex ) {
            throw new UserException.MalformedFile("Could not convert value (" + ex.getMessage() + ") from input file into a number for Data Source: " + dataSourceName);
        }

        if ( columnNames.size() != columnValues.size() ) {
            throw new UserException.BadInput("Number of columns in given header and data do not match: " + columnNames.size() + " != " + columnValues.size());
        }

        this.columnNames = columnNames;
        this.columnValues = columnValues;

        // Create our list of indices to remove:
        if ( startColumn == endColumn ) {
            // Don't add the same column more than once:
            locationColumnRemoveIndiciesInOrder = new ArrayList<>(Arrays.asList(contigColumn, startColumn));
        }
        else {
            locationColumnRemoveIndiciesInOrder = new ArrayList<>(Arrays.asList(contigColumn, startColumn, endColumn));
        }
        
        locationColumnRemoveIndiciesInOrder.sort(Collections.reverseOrder());
    }

    //==================================================================================================================
    // Override Methods:

    @Override
    public String getContig() {
        return contig;
    }

    @Override
    public int getStart() {
        return start;
    }

    @Override
    public int getEnd() {
        return end;
    }

    @Override
    public boolean equals(final Object o) {
        if ( this == o ) return true;
        if ( o == null || getClass() != o.getClass() ) return false;

        final XsvTableFeature that = (XsvTableFeature) o;

        if ( contigColumn != that.contigColumn ) return false;
        if ( startColumn != that.startColumn ) return false;
        if ( endColumn != that.endColumn ) return false;
        if ( start != that.start ) return false;
        if ( end != that.end ) return false;
        if ( locationColumnRemoveIndiciesInOrder != null ? !locationColumnRemoveIndiciesInOrder.equals(that.locationColumnRemoveIndiciesInOrder) : that.locationColumnRemoveIndiciesInOrder != null )
            return false;
        if ( dataSourceName != null ? !dataSourceName.equals(that.dataSourceName) : that.dataSourceName != null )
            return false;
        if ( contig != null ? !contig.equals(that.contig) : that.contig != null ) return false;
        if ( columnNames != null ? !columnNames.equals(that.columnNames) : that.columnNames != null ) return false;
        return columnValues != null ? columnValues.equals(that.columnValues) : that.columnValues == null;
    }

    @Override
    public int hashCode() {
        int result = contigColumn;
        result = 31 * result + startColumn;
        result = 31 * result + endColumn;
        result = 31 * result + (locationColumnRemoveIndiciesInOrder != null ? locationColumnRemoveIndiciesInOrder.hashCode() : 0);
        result = 31 * result + (dataSourceName != null ? dataSourceName.hashCode() : 0);
        result = 31 * result + (contig != null ? contig.hashCode() : 0);
        result = 31 * result + start;
        result = 31 * result + end;
        result = 31 * result + (columnNames != null ? columnNames.hashCode() : 0);
        result = 31 * result + (columnValues != null ? columnValues.hashCode() : 0);
        return result;
    }

    @Override
    public String toString() {
        return "XsvTableFeature{" +
                "contigColumn=" + contigColumn +
                ", startColumn=" + startColumn +
                ", endColumn=" + endColumn +
                ", locationColumnRemoveIndiciesInOrder=" + locationColumnRemoveIndiciesInOrder +
                ", dataSourceName='" + dataSourceName + '\'' +
                ", contig='" + contig + '\'' +
                ", start=" + start +
                ", end=" + end +
                ", columnNames=" + columnNames +
                ", columnValues=" + columnValues +
                '}';
    }

    //==================================================================================================================
    // Static Methods:

    //==================================================================================================================
    // Instance Methods:

    /**
     * @return The name of the data source from which this {@link XsvTableFeature} was derived.
     */
    public String getDataSourceName() {
        return dataSourceName;
    }

    /**
     * Get a value from this {@link XsvTableFeature}.
     * @param key The key for the desired value.
     * @return The value for the given {@code key}; {@code null} if the given {@code key} is not in this {@link XsvTableFeature}.
     */
    public String get(final String key) {
        final int position = columnNames.indexOf(key);
        if ( position < 0 ) {
            return null;
        }
        return columnValues.get(position);
    }

    /**
     * Get a value from this {@link XsvTableFeature}.
     * @param index The index for the desired value.  Must be an available index in {@link XsvTableFeature#columnValues}.
     * @return The value at the given {@code index}.
     */
    public String get(final int index) {
        ParamUtils.inRange(index, 0, columnNames.size()-1, "Index out of range: " + index);
        return columnValues.get(index);
    }

    /**
     * Get the header from this {@link XsvTableFeature}.
     * @return The header for this {@link XsvTableFeature}.
     */
    public List<String> getHeader() {
        return columnNames;
    }

    /**
     * Get the header from this {@link XsvTableFeature} without the columns that contain location information.
     * Specifically the columns specified by the following fields are not included:
     *  {@link XsvTableFeature#contigColumn}
     *  {@link XsvTableFeature#startColumn}
     *  {@link XsvTableFeature#endColumn}
     * @return The header for this {@link XsvTableFeature} without location columns.
     */
    public List<String> getHeaderWithoutLocationColumns() {
        final List<String> outList = new ArrayList<>(columnNames);
        for ( final int indx : locationColumnRemoveIndiciesInOrder ) {
            outList.remove(indx);
        }
        return outList;
    }

    /**
     * Get the data from this {@link XsvTableFeature}.
     * @return The data for this {@link XsvTableFeature}.
     */
    public List<String> getValues() {
        return columnValues;
    }

    /**
     * Get the data from this {@link XsvTableFeature} without the columns that contain location information.
     * Specifically the columns specified by the following fields are not included:
     *  {@link XsvTableFeature#contigColumn}
     *  {@link XsvTableFeature#startColumn}
     *  {@link XsvTableFeature#endColumn}
     * @return The data for this {@link XsvTableFeature} without location columns.
     */
    public List<String> getValuesWithoutLocationColumns() {
        final List<String> outList = new ArrayList<>(columnValues);
        for ( final int indx : locationColumnRemoveIndiciesInOrder ) {
            outList.remove(indx);
        }
        return outList;
    }

    /**
     * @return The number of entries in this {@link XsvTableFeature}.
     */
    public int size() {
        return columnValues.size();
    }

    //==================================================================================================================
    // Helper Data Types:

}
