package org.broadinstitute.hellbender.utils.codecs.table;

import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.Feature;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.List;

/**
 * A feature representing a single row out of a text table.
 * Note: stop positions not explicitly specified in the
 * interval column of records will be filled in with {@link Integer#MAX_VALUE}.
 */
public final class TableFeature implements Feature {
    // stores the values for the columns separated out
    private final List<String> values;

    // if we have column names, we store them here
    private final List<String> keys;

    // our location
    private final Locatable position;

    public TableFeature(Locatable position, List<String> values, List<String> keys) {
        this.values = values;
        this.keys = keys;
        this.position = position;
    }

    @Override
    @SuppressWarnings("deprecation")
    //Suppressing because we have to implement this
    @Deprecated
    public String getChr() {
        return getContig();
    }

    @Override
    public String getContig() {
        return position.getContig();
    }

    @Override
    public int getStart() {
        return position.getStart();
    }

    @Override
    public int getEnd() {
        return position.getEnd();
    }

    public int columnCount(){
        return values.size();
    }

    public String getValue(int columnPosition) {
        if (columnPosition >= columnCount()) throw new IllegalArgumentException("We only have " + columnCount() + " columns, the requested column = " + columnPosition);
        return values.get(columnPosition);
    }

    public String toString() {
        return String.format("%s\t%s", position.toString(), Utils.join("\t", values));
    }

    public String get(String columnName) {
        int position = keys.indexOf(columnName);
        if (position < 0) throw new IllegalArgumentException("We don't have a column named " + columnName);
        return values.get(position);
    }

    public Locatable getLocation() {
        return this.position;
    }

    public List<String> getAllValues() {
        return getValuesTo(values.size());
    }

    public List<String> getValuesTo(int columnPosition) {
        return values.subList(0, columnPosition);
    }

    public List<String> getHeader() {
        return keys;
    }
}
