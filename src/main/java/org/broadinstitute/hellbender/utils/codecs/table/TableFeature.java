package org.broadinstitute.hellbender.utils.codecs.table;

import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.Feature;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.List;

/**
 * Feature representing a row in a text table.
 *
 * <p>Note: stop positions not explicitly specified in the interval column will be set to {@link Integer#MAX_VALUE}.
 *
 * @see TableCodec for more information about the table format.
 */
public final class TableFeature implements Feature {
    // stores the values for the columns separated out
    private final List<String> values;

    // if we have column names, we store them here
    private final List<String> keys;

    // our location
    private final Locatable position;

    /**
     * Public constructor.
     *
     * @param position the coordinates for the feature.
     * @param values values for each of the columns.
     * @param keys column names.
     */
    public TableFeature(final Locatable position, final List<String> values, final List<String> keys) {
        this.values = values;
        this.keys = keys;
        this.position = position;
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

    /** Gets the number of columns in the feature. */
    public int columnCount(){
        return values.size();
    }

    /** Gets the value associated with the ith column (first column is 0). */
    public String getValue(int columnPosition) {
        Utils.validateArg(columnPosition >= 0, () -> "Requested a negative column: " + columnPosition);
        Utils.validateArg(columnPosition < columnCount(), () -> "We only have " + columnCount() + " columns, the requested column = " + columnPosition);
        return values.get(columnPosition);
    }

    /** Format as a tab-delimited row. */
    public String toString() {
        return String.format("%s\t%s", position.toString(), Utils.join("\t", values));
    }

    /** Gets the value associated with the key. */
    public String get(String columnName) {
        int position = keys.indexOf(columnName);
        Utils.validateArg(position >= 0, () -> "We don't have a column named " + columnName);
        return values.get(position);
    }

    /** Gets the position for the feature. */
    public Locatable getLocation() {
        return this.position;
    }

    /** Gets all the values in order.  */
    public List<String> getAllValues() {
        return getValuesTo(columnCount());
    }

    /** Gets all the values from the first column (0) to the last one specified (exclusive). */
    public List<String> getValuesTo(int columnPosition) {
        Utils.validateArg(columnPosition >= 0, () -> "Requested a negative column: " + columnPosition);
        Utils.validateArg(columnPosition <= columnCount(), () -> "We only have " + columnCount() + " columns, the requested column = " + columnPosition);
        return values.subList(0, columnPosition);
    }

    /** Gets the column names. */
    public List<String> getHeader() {
        return keys;
    }
}
