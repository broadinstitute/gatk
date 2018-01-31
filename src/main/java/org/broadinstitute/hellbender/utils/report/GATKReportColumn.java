package org.broadinstitute.hellbender.utils.report;

import org.apache.commons.lang3.math.NumberUtils;

import java.util.Arrays;
import java.util.Collection;

/**
 * column information within a GATK report table
 */
public final class GATKReportColumn {
    private final String columnName;
    private final String format;
    private final GATKReportDataType dataType;

    private GATKReportColumnFormat columnFormat;
    private GATKReportColumnFormat.Alignment alignment = GATKReportColumnFormat.Alignment.RIGHT;  // default alignment is to the right unless values added ask for a left alignment
    private int maxWidth = 0;

    /**
     * Construct the column object, specifying the column name, default value, whether or not the column should be
     * displayed, and the format string. This cannot be null.
     *
     * @param columnName   the name of the column
     * @param format       format string
     */
    public GATKReportColumn(final String columnName, final String format) {
        this.columnName = columnName;
        this.maxWidth = columnName.length();
        if ( format.equals("") ) {
            this.format = "%s";
            this.dataType = GATKReportDataType.Unknown;
        }
        else {
            this.format = format;
            this.dataType = GATKReportDataType.fromFormatString(format);
        }
    }

    /**
     * Get the display width for this column.  This allows the entire column to be displayed with the appropriate, fixed
     * width.
     *
     * @return the format string for this column
     */
    public GATKReportColumnFormat getColumnFormat() {
        if (columnFormat != null)
            return columnFormat;

        columnFormat = new GATKReportColumnFormat(maxWidth, alignment);
        return columnFormat;
    }

    private static final Collection<String> RIGHT_ALIGN_STRINGS = Arrays.asList(
            "null",
            "NA",
            String.valueOf(Double.POSITIVE_INFINITY),
            String.valueOf(Double.NEGATIVE_INFINITY),
            String.valueOf(Double.NaN));

    /**
     * Check if the value can be right aligned. Does not trim the values before checking if numeric since it assumes
     * the spaces mean that the value is already padded.
     *
     * @param value to check
     * @return true if the value is a right alignable
     */
    protected static boolean isRightAlign(final String value) {
        return value == null || RIGHT_ALIGN_STRINGS.contains(value) || NumberUtils.isCreatable(value.trim());
    }

    /**
     * Returns a string version of the values.
     *
     * @param obj The object to convert to a string
     * @return The string representation of the column
     */
    private String formatValue(final Object obj) {
        String value;
        if (obj == null) {
            value = "null";
        }
        else if ( dataType.equals(GATKReportDataType.Unknown) && (obj instanceof Double || obj instanceof Float) ) {
            value = String.format("%.8f", obj);
        }
        else
            value = String.format(format, obj);

        return value;
    }

    public GATKReportDataType getDataType() {
        return dataType;
    }

    public String getColumnName() {
        return columnName;
    }

    public String getFormat() {
        return dataType.equals(GATKReportDataType.Unknown) ? "%s" : format;
    }

    public void updateFormatting(final Object value) {
        if (value != null) {
            final String formatted = formatValue(value);
            if (!formatted.isEmpty()) {
                updateMaxWidth(formatted);
                updateFormat(formatted);
            }
        }
    }

    private void updateMaxWidth(final String formatted) {
        maxWidth = Math.max(formatted.length(), maxWidth);
    }

    private void updateFormat(final String formatted) {
        if (alignment == GATKReportColumnFormat.Alignment.RIGHT)
            alignment = isRightAlign(formatted) ? GATKReportColumnFormat.Alignment.RIGHT : GATKReportColumnFormat.Alignment.LEFT;
    }
}
