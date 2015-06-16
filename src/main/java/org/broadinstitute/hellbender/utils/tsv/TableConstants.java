package org.broadinstitute.hellbender.utils.tsv;

/**
 * Common constants for table readers and writers.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class TableConstants {

    /**
     * Column separator {@value #COLUMN_SEPARATOR_STRING}.
     */
    public static final char COLUMN_SEPARATOR = '\t';

    /**
     * Column separator as an string.
     */
    public static final String COLUMN_SEPARATOR_STRING = "" + COLUMN_SEPARATOR;

    /**
     * Comment line prefix string {@value}.
     * <p>
     * Lines that start with this prefix (spaces are not ignored), will be considered comment
     * lines (neither a header line nor data line).
     * </p>
     */
    public static final String COMMENT_PREFIX = "#";

    /**
     * Quote character {@value #QUOTE_STRING}.
     * <p>
     * Character used to quote table values that contain special characters.
     * </p>
     */
    public static final char QUOTE_CHARACTER = '\"';

    /**
     * Quote character as a string.
     */
    public static final String QUOTE_STRING = "" + QUOTE_CHARACTER;

    /**
     * Escape character {@value #ESCAPE_STRING}.
     * <p>
     * Within {@value #QUOTE_STRING quotes},
     * the user must prepend this character when including {@value #QUOTE_STRING quotes}
     * or the {@value #ESCAPE_STRING escape}
     * character itself as part of the string.
     * </p>
     */
    public static final char ESCAPE_CHARACTER = '\\';

    /**
     * Escape character as a string.
     */
    public static final String ESCAPE_STRING = "" + ESCAPE_CHARACTER;

    /**
     * Declared to make instantiation impossible.
     */
    private TableConstants() {
        throw new UnsupportedOperationException();
    }
}
