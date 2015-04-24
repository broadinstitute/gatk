package org.broadinstitute.hellbender.utils.text.parsers;

import java.io.File;
import java.io.InputStream;

/**
 * Parser for tab-delimited files
 *
 * @author Kathleen Tibbetts
 */
public final class TabbedInputParser extends BasicInputParser {

    /**
     * Constructor
     *
     * @param stream  The input stream(s) to parse
     */
    public TabbedInputParser(boolean treatGroupedDelimitersAsOne, InputStream... stream) {
        super(treatGroupedDelimitersAsOne, stream);
    }

    /**
     * Constructor
     *
     * @param file  The file(s) to parse
     */
    public TabbedInputParser(boolean treatGroupedDelimitersAsOne, File... file) {
        super(treatGroupedDelimitersAsOne, file);
    }

    /**
     * Determines whether a given character is a delimiter
     *
     * @param b the character to evaluate
     * @return  true if <code>b</code> is a delimiter; otherwise false
     */
    @Override
    protected boolean isDelimiter(final byte b) {
        return b == '\t';
    }
}
