package org.broadinstitute.hellbender.utils.codecs.table;

import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.readers.LineIterator;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Reads tab deliminated tabular text files
 *
 * <p>
 *     <ul>
 *     <li>Header: must begin with line HEADER or track (for IGV), followed by any number of column names,
 *     separated by whitespace.</li>
 *     <li>Header: Custom header delimiters can be provided, with a null header line being interpreted as having a non-delimeted
 *     header which consists of one line.</li>
 *     <li>Comment lines starting with # are ignored</li>
 *     <li>Each non-header and non-comment line is split into parts by whitespace,
 *     and these parts are assigned as a map to their corresponding column name in the header.
 *     Note that the first element (corresponding to the HEADER column) must be a valid genome loc
 *     such as 1, 1:1 or 1:1-10, which is the position of the Table element on the genome.  TableCodec
 *     requires that there be one value for each column in the header, and no more, on all lines.</li>
 *     </ul>
 * </p>
 *
 * </p>
 *
 * <h2>File format example 1</h2>
 * <pre>
 *     HEADER a b c
 *     1:1  1   2   3
 *     1:2  4   5   6
 *     1:3  7   8   9
 * </pre>
 *
 * <h2>File format example 2</h2>
 * <pre>
 *     a b c
 *     1:1  1   2   3
 *     1:2  4   5   6
 *     1:3  7   8   9
 * </pre>
 */
public final class TableCodec extends AsciiFeatureCodec<TableFeature> {
    protected static final String DEFAULT_HEADER_DELIMITER = "HEADER";
    protected static final String IGV_HEADER_DELIMITER = "track";
    protected static final String COMMENT_DELIMITER = "#";

    private final String headerDelimiter;

    protected String delimiter_regex = "\\s+";

    protected List<String> header = new ArrayList<>();

    private boolean havePassedHeader = false;

    /**
     * Create a TableCodec with a configured header line delimiter
     *
     * @param headerLineDelimiter the delimeter for comment header lines, or null if the header is a single commented line-
     */
    public TableCodec(final String headerLineDelimiter) {
        super(TableFeature.class);
        if ( "".equals(headerLineDelimiter) ) {
            // Note, it is valid for headerLineDelimiter to be null, just not empty as the regex breaks in that case.
            throw new GATKException("HeaderLineDelimiter must either be a valid delimiter or null");
        }
        headerDelimiter = headerLineDelimiter;
    }

    /**
     * Create a TableCodec for IGV track data.
     */
    public TableCodec() {
        this(DEFAULT_HEADER_DELIMITER);
    }

    @Override
    public TableFeature decode(final String line) {
        if ((headerDelimiter != null && line.startsWith(headerDelimiter)) ||
                (headerDelimiter == null && !havePassedHeader) ||
                line.startsWith(COMMENT_DELIMITER) || line.startsWith(IGV_HEADER_DELIMITER)) {
            havePassedHeader = true;
            return null;
        }
        final String[] split = line.split(delimiter_regex);
        if (split.length < 1) {
            throw new IllegalArgumentException("TableCodec line = " + line + " is not a valid table format");
        }
        return new TableFeature(new SimpleInterval(split[0]), Arrays.asList(split), header);
    }

    @Override
    public List<String> readActualHeader(final LineIterator reader) {
        boolean isFirst = true;
        while (reader.hasNext()) {
            final String line = reader.peek(); // Peek to avoid reading non-header data
            if ( isFirst && ! line.startsWith(COMMENT_DELIMITER) &&  headerDelimiter != null && ! line.startsWith(headerDelimiter) ) {
                throw new UserException.MalformedFile("TableCodec file does not have a header");
            }
            isFirst &= line.startsWith(COMMENT_DELIMITER);
            if (headerDelimiter == null || line.startsWith(headerDelimiter)) {
                reader.next(); // "Commit" the peek
                if (!header.isEmpty()) {
                    throw new UserException.MalformedFile("Input table file seems to have two header lines.  The second is = " + line);
                }
                final String[] spl = line.split(delimiter_regex);
                Collections.addAll(header, spl[0].equals(headerDelimiter) ? Arrays.copyOfRange(spl, 1, spl.length) : spl);
                return header;
            } else if (line.startsWith(COMMENT_DELIMITER)) {
                reader.next(); // "Commit" the peek
            } else {
                break;
            }
        }
        return header;
    }

    @Override
    public boolean canDecode(final String path) {
        return path.toLowerCase().endsWith(".table");
    }
}
