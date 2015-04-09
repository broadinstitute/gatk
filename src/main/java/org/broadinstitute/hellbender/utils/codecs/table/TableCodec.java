package org.broadinstitute.hellbender.utils.codecs.table;

import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.readers.LineIterator;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;

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
 * <h2>File format example</h2>
 * <pre>
 *     HEADER a b c
 *     1:1  1   2   3
 *     1:2  4   5   6
 *     1:3  7   8   9
 * </pre>
 */
public final class TableCodec extends AsciiFeatureCodec<TableFeature> {
    protected final static String delimiterRegex = "\\s+";
    protected final static String headerDelimiter = "HEADER";
    protected final static String igvHeaderDelimiter = "track";
    protected final static String commentDelimiter = "#";

    protected List<String> header = new ArrayList<>();

    public TableCodec() {
        super(TableFeature.class);
    }

    @Override
    public TableFeature decode(String line) {
        if (line.startsWith(headerDelimiter) || line.startsWith(commentDelimiter) || line.startsWith(igvHeaderDelimiter)) {
            return null;
        }
        String[] split = line.split(delimiterRegex);
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
            if ( isFirst && ! line.startsWith(headerDelimiter) && ! line.startsWith(commentDelimiter)) {
                throw new UserException.MalformedFile("TableCodec file does not have a header");
            }
            isFirst &= line.startsWith(commentDelimiter);
            if (line.startsWith(headerDelimiter)) {
                reader.next(); // "Commit" the peek
                if (header.size() > 0) {
                    throw new UserException.MalformedFile("Input table file seems to have two header lines.  The second is = " + line);
                }
                final String spl[] = line.split(delimiterRegex);
                Collections.addAll(header, spl);
                return header;
            } else if (line.startsWith(commentDelimiter)) {
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
