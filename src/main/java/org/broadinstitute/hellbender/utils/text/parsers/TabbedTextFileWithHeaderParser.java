package org.broadinstitute.hellbender.utils.text.parsers;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.RuntimeIOException;

import java.io.File;
import java.util.*;

/**
 * Parse a tabbed text file in which columns are found by looking at a header line rather than by position.
 *
 * @author alecw@broadinstitute.org
 */
public final class TabbedTextFileWithHeaderParser implements Iterable<TabbedTextFileWithHeaderParser.Row> {
    public final class Row {
        private final String[] fields;
        private final String currentLine;

        Row(final String[] fields, final String source) {
            this.fields = fields;
            this.currentLine = source;
        }

        /**
         * @return Array of fields in the order they appear in the file.
         */
        public String[] getFields() {
            return fields;
        }

        public String getField(final String columnLabel) {
            final Integer key = columnLabelIndices.get(columnLabel);
            if (key == null) throw new NoSuchElementException(String.format("column %s in %s", columnLabel, parser.getFileName()));
            return fields[key];
        }

        public Integer getIntegerField(final String columnLabel) {
            if (fields[columnLabelIndices.get(columnLabel)] == null)  return null;
            return Integer.parseInt(fields[columnLabelIndices.get(columnLabel)]);
        }

        public String getCurrentLine() {
            return this.currentLine;
        }
    }

    class TheIterator implements CloseableIterator<Row> {

        @Override
        public boolean hasNext() {
            return parser.hasNext();
        }

        @Override
        public Row next() {
            final String[] fields = parser.next();
            final String source = parser.getCurrentLine();
            return new Row(fields, source);
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException();
        }

        @Override
        public void close() {
            extantIterator = null;
        }
    }

    /**
     * Map from column label to positional index.
     */
    private final Map<String, Integer> columnLabelIndices = new HashMap<>();
    private final TabbedInputParser parser;
    private TheIterator extantIterator;

    public TabbedTextFileWithHeaderParser(final TabbedInputParser parser) {
        this.parser = parser;
        if (!parser.hasNext()) {
            throw new RuntimeIOException("No header line found in file " + parser.getFileName());
        }
        final String[] columnLabels = parser.next();
        for (int i = 0; i < columnLabels.length; ++i) {
            columnLabelIndices.put(columnLabels[i], i);
        }
    }
    
    public TabbedTextFileWithHeaderParser(final File file) {
        this(new TabbedInputParser(false, file));
    }

    public TabbedTextFileWithHeaderParser(final File file, final String[] columnHeaders) {
        parser = new TabbedInputParser(false, file);
        if (!parser.hasNext()) {
            throw new RuntimeIOException("No header line found in file " + file);
        }

        for (int i = 0; i < columnHeaders.length; ++i) {
            columnLabelIndices.put(columnHeaders[i], i);
        }
    }

    /**
     * @param columnLabel
     * @return True if the given column label appears in the header.
     */
    public boolean hasColumn(final String columnLabel) {
        return columnLabelIndices.containsKey(columnLabel);
    }

    /**
     *
     * @return The set of column labels for this file in no particular order.
     */
    public Set<String> columnLabels() {
        return columnLabelIndices.keySet();
    }

    /**
     * Creates the iterator object.  It is illegal to have more than one iterator extant
     * on the same parser object.
     */
    @Override
    public CloseableIterator<Row> iterator() {
        if (extantIterator != null) {
            throw new ConcurrentModificationException("Only one iterator allowed at a time.");
        }
        extantIterator = new TheIterator();
        return extantIterator;
    }

    /**
     * Release all resources associated with the parser.  Iteration will not work after this
     * has been called.
     */
    public void close() {
        parser.close();
    }

    public int getCurrentLineNumber() {
        return parser.getCurrentLineNumber();
    }

    public Set<String> getColumnNames() {
        return Collections.unmodifiableSet(this.columnLabelIndices.keySet());
    }
}
