package org.broadinstitute.hellbender.utils.tsv;

import com.opencsv.CSVWriter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.Closeable;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;

/**
 * A simple helper class wrapper around CSVWriter that has the ingrained concept of a header line with indexed fields
 *
 * The expected use case for this class is that first {@link #addHeaderLine} is called with a list of the column names
 * which will be used to determine the number of columns per line as well as how the header is indexed. Then in order to
 * construct a new line call {@link #getNewLineBuilder} to get a line builder for each line, which then has convienent
 * methods for individually assigning column values based on the header line etc. Once a line is finished being mutated
 * one simply needs to call write() on the line to validate and finalize the line.
 *
 * Why didn't I use a tableWriter here? Who really holds the patent on the wheel anyway? Certainly not TableWriter.
 */
class SimpleCSVWriterWrapperWithHeader implements Closeable {
    private int expectedColumns;
    private Map<String, Integer> headerMap = null;
    private CSVWriter outputWriter;

    // The current incomplete line in the writer.
    private SimpleCSVWriterLineBuilder currentLine = null;

    /**
     * Creates a new table writer given the file and column names.
     *
     * @param path         the destination path.
     * @param separator    separator to use for the TSV file
     * @throws IOException              if one was raised when opening the the destination file for writing.
     */
    public SimpleCSVWriterWrapperWithHeader(final Path path, final char separator) throws IOException {
        this( new OutputStreamWriter(
                Files.newOutputStream(Utils.nonNull(path, "The path cannot be null."))),
        separator);
    }

    /**
     * Creates a new table writer given the file and column names.
     *
     * @param writer       the destination writer.
     * @param separator    separator to use for the TSV file
     * @throws IOException              if one was raised when opening the the destination file for writing.
     */
    public SimpleCSVWriterWrapperWithHeader(final Writer writer, final char separator) {
        outputWriter = new CSVWriter(writer, separator);
    }

    /**
     * Provides a header line to the CSV output file. Note that this will throw an exception if all header lines
     * are not unique as it attempts to create an index for the provided header lines for convenience when building
     * rows of the CSV.
     *
     * @param columns Ordered list of header lines to be built into the CSV
     */
    public void addHeaderLine(List<String> columns) {
        if (headerMap != null) {
            throw new GATKException("Should not be adding multiple header lines to a file");
        }
        outputWriter.writeNext(columns.toArray(new String[0]), false);
        expectedColumns = columns.size();

        // Create the mapping between header and column
        headerMap = new HashMap<>();
        for (int i = 0; i < columns.size(); i++) {
            Utils.nonNull(columns.get(i));
            if (headerMap.putIfAbsent(columns.get(i), i) != null) {
                throw new GATKException("Only allow unique column headings");
            }
        }
    }

    private void writeLine(String[] line) {
        outputWriter.writeNext(line, false);
        currentLine = null;
    }

    /**
     * Builds a new SimpleCSVWriterLineBuilder and writes out the previous line if it exists.
     *
     * @return a blank SimpleCSVWriterLineBuilder to allow for defining the next line
     */
    public SimpleCSVWriterLineBuilder getNewLineBuilder() {
        if (headerMap == null) {
            throw new GATKException("Cannot construct line without first setting the header line");
        }
        if (currentLine != null) {
            currentLine.write();
        }
        currentLine = new SimpleCSVWriterLineBuilder(this, expectedColumns);
        return currentLine;
    }

    /**
     * @param column header line to get index for
     * @return zero based index corresponding to that header string, throws an exception if the headerline doesn't exist
     */
    public Integer getIndexForColumn(String column) {
        Utils.nonNull(headerMap, "Cannot request column index if the header has not been specified");
        Integer index = headerMap.get(column);
        Utils.nonNull(index, "Requested column " + column + " does not exist in the provided header");
        return index;
    }

    @Override
    public void close() throws IOException {
        if (currentLine != null) {
            currentLine.write();
        }
        outputWriter.close();
    }

    /**
     * Helper to allow for incremental construction of a body line using either indexes or column headings
     * <p>
     * Calling build() will cause the line to be written out into the underlying CSV writer in its current state. Doing
     * so will result in a validation call where an exception will be thrown if any columns of the current line have
     * not been defined. fill() can be used to provide a default value for undefined columns.
     */
    public class SimpleCSVWriterLineBuilder {
        SimpleCSVWriterWrapperWithHeader thisBuilder;
        String[] lineToBuild;
        boolean hasBuilt = false;

        SimpleCSVWriterLineBuilder(SimpleCSVWriterWrapperWithHeader me, int lineLength) {
            thisBuilder = me;
            lineToBuild = new String[lineLength];
        }

        /**
         * @param row complete line corresponding to this row of the tsv
         */
        public SimpleCSVWriterLineBuilder setRow(final String[] row) {
            checkAlteration();
            Utils.validate(row.length == lineToBuild.length, "Provided line must have the correct number of columns");
            for (int i = 0; i < row.length; i++) {
                lineToBuild[i] = row[i];
            }
            return this;
        }

        /**
         * @param index Column index to be set
         * @param value Value to be placed into the line
         */
        public SimpleCSVWriterLineBuilder setColumn(final int index, final String value) {
            checkAlteration();
            lineToBuild[index] = value;
            return this;
        }

        /**
         * @param heading Column heading to be set
         * @param value   Value to be placed into the line
         */
        public SimpleCSVWriterLineBuilder setColumn(final String heading, final String value) {
            int index = thisBuilder.getIndexForColumn(heading);
            return setColumn(index, value);
        }

        /**
         * Fills in every empty column of the pending line with the provided value
         */
        public SimpleCSVWriterLineBuilder fill(final String filling) {
            checkAlteration();
            for (int i = 0; i < lineToBuild.length; i++) {
                if (lineToBuild[i] == null) {
                    lineToBuild[i] = filling;
                }
            }
            return this;
        }

        /**
         * Constructs the line and writes it out to the output
         */
        public void write() {
            Utils.validate(!Arrays.stream(lineToBuild).anyMatch(Objects::isNull), "Attempted to construct an incomplete line, make sure all columns are filled");
            thisBuilder.writeLine(lineToBuild);
            hasBuilt = true;
        }

        // Throw an exception if we try to alter an already written out line
        private void checkAlteration() {
            Utils.validate(!hasBuilt, "Cannot make alterations to an already written out CSV line");
        }
    }
}