/*
* Copyright (c) 2012 The Broad Institute
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.hellbender.utils.collections;

import java.io.PrintStream;

/**
 * Wrapper around the basic NestedIntegerArray class that logs all updates (ie., all calls to put())
 * to the provided output stream. For testing/debugging purposes.
 *
 * Log entries are of the following form (fields are tab-separated):
 * LABEL    OPERATION    VALUE   KEY1    KEY2    ...     KEY_N
 *
 * A header line is written before the log entries giving the dimensions of this NestedIntegerArray.
 * It has the form:
 *
 * # LABEL    SIZE_OF_FIRST_DIMENSION    SIZE_OF_SECOND_DIMENSION    ...    SIZE_OF_NTH_DIMENSION
 */
public class LoggingNestedIntegerArray<T> extends NestedIntegerArray<T> {

    private PrintStream log;
    private String logEntryLabel;

    public static final String HEADER_LINE_PREFIX = "# ";
    public enum NestedIntegerArrayOperation { GET, PUT };

    /**
     *
     * @param log output stream to which to log update operations
     * @param logEntryLabel String that should be prefixed to each log entry
     * @param dimensions
     */
    public LoggingNestedIntegerArray(PrintStream log, String logEntryLabel, final int... dimensions) {
        super(dimensions);

        if ( log == null ) {
            throw new IllegalArgumentException("Log output stream must not be null");
        }
        this.log = log;
        this.logEntryLabel = logEntryLabel != null ? logEntryLabel : "";

        // Write the header line recording the dimensions of this NestedIntegerArray:
        StringBuilder logHeaderLine = new StringBuilder();

        logHeaderLine.append(HEADER_LINE_PREFIX);
        logHeaderLine.append(this.logEntryLabel);
        for ( int dimension : dimensions ) {
            logHeaderLine.append("\t");
            logHeaderLine.append(dimension);
        }

        this.log.println(logHeaderLine.toString());
    }

    @Override
    public T get( final int... keys ) {
        StringBuilder logEntry = new StringBuilder();

        logEntry.append(logEntryLabel);
        logEntry.append("\t");
        logEntry.append(NestedIntegerArrayOperation.GET);
        logEntry.append("\t");  // empty field for the datum value

        for ( int key : keys ) {
            logEntry.append("\t");
            logEntry.append(key);
        }

        log.println(logEntry.toString());

        return super.get(keys);
    }

    @Override
    public void put( final T value, final int... keys ) {
        StringBuilder logEntry = new StringBuilder();

        logEntry.append(logEntryLabel);
        logEntry.append("\t");
        logEntry.append(NestedIntegerArrayOperation.PUT);
        logEntry.append("\t");
        logEntry.append(value);
        for ( int key : keys ) {
            logEntry.append("\t");
            logEntry.append(key);
        }

        log.println(logEntry.toString());

        super.put(value, keys);
    }
}
