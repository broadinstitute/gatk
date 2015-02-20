package org.broadinstitute.hellbender.utils.text;

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

import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.*;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

/**
 * Support for Python-like xreadlines() function as a class.  This is an iterator and iterable over
 * Strings, each corresponding a line in the file (minus newline).  Enables the very simple accessing
 * of lines in a file as:
 *
 * xReadLines reader = new xReadLines(new File(file_name));
 * List<String> lines = reader.readLines();
 * reader.close();
 *
 * or
 *
 * for ( String line : new xReadLines(new File(file_name)) {
 *   doSomeWork(line);
 * }
 *
 * For the love of god, please use this system for reading lines in a file.
 */
public class XReadLines implements Iterator<String>, Iterable<String> {
    private final BufferedReader in;      // The stream we're reading from
    private String nextLine = null;       // Return value of next call to next()
    private final boolean trimWhitespace;
    private final String commentPrefix;

    public XReadLines(final File filename) throws IOException {
        this(IOUtils.makeReaderMaybeGzipped(filename), true, null);
    }

    public XReadLines(final File filename, final boolean trimWhitespace) throws IOException {
        this(IOUtils.makeReaderMaybeGzipped(filename), trimWhitespace, null);
    }

    /**
     * Creates a new xReadLines object to read lines from filename
     *
     * @param filename file name (may be gzipped if name ends on ".gz")
     * @param trimWhitespace trim whitespace
     * @param commentPrefix prefix for comments or null if no prefix is set
     * @throws FileNotFoundException when the file is not found
     */
    public XReadLines(final File filename, final boolean trimWhitespace, final String commentPrefix) throws IOException {
        this(IOUtils.makeReaderMaybeGzipped(filename), trimWhitespace, commentPrefix);
    }

    public XReadLines(final InputStream inputStream) throws FileNotFoundException {
        this(new InputStreamReader(inputStream), true, null);
    }

    public XReadLines(final InputStream inputStream, final boolean trimWhitespace) {
        this(new InputStreamReader(inputStream), trimWhitespace, null);
    }

    /**
     * Creates a new xReadLines object to read lines from an input stream
     *
     * @param inputStream input stream
     * @param trimWhitespace trim whitespace
     * @param commentPrefix prefix for comments or null if no prefix is set
     */
    public XReadLines(final InputStream inputStream, final boolean trimWhitespace, final String commentPrefix) {
        this(new InputStreamReader(inputStream), trimWhitespace, commentPrefix);
    }


    /**
     * Creates a new xReadLines object to read lines from a reader
     *
     * @param reader reader
     */
    public XReadLines(final Reader reader) {
        this(reader, true, null);
    }

    /**
     * Creates a new xReadLines object to read lines from an reader
     *
     * @param reader reader
     * @param trimWhitespace trim whitespace
     */
    public XReadLines(final Reader reader, final boolean trimWhitespace) {
        this(reader, trimWhitespace, null);
    }

    /**
     * Creates a new xReadLines object to read lines from an bufferedReader
     *
     * @param reader file name
     * @param trimWhitespace trim whitespace
     * @param commentPrefix prefix for comments or null if no prefix is set
     */
    public XReadLines(final Reader reader, final boolean trimWhitespace, final String commentPrefix) {
        this.in = (reader instanceof BufferedReader) ? (BufferedReader)reader : new BufferedReader(reader);
        this.trimWhitespace = trimWhitespace;
        this.commentPrefix = commentPrefix;
        try {
            this.nextLine = readNextLine();
        } catch(IOException e) {
            throw new IllegalArgumentException(e);
        }
    }

    /**
     * Reads all of the lines in the file, and returns them as a list of strings
     *
     * @return all of the lines in the file.
     */
    public List<String> readLines() {
        List<String> lines = new LinkedList<String>();
        for ( String line : this ) {
            lines.add(line);
        }
        return lines;
    }

    /**
     * I'm an iterator too...
     * @return an iterator
     */
    public Iterator<String> iterator() {
        return this;
    }

    public boolean hasNext() {
        return this.nextLine != null;
    }

    /**
     * Actually reads the next line from the stream, not accessible publicly
     * @return the next line or null
     * @throws IOException if an error occurs
     */
    private String readNextLine() throws IOException {
        String nextLine;
        while ((nextLine = this.in.readLine()) != null) {
            if (this.trimWhitespace) {
                nextLine = nextLine.trim();
                if (nextLine.length() == 0)
                    continue;
            }
            if (this.commentPrefix != null)
                if (nextLine.startsWith(this.commentPrefix))
                    continue;
            break;
        }
        return nextLine;
    }

    /**
     * Returns the next line (optionally minus whitespace)
     * @return the next line
     */
    public String next() {
        try {
            String result = this.nextLine;
            this.nextLine = readNextLine();

            // If we haven't reached EOF yet
            if (this.nextLine == null) {
                in.close();             // And close on EOF
            }

            // Return the line we read last time through.
            return result;
        } catch(IOException e) {
            throw new IllegalArgumentException(e);
        }
    }

    // The file is read-only; we don't allow lines to be removed.
    public void remove() {
        throw new UnsupportedOperationException();
    }

    public void close() throws IOException {
        this.in.close();
    }
}
