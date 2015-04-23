package org.broadinstitute.hellbender.utils.text;

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
 * Please use this class for reading lines in a file.
 */
public final class XReadLines implements Iterator<String>, Iterable<String>, AutoCloseable {
    private final BufferedReader in;      // The stream we're reading from
    private String nextLine = null;       // Return value of next call to next()
    private final boolean trimWhitespace;
    private final String commentPrefix;

    /**
     * Opens the given file for reading lines.
     * The file may be a text file or a gzipped text file (the distinction is made by the file extension).
     * By default, it will trim whitespaces.
     * @throws IOException if an IO error occurs during reading the file.
     */
    public XReadLines(final File filename) throws IOException {
        this(filename, true);
    }

    /**
     * Opens the given file for reading lines and optionally trim whitespaces.
     * The file may be a text file or a gzipped text file (the distinction is made by the file extension).
     * @throws IOException if an IO error occurs during reading the file.
     */
    public XReadLines(final File filename, final boolean trimWhitespace) throws IOException {
        this(IOUtils.makeReaderMaybeGzipped(filename), trimWhitespace, null);
    }

    /**
     * Creates a new xReadLines object to read lines from an bufferedReader
     *
     * @param reader file name
     * @param trimWhitespace trim whitespace
     * @param commentPrefix prefix for comments or null if no prefix is set
     * @throws IOException if an IO error occurs during reading the file.
     */
    public XReadLines(final Reader reader, final boolean trimWhitespace, final String commentPrefix) throws IOException {
        this.in = (reader instanceof BufferedReader) ? (BufferedReader)reader : new BufferedReader(reader);
        this.trimWhitespace = trimWhitespace;
        this.commentPrefix = commentPrefix;
        this.nextLine = readNextLine();
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