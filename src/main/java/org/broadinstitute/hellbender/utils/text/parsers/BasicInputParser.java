package org.broadinstitute.hellbender.utils.text.parsers;

import htsjdk.samtools.util.BufferedLineReader;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;

import java.io.File;
import java.io.InputStream;
import java.util.*;

/**
 * TextFileParser which reads a single text file.
 *
 * @author Kathleen Tibbetts
 */
public class BasicInputParser extends AbstractInputParser {
    private BufferedLineReader reader;
    private final ArrayList<InputStream> inputs = new ArrayList<>();
    private final ArrayList<String> fileNames = new ArrayList<>();
    String currentFileName = null;
    private String currentLine = null;
    private String nextLine = null;
    private int currentLineNumber = 0;
    private int nextLineNumber = 0;

    /**
     * Constructor.  Opens up a buffered reader and reads the first line.
     *
     * @param inputStreams  the file(s) to parse, in order
     */
    public BasicInputParser(final boolean treatGroupedDelimitersAsOne, final InputStream... inputStreams) {
        if (inputStreams.length == 0) {
            throw new IllegalArgumentException("At least one input must be specified.");
        }
        this.inputs.addAll(Arrays.asList(inputStreams));
        reader = new BufferedLineReader(this.inputs.remove(0));
        this.setTreatGroupedDelimitersAsOne(treatGroupedDelimitersAsOne);
    }

    public BasicInputParser(final boolean treatGroupedDelimitersAsOne, final int wordCount, final InputStream... inputStreams) {
        this(treatGroupedDelimitersAsOne, inputStreams);
        setWordCount(wordCount);
    }

    /**
     * Constructor.  Opens up a buffered reader and reads the first line.
     *
     * @param files  the file(s) to parse, in order
     */
    public BasicInputParser(final boolean treatGroupedDelimitersAsOne, final File... files) {
        this(treatGroupedDelimitersAsOne, filesToInputStreams(files));
        for (File f : files) fileNames.add(f.getAbsolutePath());
        this.currentFileName = fileNames.remove(0);
    }

    /**
     * Constructor.  In addition to opening and priming the files, it sets the number of
     * whitespace-separated "words" per line.
     *
     * @param files      the file(s) to parse
     * @param wordCount number of whitespace-separated "words" per line
     */
    public BasicInputParser(final boolean treatGroupedDelimitersAsOne, final int wordCount, final File... files) {
        this(treatGroupedDelimitersAsOne, files);
        setWordCount(wordCount);
    }

    /**
     * Workhorse method that reads the next line from the underlying reader
     *
     * @return  String or null if there is no next line
     */
    protected byte[] readNextLine() {
        final String line = reader.readLine();
        if (nextLine != null && !isComment(nextLine.getBytes())) {
            currentLineNumber = nextLineNumber;
            currentLine = nextLine;
        }
        if (line != null) {
            nextLineNumber++;
            nextLine = line;
            return line.getBytes();
        }
        if (inputs.size() > 0) {
            advanceFile();
            return readNextLine();
        }
        return null;
    }

    protected void advanceFile() {
        currentFileName = fileNames.size() > 0 ? fileNames.remove(0) : null;
        nextLineNumber = 0;
        nextLine = null;
        reader = new BufferedLineReader(inputs.remove(0));
    }

    /**
     * Closes the underlying stream
     */
    public void close() {
        if (reader != null)  {
            reader.close();
        }
        for(final InputStream stream : inputs){
            CloserUtil.close(stream);
        }
    }

    /**
     * Gets the name of the file being parsed
     *
     * @return  the name of the file being parsed
     */
    public String getFileName() {
        return this.currentFileName != null ? this.currentFileName : "(file name unavailable)";
    }

    /**
     * Provides access to the current (just parsed) line in pre-parsed format.
     * NOTE: Because AbstractInputParser pre-fetches the next line, this method actually returns the
     * next line, not the most recent line returned by next().
     */
    public String getCurrentLine() {
        return this.currentLine;
    }

    /**
     * NOTE: Because AbstractInputParser pre-fetches the next line, this method actually returns the
     * next line, not the most recent line returned by next().
     */
    public int getCurrentLineNumber() {
        return currentLineNumber;
    }

    private static InputStream[] filesToInputStreams(final File files[]) {
        final InputStream result[] = new InputStream[files.length];
        for (int i = 0; i < files.length; i++) {
            result[i] = IOUtil.openFileForReading(files[i]);
        }
        return result;
    }
}
