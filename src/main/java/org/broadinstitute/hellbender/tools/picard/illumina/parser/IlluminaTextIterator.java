package org.broadinstitute.hellbender.tools.picard.illumina.parser;

import htsjdk.samtools.util.CloserUtil;

import org.broadinstitute.hellbender.utils.text.parsers.BasicInputParser;

import java.io.File;
import java.util.*;

import static htsjdk.samtools.util.CloserUtil.close;

/**
 * Abstract class for parsing text-based whitespace-delimited Illumina output files, organized
 * by tile.  Concrete implementations must call setFiles() in order to provide the list of files
 * to be iterated over.
 *
 * @author jburke@broadinstitute.org
 */
class IlluminaTextIterator implements Iterator<String[]> {

    // Location of illumina output files to be parsed
    private final int lane;
    private int currentTile = 0;

    // List of files of the given type, sorted by tile #
    private IlluminaFileMap files;

    private boolean treatGroupedDelimitersAsOne = true;
    private BasicInputParser parser;

    public IlluminaTextIterator(final int lane, final IlluminaFileMap files) {
        this.lane = lane;
        this.files = files;
        currentTile = files.firstKey();
    }

    public IlluminaTextIterator(final int lane, final IlluminaFileMap files,
                                final boolean treatGroupedDelimitersAsOne) {
        this(lane, files);
        this.treatGroupedDelimitersAsOne = treatGroupedDelimitersAsOne;
        currentTile = files.firstKey();
    }

    /**
     * Jump so that the next record returned will be the first one from the specified tile.
     */
    public void seekToTile(final int oneBasedTileNumber) {
        close(parser);
        currentTile = oneBasedTileNumber;
        initializeParser();
    }

    /**
     * Prepare to iterate.
     */
    private void initializeParser() {
        final List<File> fileSubset = files.getFilesStartingAt(currentTile);
        parser = new BasicInputParser(treatGroupedDelimitersAsOne, fileSubset.toArray(new File[fileSubset.size()]));
    }

    /**
     * Read the next record from the list of input files, and load into data argument.
     */
    @Override
    public String[] next() {
        if (!hasNext()) {
            throw new NoSuchElementException();
        }

        return parser.next();
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException("Remove is not supported by IlluminaTextIterator");
    }

    public boolean hasNext() {
        if (parser == null) initializeParser();
        return parser.hasNext();
    }

    protected int getLane() {
        return lane;
    }

    public String getCurrentFilename() {
        if (parser == null) initializeParser();
        return parser.getFileName();
    }

    protected void validateLane(final int lane) {
        if (lane != getLane()) {
            throw new IlluminaParserException("Lane number mismatch: " + lane + " != " + getLane());
        }
    }
}
