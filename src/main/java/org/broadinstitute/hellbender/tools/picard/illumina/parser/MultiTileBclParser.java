package org.broadinstitute.hellbender.tools.picard.illumina.parser;

import htsjdk.samtools.util.CloseableIterator;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.readers.BclQualityEvaluationStrategy;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.readers.BclReader;

import java.io.File;
import java.util.*;

import static org.broadinstitute.hellbender.tools.picard.illumina.parser.readers.BclReader.makeSeekable;

/**
 * Parse .bcl.bgzf files that contain multiple tiles in a single file.  This requires an index file that tells
 * the bgzf virtual file offset of the start of each tile in the block-compressed bcl file.
 */
public class MultiTileBclParser extends BclParser {
    private final TileIndex tileIndex;
    private MultiTileBclDataCycleFileParser cycleFileParser = null;

    public MultiTileBclParser(final File directory, final int lane, final CycleIlluminaFileMap tilesToCycleFiles,
                              final OutputMapping outputMapping, final boolean applyEamssFilter,
                              final BclQualityEvaluationStrategy bclQualityEvaluationStrategy,
                              final TileIndex tileIndex) {
        super(directory, lane, tilesToCycleFiles, outputMapping, applyEamssFilter, bclQualityEvaluationStrategy);
        this.tileIndex = tileIndex;
        this.initialize();
    }

    @Override
    public void initialize() {
        if (tileIndex != null) {
            seekToTile(currentTile);
        }
    }

    private CountLimitedIterator makeReader(final List<File> files) {
        if (tileIndex != null) {
            final BclReader bclReader = makeSeekable(files, bclQualityEvaluationStrategy, outputMapping.getOutputReadLengths());
            final int numClustersInTile = bclReader.seek(files, tileIndex, currentTile);
            return new CountLimitedIterator(bclReader, numClustersInTile);
        } else {
            return null;
        }
    }

    @Override
    protected CycleFilesParser<BclData> makeCycleFileParser(final List<File> files) {
        if (cycleFileParser == null) {
            cycleFileParser = new MultiTileBclDataCycleFileParser(files, currentTile);
        } else {
            final int numClustersInTile = cycleFileParser.getReader().seek(files, tileIndex, currentTile);
            cycleFileParser.setCurrentTile(currentTile);
            cycleFileParser.resetClusterLimit(numClustersInTile);
        }
        return cycleFileParser;
    }

    /**
     * An iterator wrapper that stops when it has return a pre-determined number of records even if the underlying
     * iterator still had more records.
     */
    static class CountLimitedIterator implements CloseableIterator<BclData> {
        public BclReader getUnderlyingIterator() {
            return underlyingIterator;
        }

        private final BclReader underlyingIterator;
        private int recordLimit;
        private int numRecordsRead = 0;

        CountLimitedIterator(final BclReader underlyingIterator, final int recordLimit) {
            this.underlyingIterator = underlyingIterator;
            this.recordLimit = recordLimit;
        }

        @Override
        public void close() {
            //underlyingIterator.close();
        }

        @Override
        public boolean hasNext() {
            return numRecordsRead < recordLimit && underlyingIterator.hasNext();
        }

        @Override
        public BclData next() {
            if (!hasNext()) throw new NoSuchElementException();
            ++numRecordsRead;
            return underlyingIterator.next();
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException();
        }
    }


    private class MultiTileBclDataCycleFileParser implements CycleFilesParser<BclData> {
        final CountLimitedIterator reader;
        int currentTile;

        public MultiTileBclDataCycleFileParser(final List<File> files, final int currentTile) {
            this.currentTile = currentTile;
            reader = makeReader(files);
        }

        @Override
        public void close() {
            reader.close();
        }

        @Override
        public BclData next() {
            if (!hasNext()) {
                throw new NoSuchElementException();
            }
            return reader.next();
        }

        @Override
        public boolean hasNext() {
            try {
                return reader.hasNext();
            } catch (final NullPointerException npe) {
                return false;
            }
        }

        public BclReader getReader() {
            return reader.getUnderlyingIterator();
        }

        public void resetClusterLimit(final int numClustersInTile) {
            reader.recordLimit = numClustersInTile;
            reader.numRecordsRead = 0;
        }

        public void setCurrentTile(final int currentTile) {
            this.currentTile = currentTile;
        }
    }
}
