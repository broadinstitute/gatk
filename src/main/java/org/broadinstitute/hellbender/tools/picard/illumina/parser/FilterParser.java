package org.broadinstitute.hellbender.tools.picard.illumina.parser;

import htsjdk.samtools.util.CloseableIterator;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.readers.FilterFileReader;

import java.io.File;
import java.util.*;

import static htsjdk.samtools.util.CollectionUtil.makeSet;
import static java.util.Collections.unmodifiableSet;
import static org.broadinstitute.hellbender.tools.picard.illumina.parser.IlluminaDataType.PF;

/**
 * Sequentially parses filter files for the given tiles.  One tile is processed at a time.  IlluminaDataProvider should
 * be the ONLY client class for this class except for test classes.  For more information on the filterFile format
 * and reading it, see FilterFileReader.
 */
final class FilterParser extends PerTileParser<PfData> {
    private static Set<IlluminaDataType> supportedTypes = unmodifiableSet(makeSet(PF));

    public FilterParser(final IlluminaFileMap tilesToFiles) {
        super(tilesToFiles);
    }

    public FilterParser(final IlluminaFileMap tilesToFiles, final int startingTile) {
        super(tilesToFiles, startingTile);
    }

    /**
     * Wrap a filterFile reader in a closeable iterator and return it
     */
    @Override
    protected CloseableIterator<PfData> makeTileIterator(final File iterator) {
        return new CloseableIterator<PfData>() {
            private FilterFileReader reader = new FilterFileReader(iterator);

            public void close() {
                reader = null;
            }

            public boolean hasNext() {
                return reader.hasNext();
            }

            public PfData next() {
                final boolean nextValue = reader.next();
                return new PfData() {
                    public boolean isPf() {
                        return nextValue;
                    }
                };
            }

            public void remove() {
                throw new UnsupportedOperationException();
            }
        };
    }

    public Set<IlluminaDataType> supportedTypes() {
        return supportedTypes;
    }
}
