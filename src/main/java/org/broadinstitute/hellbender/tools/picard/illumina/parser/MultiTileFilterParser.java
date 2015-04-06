package org.broadinstitute.hellbender.tools.picard.illumina.parser;

import org.broadinstitute.hellbender.tools.picard.illumina.parser.readers.FilterFileReader;

import java.io.File;
import java.util.*;

import static java.util.Collections.singleton;
import static org.broadinstitute.hellbender.tools.picard.illumina.parser.IlluminaDataType.PF;

/**
 * Read filter file that contains multiple tiles in a single file.  A tile index is needed to parse this
 * file so that {tile number, cluster number} can be converted into absolute record number in file.
 */
public class MultiTileFilterParser extends MultiTileParser<PfData> {
    private final FilterFileReader reader;

    public MultiTileFilterParser(final TileIndex tileIndex, final List<Integer> requestedTiles, final File filterFile) {
        super(tileIndex, requestedTiles, singleton(PF));
        reader = new FilterFileReader(filterFile);
    }

    @Override
    PfData readNext() {
        final boolean nextVal = reader.next();
        return new PfData() {
            @Override
            public boolean isPf() {
                return nextVal;
            }
        };
    }

    @Override
    void skipRecords(final int numToSkip) {
        reader.skipRecords(numToSkip);
    }

    @Override
    public void close() {
        //no-op
    }
}
