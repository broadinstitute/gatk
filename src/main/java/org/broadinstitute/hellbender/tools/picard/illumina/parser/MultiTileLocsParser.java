package org.broadinstitute.hellbender.tools.picard.illumina.parser;

import org.broadinstitute.hellbender.tools.picard.illumina.parser.readers.AbstractIlluminaPositionFileReader;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.readers.LocsFileReader;

import java.io.File;
import java.util.*;

import static java.util.Collections.singleton;
import static org.broadinstitute.hellbender.tools.picard.illumina.parser.IlluminaDataType.Position;

/**
 * Read locs file that contains multiple tiles in a single file.  A tile index is needed to parse this
 * file so that {tile number, cluster number} can be converted into absolute record number in file.
 */
public class MultiTileLocsParser extends MultiTileParser<PositionalData> {
    private final LocsFileReader reader;
    private final int lane;

    public MultiTileLocsParser(final TileIndex tileIndex, final List<Integer> requestedTiles, final File locsFile, final int lane) {
        super(tileIndex, requestedTiles, singleton(Position));
        final int tileNumber;
        if (requestedTiles.size() == 1) tileNumber = requestedTiles.get(0);
        else tileNumber = -1;
        this.reader = new LocsFileReader(locsFile, lane, tileNumber);
        this.lane = lane;
    }

    @Override
    PositionalData readNext() {
        final int tile = getTileOfNextCluster();
        final AbstractIlluminaPositionFileReader.PositionInfo nextVal = reader.next();
        return new PositionalData() {
            public int getXCoordinate() {
                return nextVal.xQseqCoord;
            }

            public int getYCoordinate() {
                return nextVal.yQseqCoord;
            }
        };
    }

    @Override
    void skipRecords(final int numToSkip) {
        reader.skipRecords(numToSkip);
    }

    @Override
    public void close() {
        reader.close();
    }
}
