package org.broadinstitute.hellbender.tools.picard.illumina.parser.fakers;

import java.nio.ByteBuffer;

/**
 * Created by jcarey on 3/13/14.
 */
public final class MultiTileBclFileFaker extends FileFaker {

    @Override
    protected void fakeFile(final ByteBuffer buffer) {
        buffer.putInt(1);
        for (final Integer tile : tiles) {
            long perTileSize = size;
            while (perTileSize > 0) {
                //fill the file with no calls
                buffer.put((byte) 0);
                perTileSize--;
            }
        }
    }

    @Override
    protected boolean addLeadingZeros() {
        return false;
    }

    @Override
    protected int bufferSize() {
        return ((size - Integer.SIZE) * tiles.size()) + Integer.SIZE;
    }
}
