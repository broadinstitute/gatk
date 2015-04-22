package org.broadinstitute.hellbender.tools.picard.illumina.parser.fakers;

import java.nio.ByteBuffer;

/**
 * Created by jcarey on 3/13/14.
 */
public final class MultiTileLocsFileFaker extends FileFaker {

    @Override
    protected void fakeFile(final ByteBuffer buffer) {
        buffer.putInt(1);
        buffer.putFloat(1.0f);
        buffer.putInt(1);
        for (int count = 0; count < tiles.size(); count++) {
            buffer.putFloat(5.0f + (count * 0.5f));
            buffer.putFloat(5.0f + (count * 0.5f));
        }
    }

    @Override
    protected boolean addLeadingZeros() {
        return false;
    }

    @Override
    protected int bufferSize() {
        return (Integer.SIZE * 2) + (Float.SIZE * tiles.size()) + Float.SIZE;
    }
}
