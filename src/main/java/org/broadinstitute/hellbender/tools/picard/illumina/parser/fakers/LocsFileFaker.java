package org.broadinstitute.hellbender.tools.picard.illumina.parser.fakers;

import java.nio.ByteBuffer;

/**
 * Created by jcarey on 3/13/14.
 */
public final class LocsFileFaker extends FileFaker {

    @Override
    protected void fakeFile(final ByteBuffer buffer) {
        buffer.putInt(1);
        buffer.putFloat(1.0f);
        buffer.putInt(1);
        buffer.putFloat(5.0f);
        buffer.putFloat(5.0f);
    }

    @Override
    protected boolean addLeadingZeros() {
        return false;
    }

    @Override
    protected int bufferSize() {
        return (Integer.SIZE * 2) + (Float.SIZE * 3);
    }
}
