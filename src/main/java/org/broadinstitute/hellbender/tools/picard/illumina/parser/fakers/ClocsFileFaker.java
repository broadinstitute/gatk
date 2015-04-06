package org.broadinstitute.hellbender.tools.picard.illumina.parser.fakers;

import java.nio.ByteBuffer;

public class ClocsFileFaker extends FileFaker {

    @Override
    protected void fakeFile(final ByteBuffer buffer) {
        buffer.put((byte) 1);
        buffer.putInt(1);
        buffer.put((byte) (0xff & 1));
        buffer.put((byte) (0xff & 5));
        buffer.put((byte) (0xff & 5));
    }

    @Override
    protected boolean addLeadingZeros() {
        return false;
    }

    @Override
    protected int bufferSize() {
        return Integer.SIZE + (Byte.SIZE * 4);
    }
}
