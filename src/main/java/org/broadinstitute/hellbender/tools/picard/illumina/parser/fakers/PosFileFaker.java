package org.broadinstitute.hellbender.tools.picard.illumina.parser.fakers;

import java.nio.ByteBuffer;

/**
 * Created by jcarey on 3/13/14.
 */
public class PosFileFaker extends FileFaker {
    private final String posFileString = "102.0\t303.3\n";

    @Override
    protected void fakeFile(final ByteBuffer buffer) {
        buffer.put(posFileString.getBytes());
    }

    @Override
    protected boolean addLeadingZeros() {
        return false;
    }

    @Override
    protected int bufferSize() {
        return posFileString.getBytes().length;
    }
}
