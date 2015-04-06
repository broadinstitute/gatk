package org.broadinstitute.hellbender.tools.picard.illumina.parser.fakers;

import java.nio.ByteBuffer;

/**
 * Created by jcarey on 3/13/14.
 */
public class BarcodeFileFaker extends FileFaker {
    private final String barcodeString = "1\tn\t \n";

    @Override
    protected void fakeFile(final ByteBuffer buffer) {
        buffer.put(barcodeString.getBytes());
    }

    @Override
    protected boolean addLeadingZeros() {
        return false;
    }

    @Override
    protected int bufferSize() {
        return barcodeString.getBytes().length;
    }
}
