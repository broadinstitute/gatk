package org.broadinstitute.hellbender.utils.io;

import org.apache.commons.io.output.ThresholdingOutputStream;

import java.io.IOException;

/**
 * An output stream which stops at the threshold
 * instead of potentially triggering early.
 */
public abstract class HardThresholdingOutputStream extends ThresholdingOutputStream {
    protected HardThresholdingOutputStream(int threshold) {
        super(threshold);
    }

    @Override
    public void write(byte[] b) throws IOException {
        write(b, 0, b.length);
    }

    @Override
    public void write(byte[] b, int off, int len) throws IOException {
        int remaining = this.getThreshold() - (int)this.getByteCount();
        if (!isThresholdExceeded() && len > remaining) {
            super.write(b, off, remaining);
            super.write(b, off + remaining, len - remaining);
        } else {
            super.write(b, off, len);
        }
    }
}
