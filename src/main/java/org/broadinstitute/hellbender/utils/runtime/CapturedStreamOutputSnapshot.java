package org.broadinstitute.hellbender.utils.runtime;

import org.broadinstitute.hellbender.utils.Utils;

import java.io.*;

/**
 * Stream output captured from a streaming stream.
 */
public final class CapturedStreamOutputSnapshot extends CapturedStreamOutput {

    /**
     * @param settings       Settings that define what to capture.
     * @param processStream  Stream to capture output.
     * @param standardStream Stream to write debug output.
     */
    public CapturedStreamOutputSnapshot(final OutputStreamSettings settings, final InputStream processStream, final PrintStream standardStream) {
        super(settings, processStream, standardStream);

        // strictly speaking, we could support these, but their use in streaming applications seems questionable
        Utils.validate(!settings.getStreamLocations().contains(StreamLocation.Standard), "Stream snapshots can't go to standard out");
        Utils.validate(!settings.getStreamLocations().contains(StreamLocation.File), "Stream snapshots can't go to a file");
    }

    /**
     * Drain the input stream to keep the process from backing up until it's empty.
     *
     * @throws IOException When unable to read or write.
     */
    @Override
    public void read() throws IOException {
        try {
            byte[] buf = new byte[CapturedStreamOutput.STREAM_BLOCK_TRANSFER_SIZE];
            // keep reading until there is no data more available
            do {
                final int readCount = processStream.read(buf);
                if (readCount == -1) {
                    return;
                }
                else if (readCount >= 0) {
                    bufferStream.write(buf, 0, readCount);
                }
            } while (processStream.available() > 0);
        } finally {
            bufferStream.flush();
        }
    }
}
