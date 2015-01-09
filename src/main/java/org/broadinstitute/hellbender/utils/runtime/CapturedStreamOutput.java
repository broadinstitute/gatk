/*
* Copyright (c) 2012 The Broad Institute
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.hellbender.utils.runtime;

import org.apache.commons.io.IOUtils;
import org.apache.commons.io.output.NullOutputStream;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.io.HardThresholdingOutputStream;

import java.io.*;
import java.util.EnumMap;

/**
 * Stream output captured from a stream.
 */
public class CapturedStreamOutput extends StreamOutput {
    private final InputStream processStream;
    private final EnumMap<StreamLocation, OutputStream> outputStreams = new EnumMap<>(StreamLocation.class);

    /**
     * The byte stream to capture content or null if no output string content was requested.
     */
    private final ByteArrayOutputStream bufferStream;

    /**
     * True if the buffer is truncated.
     */
    private boolean bufferTruncated = false;

    /**
     * @param settings       Settings that define what to capture.
     * @param processStream  Stream to capture output.
     * @param standardStream Stream to write debug output.
     */
    public CapturedStreamOutput(OutputStreamSettings settings, InputStream processStream, PrintStream standardStream) {
        this.processStream = processStream;
        int bufferSize = settings.getBufferSize();
        this.bufferStream = (bufferSize < 0) ? new ByteArrayOutputStream() : new ByteArrayOutputStream(bufferSize);

        for (StreamLocation location : settings.getStreamLocations()) {
            OutputStream outputStream;
            switch (location) {
                case Buffer:
                    if (bufferSize < 0) {
                        outputStream = this.bufferStream;
                    } else {
                        outputStream = new HardThresholdingOutputStream(bufferSize) {
                            @Override
                            protected OutputStream getStream() throws IOException {
                                return bufferTruncated ? NullOutputStream.NULL_OUTPUT_STREAM : bufferStream;
                            }

                            @Override
                            protected void thresholdReached() throws IOException {
                                bufferTruncated = true;
                            }
                        };
                    }
                    break;
                case File:
                    try {
                        outputStream = new FileOutputStream(settings.getOutputFile(), settings.isAppendFile());
                    } catch (IOException e) {
                        throw new UserException.BadInput(e.getMessage());
                    }
                    break;
                case Standard:
                    outputStream = standardStream;
                    break;
                default:
                    throw new GATKException("Unexpected stream location: " + location);
            }
            this.outputStreams.put(location, outputStream);
        }
    }

    @Override
    public byte[] getBufferBytes() {
        return bufferStream.toByteArray();
    }

    @Override
    public boolean isBufferTruncated() {
        return bufferTruncated;
    }

    /**
     * Drain the input stream to keep the process from backing up until it's empty.
     * File streams will be closed automatically when this method returns.
     *
     * @throws java.io.IOException When unable to read or write.
     */
    public void readAndClose() throws IOException {
        try {
            byte[] buf = new byte[4096];
            int readCount;
            while ((readCount = processStream.read(buf)) >= 0)
                for (OutputStream outputStream : this.outputStreams.values()) {
                    outputStream.write(buf, 0, readCount);
                }
        } finally {
            for (StreamLocation location : this.outputStreams.keySet()) {
                OutputStream outputStream = this.outputStreams.get(location);
                outputStream.flush();
                if (location != StreamLocation.Standard)
                    IOUtils.closeQuietly(outputStream);
            }
        }
    }
}
