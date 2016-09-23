package org.broadinstitute.hellbender.utils.nio;

import htsjdk.samtools.seekablestream.SeekableStream;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.channels.SeekableByteChannel;

/**
 * A SeekableByteChannel in SeekableStream clothes.
 * This is key to using NIO Paths with htsjdk (enabling it to work with gcloud-java-nio).
 */
public class ChannelAsSeekableStream extends SeekableStream {

    private final SeekableByteChannel chan;
    private final String source;
    private ByteBuffer oneByte;

    public ChannelAsSeekableStream(SeekableByteChannel chan) {
        this.chan = chan;
        this.source = null;
    }

    public ChannelAsSeekableStream(SeekableByteChannel chan, String source) {
        this.chan = chan;
        this.source = source;
    }

    @Override
    public long length() {
        try {
            return chan.size();
        } catch (IOException x) {
            throw new RuntimeException(x);
        }
    }

    @Override
    public long position() throws IOException {
        return chan.position();
    }

    @Override
    public void seek(long position) throws IOException {
        chan.position(position);
    }

    /**
     * Reads the next byte of data from the input stream. The value byte is
     * returned as an <code>int</code> in the range <code>0</code> to
     * <code>255</code>. If no byte is available because the end of the stream
     * has been reached, the value <code>-1</code> is returned. This method
     * blocks until input data is available, the end of the stream is detected,
     * or an exception is thrown.
     *
     * @return the next byte of data, or <code>-1</code> if the end of the
     * stream is reached.
     * @throws IOException if an I/O error occurs.
     */
    @Override
    public int read() throws IOException {
        ByteBuffer buf = oneByte;
        if (null == buf) {
            buf = ByteBuffer.allocate(1);
            oneByte = buf;
        }
        chan.read(buf);
        if (buf.remaining() == 0) {
            return -1;
        }
        return buf.get(0);
    }

    @Override
    public int read(byte[] buffer, int offset, int length) throws IOException {
        ByteBuffer buf = ByteBuffer.wrap(buffer, offset, length);
        return chan.read(buf);
    }

    @Override
    public void close() throws IOException {
        chan.close();
    }

    @Override
    public boolean eof() throws IOException {
        return chan.position() >= chan.size();
    }

    /**
     * @return String representation of source (e.g. URL, file path, etc.), or null if not available.
     */
    @Override
    public String getSource() {
        return source;
    }
}
