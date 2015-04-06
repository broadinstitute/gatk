package org.broadinstitute.hellbender.tools.picard.illumina.parser.readers;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.FileChannel;
import java.util.*;

import static htsjdk.samtools.util.CloserUtil.close;
import static htsjdk.samtools.util.IOUtil.assertFileIsReadable;
import static java.nio.ByteBuffer.allocate;
import static java.nio.ByteBuffer.wrap;
import static java.nio.ByteOrder.LITTLE_ENDIAN;
import static java.nio.channels.FileChannel.MapMode;
import static java.nio.channels.FileChannel.MapMode.READ_ONLY;

/**
 * MMapBackedIteratorFactory a file reader that takes a header size and a binary file, maps the file to
 * a read-only byte buffer and provides methods to retrieve the header as it's own bytebuffer and create
 * iterators of different data types over the values of file (starting after the end of the header).
 * Values provided by the MMappedBinaryFileReader are read as if they are little endian.
 * <p>
 * Note (read to end):
 * This class IS thread-safe and immutable though the iterator and ByteBuffers it produces are NOT.
 * The values read are assumed to be signed, NO promoting/sign conversion happens in this class.
 */
public class MMapBackedIteratorFactory {
    private static int BYTE_SIZE = 1;
    private static int INT_SIZE = 4;
    private static int FLOAT_SIZE = 4;
    private static int LONG_SIZE = 8;

    public static BinaryFileIterator<Integer> getIntegerIterator(final int headerSize, final File binaryFile) {
        checkFactoryVars(headerSize, binaryFile);
        final ByteBuffer buf = getBuffer(binaryFile);
        final byte[] header = getHeader(buf, headerSize);

        return new IntegerMMapIterator(header, binaryFile, buf);
    }

    public static BinaryFileIterator<Byte> getByteIterator(final int headerSize, final File binaryFile) {
        checkFactoryVars(headerSize, binaryFile);
        final ByteBuffer buf = getBuffer(binaryFile);
        final byte[] header = getHeader(buf, headerSize);

        return new ByteMMapIterator(header, binaryFile, buf);
    }

    public static BinaryFileIterator<Float> getFloatIterator(final int headerSize, final File binaryFile) {
        checkFactoryVars(headerSize, binaryFile);
        final ByteBuffer buf = getBuffer(binaryFile);
        final byte[] header = getHeader(buf, headerSize);

        return new FloatMMapIterator(header, binaryFile, buf);
    }

    public static BinaryFileIterator<Long> getLongIterator(final int headerSize, final File binaryFile) {
        checkFactoryVars(headerSize, binaryFile);
        final ByteBuffer buf = getBuffer(binaryFile);
        final byte[] header = getHeader(buf, headerSize);

        return new LongMMapIterator(header, binaryFile, buf);
    }

    public static BinaryFileIterator<ByteBuffer> getByteBufferIterator(final int headerSize, final int elementSize, final File binaryFile) {
        checkFactoryVars(headerSize, binaryFile);
        final ByteBuffer buf = getBuffer(binaryFile);
        final byte[] header = getHeader(buf, headerSize);

        return new ByteBufferMMapIterator(header, binaryFile, elementSize, buf);
    }

    private static void checkFactoryVars(final int headerSize, final File binaryFile) {
        assertFileIsReadable(binaryFile);

        if (headerSize < 0) {
            throw new IlluminaReaderException("Header size cannot be negative.  HeaderSize(" + headerSize + ") for file " + binaryFile.getAbsolutePath());
        }

        if (headerSize > binaryFile.length()) {
            throw new IlluminaReaderException("Header size(" + headerSize + ") is greater than file size(" + binaryFile.length() + ") for file " + binaryFile.getAbsolutePath());
        }
    }

    private static ByteBuffer getBuffer(final File binaryFile) {
        final ByteBuffer buf;
        try {
            final FileInputStream is = new FileInputStream(binaryFile);
            final FileChannel channel = is.getChannel();
            final long fileSize = channel.size();
            buf = channel.map(READ_ONLY, 0, fileSize);
            buf.order(LITTLE_ENDIAN);
            close(channel);
            close(is);
        } catch (IOException e) {
            throw new RuntimeIOException("IOException opening cluster binary file " + binaryFile, e);
        }

        return buf;
    }

    private static byte[] getHeader(final ByteBuffer buf, final int headerSize) {
        final byte[] headerBytes = new byte[headerSize];
        if (headerSize > 0) {
            buf.get(headerBytes);
        }
        return headerBytes;
    }

    /**
     * A simple iterator that uses a reference to the enclosing ByteBuffer and a member position
     * value to iterate over values in the buffer, starting after headerSize bytes
     */
    static abstract class MMapBackedIterator<TYPE> extends BinaryFileIterator<TYPE> {
        protected final ByteBuffer buffer;

        protected MMapBackedIterator(final byte[] header, final File file, final int elementSize, final ByteBuffer buffer) {
            super(header, file, elementSize);
            this.buffer = buffer;
        }

        public boolean hasNext() {
            return buffer.limit() - buffer.position() >= elementSize;
        }

        public void skipElements(final int numElements) {
            buffer.position(buffer.position() + (numElements * elementSize));
        }

        /**
         * The method that actually retrieves the data from the enclosing buffer
         */
        protected abstract TYPE getElement();

        public Iterator<TYPE> iterator() {
            return this;
        }
    }

    private static class IntegerMMapIterator extends MMapBackedIterator<Integer> {
        public IntegerMMapIterator(final byte[] header, final File file, final ByteBuffer buf) {
            super(header, file, INT_SIZE, buf);
        }

        @Override
        protected Integer getElement() {
            return buffer.getInt();
        }
    }

    private static class ByteMMapIterator extends MMapBackedIterator<Byte> {
        public ByteMMapIterator(final byte[] header, final File file, final ByteBuffer buf) {
            super(header, file, BYTE_SIZE, buf);
        }

        @Override
        protected Byte getElement() {
            return buffer.get();
        }
    }

    private static class FloatMMapIterator extends MMapBackedIterator<Float> {
        public FloatMMapIterator(final byte[] header, final File file, final ByteBuffer buf) {
            super(header, file, FLOAT_SIZE, buf);
        }

        @Override
        protected Float getElement() {
            return buffer.getFloat();
        }
    }

    private static class LongMMapIterator extends MMapBackedIterator<Long> {
        public LongMMapIterator(final byte[] header, final File file, final ByteBuffer buf) {
            super(header, file, LONG_SIZE, buf);
        }

        @Override
        protected Long getElement() {
            return buffer.getLong();
        }
    }

    //TODO: Add test
    //TODO: Make a note that if you want to multithread over this then you have to copy the contents
    private static class ByteBufferMMapIterator extends MMapBackedIterator<ByteBuffer> {
        private byte[] localBacking;
        private ByteBuffer localBuffer;

        public ByteBufferMMapIterator(final byte[] header, final File file, final int elementBufferSize, final ByteBuffer buf) {
            super(header, file, elementBufferSize, buf);
            this.localBacking = new byte[elementBufferSize];
            this.localBuffer = wrap(localBacking);
            this.localBuffer.order(LITTLE_ENDIAN);
        }

        @Override
        protected ByteBuffer getElement() {
            localBuffer.position(0);
            buffer.get(this.localBacking);
            localBuffer.position(0);
            return localBuffer;
        }
    }
}

abstract class BinaryFileIterator<TYPE> implements Iterator<TYPE>, Iterable<TYPE> {
    protected final File file;
    protected final long fileSize;
    protected final int elementSize;
    private final byte[] header;

    public BinaryFileIterator(final byte[] header, final File file, final int elementSize) {
        this.header = header;
        this.file = file;
        this.fileSize = file.length();
        this.elementSize = elementSize;
    }

    /**
     * Return the bytes found in the first headerSize bytes of the file, wrapped as a
     * ByteBuffer
     */
    public ByteBuffer getHeaderBytes() {
        final ByteBuffer bb = allocate(header.length);
        bb.order(LITTLE_ENDIAN);
        bb.put(header);
        bb.position(0);
        return bb;
    }

    public void assertTotalElementsEqual(final long numElements) {
        if (getElementsInFile() != numElements) {
            throw new IlluminaReaderException("Expected " + numElements + " elements in file but found " + getElementsInFile()
                    + " elements! File(" + file.getAbsolutePath() + ")");
        }

        if (getExtraBytes() != 0) {
            throw new IlluminaReaderException("Malformed file, expected " + (header.length + numElements * elementSize)
                    + " bytes in file, found " + fileSize + " bytes for file(" + file.getAbsolutePath() + ")");
        }
    }

    public int getElementSize() {
        return elementSize;
    }

    public long getExtraBytes() {
        return fileSize - header.length - (getElementsInFile() * elementSize);
    }

    public long getElementsInFile() {
        return (fileSize - header.length) / elementSize;
    }

    public File getFile() {
        return file;
    }

    public TYPE next() {
        if (!hasNext()) {
            throw new NoSuchElementException();
        }
        return getElement();
    }

    public void remove() {
        throw new UnsupportedOperationException();
    }

    public Iterator<TYPE> iterator() {
        return this;
    }

    /**
     * The method that actually retrieves the data from the enclosing buffer
     */
    protected abstract TYPE getElement();

    public abstract void skipElements(final int numElementsToSkip);

    public abstract boolean hasNext();
}



