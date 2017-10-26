package org.broadinstitute.hellbender.utils.nio;

import com.google.cloud.storage.StorageException;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.channels.ReadableByteChannel;
import java.nio.channels.SeekableByteChannel;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import org.testng.Assert;
import org.testng.annotations.Test;

/**
 * Tests for SeekableByteChannelPrefetcher
 */
public class SeekableByteChannelPrefetcherTest {
    // A file big enough to try seeks on.
    private final String input = "src/test/resources/exampleFASTA.fasta";

    public class SeekableByteChannelAdapter implements SeekableByteChannel {
        protected final SeekableByteChannel inner;

        public SeekableByteChannelAdapter(SeekableByteChannel inner) {
            this.inner = inner;
        }

        /**
         * Reads a sequence of bytes from this channel into the given buffer.
         *
         * <p> Bytes are read starting at this channel's current position, and
         * then the position is updated with the number of bytes actually read.
         * Otherwise this method behaves exactly as specified in the {@link
         * ReadableByteChannel} interface.
         */
        @Override
        public int read(ByteBuffer dst) throws IOException {

            return inner.read(dst);
        }

        /**
         * Writes a sequence of bytes to this channel from the given buffer.
         *
         * <p> Bytes are written starting at this channel's current position, unless
         * the channel is connected to an entity such as a file that is opened with
         * the {@link StandardOpenOption#APPEND APPEND} option, in
         * which case the position is first advanced to the end. The entity to which
         * the channel is connected is grown, if necessary, to accommodate the
         * written bytes, and then the position is updated with the number of bytes
         * actually written. Otherwise this method behaves exactly as specified by
         * the {@link WritableByteChannel} interface.
         */
        @Override
        public int write(ByteBuffer src) throws IOException {
            return inner.write(src);
        }

        /**
         * Returns this channel's position.
         *
         * @return This channel's position, a non-negative integer counting the number of bytes from the
         * beginning of the entity to the current position
         * @throws ClosedChannelException If this channel is closed
         * @throws IOException If some other I/O error occurs
         */
        @Override
        public long position() throws IOException {
            return inner.position();
        }

        /**
         * Sets this channel's position.
         *
         * <p> Setting the position to a value that is greater than the current size
         * is legal but does not change the size of the entity.  A later attempt to
         * read bytes at such a position will immediately return an end-of-file
         * indication.  A later attempt to write bytes at such a position will cause
         * the entity to grow to accommodate the new bytes; the values of any bytes
         * between the previous end-of-file and the newly-written bytes are
         * unspecified.
         *
         * <p> Setting the channel's position is not recommended when connected to
         * an entity, typically a file, that is opened with the {@link
         * StandardOpenOption#APPEND APPEND} option. When opened for
         * append, the position is first advanced to the end before writing.
         *
         * @param newPosition The new position, a non-negative integer counting the number of bytes from
         * the beginning of the entity
         * @return This channel
         * @throws ClosedChannelException If this channel is closed
         * @throws IllegalArgumentException If the new position is negative
         * @throws IOException If some other I/O error occurs
         */
        @Override
        public SeekableByteChannel position(long newPosition) throws IOException {
            inner.position(newPosition);
            return this;
        }

        /**
         * Returns the current size of entity to which this channel is connected.
         *
         * @return The current size, measured in bytes
         * @throws ClosedChannelException If this channel is closed
         * @throws IOException If some other I/O error occurs
         */
        @Override
        public long size() throws IOException {
            return inner.size();
        }

        /**
         * Truncates the entity, to which this channel is connected, to the given
         * size.
         *
         * <p> If the given size is less than the current size then the entity is
         * truncated, discarding any bytes beyond the new end. If the given size is
         * greater than or equal to the current size then the entity is not modified.
         * In either case, if the current position is greater than the given size
         * then it is set to that size.
         *
         * <p> An implementation of this interface may prohibit truncation when
         * connected to an entity, typically a file, opened with the {@link
         * StandardOpenOption#APPEND APPEND} option.
         *
         * @param size The new size, a non-negative byte count
         * @return This channel
         * @throws NonWritableChannelException If this channel was not opened for writing
         * @throws ClosedChannelException If this channel is closed
         * @throws IllegalArgumentException If the new size is negative
         * @throws IOException If some other I/O error occurs
         */
        @Override
        public SeekableByteChannel truncate(long size) throws IOException {
            inner.truncate(size);
            return this;
        }

        /**
         * Tells whether or not this channel is open.
         *
         * @return <tt>true</tt> if, and only if, this channel is open
         */
        @Override
        public boolean isOpen() {
            return inner.isOpen();
        }

        /**
         * Closes this channel.
         *
         * <p> After a channel is closed, any further attempt to invoke I/O
         * operations upon it will cause a {@link ClosedChannelException} to be
         * thrown.
         *
         * <p> If this channel is already closed then invoking this method has no
         * effect.
         *
         * <p> This method may be invoked at any time.  If some other thread has
         * already invoked it, however, then another invocation will block until
         * the first invocation is complete, after which it will return without
         * effect. </p>
         *
         * @throws IOException If an I/O error occurs
         */
        @Override
        public void close() throws IOException {
            inner.close();
        }
    }

    public class RandomlyErroringSeekableByteChannel extends SeekableByteChannelAdapter {
        private int reads = 0;
        private final int pctError;
        public int errors = 0;

        public RandomlyErroringSeekableByteChannel(SeekableByteChannel inner, int pctError) {
            super(inner);
            this.pctError = pctError;
        }

        @Override
        public int read(ByteBuffer dst) throws IOException {
            if (reads++ % 100 < pctError) {
                errors++;
                throw new StorageException(403, "pretend-forbidden");
            }
            return super.read(dst);
        }
    }

    @Test
    public void testRead() throws Exception {
        SeekableByteChannel chan1 = Files.newByteChannel(Paths.get(input));
        SeekableByteChannel chan2 = new SeekableByteChannelPrefetcher(Files.newByteChannel(Paths.get(input)), 1024);

        testReading(chan1, chan2, 0);
        testReading(chan1, chan2, 128);
        testReading(chan1, chan2, 1024);
        testReading(chan1, chan2, 1500);
        testReading(chan1, chan2, 2048);
        testReading(chan1, chan2, 3000);
        testReading(chan1, chan2, 6000);
    }

    @Test
    public void testRetries() throws Exception {
        SeekableByteChannel chan1 = Files.newByteChannel(Paths.get(input));
        final RandomlyErroringSeekableByteChannel clumsy = new RandomlyErroringSeekableByteChannel(
            Files.newByteChannel(Paths.get(input)), 2);
        SeekableByteChannel chan2 = new SeekableByteChannelPrefetcher(
            clumsy,
            1024);

        testReading(chan1, chan2, 0);
        testReading(chan1, chan2, 128);
        Assert.assertTrue(clumsy.errors > 0);
    }

    @Test(expectedExceptions = RuntimeException.class)
    public void testRetriesEventuallyGiveUp() throws Exception {
        SeekableByteChannel chan1 = Files.newByteChannel(Paths.get(input));
        final RandomlyErroringSeekableByteChannel clumsy = new RandomlyErroringSeekableByteChannel(
            Files.newByteChannel(Paths.get(input)), 100);
        SeekableByteChannel chan2 = new SeekableByteChannelPrefetcher(
            clumsy,
            1024);

        testReading(chan1, chan2, 128);
        Assert.assertTrue(clumsy.errors > 0);
    }

    @Test
    public void testSeek() throws Exception {
        SeekableByteChannel chan1 = Files.newByteChannel(Paths.get(input));
        SeekableByteChannel chan2 = new SeekableByteChannelPrefetcher(Files.newByteChannel(Paths.get(input)), 1024);

        testSeeking(chan1, chan2, 1024);
        testSeeking(chan1, chan2, 1500);
        testSeeking(chan1, chan2, 128);
        testSeeking(chan1, chan2, 256);
        testSeeking(chan1, chan2, 128);
        // yes, testReading - let's make sure that reading more than one block still works
        // even after a seek.
        testReading(chan1, chan2, 1500);
        testSeeking(chan1, chan2, 2048);
        testSeeking(chan1, chan2, 0);
        testSeeking(chan1, chan2, 3000);
        testSeeking(chan1, chan2, 6000);
        testSeeking(chan1, chan2, (int)chan1.size()-127);
        testSeeking(chan1, chan2, (int)chan1.size()-128);
        testSeeking(chan1, chan2, (int)chan1.size()-129);
    }

    @Test
    public void testPartialBuffers() throws Exception {
        SeekableByteChannel chan1 = Files.newByteChannel(Paths.get(input));
        SeekableByteChannel chan2 = new SeekableByteChannelPrefetcher(
          Files.newByteChannel(Paths.get(input)), 1024);
        // get a partial buffer
        testSeeking(chan1, chan2, (int) chan1.size() - 127);
        // make sure normal reads can use the full buffer
        for (int i = 0; i < 2; i++) {
            testSeeking(chan1, chan2, i * 1024);
        }
        // get a partial buffer, replacing one of the full ones
        testSeeking(chan1, chan2, (int) chan1.size() - 127);
        // make sure the buffers are still OK
        for (int i = 0; i < 2; i++) {
            testSeeking(chan1, chan2, i * 1024);
        }
  }

    @Test
    public void testEOF() throws Exception {
        SeekableByteChannel chan1 = Files.newByteChannel(Paths.get(input));
        SeekableByteChannel chan2 = new SeekableByteChannelPrefetcher(
            Files.newByteChannel(Paths.get(input)), 1024);
        // read the final 128 bytes, exactly.
        testSeeking(chan1, chan2, (int) chan1.size() - 128);
        // read truncated because we're asking for beyond EOF
        testSeeking(chan1, chan2, (int) chan1.size() - 64);
        // read starting past EOF
        testSeeking(chan1, chan2, (int) chan1.size() + 128);
        // read more than a whole block past EOF
        testSeeking(chan1, chan2, (int) chan1.size() + 1024 * 2);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testDoubleWrapping() throws Exception {
        SeekableByteChannel chan1 = new SeekableByteChannelPrefetcher(
            Files.newByteChannel(Paths.get(input)), 1024);
        new SeekableByteChannelPrefetcher(chan1, 1024);
    }

    @Test
    public void testCloseWhilePrefetching() throws Exception {
        SeekableByteChannel chan = new SeekableByteChannelPrefetcher(
            Files.newByteChannel(Paths.get(input)), 10*1024*1024);
        // read just 1 byte, get the prefetching going
        ByteBuffer one = ByteBuffer.allocate(1);
        readFully(chan, one);
        // closing must not throw an exception, even if the prefetching
        // thread is active.
        chan.close();
    }

    private void testReading(SeekableByteChannel chan1, SeekableByteChannel chan2, int howMuch) throws IOException {
        ByteBuffer one = ByteBuffer.allocate(howMuch);
        ByteBuffer two = ByteBuffer.allocate(howMuch);

        readFully(chan1, one);
        readFully(chan2, two);

        Assert.assertEquals(one.position(), two.position());
        Assert.assertEquals(one.array(), two.array());
    }

    private void testSeeking(SeekableByteChannel chan1, SeekableByteChannel chan2, int position) throws IOException {
        ByteBuffer one = ByteBuffer.allocate(128);
        ByteBuffer two = ByteBuffer.allocate(128);

        chan1.position(position);
        chan2.position(position);

        readFully(chan1, one);
        readFully(chan2, two);

        Assert.assertEquals(one.position(), two.position());
        Assert.assertEquals(one.array(), two.array());
    }

    private void readFully(ReadableByteChannel chan, ByteBuffer buf) throws IOException {
        // the countdown isn't strictly necessary but it protects us against infinite loops
        // for some potential bugs in the channel implementation.
        int countdown = buf.capacity();
        while (chan.read(buf) > 0  && countdown-- > 0) {}
    }

}
