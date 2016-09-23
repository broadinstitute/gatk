package org.broadinstitute.hellbender.utils.nio;

import org.testng.annotations.Test;
import org.testng.Assert;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.channels.ReadableByteChannel;
import java.nio.channels.SeekableByteChannel;
import java.nio.file.Files;
import java.nio.file.Paths;

/**
 *
 */
public class SeekableByteChannelPrefetcherTest {
    // A file big enough to try seeks on.
    private final String input = "src/test/resources/exampleFASTA.fasta";

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
