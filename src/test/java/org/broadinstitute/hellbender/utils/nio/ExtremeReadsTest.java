package org.broadinstitute.hellbender.utils.nio;

import com.google.common.base.Stopwatch;
import java.nio.ByteBuffer;
import java.nio.channels.SeekableByteChannel;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.TimeUnit;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.Test;
import org.testng.Assert;

/**
 * Stress test for reading lots of data from the cloud using a very small prefetch buffer.
 * Do not run this too often.
 */
public final class ExtremeReadsTest extends BaseTest implements Runnable {

    String fname = GCS_GATK_TEST_RESOURCES + "large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam";

    static int THREAD_COUNT = 1000;
    static int CHANNELS_PER_THREAD = 1000;
    static volatile int errors = 0;

    /**
     * Read a bunch of bytes. Part of the manyParallelReads test.
     */
    @Override
    public void run() {
        try {
            Path path = IOUtils.getPath(fname);
            ArrayList<SeekableByteChannel> chans = new ArrayList<SeekableByteChannel>();
            for (int i=0; i<CHANNELS_PER_THREAD; i++) {
                SeekableByteChannelPrefetcher chan = new SeekableByteChannelPrefetcher(
                    Files.newByteChannel(path), 2 * 1024 * 1024);
                // skip the first half
                chan.position(chan.position()/2);
                chans.add(chan);
            }
            long size = chans.get(0).size();
            ByteBuffer buf = ByteBuffer.allocate(1024*1024 - 5);
            while (!chans.isEmpty()) {
                SeekableByteChannel chan = chans.remove(0);
                buf.clear();
                int read = chan.read(buf);
                if (read>=0) {
                    chans.add(chan);
                    continue;
                }
                // EOF
                long position = chan.position();
                if (size != position) {
                    System.out.println("Done at wrong position! " + position + " != " + size);
                    errors++;
                }
            }
        } catch (Exception x) {
            errors++;
            System.out.println("Caught: " + x.getMessage());
            x.printStackTrace();
            System.out.println();
        }
    }

    /**
     * This test takes about a half hour and reads a fair amount of data.
     * It definitely shouldn't be part of the normal test suite (that's why it's disabled)
     * but it's kept here so we can manually run it should we need to investigate mysterious
     * disconnects again.
     **/
    @Test(groups={"bucket"}, enabled=false)
    public void manyParallelReads() throws InterruptedException {
        List<Thread> threads = new ArrayList<>();
        Stopwatch sw = Stopwatch.createStarted();
        errors = 0;
        int count = THREAD_COUNT;
        for (int i=0; i<count; i++) {
            Thread t = new Thread(this);
            threads.add(t);
            t.start();
        }
        long parallel_reads = THREAD_COUNT * CHANNELS_PER_THREAD;
        System.out.println(parallel_reads + " parallel reads via " + THREAD_COUNT + " threads (this will take a while).");
        for (Thread t : threads) {
            t.join();
        }
        sw.stop();
        System.out.println("All done. Elapsed: " + sw.elapsed(TimeUnit.MINUTES) + " min.");
        System.out.println("There were " + errors + " error(s).");
        Assert.assertEquals(errors, 0);
    }
}
