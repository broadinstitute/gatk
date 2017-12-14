package org.broadinstitute.hellbender.utils.nio;

import com.google.common.base.Stopwatch;
import java.nio.ByteBuffer;
import java.nio.channels.SeekableByteChannel;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadFactory;
import java.util.concurrent.TimeUnit;
import org.apache.logging.log4j.LogManager;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.annotations.Test;
import org.testng.Assert;
import org.apache.logging.log4j.Logger;

/**
 * Stress test for reading lots of data from the cloud using a very small prefetch buffer.
 * Do not run this too often.
 */
public final class ExtremeReadsTest extends GATKBaseTest {

    static final String fname = GCS_GATK_TEST_RESOURCES + "large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam";
    static final int THREAD_COUNT = 1000;
    static final int CHANNELS_PER_THREAD = 1000;
    private static Logger logger = LogManager.getLogger(ExtremeReadsTest.class);

    static volatile int errors = 0;

    private static class Runner implements Runnable {
        /**
         * Read a bunch of bytes. Part of the manyParallelReads test.
         */
        @Override
        public void run() {
            try {
                Path path = IOUtils.getPath(ExtremeReadsTest.fname);
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
                        logger.info("Done at wrong position! " + position + " != " + size);
                        ExtremeReadsTest.errors++;
                    }
                }
            } catch (Exception x) {
                ExtremeReadsTest.errors++;
                logger.info("Caught: " + x.getMessage());
                x.printStackTrace();
            }
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
        final ExecutorService executor = Executors.newFixedThreadPool(THREAD_COUNT,
            new ThreadFactory() {
                public Thread newThread(Runnable r) {
                    Thread t = Executors.defaultThreadFactory().newThread(r);
                    t.setDaemon(true);
                    return t;
                }
            });
        Stopwatch sw = Stopwatch.createStarted();
        errors = 0;
        for (int i=0; i<THREAD_COUNT; i++) {
            executor.execute(new Runner());
        }
        long parallel_reads = THREAD_COUNT * CHANNELS_PER_THREAD;
        logger.info(parallel_reads + " parallel reads via " + THREAD_COUNT + " threads (this will take a while).");
        executor.shutdown();
        executor.awaitTermination(1, TimeUnit.DAYS);
        sw.stop();
        logger.info("All done. Elapsed: " + sw.elapsed(TimeUnit.MINUTES) + " min.");
        logger.info("There were " + errors + " error(s).");
        Assert.assertEquals(errors, 0);
    }
}
