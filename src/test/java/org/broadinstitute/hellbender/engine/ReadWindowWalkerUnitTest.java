package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;

public class ReadWindowWalkerUnitTest extends BaseTest {

    @DataProvider(name = "ShardIntervalTestData")
    public Object[][] shardIntervalTestData() {
        // Doesn't matter which bam we use to back the reads source for the purposes of these tests.
        final ReadsDataSource readsSource = new ReadsDataSource(new File(publicTestDir + "org/broadinstitute/hellbender/engine/reads_data_source_test1.bam"));
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary(Arrays.asList(new SAMSequenceRecord("1", 16000)));

        return new Object[][] {
                {
                        new SimpleInterval("1", 1, 100), 100, 0, readsSource, dictionary,
                        Arrays.asList(new ReadWindow(new SimpleInterval("1", 1, 100), new SimpleInterval("1", 1, 100), readsSource))
                },

                {
                        new SimpleInterval("1", 1, 100), 100, 10, readsSource, dictionary,
                        Arrays.asList(new ReadWindow(new SimpleInterval("1", 1, 100), new SimpleInterval("1", 1, 110), readsSource))
                },

                {
                        new SimpleInterval("1", 1, 100), 50, 0, readsSource, dictionary,
                        Arrays.asList(new ReadWindow(new SimpleInterval("1", 1, 50), new SimpleInterval("1", 1, 50), readsSource),
                                      new ReadWindow(new SimpleInterval("1", 51, 100), new SimpleInterval("1", 51, 100), readsSource))
                },

                {
                        new SimpleInterval("1", 1, 100), 50, 10, readsSource, dictionary,
                        Arrays.asList(new ReadWindow(new SimpleInterval("1", 1, 50), new SimpleInterval("1", 1, 60), readsSource),
                                      new ReadWindow(new SimpleInterval("1", 51, 100), new SimpleInterval("1", 41, 110), readsSource))
                },

                {
                        new SimpleInterval("1", 51, 200), 50, 0, readsSource, dictionary,
                        Arrays.asList(new ReadWindow(new SimpleInterval("1", 51, 100), new SimpleInterval("1", 51, 100), readsSource),
                                      new ReadWindow(new SimpleInterval("1", 101, 150), new SimpleInterval("1", 101, 150), readsSource),
                                      new ReadWindow(new SimpleInterval("1", 151, 200), new SimpleInterval("1", 151, 200), readsSource))
                },

                {
                        new SimpleInterval("1", 51, 200), 50, 10, readsSource, dictionary,
                        Arrays.asList(new ReadWindow(new SimpleInterval("1", 51, 100), new SimpleInterval("1", 41, 110), readsSource),
                                      new ReadWindow(new SimpleInterval("1", 101, 150), new SimpleInterval("1", 91, 160), readsSource),
                                      new ReadWindow(new SimpleInterval("1", 151, 200), new SimpleInterval("1", 141, 210), readsSource))
                },

                {
                        new SimpleInterval("1", 51, 200), 70, 0, readsSource, dictionary,
                        Arrays.asList(new ReadWindow(new SimpleInterval("1", 51, 120), new SimpleInterval("1", 51, 120), readsSource),
                                      new ReadWindow(new SimpleInterval("1", 121, 190), new SimpleInterval("1", 121, 190), readsSource),
                                      new ReadWindow(new SimpleInterval("1", 191, 200), new SimpleInterval("1", 191, 200), readsSource))
                },

                {
                        new SimpleInterval("1", 51, 200), 70, 10, readsSource, dictionary,
                        Arrays.asList(new ReadWindow(new SimpleInterval("1", 51, 120), new SimpleInterval("1", 41, 130), readsSource),
                                      new ReadWindow(new SimpleInterval("1", 121, 190), new SimpleInterval("1", 111, 200), readsSource),
                                      new ReadWindow(new SimpleInterval("1", 191, 200), new SimpleInterval("1", 181, 210), readsSource))
                },

                {
                        new SimpleInterval("1", 15999, 16000), 100, 10, readsSource, dictionary,
                        Arrays.asList(new ReadWindow(new SimpleInterval("1", 15999, 16000), new SimpleInterval("1", 15989, 16000), readsSource))
                },
        };
    }

    @Test(dataProvider = "ShardIntervalTestData")
    public void testShardInterval( final SimpleInterval originalInterval, final int windowSize, final int windowPadding, final ReadsDataSource readsSource, final SAMSequenceDictionary dictionary, List<ReadWindow> expectedWindows  ) {
        final List<ReadWindow> windows = ReadWindowWalker.shardInterval(originalInterval, windowSize, windowPadding, readsSource, dictionary);

        Assert.assertEquals(windows.size(), expectedWindows.size(), "Wrong number of windows created by ReadWindowWalker.shardInterval()");
        for ( int i = 0; i < windows.size(); ++i ) {
            final ReadWindow window = windows.get(i);
            Assert.assertEquals(window.getInterval(), expectedWindows.get(i).getInterval(), "Window has wrong size");
            Assert.assertEquals(window.getPaddedInterval(), expectedWindows.get(i).getPaddedInterval(), "Window has wrong amount of padding");
        }
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testShardIntervalThrowsOnInvalidContig() {
        final ReadsDataSource readsSource = new ReadsDataSource(new File(publicTestDir + "org/broadinstitute/hellbender/engine/reads_data_source_test1.bam"));
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary(Arrays.asList(new SAMSequenceRecord("1", 16000)));

        // contig not in dictionary
        ReadWindowWalker.shardInterval(new SimpleInterval("2", 1, 10), 100, 0, readsSource, dictionary);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testShardIntervalThrowsOnOutOfBoundsInterval() {
        final ReadsDataSource readsSource = new ReadsDataSource(new File(publicTestDir + "org/broadinstitute/hellbender/engine/reads_data_source_test1.bam"));
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary(Arrays.asList(new SAMSequenceRecord("1", 16000)));

        // interval not within contig bounds
        ReadWindowWalker.shardInterval(new SimpleInterval("2", 1, 17000), 100, 0, readsSource, dictionary);
    }
}
