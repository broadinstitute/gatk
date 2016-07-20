package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.downsampling.ReadsDownsampler;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class ReadShardUnitTest extends BaseTest {

    @DataProvider(name = "ShardIterationTestData")
    public Object[][] shardIterationTestData() {
        final ReadsDataSource readsSource = new ReadsDataSource(new File(publicTestDir + "org/broadinstitute/hellbender/engine/reads_data_source_test1.bam"));

        final ReadFilter keepReadBOnly = read -> read.getName().equals("b");
        final ReadShard filteredShard = new ReadShard(new SimpleInterval("1", 200, 210), new SimpleInterval("1", 200, 210), readsSource);
        filteredShard.setReadFilter(keepReadBOnly);

        final ReadsDownsampler readsBAndCOnlyDownsampler = new KeepReadsBAndCOnlyDownsampler();
        final ReadShard downsampledShard = new ReadShard(new SimpleInterval("1", 1, 5000), new SimpleInterval("1", 1, 5000), readsSource);
        downsampledShard.setDownsampler(readsBAndCOnlyDownsampler);

        return new Object[][] {
                { new ReadShard(new SimpleInterval("1", 200, 210), new SimpleInterval("1", 200, 210), readsSource), Arrays.asList("a", "b", "c") },
                { new ReadShard(new SimpleInterval("1", 200, 209), new SimpleInterval("1", 200, 209), readsSource), Arrays.asList("a", "b") },
                { new ReadShard(new SimpleInterval("1", 200, 204), new SimpleInterval("1", 200, 204), readsSource), Arrays.asList("a") },
                { new ReadShard(new SimpleInterval("1", 200, 204), new SimpleInterval("1", 190, 210), readsSource), Arrays.asList("a", "b", "c") },
                { new ReadShard(new SimpleInterval("1", 200, 204), new SimpleInterval("1", 200, 205), readsSource), Arrays.asList("a", "b") },
                { new ReadShard(new SimpleInterval("1", 400, 500), new SimpleInterval("1", 400, 500), readsSource), Collections.<String>emptyList() },
                { filteredShard, Arrays.asList("b") },
                { downsampledShard, Arrays.asList("b", "c")}
        };
    }

    @Test(dataProvider = "ShardIterationTestData")
    public void testShardIteration( final ReadShard shard, final List<String> expectedReadNames ) {
        final List<String> actualReadNames = new ArrayList<>();

        for ( final GATKRead read : shard ) {
            actualReadNames.add(read.getName());
        }

        Assert.assertEquals(actualReadNames.size(), expectedReadNames.size(), "Wrong number of reads returned");
        Assert.assertEquals(actualReadNames, expectedReadNames, "Wrong reads returned");
    }

    @Test(dataProvider = "ShardIterationTestData")
    public void testShardLoadReads( final ReadShard shard, final List<String> expectedReadNames ) {
        final List<String> actualReadNames = new ArrayList<>();

        for ( final GATKRead read : shard.loadAllReads() ) {
            actualReadNames.add(read.getName());
        }

        Assert.assertEquals(actualReadNames.size(), expectedReadNames.size(), "Wrong number of reads returned");
        Assert.assertEquals(actualReadNames, expectedReadNames, "Wrong reads returned");
    }

    @DataProvider(name = "DivideIntervalIntoShardsTestData")
    public Object[][] divideIntervalIntoShardsTestData() {
        // Doesn't matter which bam we use to back the reads source for the purposes of these tests.
        final ReadsDataSource readsSource = new ReadsDataSource(new File(publicTestDir + "org/broadinstitute/hellbender/engine/reads_data_source_test1.bam"));
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary(Arrays.asList(new SAMSequenceRecord("1", 16000)));

        return new Object[][] {
                // shardSize == 100, shardStep == 100, shardPadding == 0, start of contig
                {
                        new SimpleInterval("1", 1, 100), 100, 100, 0, readsSource, dictionary,
                        Arrays.asList(new ReadShard(new SimpleInterval("1", 1, 100), new SimpleInterval("1", 1, 100), readsSource))
                },

                // shardSize == 100, shardStep == 100, shardPadding == 10, start of contig
                {
                        new SimpleInterval("1", 1, 100), 100, 100, 10, readsSource, dictionary,
                        Arrays.asList(new ReadShard(new SimpleInterval("1", 1, 100), new SimpleInterval("1", 1, 110), readsSource))
                },

                // shardSize == 50, shardStep == 50, shardPadding == 0, start of contig
                {
                        new SimpleInterval("1", 1, 100), 50, 50, 0, readsSource, dictionary,
                        Arrays.asList(new ReadShard(new SimpleInterval("1", 1, 50), new SimpleInterval("1", 1, 50), readsSource),
                                      new ReadShard(new SimpleInterval("1", 51, 100), new SimpleInterval("1", 51, 100), readsSource))
                },

                // shardSize == 50, shardStep == 50, shardPadding == 10, start of contig
                {
                        new SimpleInterval("1", 1, 100), 50, 50, 10, readsSource, dictionary,
                        Arrays.asList(new ReadShard(new SimpleInterval("1", 1, 50), new SimpleInterval("1", 1, 60), readsSource),
                                      new ReadShard(new SimpleInterval("1", 51, 100), new SimpleInterval("1", 41, 110), readsSource))
                },

                // shardSize == 50, shardStep == 50, shardPadding == 0
                {
                        new SimpleInterval("1", 51, 200), 50, 50, 0, readsSource, dictionary,
                        Arrays.asList(new ReadShard(new SimpleInterval("1", 51, 100), new SimpleInterval("1", 51, 100), readsSource),
                                      new ReadShard(new SimpleInterval("1", 101, 150), new SimpleInterval("1", 101, 150), readsSource),
                                      new ReadShard(new SimpleInterval("1", 151, 200), new SimpleInterval("1", 151, 200), readsSource))
                },

                // shardSize == 50, shardStep == 50, shardPadding == 10
                {
                        new SimpleInterval("1", 51, 200), 50, 50, 10, readsSource, dictionary,
                        Arrays.asList(new ReadShard(new SimpleInterval("1", 51, 100), new SimpleInterval("1", 41, 110), readsSource),
                                      new ReadShard(new SimpleInterval("1", 101, 150), new SimpleInterval("1", 91, 160), readsSource),
                                      new ReadShard(new SimpleInterval("1", 151, 200), new SimpleInterval("1", 141, 210), readsSource))
                },

                // shardSize == 70, shardStep == 70, shardPadding == 0
                {
                        new SimpleInterval("1", 51, 200), 70, 70, 0, readsSource, dictionary,
                        Arrays.asList(new ReadShard(new SimpleInterval("1", 51, 120), new SimpleInterval("1", 51, 120), readsSource),
                                      new ReadShard(new SimpleInterval("1", 121, 190), new SimpleInterval("1", 121, 190), readsSource),
                                      new ReadShard(new SimpleInterval("1", 191, 200), new SimpleInterval("1", 191, 200), readsSource))
                },

                // shardSize == 70, shardStep == 70, shardPadding == 10
                {
                        new SimpleInterval("1", 51, 200), 70, 70, 10, readsSource, dictionary,
                        Arrays.asList(new ReadShard(new SimpleInterval("1", 51, 120), new SimpleInterval("1", 41, 130), readsSource),
                                      new ReadShard(new SimpleInterval("1", 121, 190), new SimpleInterval("1", 111, 200), readsSource),
                                      new ReadShard(new SimpleInterval("1", 191, 200), new SimpleInterval("1", 181, 210), readsSource))
                },

                // shardSize == 100, shardStep == 100, shardPadding == 10, end of contig
                {
                        new SimpleInterval("1", 15999, 16000), 100, 100, 10, readsSource, dictionary,
                        Arrays.asList(new ReadShard(new SimpleInterval("1", 15999, 16000), new SimpleInterval("1", 15989, 16000), readsSource))
                },

                // shardSize == 100, shardStep == 100, shardPadding == 0
                {
                        new SimpleInterval("1", 51, 300), 100, 100, 0, readsSource, dictionary,
                        Arrays.asList(new ReadShard(new SimpleInterval("1", 51, 150), new SimpleInterval("1", 51, 150), readsSource),
                                new ReadShard(new SimpleInterval("1", 151, 250), new SimpleInterval("1", 151, 250), readsSource),
                                new ReadShard(new SimpleInterval("1", 251, 300), new SimpleInterval("1", 251, 300), readsSource))
                },

                // shardSize == 100, shardStep == 50, shardPadding == 0
                {
                        new SimpleInterval("1", 51, 300), 100, 50, 0, readsSource, dictionary,
                        Arrays.asList(new ReadShard(new SimpleInterval("1", 51, 150), new SimpleInterval("1", 51, 150), readsSource),
                                new ReadShard(new SimpleInterval("1", 101, 200), new SimpleInterval("1", 101, 200), readsSource),
                                new ReadShard(new SimpleInterval("1", 151, 250), new SimpleInterval("1", 151, 250), readsSource),
                                new ReadShard(new SimpleInterval("1", 201, 300), new SimpleInterval("1", 201, 300), readsSource),
                                new ReadShard(new SimpleInterval("1", 251, 300), new SimpleInterval("1", 251, 300), readsSource))
                },

                // shardSize == 100, shardStep == 50, shardPadding == 10
                {
                        new SimpleInterval("1", 51, 300), 100, 50, 10, readsSource, dictionary,
                        Arrays.asList(new ReadShard(new SimpleInterval("1", 51, 150), new SimpleInterval("1", 41, 160), readsSource),
                                new ReadShard(new SimpleInterval("1", 101, 200), new SimpleInterval("1", 91, 210), readsSource),
                                new ReadShard(new SimpleInterval("1", 151, 250), new SimpleInterval("1", 141, 260), readsSource),
                                new ReadShard(new SimpleInterval("1", 201, 300), new SimpleInterval("1", 191, 310), readsSource),
                                new ReadShard(new SimpleInterval("1", 251, 300), new SimpleInterval("1", 241, 310), readsSource))
                },

                // shardSize == 100, shardStep == 25 shardPadding == 0
                {
                        new SimpleInterval("1", 51, 300), 100, 25, 0, readsSource, dictionary,
                        Arrays.asList(new ReadShard(new SimpleInterval("1", 51, 150), new SimpleInterval("1", 51, 150), readsSource),
                                new ReadShard(new SimpleInterval("1", 76, 175), new SimpleInterval("1", 76, 175), readsSource),
                                new ReadShard(new SimpleInterval("1", 101, 200), new SimpleInterval("1", 101, 200), readsSource),
                                new ReadShard(new SimpleInterval("1", 126, 225), new SimpleInterval("1", 126, 225), readsSource),
                                new ReadShard(new SimpleInterval("1", 151, 250), new SimpleInterval("1", 151, 250), readsSource),
                                new ReadShard(new SimpleInterval("1", 176, 275), new SimpleInterval("1", 176, 275), readsSource),
                                new ReadShard(new SimpleInterval("1", 201, 300), new SimpleInterval("1", 201, 300), readsSource),
                                new ReadShard(new SimpleInterval("1", 226, 300), new SimpleInterval("1", 226, 300), readsSource),
                                new ReadShard(new SimpleInterval("1", 251, 300), new SimpleInterval("1", 251, 300), readsSource),
                                new ReadShard(new SimpleInterval("1", 276, 300), new SimpleInterval("1", 276, 300), readsSource))
                },

                // shardSize == 100, shardStep == 25, shardPadding == 10
                {
                        new SimpleInterval("1", 51, 300), 100, 25, 10, readsSource, dictionary,
                        Arrays.asList(new ReadShard(new SimpleInterval("1", 51, 150), new SimpleInterval("1", 41, 160), readsSource),
                                new ReadShard(new SimpleInterval("1", 76, 175), new SimpleInterval("1", 66, 185), readsSource),
                                new ReadShard(new SimpleInterval("1", 101, 200), new SimpleInterval("1", 91, 210), readsSource),
                                new ReadShard(new SimpleInterval("1", 126, 225), new SimpleInterval("1", 116, 235), readsSource),
                                new ReadShard(new SimpleInterval("1", 151, 250), new SimpleInterval("1", 141, 260), readsSource),
                                new ReadShard(new SimpleInterval("1", 176, 275), new SimpleInterval("1", 166, 285), readsSource),
                                new ReadShard(new SimpleInterval("1", 201, 300), new SimpleInterval("1", 191, 310), readsSource),
                                new ReadShard(new SimpleInterval("1", 226, 300), new SimpleInterval("1", 216, 310), readsSource),
                                new ReadShard(new SimpleInterval("1", 251, 300), new SimpleInterval("1", 241, 310), readsSource),
                                new ReadShard(new SimpleInterval("1", 276, 300), new SimpleInterval("1", 266, 310), readsSource))
                }
        };
    }

    @Test(dataProvider = "DivideIntervalIntoShardsTestData")
    public void testDivideIntervalIntoShards( final SimpleInterval originalInterval, final int shardSize, final int shardStep, final int shardPadding, final ReadsDataSource readsSource, final SAMSequenceDictionary dictionary, List<ReadShard> expectedShards  ) {
        // If we're invoked with shardSize == shardStep, invoke the version of divideIntervalIntoShards() that does
        // not take a shardStep parameter, in order to give it test coverage.
        final List<ReadShard> shards = (shardSize == shardStep) ?
                ReadShard.divideIntervalIntoShards(originalInterval, shardSize, shardPadding, readsSource, dictionary) :
                ReadShard.divideIntervalIntoShards(originalInterval, shardSize, shardStep, shardPadding, readsSource, dictionary);

        Assert.assertEquals(shards.size(), expectedShards.size(), "Wrong number of shards");
        for ( int i = 0; i < shards.size(); ++i ) {
            final ReadShard shard = shards.get(i);
            Assert.assertEquals(shard.getInterval(), expectedShards.get(i).getInterval(), "Shard has wrong size");
            Assert.assertEquals(shard.getPaddedInterval(), expectedShards.get(i).getPaddedInterval(), "Shard has wrong amount of padding");
        }
    }

    @DataProvider(name = "DivideIntervalIntoShardsInvalidTestData")
    public Object[][] divideIntervalIntoShardsInvalidTestData() {
        final ReadsDataSource readsSource = new ReadsDataSource(new File(publicTestDir + "org/broadinstitute/hellbender/engine/reads_data_source_test1.bam"));
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary(Arrays.asList(new SAMSequenceRecord("1", 16000)));

        return new Object[][] {
                // contig not in dictionary
                { new SimpleInterval("2", 1, 10), 1, 1, 0, readsSource, dictionary },
                // interval not within contig bounds
                { new SimpleInterval("1", 1, 17000), 1, 1, 0, readsSource, dictionary },
                // shardSize == 0
                { new SimpleInterval("1", 1, 10), 0, 1, 0, readsSource, dictionary },
                // shardSize < 0
                { new SimpleInterval("1", 1, 10), -1, 1, 0, readsSource, dictionary },
                // shardStep == 0
                { new SimpleInterval("1", 1, 10), 1, 0, 0, readsSource, dictionary },
                // shardStep < 0
                { new SimpleInterval("1", 1, 10), 1, -1, 0, readsSource, dictionary },
                // shardPadding < 0
                { new SimpleInterval("1", 1, 10), 1, 1, -1, readsSource, dictionary },
                // null interval
                { null, 1, 1, 0, readsSource, dictionary },
                // null data source
                { new SimpleInterval("1", 1, 10), 1, 1, 0, null, dictionary },
                // null dictionary
                { new SimpleInterval("1", 1, 10), 1, 1, 0, readsSource, null }
        };
    }

    @Test(dataProvider = "DivideIntervalIntoShardsInvalidTestData", expectedExceptions = IllegalArgumentException.class)
    public void testDivideIntervalIntoShardInvalidArgument( final SimpleInterval originalInterval, final int shardSize, final int shardStep, final int shardPadding, final ReadsDataSource readsSource, final SAMSequenceDictionary dictionary ) {
        ReadShard.divideIntervalIntoShards(originalInterval, shardSize, shardStep, shardPadding, readsSource, dictionary);
    }

    // Toy downsampler that keeps only reads with names "b" or "c"
    private static class KeepReadsBAndCOnlyDownsampler extends ReadsDownsampler {
        private List<GATKRead> finalizedReads = new ArrayList<>();

        @Override
        public boolean requiresCoordinateSortOrder() {
            return false;
        }

        @Override
        public void signalNoMoreReadsBefore( GATKRead read ) {
            // no-op
        }

        @Override
        public void submit( GATKRead item ) {
            if ( item.getName() != null && (item.getName().equals("b") || item.getName().equals("c")) ) {
                finalizedReads.add(item);
            }
            else {
                incrementNumberOfDiscardedItems(1);
            }
        }

        @Override
        public boolean hasFinalizedItems() {
            return ! finalizedReads.isEmpty();
        }

        @Override
        public List<GATKRead> consumeFinalizedItems() {
            final List<GATKRead> toReturn = finalizedReads;
            finalizedReads = new ArrayList<>();
            return toReturn;
        }

        @Override
        public boolean hasPendingItems() {
            return false;
        }

        @Override
        public GATKRead peekFinalized() {
            return hasFinalizedItems() ? finalizedReads.get(0) : null;
        }

        @Override
        public GATKRead peekPending() {
            return null;
        }

        @Override
        public int size() {
            return finalizedReads.size();
        }

        @Override
        public void signalEndOfInput() {
            // no-op
        }

        @Override
        public void clearItems() {
            finalizedReads.clear();
        }
    };
}
