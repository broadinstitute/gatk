package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;

public class ShardUnitTest extends GATKBaseTest {

    @DataProvider(name = "DivideIntervalIntoShardsTestData")
    public Object[][] divideIntervalIntoShardsTestData() {
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary(Arrays.asList(new SAMSequenceRecord("1", 16000)));

        return new Object[][] {
                // shardSize == 100, shardStep == 100, shardPadding == 0, start of contig
                {
                        new SimpleInterval("1", 1, 100), 100, 100, 0, dictionary,
                        Arrays.asList(new ShardBoundary(new SimpleInterval("1", 1, 100), new SimpleInterval("1", 1, 100)))
                },

                // shardSize == 100, shardStep == 100, shardPadding == 10, start of contig
                {
                        new SimpleInterval("1", 1, 100), 100, 100, 10, dictionary,
                        Arrays.asList(new ShardBoundary(new SimpleInterval("1", 1, 100), new SimpleInterval("1", 1, 110)))
                },

                // shardSize == 50, shardStep == 50, shardPadding == 0, start of contig
                {
                        new SimpleInterval("1", 1, 100), 50, 50, 0, dictionary,
                        Arrays.asList(new ShardBoundary(new SimpleInterval("1", 1, 50), new SimpleInterval("1", 1, 50)),
                                new ShardBoundary(new SimpleInterval("1", 51, 100), new SimpleInterval("1", 51, 100)))
                },

                // shardSize == 50, shardStep == 50, shardPadding == 10, start of contig
                {
                        new SimpleInterval("1", 1, 100), 50, 50, 10, dictionary,
                        Arrays.asList(new ShardBoundary(new SimpleInterval("1", 1, 50), new SimpleInterval("1", 1, 60)),
                                new ShardBoundary(new SimpleInterval("1", 51, 100), new SimpleInterval("1", 41, 110)))
                },

                // shardSize == 50, shardStep == 50, shardPadding == 0
                {
                        new SimpleInterval("1", 51, 200), 50, 50, 0, dictionary,
                        Arrays.asList(new ShardBoundary(new SimpleInterval("1", 51, 100), new SimpleInterval("1", 51, 100)),
                                new ShardBoundary(new SimpleInterval("1", 101, 150), new SimpleInterval("1", 101, 150)),
                                new ShardBoundary(new SimpleInterval("1", 151, 200), new SimpleInterval("1", 151, 200)))
                },

                // shardSize == 50, shardStep == 50, shardPadding == 10
                {
                        new SimpleInterval("1", 51, 200), 50, 50, 10, dictionary,
                        Arrays.asList(new ShardBoundary(new SimpleInterval("1", 51, 100), new SimpleInterval("1", 41, 110)),
                                new ShardBoundary(new SimpleInterval("1", 101, 150), new SimpleInterval("1", 91, 160)),
                                new ShardBoundary(new SimpleInterval("1", 151, 200), new SimpleInterval("1", 141, 210)))
                },

                // shardSize == 70, shardStep == 70, shardPadding == 0
                {
                        new SimpleInterval("1", 51, 200), 70, 70, 0, dictionary,
                        Arrays.asList(new ShardBoundary(new SimpleInterval("1", 51, 120), new SimpleInterval("1", 51, 120)),
                                new ShardBoundary(new SimpleInterval("1", 121, 190), new SimpleInterval("1", 121, 190)),
                                new ShardBoundary(new SimpleInterval("1", 191, 200), new SimpleInterval("1", 191, 200)))
                },

                // shardSize == 70, shardStep == 70, shardPadding == 10
                {
                        new SimpleInterval("1", 51, 200), 70, 70, 10, dictionary,
                        Arrays.asList(new ShardBoundary(new SimpleInterval("1", 51, 120), new SimpleInterval("1", 41, 130)),
                                new ShardBoundary(new SimpleInterval("1", 121, 190), new SimpleInterval("1", 111, 200)),
                                new ShardBoundary(new SimpleInterval("1", 191, 200), new SimpleInterval("1", 181, 210)))
                },

                // shardSize == 100, shardStep == 100, shardPadding == 10, end of contig
                {
                        new SimpleInterval("1", 15999, 16000), 100, 100, 10, dictionary,
                        Arrays.asList(new ShardBoundary(new SimpleInterval("1", 15999, 16000), new SimpleInterval("1", 15989, 16000)))
                },

                // shardSize == 100, shardStep == 100, shardPadding == 0
                {
                        new SimpleInterval("1", 51, 300), 100, 100, 0, dictionary,
                        Arrays.asList(new ShardBoundary(new SimpleInterval("1", 51, 150), new SimpleInterval("1", 51, 150)),
                                new ShardBoundary(new SimpleInterval("1", 151, 250), new SimpleInterval("1", 151, 250)),
                                new ShardBoundary(new SimpleInterval("1", 251, 300), new SimpleInterval("1", 251, 300)))
                },

                // shardSize == 100, shardStep == 50, shardPadding == 0
                {
                        new SimpleInterval("1", 51, 300), 100, 50, 0, dictionary,
                        Arrays.asList(new ShardBoundary(new SimpleInterval("1", 51, 150), new SimpleInterval("1", 51, 150)),
                                new ShardBoundary(new SimpleInterval("1", 101, 200), new SimpleInterval("1", 101, 200)),
                                new ShardBoundary(new SimpleInterval("1", 151, 250), new SimpleInterval("1", 151, 250)),
                                new ShardBoundary(new SimpleInterval("1", 201, 300), new SimpleInterval("1", 201, 300)),
                                new ShardBoundary(new SimpleInterval("1", 251, 300), new SimpleInterval("1", 251, 300)))
                },

                // shardSize == 100, shardStep == 50, shardPadding == 10
                {
                        new SimpleInterval("1", 51, 300), 100, 50, 10, dictionary,
                        Arrays.asList(new ShardBoundary(new SimpleInterval("1", 51, 150), new SimpleInterval("1", 41, 160)),
                                new ShardBoundary(new SimpleInterval("1", 101, 200), new SimpleInterval("1", 91, 210)),
                                new ShardBoundary(new SimpleInterval("1", 151, 250), new SimpleInterval("1", 141, 260)),
                                new ShardBoundary(new SimpleInterval("1", 201, 300), new SimpleInterval("1", 191, 310)),
                                new ShardBoundary(new SimpleInterval("1", 251, 300), new SimpleInterval("1", 241, 310)))
                },

                // shardSize == 100, shardStep == 25 shardPadding == 0
                {
                        new SimpleInterval("1", 51, 300), 100, 25, 0, dictionary,
                        Arrays.asList(new ShardBoundary(new SimpleInterval("1", 51, 150), new SimpleInterval("1", 51, 150)),
                                new ShardBoundary(new SimpleInterval("1", 76, 175), new SimpleInterval("1", 76, 175)),
                                new ShardBoundary(new SimpleInterval("1", 101, 200), new SimpleInterval("1", 101, 200)),
                                new ShardBoundary(new SimpleInterval("1", 126, 225), new SimpleInterval("1", 126, 225)),
                                new ShardBoundary(new SimpleInterval("1", 151, 250), new SimpleInterval("1", 151, 250)),
                                new ShardBoundary(new SimpleInterval("1", 176, 275), new SimpleInterval("1", 176, 275)),
                                new ShardBoundary(new SimpleInterval("1", 201, 300), new SimpleInterval("1", 201, 300)),
                                new ShardBoundary(new SimpleInterval("1", 226, 300), new SimpleInterval("1", 226, 300)),
                                new ShardBoundary(new SimpleInterval("1", 251, 300), new SimpleInterval("1", 251, 300)),
                                new ShardBoundary(new SimpleInterval("1", 276, 300), new SimpleInterval("1", 276, 300)))
                },

                // shardSize == 100, shardStep == 25, shardPadding == 10
                {
                        new SimpleInterval("1", 51, 300), 100, 25, 10, dictionary,
                        Arrays.asList(new ShardBoundary(new SimpleInterval("1", 51, 150), new SimpleInterval("1", 41, 160)),
                                new ShardBoundary(new SimpleInterval("1", 76, 175), new SimpleInterval("1", 66, 185)),
                                new ShardBoundary(new SimpleInterval("1", 101, 200), new SimpleInterval("1", 91, 210)),
                                new ShardBoundary(new SimpleInterval("1", 126, 225), new SimpleInterval("1", 116, 235)),
                                new ShardBoundary(new SimpleInterval("1", 151, 250), new SimpleInterval("1", 141, 260)),
                                new ShardBoundary(new SimpleInterval("1", 176, 275), new SimpleInterval("1", 166, 285)),
                                new ShardBoundary(new SimpleInterval("1", 201, 300), new SimpleInterval("1", 191, 310)),
                                new ShardBoundary(new SimpleInterval("1", 226, 300), new SimpleInterval("1", 216, 310)),
                                new ShardBoundary(new SimpleInterval("1", 251, 300), new SimpleInterval("1", 241, 310)),
                                new ShardBoundary(new SimpleInterval("1", 276, 300), new SimpleInterval("1", 266, 310)))
                }
        };
    }

    @Test(dataProvider = "DivideIntervalIntoShardsTestData")
    public void testDivideIntervalIntoShards( final SimpleInterval originalInterval, final int shardSize, final int shardStep, final int shardPadding, final SAMSequenceDictionary dictionary, List<ShardBoundary> expectedShards  ) {
        final List<ShardBoundary> shards = Shard.divideIntervalIntoShards(originalInterval, shardSize, shardStep, shardPadding, dictionary);

        Assert.assertEquals(shards.size(), expectedShards.size(), "Wrong number of shards");
        for ( int i = 0; i < shards.size(); ++i ) {
            final ShardBoundary shard = shards.get(i);
            Assert.assertEquals(shard.getInterval(), expectedShards.get(i).getInterval(), "Shard has wrong size");
            Assert.assertEquals(shard.getPaddedInterval(), expectedShards.get(i).getPaddedInterval(), "Shard has wrong amount of padding");
        }
    }

    @DataProvider(name = "DivideIntervalIntoShardsInvalidTestData")
    public Object[][] divideIntervalIntoShardsInvalidTestData() {
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary(Arrays.asList(new SAMSequenceRecord("1", 16000)));

        return new Object[][] {
                // contig not in dictionary
                { new SimpleInterval("2", 1, 10), 1, 1, 0, dictionary },
                // interval not within contig bounds
                { new SimpleInterval("1", 1, 17000), 1, 1, 0, dictionary },
                // shardSize == 0
                { new SimpleInterval("1", 1, 10), 0, 1, 0, dictionary },
                // shardSize < 0
                { new SimpleInterval("1", 1, 10), -1, 1, 0, dictionary },
                // shardStep == 0
                { new SimpleInterval("1", 1, 10), 1, 0, 0, dictionary },
                // shardStep < 0
                { new SimpleInterval("1", 1, 10), 1, -1, 0, dictionary },
                // shardPadding < 0
                { new SimpleInterval("1", 1, 10), 1, 1, -1, dictionary },
                // null interval
                { null, 1, 1, 0, dictionary },
                // null dictionary
                { new SimpleInterval("1", 1, 10), 1, 1, 0, null }
        };
    }

    @Test(dataProvider = "DivideIntervalIntoShardsInvalidTestData", expectedExceptions = IllegalArgumentException.class)
    public void testDivideIntervalIntoShardInvalidArgument( final SimpleInterval originalInterval, final int shardSize, final int shardStep, final int shardPadding, final SAMSequenceDictionary dictionary ) {
        Shard.divideIntervalIntoShards(originalInterval, shardSize, shardStep, shardPadding, dictionary);
    }

}
