package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.downsampling.ReadsDownsampler;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

public class MultiIntervalLocalReadShardUnitTest extends GATKBaseTest {

    @DataProvider
    public Object[][] shardBoundsTestData() {
        final ReadsDataSource readsSource = new ReadsDataSource(IOUtils.getPath(publicTestDir + "org/broadinstitute/hellbender/engine/reads_data_source_test1.bam"));

        return new Object[][] {
                // Shard, expected shard intervals, expected padded shard intervals
                { new MultiIntervalLocalReadShard(Arrays.asList(new SimpleInterval("1", 200, 250)), 0, readsSource), Arrays.asList(new SimpleInterval("1", 200, 250)), Arrays.asList(new SimpleInterval("1", 200, 250)) },
                { new MultiIntervalLocalReadShard(Arrays.asList(new SimpleInterval("1", 200, 250)), 25, readsSource), Arrays.asList(new SimpleInterval("1", 200, 250)), Arrays.asList(new SimpleInterval("1", 175, 275)) },
                { new MultiIntervalLocalReadShard(Arrays.asList(new SimpleInterval("1", 200, 250)), 50, readsSource), Arrays.asList(new SimpleInterval("1", 200, 250)), Arrays.asList(new SimpleInterval("1", 150, 300)) },
                { new MultiIntervalLocalReadShard(Arrays.asList(new SimpleInterval("1", 200, 250)), 500, readsSource), Arrays.asList(new SimpleInterval("1", 200, 250)), Arrays.asList(new SimpleInterval("1", 1, 750)) },

                { new MultiIntervalLocalReadShard(Arrays.asList(new SimpleInterval("1", 100, 200), new SimpleInterval("1", 250, 350)), 0, readsSource), Arrays.asList(new SimpleInterval("1", 100, 200), new SimpleInterval("1", 250, 350)), Arrays.asList(new SimpleInterval("1", 100, 200), new SimpleInterval("1", 250, 350)) },
                { new MultiIntervalLocalReadShard(Arrays.asList(new SimpleInterval("1", 100, 200), new SimpleInterval("1", 250, 350)), 20, readsSource), Arrays.asList(new SimpleInterval("1", 100, 200), new SimpleInterval("1", 250, 350)), Arrays.asList(new SimpleInterval("1", 80, 220), new SimpleInterval("1", 230, 370)) },
                { new MultiIntervalLocalReadShard(Arrays.asList(new SimpleInterval("1", 100, 200), new SimpleInterval("1", 250, 350)), 50, readsSource), Arrays.asList(new SimpleInterval("1", 100, 200), new SimpleInterval("1", 250, 350)), Arrays.asList(new SimpleInterval("1", 50, 400)) },
                { new MultiIntervalLocalReadShard(Arrays.asList(new SimpleInterval("1", 100, 200), new SimpleInterval("1", 250, 350)), 100, readsSource), Arrays.asList(new SimpleInterval("1", 100, 200), new SimpleInterval("1", 250, 350)), Arrays.asList(new SimpleInterval("1", 1, 450)) },

                { new MultiIntervalLocalReadShard(Arrays.asList(new SimpleInterval("1", 100, 200), new SimpleInterval("1", 250, 350), new SimpleInterval("2", 1, 10)), 0, readsSource), Arrays.asList(new SimpleInterval("1", 100, 200), new SimpleInterval("1", 250, 350), new SimpleInterval("2", 1, 10)), Arrays.asList(new SimpleInterval("1", 100, 200), new SimpleInterval("1", 250, 350), new SimpleInterval("2", 1, 10)) },
                { new MultiIntervalLocalReadShard(Arrays.asList(new SimpleInterval("1", 100, 200), new SimpleInterval("1", 250, 350), new SimpleInterval("2", 1, 10)), 50, readsSource), Arrays.asList(new SimpleInterval("1", 100, 200), new SimpleInterval("1", 250, 350), new SimpleInterval("2", 1, 10)), Arrays.asList(new SimpleInterval("1", 50, 400), new SimpleInterval("2", 1, 60)) },
        };
    }

    @Test(dataProvider = "shardBoundsTestData")
    public void testShardBounds( final MultiIntervalLocalReadShard shard, final List<SimpleInterval> expectedShardIntervals, final List<SimpleInterval> expectedPaddedShardIntervals) {
        Assert.assertEquals(shard.getIntervals(), expectedShardIntervals, "shard has wrong intervals");
        Assert.assertEquals(shard.getPaddedIntervals(), expectedPaddedShardIntervals, "shard has wrong padded intervals");
    }

    private static final class KeepNamedReadsDownsampler extends ReadsDownsampler {
        private List<GATKRead> finalizedReads = new ArrayList<>();
        private final List<String> readNamesToKeep;

        public KeepNamedReadsDownsampler(final List<String> readNamesToKeep) {
            this.readNamesToKeep = readNamesToKeep;
        }

        @Override
        public boolean requiresCoordinateSortOrder() {
            return false;
        }

        @Override
        public void signalNoMoreReadsBefore( GATKRead read ) {
            // no-op
        }

        @Override
        public void submit( GATKRead read ) {
            if ( read.getName() != null && readNamesToKeep.contains(read.getName()) ) {
                finalizedReads.add(read);
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
    }


    @DataProvider
    public Object[][] shardIterationTestData() {
        final ReadsDataSource readsSource = new ReadsDataSource(IOUtils.getPath(publicTestDir + "org/broadinstitute/hellbender/engine/reads_data_source_test1.bam"));

        final ReadFilter keepReadBOnly = new ReadFilter() {
            private static final long serialVersionUID = 1l;
            @Override
            public boolean test( GATKRead read ) { return read.getName().equals("b"); };
        };
        final MultiIntervalLocalReadShard filteredShard = new MultiIntervalLocalReadShard(Arrays.asList(new SimpleInterval("1", 200, 210)), 0, readsSource);
        filteredShard.setReadFilter(keepReadBOnly);

        final ReadsDownsampler readsBAndCOnlyDownsampler = new KeepNamedReadsDownsampler(Arrays.asList("b", "c"));
        final MultiIntervalLocalReadShard downsampledShard = new MultiIntervalLocalReadShard(Arrays.asList(new SimpleInterval("1", 1, 5000)), 0, readsSource);
        downsampledShard.setDownsampler(readsBAndCOnlyDownsampler);

        final ReadsDownsampler allReadsExceptDAndFDownsampler = new KeepNamedReadsDownsampler(Arrays.asList("a", "b", "c", "e", "g", "h", "i", "j", "k"));

        final MultiIntervalLocalReadShard multiIntervalDownsampledShard = new MultiIntervalLocalReadShard(Arrays.asList(new SimpleInterval("1", 205, 1500), new SimpleInterval("2", 1, 600)), 0, readsSource);
        multiIntervalDownsampledShard.setDownsampler(allReadsExceptDAndFDownsampler);

        final MultiIntervalLocalReadShard multiIntervalDownsampledShardWithPadding = new MultiIntervalLocalReadShard(Arrays.asList(new SimpleInterval("1", 205, 1500), new SimpleInterval("2", 1, 600)), 50, readsSource);
        multiIntervalDownsampledShardWithPadding.setDownsampler(allReadsExceptDAndFDownsampler);

        final ReadTransformer renameToA = read -> {
            read.setName("a");
            return read;
        };
        final ReadTransformer renameToB = read -> {
            read.setName("b");
            return read;
        };
        final MultiIntervalLocalReadShard transformedInOrderToPassFilterShard = new MultiIntervalLocalReadShard(Arrays.asList(new SimpleInterval("1", 200, 210)), 0, readsSource);
        transformedInOrderToPassFilterShard.setPreReadFilterTransformer(renameToB);
        transformedInOrderToPassFilterShard.setReadFilter(keepReadBOnly);

        final MultiIntervalLocalReadShard transformedAfterFilterShard = new MultiIntervalLocalReadShard(Arrays.asList(new SimpleInterval("1", 200, 210)), 0, readsSource);
        transformedAfterFilterShard.setPostReadFilterTransformer(renameToA);
        transformedAfterFilterShard.setReadFilter(keepReadBOnly);

        return new Object[][] {
                // shard, expected read names after iteration
                {new MultiIntervalLocalReadShard(Arrays.asList(new SimpleInterval("1", 200, 210)), 0, readsSource), Arrays.asList("a", "b", "c") },
                {new MultiIntervalLocalReadShard(Arrays.asList(new SimpleInterval("1", 200, 209)), 0, readsSource), Arrays.asList("a", "b") },
                {new MultiIntervalLocalReadShard(Arrays.asList(new SimpleInterval("1", 200, 204)), 0, readsSource), Arrays.asList("a") },
                {new MultiIntervalLocalReadShard(Arrays.asList(new SimpleInterval("1", 400, 500)), 0, readsSource), Collections.<String>emptyList() },
                {new MultiIntervalLocalReadShard(Arrays.asList(new SimpleInterval("1", 205, 950), new SimpleInterval("1", 1200, 1500)), 0, readsSource), Arrays.asList("a", "b", "c")},
                {new MultiIntervalLocalReadShard(Arrays.asList(new SimpleInterval("1", 205, 950), new SimpleInterval("1", 1200, 1500)), 50, readsSource), Arrays.asList("a", "b", "c", "d", "e")},
                {new MultiIntervalLocalReadShard(Arrays.asList(new SimpleInterval("1", 205, 208), new SimpleInterval("1", 1000, 1070)), 0, readsSource), Arrays.asList("a", "b", "d")},
                {new MultiIntervalLocalReadShard(Arrays.asList(new SimpleInterval("1", 205, 208), new SimpleInterval("1", 1000, 1070)), 2, readsSource), Arrays.asList("a", "b", "c", "d")},
                {new MultiIntervalLocalReadShard(Arrays.asList(new SimpleInterval("1", 205, 208), new SimpleInterval("1", 1000, 1070)), 30, readsSource), Arrays.asList("a", "b", "c", "d", "e")},
                {filteredShard, Arrays.asList("b") },
                {downsampledShard, Arrays.asList("b", "c")},
                {transformedInOrderToPassFilterShard, Arrays.asList("b", "b", "b")},
                {transformedAfterFilterShard, Arrays.asList("a")},   
                {multiIntervalDownsampledShard, Arrays.asList("a", "b", "c", "e", "g")},
                {multiIntervalDownsampledShardWithPadding, Arrays.asList("a", "b", "c", "e", "g", "h")}
        };
    }

    @Test(dataProvider = "shardIterationTestData")
    public void testShardIteration(final MultiIntervalLocalReadShard shard, final List<String> expectedReadNames) {
        final List<String> actualReadNames = new ArrayList<>();

        for ( final GATKRead read : shard ) {
            actualReadNames.add(read.getName());
        }

        Assert.assertEquals(actualReadNames, expectedReadNames, "Wrong reads returned");
    }
}
