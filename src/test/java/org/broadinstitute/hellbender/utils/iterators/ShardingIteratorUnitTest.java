package org.broadinstitute.hellbender.utils.iterators;

import com.google.common.collect.Iterators;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.Shard;
import org.broadinstitute.hellbender.engine.ShardBoundary;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class ShardingIteratorUnitTest extends GATKBaseTest {

    private static final SAMFileHeader TEST_HEADER = ArtificialReadUtils.createArtificialSamHeader();

    private static final List<GATKRead> oneBaseReadsUntilPosition(final int refIndex, final int endPos) {
        final List<GATKRead> reads = new ArrayList<>(endPos);
        // create a one-base read per position
        for (int i = 1; i <= endPos; i++) {
            reads.add(ArtificialReadUtils.createArtificialRead(TEST_HEADER, "read" + i, refIndex, i, 1));
        }
        return reads;
    }

    private static void testNextElement(final ShardingIterator<GATKRead> it,
            final SimpleInterval expectedInterval, final SimpleInterval expectedPaddedInterval,
            final int expectedCount) {
        // assert that there is a next one
        Assert.assertTrue(it.hasNext());
        final Shard<GATKRead> shard = it.next();

        // asserts that the shard is correct
        Assert.assertEquals(shard.getInterval(), expectedInterval);
        Assert.assertEquals(shard.getPaddedInterval(), expectedPaddedInterval);
        Assert.assertEquals(Utils.stream(shard.iterator()).count(), expectedCount);
        // assert that all elements overlaps
        assertShardElementsOverlaps(shard);
    }

    private static final void assertShardElementsOverlaps(final Shard<? extends Locatable> shard) {
        for (final Locatable locatable : shard) {
            Assert.assertTrue(locatable.overlaps(shard.getPaddedInterval()), String.format("%s does not overlap shard padded interval (%s)", locatable, shard.getPaddedInterval()));
        }
    }

    @DataProvider
    public Object[][] shardIteratorData() {
        // test data:
        // 1. start interval
        // 2. end interval
        // 3. padding
        return new Object[][] {
                // no padding interval
                {100, 200, 0},
                // padding interval
                {100, 200, 20}
        };
    }

    @Test(dataProvider = "shardIteratorData")
    public void testSingleShardIterator(final int start, final int end, final int pad) {
        // create a simple interval
        final SimpleInterval interval = new SimpleInterval(TEST_HEADER.getSequence(0).getSequenceName(), start, end);
        final SimpleInterval paddedInterval = interval.expandWithinContig(pad, TEST_HEADER.getSequenceDictionary());
        // create a list of reads larger than the interval
        final ShardingIterator<GATKRead> it = new ShardingIterator<>(
                // always from the beginning and with at least one read more
                oneBaseReadsUntilPosition(0, paddedInterval.getEnd() + 1).iterator(),
                Collections.singletonList(new ShardBoundary(interval, paddedInterval)),
                TEST_HEADER.getSequenceDictionary());

        // test the next element
        testNextElement(it, interval, paddedInterval, paddedInterval.getLengthOnReference());
        // assert that it is exhausted
        Assert.assertFalse(it.hasNext());
    }

    @Test(dataProvider = "shardIteratorData")
    public void testDataOnlyAfterShards(final int start, final int end, final int pad) {
        // create a simple interval
        final SimpleInterval interval = new SimpleInterval(TEST_HEADER.getSequence(0).getSequenceName(), start, end);
        final SimpleInterval paddedInterval = interval.expandWithinContig(pad, TEST_HEADER.getSequenceDictionary());
        // create a list of reads larger than the interval
        final ShardingIterator<GATKRead> it = new ShardingIterator<>(
                // create on a different contig (after)
                oneBaseReadsUntilPosition(1, 10).iterator(),
                Collections.singletonList(new ShardBoundary(interval, paddedInterval)),
                TEST_HEADER.getSequenceDictionary());

        // test the next element
        testNextElement(it, interval, paddedInterval, 0);
        // assert that it is exhausted
        Assert.assertFalse(it.hasNext());
    }

    @Test(dataProvider = "shardIteratorData")
    public void testDataOnlyBeforeShards(final int start, final int end, final int pad) {
        // create a simple interval
        final SimpleInterval interval = new SimpleInterval(TEST_HEADER.getSequence(1).getSequenceName(), start, end);
        final SimpleInterval paddedInterval = interval.expandWithinContig(pad, TEST_HEADER.getSequenceDictionary());
        // create a list of reads larger than the interval
        final ShardingIterator<GATKRead> it = new ShardingIterator<>(
                // create on a previous contig
                oneBaseReadsUntilPosition(0, 10).iterator(),
                Collections.singletonList(new ShardBoundary(interval, paddedInterval)),
                TEST_HEADER.getSequenceDictionary());

        // test the next element
        testNextElement(it, interval, paddedInterval, 0);
        // assert that it is exhausted
        Assert.assertFalse(it.hasNext());
    }

    @Test(dataProvider = "shardIteratorData")
    public void testNonOverlappingShards(final int start, final int end, final int pad) {
        // create boundaries with intervals
        final List<GATKRead> reads = new ArrayList<>();
        final List<ShardBoundary> shards = new ArrayList<>();
        for (int i = 0; i < 2; i++) {
            final SimpleInterval interval = new SimpleInterval(TEST_HEADER.getSequence(i).getSequenceName(), start, end);
            final SimpleInterval paddedInterval = interval.expandWithinContig(pad, TEST_HEADER.getSequenceDictionary());
            reads.addAll(oneBaseReadsUntilPosition(i, paddedInterval.getEnd() + 1));
            shards.add(new ShardBoundary(interval, paddedInterval));
        }

        // create shard iterator
        // create a list of reads larger than the interval
        final ShardingIterator<GATKRead> it = new ShardingIterator<>(
                reads.iterator(),
                shards,
                TEST_HEADER.getSequenceDictionary());

        for (final ShardBoundary s: shards) {
            testNextElement(it, s.getInterval(), s.getPaddedInterval(), s.getPaddedInterval().getLengthOnReference());
        }

        // assert that the iterator is exhausted
        Assert.assertFalse(it.hasNext());
    }

    @Test(dataProvider = "shardIteratorData")
    public void testOverlappingShards(final int start, final int end, final int pad) {
        // first create the
        // create a simple interval
        final SimpleInterval interval = new SimpleInterval(TEST_HEADER.getSequence(0).getSequenceName(), start, end);
        final List<GATKRead> reads = oneBaseReadsUntilPosition(0, interval.getEnd() + pad + 1);
        final List<ShardBoundary> shards = Shard.divideIntervalIntoShards(interval, end/20, end/10, pad, TEST_HEADER.getSequenceDictionary());

        // create shard iterator
        // create a list of reads larger than the interval
        final ShardingIterator<GATKRead> it = new ShardingIterator<>(
                reads.iterator(),
                shards,
                TEST_HEADER.getSequenceDictionary());


        for (final ShardBoundary s: shards) {
            testNextElement(it, s.getInterval(), s.getPaddedInterval(), s.getPaddedInterval().getLengthOnReference());
        }

        // assert that the iterator is exhausted
        Assert.assertFalse(it.hasNext());
    }

    @Test
    public void testDifferentSizeReads() {
        final SimpleInterval region = new SimpleInterval(TEST_HEADER.getSequence(0).getSequenceName(), 1, 100);
        final List<ShardBoundary> shards = Shard.divideIntervalIntoShards(region, 50, 50, 0, TEST_HEADER.getSequenceDictionary());

        final GATKRead regionRead = ArtificialReadUtils.createArtificialRead(TEST_HEADER, "regionRead", 0, region.getStart(), region.getLengthOnReference());
        final GATKRead firstWindowMinusOneRead = ArtificialReadUtils.createArtificialRead(TEST_HEADER, "firstRead", 0, shards.get(0).getStart() + 1, shards.get(0).getLengthOnReference() - 1);

        final List<GATKRead> reads = Arrays.asList(regionRead, firstWindowMinusOneRead);

        final ShardingIterator<GATKRead> it = new ShardingIterator<>(
                reads.iterator(),
                shards,
                TEST_HEADER.getSequenceDictionary());

        // the first window should contain both reads
        final Shard<GATKRead> first = it.next();
        assertShardElementsOverlaps(first);
        Assert.assertEquals(first.iterator(), reads.iterator());


        // the second window should contain only the overlapping read
        final Shard<GATKRead> second = it.next();
        assertShardElementsOverlaps(second);
        Assert.assertEquals(second.iterator(), Iterators.singletonIterator(regionRead));

        // no more if the test is correctly implemented
        Assert.assertFalse(it.hasNext());
    }
}