package org.broadinstitute.hellbender.utils.downsampling;

import com.google.common.collect.ImmutableMap;
import htsjdk.samtools.SAMFileHeader;
import org.apache.commons.lang.mutable.MutableInt;
import org.broadinstitute.hellbender.tools.exome.HashedListTargetCollection;
import org.broadinstitute.hellbender.tools.exome.TargetCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.*;
import org.broadinstitute.hellbender.utils.test.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;

public class PositionalDownsamplerUnitTest extends GATKBaseTest {

    private final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();

    private static final SimpleInterval DUMMY_INTERVAL_FOR_READS_WITH_NO_POSITION =
            new SimpleInterval(ReadConstants.UNSET_CONTIG, 1, 1);
    private static final SimpleInterval chr1_1_1 = new SimpleInterval("1", 1, 1);
    private static final SimpleInterval chr1_1_5 = new SimpleInterval("1", 1, 5);
    private static final SimpleInterval chr1_2_2 = new SimpleInterval("1", 2, 2);
    private static final SimpleInterval chr2_1_1 = new SimpleInterval("2", 1, 1);
    private static final SimpleInterval chr2_1_5 = new SimpleInterval("2", 1, 5);
    private static final SimpleInterval chr2_2_2 = new SimpleInterval("2", 2, 2);

    @DataProvider(name = "PositionalDownsamplerTestData")
    public Object[][] positionalDownsamplerTestData() {
        return new Object[][] {
                // List<List<>> of reads, downsampling coverage, stride length, and Map of stride intervals to number of reads contained

                // all mapped reads
                { Arrays.asList(createStackOfMappedReads(10, "1", 1)), 1, 1, ImmutableMap.of(chr1_1_1, 1)},
                { Arrays.asList(createStackOfMappedReads(10, "1", 1)), 1, 5, ImmutableMap.of(chr1_1_5, 1)},
                { Arrays.asList(createStackOfMappedReads(10, "1", 1)), 5, 1, ImmutableMap.of(chr1_1_1, 5)},
                { Arrays.asList(createStackOfMappedReads(10, "1", 1)), 5, 5, ImmutableMap.of(chr1_1_5, 5)},
                { Arrays.asList(createStackOfMappedReads(10, "1", 1)), 10, 1, ImmutableMap.of(chr1_1_1, 10)},
                { Arrays.asList(createStackOfMappedReads(10, "1", 1)), 15, 1, ImmutableMap.of(chr1_1_1, 10)},
                { Arrays.asList(createStackOfMappedReads(10, "1", 1), createStackOfMappedReads(8, "1", 2)), 5, 1, ImmutableMap.of(chr1_1_1, 5, chr1_2_2, 5)},
                { Arrays.asList(createStackOfMappedReads(1, "1", 1), createStackOfMappedReads(1, "1", 4), createStackOfMappedReads(1, "1", 6), createStackOfMappedReads(1, "1", 11)), 1, 5, ImmutableMap.of(new SimpleInterval("1", 1, 4), 1, new SimpleInterval("1", 6, 10), 1, new SimpleInterval("1", 11, 11), 1)},
                { Arrays.asList(createStackOfMappedReads(10, "1", 1), createStackOfMappedReads(8, "1", 2)), 5, 5, ImmutableMap.of(chr1_1_5, 5)},
                { Arrays.asList(createStackOfMappedReads(10, "1", 1), createStackOfMappedReads(8, "1", 2)), 14, 5, ImmutableMap.of(chr1_1_5, 14)},
                { Arrays.asList(createStackOfMappedReads(10, "1", 1), createStackOfMappedReads(5, "1", 2)), 5, 1, ImmutableMap.of(chr1_1_1, 5, chr1_2_2, 5) },
                { Arrays.asList(createStackOfMappedReads(10, "1", 1), createStackOfMappedReads(3, "1", 2)), 5, 1, ImmutableMap.of(chr1_1_1, 5, chr1_2_2, 3)},
                { Arrays.asList(createStackOfMappedReads(10, "1", 1), createStackOfMappedReads(10, "2", 1)), 5, 1, ImmutableMap.of(chr1_1_1, 5, chr2_1_1, 5)},
                { Arrays.asList(createStackOfMappedReads(10, "1", 1), createStackOfMappedReads(10, "2", 1)), 5, 5, ImmutableMap.of(chr1_1_1, 5, chr2_1_1, 5)},
                { Arrays.asList(createStackOfMappedReads(10, "1", 1), createStackOfMappedReads(10, "2", 1), createStackOfMappedReads(10, "2", 2)), 5, 1, ImmutableMap.of(chr1_1_1, 5, chr2_1_1, 5, chr2_2_2, 5)},
                { Arrays.asList(createStackOfMappedReads(10, "1", 1), createStackOfMappedReads(10, "2", 1), createStackOfMappedReads(10, "2", 2)), 10, 1, ImmutableMap.of(chr1_1_1, 10, chr2_1_1, 10, chr2_2_2, 10) },
                { Arrays.asList(createStackOfMappedReads(10, "1", 1), createStackOfMappedReads(10, "2", 1), createStackOfMappedReads(10, "2", 2)), 12, 5, ImmutableMap.of(chr1_1_5, 10, chr2_1_5, 12) },
                { Arrays.asList(createStackOfMappedReads(3, "1", 1), createStackOfMappedReads(10, "2", 1), createStackOfMappedReads(4, "2", 2)), 5, 1, ImmutableMap.of(chr1_1_1, 3, chr2_1_1, 5, chr2_2_2, 4)},
                { Arrays.asList(createStackOfMappedReads(3, "1", 1), createStackOfMappedReads(10, "2", 1), createStackOfMappedReads(4, "2", 2)), 1, 1, ImmutableMap.of(chr1_1_1, 1, chr2_1_1, 1, chr2_2_2, 1)},

                // Unmapped reads with no positions should not get downsampled *and* they should not interrupt a stride
                { Arrays.asList(createStackOfMappedReads(10, "1", 1), createStackOfMappedReads(10, "2", 1), createStackOfMappedReads(10, "2", 2), createStackOfUnmappedReads(100)), 5, 1, ImmutableMap.of(chr1_1_1, 5, chr2_1_1, 5, chr2_2_2, 5, DUMMY_INTERVAL_FOR_READS_WITH_NO_POSITION, 100)},
                { Arrays.asList(createStackOfMappedReads(3, "1", 1), createStackOfUnmappedReads(100), createStackOfMappedReads(3, "1", 2)), 3, 5, ImmutableMap.of(chr1_1_5, 3, DUMMY_INTERVAL_FOR_READS_WITH_NO_POSITION, 100)},

                // Unmapped reads with assigned positions should get downsampled
                { Arrays.asList(createStackOfMappedReads(10, "1", 1), createStackOfUnmappedReadsWithPosition(10, "1", 1)), 5, 1, ImmutableMap.of(chr1_1_1, 5)},
                { Arrays.asList(createStackOfMappedReads(10, "1", 1), createStackOfUnmappedReadsWithPosition(10, "1", 1), createStackOfMappedReads(5, "2", 1), createStackOfUnmappedReadsWithPosition(5, "2", 1)), 5, 1, ImmutableMap.of(chr1_1_1, 5, chr2_1_1, 5) },
                { Arrays.asList(createStackOfMappedReads(3, "1", 1), createStackOfUnmappedReadsWithPosition(100, "1",2), createStackOfMappedReads(3, "1", 3)), 3, 5, ImmutableMap.of(chr1_1_5, 3)},

                // Mix of mapped, pure unmapped, and mapped with assigned position
                { Arrays.asList(createStackOfMappedReads(10, "1", 1), createStackOfUnmappedReadsWithPosition(10, "1", 1), createStackOfMappedReads(5, "2", 1), createStackOfUnmappedReadsWithPosition(5, "2", 1), createStackOfUnmappedReads(100)), 5, 1, ImmutableMap.of(chr1_1_1, 5, chr2_1_1, 5, DUMMY_INTERVAL_FOR_READS_WITH_NO_POSITION, 100) },

                // Unmapped only
                { Arrays.asList(createStackOfUnmappedReads(10)), 1, 1, ImmutableMap.of(DUMMY_INTERVAL_FOR_READS_WITH_NO_POSITION, 10)},

                // No reads
                { Collections.emptyList(), 1, 1, Collections.emptyMap()}
        };
    }

    @Test(dataProvider = "PositionalDownsamplerTestData")
    public void testPositionalDownsampler( final List<List<GATKRead>> reads, final int targetCoverage, final int strideLength,
                                           final Map<SimpleInterval, Integer> expectedCountsByStride ) {
        final boolean biasToHigherMappingQuality = false;
        final List<GATKRead> allReads = new ArrayList<>();
        for ( final List<GATKRead> readStack : reads ) {
            allReads.addAll(readStack);
        }
        final int expectedDownsampledReads = expectedCountsByStride.values().stream().mapToInt(Integer::intValue).sum();

        final ReadsDownsampler downsampler = new PositionalDownsampler(targetCoverage, strideLength, biasToHigherMappingQuality, Integer.MAX_VALUE);
        Assert.assertTrue(downsampler.requiresCoordinateSortOrder());
        Assert.assertEquals(downsampler.size(), 0);

        Utils.resetRandomGenerator();
        downsampler.submit(allReads);
        Assert.assertEquals(downsampler.size(), expectedDownsampledReads);

        // there will be finalized reads if either 1) any read had no assigned position or 2) at least one stride was completed,
        // triggering finalization of its reads
        final boolean expectFinalizedItemsAfterSubmission =
                expectedCountsByStride.getOrDefault(DUMMY_INTERVAL_FOR_READS_WITH_NO_POSITION, 0) > 0
                || expectedCountsByStride.keySet().size() > 1;
        Assert.assertEquals(downsampler.hasFinalizedItems(), expectFinalizedItemsAfterSubmission);
        Assert.assertEquals(downsampler.peekFinalized() != null, expectFinalizedItemsAfterSubmission);

        // we expect pending items if there are any reads with assigned position in the reservoir.  The only way
        // to flush reads from the reservoir without invoking signalEndOfInput() or consumeFinalizedItems() is to encounter
        // a new read with assigned position at least one stride length ahead.  Since this new read enters the reservoir
        // and becomes a pending item, the only way not to have pending items at this point is if every read had no
        // assigned position
        final boolean expectPendingItemsAfterSubmission = !allReads.stream().allMatch(ReadUtils::readHasNoAssignedPosition);
        Assert.assertEquals(downsampler.hasPendingItems(), expectPendingItemsAfterSubmission);
        Assert.assertEquals(downsampler.peekPending() != null, expectPendingItemsAfterSubmission);

        downsampler.signalEndOfInput();
        Assert.assertEquals(downsampler.size(), expectedDownsampledReads);

        final boolean expectFinalizedItemsAfterEndOfInput = expectedDownsampledReads > 0;
        Assert.assertEquals(downsampler.hasFinalizedItems(), expectFinalizedItemsAfterEndOfInput);
        Assert.assertEquals(downsampler.peekFinalized() != null, expectFinalizedItemsAfterEndOfInput);

        Assert.assertFalse(downsampler.hasPendingItems());
        Assert.assertNull(downsampler.peekPending());

        List<GATKRead> downsampledReads = downsampler.consumeFinalizedItems();
        Assert.assertFalse(downsampler.hasFinalizedItems());
        Assert.assertFalse(downsampler.hasPendingItems());
        Assert.assertNull(downsampler.peekFinalized());
        Assert.assertNull(downsampler.peekPending());

        verifySortedness(downsampledReads);

        final Map<SimpleInterval, MutableInt> strideReadCounts = getStrideReadCounts(expectedCountsByStride.keySet(), downsampledReads);
        for (final SimpleInterval stride : expectedCountsByStride.keySet()) {
            Assert.assertEquals(strideReadCounts.get(stride).toInteger(), expectedCountsByStride.get(stride));
        }

        Assert.assertEquals(downsampler.getNumberOfDiscardedItems(), allReads.size() - downsampledReads.size());

        downsampler.resetStats();
        Assert.assertEquals(downsampler.getNumberOfDiscardedItems(), 0);

        // Now test with a PositionalDownsampler wrapped in an iterator, and make sure we get the same results.
        // It's crucial to reset the random number generator again in order to match the selections made by the
        // first downsampling pass.
        Utils.resetRandomGenerator();
        final ReadsDownsamplingIterator downsamplingIter = new ReadsDownsamplingIterator(allReads.iterator(), new PositionalDownsampler(targetCoverage, strideLength, biasToHigherMappingQuality, Integer.MAX_VALUE));
        final List<GATKRead> downsampledReadsFromIter = new ArrayList<>();
        for ( final GATKRead downsampledRead : downsamplingIter ) {
            downsampledReadsFromIter.add(downsampledRead);
        }

        Assert.assertEquals(downsampledReadsFromIter, downsampledReads, "Results from PositionalDownsampler wrapped in a ReadsDownsamplingIterator do not match results from standalone PositionalDownsampler");
    }

    @Test
    public void testSignalNoMoreReadsBefore() {
        final List<GATKRead> reads = createStackOfMappedReads(10, "1", 1);
        final ReadsDownsampler downsampler = new PositionalDownsampler(5);
        downsampler.submit(reads);

        Assert.assertTrue(downsampler.hasPendingItems());
        Assert.assertFalse(downsampler.hasFinalizedItems());

        final GATKRead nextRead = ArtificialReadUtils.createArtificialRead(header, "foo", "1", 2, new byte[]{'A'}, new byte[]{30});
        downsampler.signalNoMoreReadsBefore(nextRead);

        // Calling signalNoMoreReadsBefore() allows the downsampler to finalize the reads we gave it earlier than usual
        Assert.assertFalse(downsampler.hasPendingItems());
        Assert.assertTrue(downsampler.hasFinalizedItems());
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testPositionalDownsamplerZeroTargetCoverage() {
        final PositionalDownsampler downsampler = new PositionalDownsampler(0);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testPositionalDownsamplerNegativeTargetCoverage() {
        final PositionalDownsampler downsampler = new PositionalDownsampler(-1);
    }

    private List<GATKRead> createStackOfMappedReads( final int numReads, final String contig, final int startPosition ) {
        final List<GATKRead> reads = new ArrayList<>();

        for ( int i = 1; i <= numReads; ++i ) {
            reads.add(ArtificialReadUtils.createArtificialRead(header, "foo", contig, startPosition, new byte[]{'A'}, new byte[]{30}));
        }

        return reads;
    }

    private List<GATKRead> createStackOfUnmappedReads( final int numReads ) {
        final List<GATKRead> reads = new ArrayList<>();

        for ( int i = 1; i <= numReads; ++i ) {
            reads.add(ArtificialReadUtils.createArtificialUnmappedRead(header, new byte[]{'A'}, new byte[]{30}));
        }

        return reads;
    }

    private List<GATKRead> createStackOfUnmappedReadsWithPosition( final int numReads, final String contig, final int startPosition ) {
        final List<GATKRead> reads = new ArrayList<>();

        for ( int i = 1; i <= numReads; ++i ) {
            reads.add(ArtificialReadUtils.createArtificialUnmappedReadWithAssignedPosition(header, contig, startPosition, new byte[]{'A'}, new byte[]{30}));
        }

        return reads;
    }

    private Map<SimpleInterval, MutableInt> getStrideReadCounts(final Collection<SimpleInterval> strides, final List<GATKRead> downsampledReads) {
        final TargetCollection<SimpleInterval> intervals = new HashedListTargetCollection<>(strides.stream().collect(Collectors.toList()));
        final Map<SimpleInterval, MutableInt> result = strides.stream()
                .collect(Collectors.toMap(interval -> interval, interval -> new MutableInt(0)));

        for (final GATKRead read : downsampledReads) {
            if (ReadUtils.readHasNoAssignedPosition(read)) {
                result.get(DUMMY_INTERVAL_FOR_READS_WITH_NO_POSITION).increment();
            } else {
                final SimpleInterval interval = intervals.target(new SimpleInterval(read.getAssignedContig(), read.getAssignedStart(), read.getAssignedStart()));
                if (interval == null) {
                    Assert.fail("The given stride intervals do not cover one of the downsampled reads even though it has an assigned position.  Double-check the test code.");
                } else {
                    result.get(interval).increment();
                }
            }
        }
        return result;
    }

    private void verifySortedness(List<GATKRead> downsampledReads) {
        final List<GATKRead> readsWithAssignedPosition = downsampledReads.stream()
                .filter(read -> !ReadUtils.readHasNoAssignedPosition(read))
                .collect(Collectors.toList());
        for (int n = 0; n < readsWithAssignedPosition.size() - 1; n++) {
            final GATKRead read1 = readsWithAssignedPosition.get(n);
            final GATKRead read2 = readsWithAssignedPosition.get(n + 1);
            Assert.assertTrue(ReadCoordinateComparator.compareCoordinates(read1, read2, header) <= 0);
        }
    }
}
