package org.broadinstitute.hellbender.utils.downsampling;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadCoordinateComparator;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

public class PositionalDownsamplerUnitTest extends GATKBaseTest {

    private final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();

    @DataProvider(name = "PositionalDownsamplerTestData")
    public Object[][] positionalDownsamplerTestData() {
        return new Object[][] {
                // reads, downsampling target coverage, expected size and number of position-based stacks after downsampling, and whether the stream ends with a stack of pure unmapped reads
                { Arrays.asList(createStackOfMappedReads(10, "1", 1)), 1, Arrays.asList(1), false },
                { Arrays.asList(createStackOfMappedReads(10, "1", 1)), 5, Arrays.asList(5), false },
                { Arrays.asList(createStackOfMappedReads(10, "1", 1)), 10, Arrays.asList(10), false },
                { Arrays.asList(createStackOfMappedReads(10, "1", 1)), 15, Arrays.asList(10), false },
                { Arrays.asList(createStackOfMappedReads(10, "1", 1), createStackOfMappedReads(8, "1", 2)), 5, Arrays.asList(5, 5), false },
                { Arrays.asList(createStackOfMappedReads(10, "1", 1), createStackOfMappedReads(5, "1", 2)), 5, Arrays.asList(5, 5), false },
                { Arrays.asList(createStackOfMappedReads(10, "1", 1), createStackOfMappedReads(3, "1", 2)), 5, Arrays.asList(5, 3), false },
                { Arrays.asList(createStackOfMappedReads(10, "1", 1), createStackOfMappedReads(10, "2", 1)), 5, Arrays.asList(5, 5), false },
                { Arrays.asList(createStackOfMappedReads(10, "1", 1), createStackOfMappedReads(10, "2", 1), createStackOfMappedReads(10, "2", 2)), 5, Arrays.asList(5, 5, 5), false },
                { Arrays.asList(createStackOfMappedReads(10, "1", 1), createStackOfMappedReads(10, "2", 1), createStackOfMappedReads(10, "2", 2)), 10, Arrays.asList(10, 10, 10), false },
                { Arrays.asList(createStackOfMappedReads(3, "1", 1), createStackOfMappedReads(10, "2", 1), createStackOfMappedReads(4, "2", 2)), 5, Arrays.asList(3, 5, 4), false },
                { Arrays.asList(createStackOfMappedReads(3, "1", 1), createStackOfMappedReads(10, "2", 1), createStackOfMappedReads(4, "2", 2)), 1, Arrays.asList(1, 1, 1), false },
                // Unmapped reads with no positions should not get downsampled
                { Arrays.asList(createStackOfMappedReads(10, "1", 1), createStackOfMappedReads(10, "2", 1), createStackOfMappedReads(10, "2", 2), createStackOfUnmappedReads(100)), 5, Arrays.asList(5, 5, 5, 100), true },
                // Unmapped reads with assigned positions should get downsampled
                { Arrays.asList(createStackOfMappedReads(10, "1", 1), createStackOfUnmappedReadsWithPosition(10, "1", 1)), 5, Arrays.asList(5), false },
                { Arrays.asList(createStackOfMappedReads(10, "1", 1), createStackOfUnmappedReadsWithPosition(10, "1", 1), createStackOfMappedReads(5, "2", 1), createStackOfUnmappedReadsWithPosition(5, "2", 1)), 5, Arrays.asList(5, 5), false },
                // Mix of mapped, pure unmapped, and mapped with assigned position
                { Arrays.asList(createStackOfMappedReads(10, "1", 1), createStackOfUnmappedReadsWithPosition(10, "1", 1), createStackOfMappedReads(5, "2", 1), createStackOfUnmappedReadsWithPosition(5, "2", 1), createStackOfUnmappedReads(100)), 5, Arrays.asList(5, 5, 100), true },
                // Unmapped only
                { Arrays.asList(createStackOfUnmappedReads(10)), 1, Arrays.asList(10), true },
                // No reads
                { Collections.emptyList(), 1, Collections.emptyList(), false }
        };
    }

    @Test(dataProvider = "PositionalDownsamplerTestData")
    public void testPositionalDownsampler( final List<List<GATKRead>> reads, final int targetCoverage, final List<Integer> expectedStackSizes, final boolean lastStackIsPureUnmapped ) {
        final List<GATKRead> allReads = new ArrayList<>();
        for ( final List<GATKRead> readStack : reads ) {
            allReads.addAll(readStack);
        }
        final int expectedStackSizesSum = expectedStackSizes.stream().mapToInt(i -> i).sum();

        final ReadsDownsampler downsampler = new PositionalDownsampler(targetCoverage, header);
        Assert.assertTrue(downsampler.requiresCoordinateSortOrder());
        Assert.assertEquals(downsampler.size(), 0);

        Utils.resetRandomGenerator();
        downsampler.submit(allReads);

        Assert.assertEquals(downsampler.size(), expectedStackSizesSum);

        final boolean expectFinalizedItemsAfterSubmission = expectedStackSizes.size() > 1 || (expectedStackSizes.size() == 1 && lastStackIsPureUnmapped);
        Assert.assertEquals(downsampler.hasFinalizedItems(), expectFinalizedItemsAfterSubmission);
        if ( expectFinalizedItemsAfterSubmission ) {
            Assert.assertNotNull(downsampler.peekFinalized());
        }
        else {
            Assert.assertNull(downsampler.peekFinalized());
        }

        final boolean expectPendingItemsAfterSubmission = ! expectedStackSizes.isEmpty() && ! lastStackIsPureUnmapped;
        Assert.assertEquals(downsampler.hasPendingItems(), expectPendingItemsAfterSubmission);
        if ( expectPendingItemsAfterSubmission ) {
            Assert.assertNotNull(downsampler.peekPending());
        }
        else {
            Assert.assertNull(downsampler.peekPending());
        }

        downsampler.signalEndOfInput();

        Assert.assertEquals(downsampler.size(), expectedStackSizesSum);

        final boolean expectFinalizedItemsAfterEndOfInput = ! expectedStackSizes.isEmpty();
        Assert.assertEquals(downsampler.hasFinalizedItems(), expectFinalizedItemsAfterEndOfInput);
        if ( expectFinalizedItemsAfterEndOfInput ) {
            Assert.assertNotNull(downsampler.peekFinalized());
        }
        else {
            Assert.assertNull(downsampler.peekFinalized());
        }

        Assert.assertFalse(downsampler.hasPendingItems());
        Assert.assertNull(downsampler.peekPending());

        List<GATKRead> downsampledReads = downsampler.consumeFinalizedItems();
        Assert.assertFalse(downsampler.hasFinalizedItems());
        Assert.assertFalse(downsampler.hasPendingItems());
        Assert.assertNull(downsampler.peekFinalized());
        Assert.assertNull(downsampler.peekPending());

        if ( expectedStackSizes.size() == 0 ) {
            Assert.assertTrue(downsampledReads.isEmpty());
        }
        else {
            List<Integer> downsampledStackSizes = getDownsampledStackSizesAndVerifySortedness(downsampledReads);

            Assert.assertEquals(downsampledStackSizes.size(), expectedStackSizes.size());
            Assert.assertEquals(downsampledStackSizes, expectedStackSizes);

            int numReadsActuallyEliminated = allReads.size() - downsampledReads.size();
            int numReadsReportedEliminated = downsampler.getNumberOfDiscardedItems();
            Assert.assertEquals(numReadsActuallyEliminated, numReadsReportedEliminated);
        }

        downsampler.resetStats();
        Assert.assertEquals(downsampler.getNumberOfDiscardedItems(), 0);

        // Now test with a PositionalDownsampler wrapped in an iterator, and make sure we get the same results.
        // It's crucial to reset the random number generator again in order to match the selections made by the
        // first downsampling pass.
        Utils.resetRandomGenerator();
        final ReadsDownsamplingIterator downsamplingIter = new ReadsDownsamplingIterator(allReads.iterator(), new PositionalDownsampler(targetCoverage, header));
        final List<GATKRead> downsampledReadsFromIter = new ArrayList<>();
        for ( final GATKRead downsampledRead : downsamplingIter ) {
            downsampledReadsFromIter.add(downsampledRead);
        }

        Assert.assertEquals(downsampledReadsFromIter, downsampledReads, "Results from PositionalDownsampler wrapped in a ReadsDownsamplingIterator do not match results from standalone PositionalDownsampler");
    }

    @Test
    public void testSignalNoMoreReadsBefore() {
        final List<GATKRead> reads = createStackOfMappedReads(10, "1", 1);
        final ReadsDownsampler downsampler = new PositionalDownsampler(5, header);
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
        final PositionalDownsampler downsampler = new PositionalDownsampler(0, header);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testPositionalDownsamplerNegativeTargetCoverage() {
        final PositionalDownsampler downsampler = new PositionalDownsampler(-1, header);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testPositionalDownsamplerInvalidheader() {
        final PositionalDownsampler downsampler = new PositionalDownsampler(1, null);
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

    private List<Integer> getDownsampledStackSizesAndVerifySortedness( final List<GATKRead> downsampledReads ) {
        List<Integer> stackSizes = new ArrayList<>();

        if ( downsampledReads.isEmpty() ) {
            return stackSizes;
        }

        Iterator<GATKRead> iter = downsampledReads.iterator();
        Assert.assertTrue(iter.hasNext());

        GATKRead previousRead = iter.next();
        int currentStackSize = 1;

        while ( iter.hasNext() ) {
            GATKRead currentRead = iter.next();

            final int positionComparison = ReadCoordinateComparator.compareCoordinates(previousRead, currentRead, header);

            if ( positionComparison < 0  ) {
                stackSizes.add(currentStackSize);
                currentStackSize = 1;
            }
            else if ( positionComparison > 0 ) {
                Assert.fail(String.format("Reads are out of order: %s %s", previousRead, currentRead));
            }
            else {
                ++currentStackSize;
            }

            previousRead = currentRead;
        }

        stackSizes.add(currentStackSize);
        return stackSizes;
    }
}
