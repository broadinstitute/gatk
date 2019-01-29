package org.broadinstitute.hellbender.utils.downsampling;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

public final class ReservoirDownsamplerUnitTest extends GATKBaseTest {

    private static class ReservoirDownsamplerTest extends TestDataProvider {
        int reservoirSize;
        int totalReads;
        int expectedNumReadsAfterDownsampling;
        int expectedNumDiscardedItems;

        public ReservoirDownsamplerTest(final int reservoirSize, final int totalReads ) {
            super(ReservoirDownsamplerTest.class);

            this.reservoirSize = reservoirSize;
            this.totalReads = totalReads;

            expectedNumReadsAfterDownsampling = Math.min(reservoirSize, totalReads);
            expectedNumDiscardedItems = totalReads <= reservoirSize ? 0 : totalReads - reservoirSize;

            setName(String.format("%s: reservoirSize=%d totalReads=%d expectedNumReadsAfterDownsampling=%d expectedNumDiscardedItems=%d",
                    getClass().getSimpleName(), reservoirSize, totalReads, expectedNumReadsAfterDownsampling, expectedNumDiscardedItems));
        }

        public Collection<GATKRead> createReads() {
            final Collection<GATKRead> reads = new ArrayList<>(totalReads);

            final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(1, 1, 1000000);
            reads.addAll(ArtificialReadUtils.createIdenticalArtificialReads(totalReads, header, "foo", 0, 1, 100));

            return reads;
        }
    }

    @DataProvider(name = "ReservoirDownsamplerTestDataProvider")
    public Object[][] createReservoirDownsamplerTestData() {
        for ( int reservoirSize = 1; reservoirSize <= 10000; reservoirSize *= 10 ) {
            new ReservoirDownsamplerTest(reservoirSize, 0);
            for ( int totalReads = 1; totalReads <= 10000; totalReads *= 10 ) {
                new ReservoirDownsamplerTest(reservoirSize, totalReads);
            }
        }

        return ReservoirDownsamplerTest.getTests(ReservoirDownsamplerTest.class);
    }

    @Test(dataProvider = "ReservoirDownsamplerTestDataProvider")
    public void testReservoirDownsampler(final ReservoirDownsamplerTest test ) {
        logger.warn("Running test: " + test);

        Utils.resetRandomGenerator();
        final ReadsDownsampler downsampler = new ReservoirDownsampler(test.reservoirSize);
        downsampler.submit(test.createReads());

        // after submit, but before signalEndOfInput, all reads are pending, none are finalized
        if ( test.totalReads > 0 ) {
            Assert.assertFalse(downsampler.hasFinalizedItems());
            Assert.assertNull(downsampler.peekFinalized());
            Assert.assertTrue(downsampler.hasPendingItems());
            Assert.assertNotNull(downsampler.peekPending());
        }
        else {
            Assert.assertFalse(downsampler.hasFinalizedItems() || downsampler.hasPendingItems());
            Assert.assertTrue(downsampler.peekFinalized() == null && downsampler.peekPending() == null);
        }

        // after signalEndOfInput, no reads are pending, all are finalized
        downsampler.signalEndOfInput();

        if ( test.totalReads > 0 ) {
            Assert.assertTrue(downsampler.hasFinalizedItems());
            Assert.assertNotNull(downsampler.peekFinalized());
            Assert.assertFalse(downsampler.hasPendingItems());
            Assert.assertNull(downsampler.peekPending());
        }
        else {
            Assert.assertFalse(downsampler.hasFinalizedItems() || downsampler.hasPendingItems());
            Assert.assertTrue(downsampler.peekFinalized() == null && downsampler.peekPending() == null);
        }

        Assert.assertEquals(downsampler.size(), test.expectedNumReadsAfterDownsampling);
        final List<GATKRead> downsampledReads = downsampler.consumeFinalizedItems();
        Assert.assertFalse(downsampler.hasFinalizedItems() || downsampler.hasPendingItems());
        Assert.assertTrue(downsampler.peekFinalized() == null && downsampler.peekPending() == null);

        Assert.assertEquals(downsampledReads.size(), test.expectedNumReadsAfterDownsampling);

        Assert.assertEquals(downsampler.getNumberOfDiscardedItems(), test.expectedNumDiscardedItems);
        Assert.assertEquals(test.totalReads - downsampledReads.size(), test.expectedNumDiscardedItems);

        downsampler.resetStats();
        Assert.assertEquals(downsampler.getNumberOfDiscardedItems(), 0);

        // use the same downsampling parameters, but this time consume the reads through a
        // ReadsDownsamplingIterator, and validate that we get the same results as using downsampler directly
        Utils.resetRandomGenerator();
        final ReadsDownsampler downsamplerForIterator = new ReservoirDownsampler(test.reservoirSize);
        final ReadsDownsamplingIterator downsamplingIterator = new ReadsDownsamplingIterator(
                test.createReads().iterator(),
                downsamplerForIterator);
        final List<GATKRead> downsampledReadsFromIterator = new ArrayList<>(test.reservoirSize);
        downsamplingIterator.forEach(downsampledReadsFromIterator::add);

        Assert.assertEquals(downsamplerForIterator.getNumberOfDiscardedItems(), test.expectedNumDiscardedItems);
        Assert.assertEquals(downsampledReadsFromIterator, downsampledReads);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInvalidCtorArg() throws Exception {
        new ReservoirDownsampler(0);
    }

    @Test
    public void testRequiredOrder() throws Exception {
        Assert.assertFalse(new ReservoirDownsampler(1).requiresCoordinateSortOrder());
    }

    @Test
    public void testSignalNoMoreReadsBefore() throws Exception {
        ReservoirDownsampler rd = new ReservoirDownsampler(1, true);
        final GATKRead r1= ArtificialReadUtils.createArtificialRead("100M");
        final GATKRead r2= ArtificialReadUtils.createArtificialRead("101M");
        rd.submit(r1);
        rd.signalNoMoreReadsBefore(r2);//no op
        rd.submit(r2);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNoNullSignalNoMoreReadsBefore() throws Exception {
        ReadsDownsampler rd = new ReservoirDownsampler(1, true);
        rd.signalNoMoreReadsBefore(null);
    }

    // Calling consumeFinalizeItems() before end of input on a non-empty downsampler should
    // return an empty List and not change the state of the downsampler:
    @Test
    public void testConsumeFinalizedItemsBeforeEndOfInput() {
        final ReservoirDownsampler downsampler = new ReservoirDownsampler(10, true);
        final GATKRead read = ArtificialReadUtils.createArtificialRead("100M");
        downsampler.submit(read);

        Assert.assertFalse(downsampler.hasFinalizedItems());
        Assert.assertTrue(downsampler.hasPendingItems());

        final List<GATKRead> returnedReads = downsampler.consumeFinalizedItems();
        Assert.assertTrue(returnedReads.isEmpty());

        // The downsampler should still have pending reads after the call to consumeFinalizedItems()
        // (ie., the downsampling process should still be ongoing, and the state of the downsampler
        // should have been preserved):
        Assert.assertFalse(downsampler.hasFinalizedItems());
        Assert.assertTrue(downsampler.hasPendingItems());
    }

    // Calling consumeFinalizedItems() after end of input on an empty downsampler should
    // return an empty List and reset the end of input flag in the downsampler:
    @Test
    public void testConsumeFinalizedItemsAfterEndOfInputOnEmptyDownsampler() {
        final ReservoirDownsampler downsampler = new ReservoirDownsampler(10, true);

        Assert.assertFalse(downsampler.hasFinalizedItems());
        Assert.assertFalse(downsampler.hasPendingItems());

        downsampler.signalEndOfInput();

        final List<GATKRead> returnedReads = downsampler.consumeFinalizedItems();
        Assert.assertTrue(returnedReads.isEmpty());

        Assert.assertFalse(downsampler.hasFinalizedItems());
        Assert.assertFalse(downsampler.hasPendingItems());

        // Submitting an additional read after both signalEndOfInput() and consumeFinalizedItems()
        // should succeed (not throw), since the call to consumeFinalizedItems() should reset the
        // end of input flag:
        final GATKRead read = ArtificialReadUtils.createArtificialRead("100M");
        downsampler.submit(read);

        Assert.assertFalse(downsampler.hasFinalizedItems());
        Assert.assertTrue(downsampler.hasPendingItems());
    }

    // Calling submit() directly after signalEndOfInput() should throw:
    @Test(expectedExceptions = IllegalStateException.class)
    public void testSubmitAfterEndOfInputThrows() {
        final ReservoirDownsampler downsampler = new ReservoirDownsampler(10, true);
        downsampler.signalEndOfInput();
        final GATKRead read = ArtificialReadUtils.createArtificialRead("100M");
        downsampler.submit(read);
    }
}
