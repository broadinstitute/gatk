package org.broadinstitute.hellbender.utils.downsampling;

import htsjdk.samtools.SAMFileHeader;
import org.apache.commons.math3.distribution.BinomialDistribution;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

public final class FractionalDownsamplerUnitTest extends GATKBaseTest {

    private static class FractionalDownsamplerTest extends TestDataProvider {
        double fraction;
        int totalReads;
        int expectedMinNumReadsAfterDownsampling;
        int expectedMaxNumReadsAfterDownsampling;
        int expectedMinDiscardedItems;
        int expectedMaxDiscardedItems;

        private static final double EXPECTED_ACCURACY = 0.05; // should be accurate to within +/- this percent

        public FractionalDownsamplerTest(double fraction, int totalReads) {
            super(FractionalDownsamplerTest.class);

            this.fraction = fraction;
            this.totalReads = totalReads;

            calculateExpectations();

            setName(String.format("%s: fraction=%.2f totalReads=%d expectedMinNumReadsAfterDownsampling=%d expectedMaxNumReadsAfterDownsampling=%d",
                    getClass().getSimpleName(), fraction, totalReads, expectedMinNumReadsAfterDownsampling, expectedMaxNumReadsAfterDownsampling));
        }

        private void calculateExpectations() {
            // Require an exact match in the 0% and 100% cases
            if (fraction == 0.0) {
                expectedMinNumReadsAfterDownsampling = expectedMaxNumReadsAfterDownsampling = 0;
                expectedMinDiscardedItems = expectedMaxDiscardedItems = totalReads;
            } else if (fraction == 1.0) {
                expectedMinNumReadsAfterDownsampling = expectedMaxNumReadsAfterDownsampling = totalReads;
                expectedMinDiscardedItems = expectedMaxDiscardedItems = 0;
            } else {
                expectedMinNumReadsAfterDownsampling = Math.max((int) ((fraction - EXPECTED_ACCURACY) * totalReads), 0);
                expectedMaxNumReadsAfterDownsampling = Math.min((int) ((fraction + EXPECTED_ACCURACY) * totalReads), totalReads);
                expectedMinDiscardedItems = totalReads - expectedMaxNumReadsAfterDownsampling;
                expectedMaxDiscardedItems = totalReads - expectedMinNumReadsAfterDownsampling;
            }
        }

        public Collection<GATKRead> createReads() {
            Collection<GATKRead> reads = new ArrayList<>(totalReads);

            SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(1, 1, 1000000);
            reads.addAll(ArtificialReadUtils.createIdenticalArtificialReads(totalReads, header, "foo", 0, 1, 100));

            return reads;
        }
    }

    @DataProvider(name = "FractionalDownsamplerTestDataProvider")
    public Object[][] createFractionalDownsamplerTestData() {
        for (double fraction : Arrays.asList(0.0, 0.25, 0.5, 0.75, 1.0)) {
            for (int totalReads : Arrays.asList(0, 1000, 10000)) {
                new FractionalDownsamplerTest(fraction, totalReads);
            }
        }

        return FractionalDownsamplerTest.getTests(FractionalDownsamplerTest.class);
    }

    @Test(dataProvider = "FractionalDownsamplerTestDataProvider")
    public void runFractionalDownsamplerTest(FractionalDownsamplerTest test) {
        logger.warn("Running test: " + test);

        Utils.resetRandomGenerator();

        final ReadsDownsampler downsampler = new FractionalDownsampler(test.fraction);

        downsampler.submit(test.createReads());

        if (test.totalReads > 0) {
            if (test.fraction > FractionalDownsamplerTest.EXPECTED_ACCURACY) {
                Assert.assertTrue(downsampler.hasFinalizedItems());
                Assert.assertTrue(downsampler.peekFinalized() != null);
            }
            Assert.assertFalse(downsampler.hasPendingItems());
            Assert.assertTrue(downsampler.peekPending() == null);
        } else {
            Assert.assertFalse(downsampler.hasFinalizedItems() || downsampler.hasPendingItems());
            Assert.assertTrue(downsampler.peekFinalized() == null && downsampler.peekPending() == null);
        }

        downsampler.signalEndOfInput();

        if (test.totalReads > 0) {
            if (test.fraction > FractionalDownsamplerTest.EXPECTED_ACCURACY) {
                Assert.assertTrue(downsampler.hasFinalizedItems());
                Assert.assertTrue(downsampler.peekFinalized() != null);
            }
            Assert.assertFalse(downsampler.hasPendingItems());
            Assert.assertTrue(downsampler.peekPending() == null);
        } else {
            Assert.assertFalse(downsampler.hasFinalizedItems() || downsampler.hasPendingItems());
            Assert.assertTrue(downsampler.peekFinalized() == null && downsampler.peekPending() == null);
        }

        List<GATKRead> downsampledReads = downsampler.consumeFinalizedItems();
        Assert.assertFalse(downsampler.hasFinalizedItems() || downsampler.hasPendingItems());
        Assert.assertTrue(downsampler.peekFinalized() == null && downsampler.peekPending() == null);

        Assert.assertTrue(downsampledReads.size() >= test.expectedMinNumReadsAfterDownsampling &&
                downsampledReads.size() <= test.expectedMaxNumReadsAfterDownsampling);

        Assert.assertTrue(downsampler.getNumberOfDiscardedItems() >= test.expectedMinDiscardedItems &&
                downsampler.getNumberOfDiscardedItems() <= test.expectedMaxDiscardedItems);

        Assert.assertEquals(downsampler.getNumberOfDiscardedItems(), test.totalReads - downsampledReads.size());

        downsampler.resetStats();
        Assert.assertEquals(downsampler.getNumberOfDiscardedItems(), 0);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testBadArguments_negative() throws Exception {
        new FractionalDownsampler(-0.1);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testBadArguments_largerThan1() throws Exception {
        new FractionalDownsampler(1.1);
    }

    @Test
    public void testClear() throws Exception {
        final double f = 0.01;
        final int N = 100000;
        final ReadsDownsampler fd = new FractionalDownsampler(f);
        for (int i = 0; i < N; i++) {
            fd.submit(ArtificialReadUtils.createArtificialRead("100M"));
        }
        final BinomialDistribution bin = new BinomialDistribution(N, f);
        final double errorRate = 1.0 / 1000;   //we can fail this often
        Assert.assertTrue(fd.size() <= bin.inverseCumulativeProbability(1 - errorRate / 2));
        Assert.assertTrue(fd.size() >= bin.inverseCumulativeProbability(errorRate / 2));

        Assert.assertTrue(fd.hasFinalizedItems());
        Assert.assertFalse(fd.hasPendingItems());
        Assert.assertFalse(fd.requiresCoordinateSortOrder());
        Assert.assertNotNull(fd.peekFinalized());
        Assert.assertNull(fd.peekPending());
        fd.clearItems();
        Assert.assertEquals(fd.size(), 0);
        Assert.assertFalse(fd.hasFinalizedItems());
        Assert.assertFalse(fd.hasPendingItems());
        Assert.assertFalse(fd.requiresCoordinateSortOrder());
        Assert.assertNull(fd.peekFinalized());
        Assert.assertNull(fd.peekPending());

    }

    @Test
    public void testSignalNoMoreReadsBefore() throws Exception {
        FractionalDownsampler rd = new FractionalDownsampler(0.1);
        final GATKRead r1 = ArtificialReadUtils.createArtificialRead("100M");
        final GATKRead r2 = ArtificialReadUtils.createArtificialRead("101M");
        rd.submit(r1);
        rd.signalNoMoreReadsBefore(r2);//no op
        rd.submit(r2);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNoNullSignalNoMoreReadsBefore() throws Exception {
        ReadsDownsampler rd = new FractionalDownsampler(0.1);
        rd.signalNoMoreReadsBefore(null);
    }
}