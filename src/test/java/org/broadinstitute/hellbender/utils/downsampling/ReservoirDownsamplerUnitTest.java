/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.engine.downsampling;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.downsampling.ReadsDownsampler;
import org.broadinstitute.gatk.utils.downsampling.ReservoirDownsampler;
import org.broadinstitute.gatk.utils.sam.ArtificialSAMUtils;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import org.testng.Assert;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

public class ReservoirDownsamplerUnitTest extends BaseTest {

    private static class ReservoirDownsamplerTest extends TestDataProvider {
        int reservoirSize;
        int totalReads;
        int expectedNumReadsAfterDownsampling;
        int expectedNumDiscardedItems;

        public ReservoirDownsamplerTest( int reservoirSize, int totalReads ) {
            super(ReservoirDownsamplerTest.class);

            this.reservoirSize = reservoirSize;
            this.totalReads = totalReads;

            expectedNumReadsAfterDownsampling = Math.min(reservoirSize, totalReads);
            expectedNumDiscardedItems = totalReads <= reservoirSize ? 0 : totalReads - reservoirSize;

            setName(String.format("%s: reservoirSize=%d totalReads=%d expectedNumReadsAfterDownsampling=%d expectedNumDiscardedItems=%d",
                    getClass().getSimpleName(), reservoirSize, totalReads, expectedNumReadsAfterDownsampling, expectedNumDiscardedItems));
        }

        public Collection<SAMRecord> createReads() {
            Collection<SAMRecord> reads = new ArrayList<SAMRecord>(totalReads);

            SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000000);
            reads.addAll(ArtificialSAMUtils.createStackOfIdenticalArtificialReads(totalReads, header, "foo", 0, 1, 100));

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
    public void testReservoirDownsampler( ReservoirDownsamplerTest test ) {
        logger.warn("Running test: " + test);

        Utils.resetRandomGenerator();

        ReadsDownsampler<SAMRecord> downsampler = new ReservoirDownsampler<SAMRecord>(test.reservoirSize);

        downsampler.submit(test.createReads());

        if ( test.totalReads > 0 ) {
            Assert.assertTrue(downsampler.hasFinalizedItems());
            Assert.assertTrue(downsampler.peekFinalized() != null);
            Assert.assertFalse(downsampler.hasPendingItems());
            Assert.assertTrue(downsampler.peekPending() == null);
        }
        else {
            Assert.assertFalse(downsampler.hasFinalizedItems() || downsampler.hasPendingItems());
            Assert.assertTrue(downsampler.peekFinalized() == null && downsampler.peekPending() == null);
        }

        downsampler.signalEndOfInput();

        if ( test.totalReads > 0 ) {
            Assert.assertTrue(downsampler.hasFinalizedItems());
            Assert.assertTrue(downsampler.peekFinalized() != null);
            Assert.assertFalse(downsampler.hasPendingItems());
            Assert.assertTrue(downsampler.peekPending() == null);
        }
        else {
            Assert.assertFalse(downsampler.hasFinalizedItems() || downsampler.hasPendingItems());
            Assert.assertTrue(downsampler.peekFinalized() == null && downsampler.peekPending() == null);
        }

        Assert.assertEquals(downsampler.size(), test.expectedNumReadsAfterDownsampling);
        List<SAMRecord> downsampledReads = downsampler.consumeFinalizedItems();
        Assert.assertFalse(downsampler.hasFinalizedItems() || downsampler.hasPendingItems());
        Assert.assertTrue(downsampler.peekFinalized() == null && downsampler.peekPending() == null);

        Assert.assertEquals(downsampledReads.size(), test.expectedNumReadsAfterDownsampling);

        Assert.assertEquals(downsampler.getNumberOfDiscardedItems(), test.expectedNumDiscardedItems);
        Assert.assertEquals(test.totalReads - downsampledReads.size(), test.expectedNumDiscardedItems);

        downsampler.resetStats();
        Assert.assertEquals(downsampler.getNumberOfDiscardedItems(), 0);
    }
}
