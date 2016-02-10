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
import org.broadinstitute.gatk.utils.downsampling.SimplePositionalDownsampler;
import org.broadinstitute.gatk.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import org.testng.Assert;

import java.util.*;

public class SimplePositionalDownsamplerUnitTest extends BaseTest {

    private static class SimplePositionalDownsamplerTest extends TestDataProvider {
        int targetCoverage;
        int numStacks;
        List<Integer> stackSizes;
        List<Integer> expectedStackSizes;
        boolean multipleContigs;
        int totalInitialReads;

        public SimplePositionalDownsamplerTest( int targetCoverage, List<Integer> stackSizes, boolean multipleContigs ) {
            super(SimplePositionalDownsamplerTest.class);

            this.targetCoverage = targetCoverage;
            this.numStacks = stackSizes.size();
            this.stackSizes = stackSizes;
            this.multipleContigs = multipleContigs;

            calculateExpectedDownsampledStackSizes();

            totalInitialReads = 0;
            for ( Integer stackSize : stackSizes ) {
                totalInitialReads += stackSize;
            }

            setName(String.format("%s: targetCoverage=%d numStacks=%d stackSizes=%s expectedSizes=%s multipleContigs=%b",
                    getClass().getSimpleName(), targetCoverage, numStacks, stackSizes, expectedStackSizes, multipleContigs));
        }

        public Collection<SAMRecord> createReads() {
            Collection<SAMRecord> reads = new ArrayList<SAMRecord>();
            SAMFileHeader header = multipleContigs ?
                                   ArtificialSAMUtils.createArtificialSamHeader(2, 1, 1000000) :
                                   ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000000);

            int refIndex = 0;
            int alignmentStart = 1;
            int readLength = 100;

            for ( int i = 0; i < numStacks; i++ ) {
                if ( multipleContigs && refIndex == 0 && i >= numStacks / 2 ) {
                    refIndex++;
                }

                reads.addAll(ArtificialSAMUtils.createStackOfIdenticalArtificialReads(stackSizes.get(i), header, "foo",
                                                                                      refIndex, alignmentStart, readLength));

                alignmentStart += 10;
            }

            return reads;
        }

        private void calculateExpectedDownsampledStackSizes() {
            expectedStackSizes = new ArrayList<Integer>(numStacks);

            for ( Integer stackSize : stackSizes ) {
                int expectedSize = targetCoverage >= stackSize ? stackSize : targetCoverage;
                expectedStackSizes.add(expectedSize);
            }
        }
    }

    @DataProvider(name = "SimplePositionalDownsamplerTestDataProvider")
    public Object[][] createSimplePositionalDownsamplerTestData() {
        Utils.resetRandomGenerator();

        for ( int targetCoverage = 1; targetCoverage <= 10000; targetCoverage *= 10 ) {
            for ( int contigs = 1; contigs <= 2; contigs++ ) {
                for ( int numStacks = 0; numStacks <= 10; numStacks++ ) {
                    List<Integer> stackSizes = new ArrayList<Integer>(numStacks);
                    for ( int stack = 1; stack <= numStacks; stack++ ) {
                        stackSizes.add(Utils.getRandomGenerator().nextInt(targetCoverage * 2) + 1);
                    }
                    new SimplePositionalDownsamplerTest(targetCoverage, stackSizes, contigs > 1);
                }
            }
        }

        return SimplePositionalDownsamplerTest.getTests(SimplePositionalDownsamplerTest.class);
    }

    @Test( dataProvider = "SimplePositionalDownsamplerTestDataProvider" )
    public void testSimplePostionalDownsampler( SimplePositionalDownsamplerTest test ) {
        logger.warn("Running test: " + test);

        Utils.resetRandomGenerator();

        ReadsDownsampler<SAMRecord> downsampler = new SimplePositionalDownsampler<SAMRecord>(test.targetCoverage);

        downsampler.submit(test.createReads());

        if ( test.numStacks > 1 ) {
            Assert.assertTrue(downsampler.hasFinalizedItems());
            Assert.assertTrue(downsampler.peekFinalized() != null);
            Assert.assertTrue(downsampler.hasPendingItems());
            Assert.assertTrue(downsampler.peekPending() != null);
        }
        else if ( test.numStacks == 1 ) {
            Assert.assertFalse(downsampler.hasFinalizedItems());
            Assert.assertTrue(downsampler.peekFinalized() == null);
            Assert.assertTrue(downsampler.hasPendingItems());
            Assert.assertTrue(downsampler.peekPending() != null);
        }
        else {
            Assert.assertFalse(downsampler.hasFinalizedItems() || downsampler.hasPendingItems());
            Assert.assertTrue(downsampler.peekFinalized() == null && downsampler.peekPending() == null);
        }

        downsampler.signalEndOfInput();

        if ( test.numStacks > 0 ) {
            Assert.assertTrue(downsampler.hasFinalizedItems());
            Assert.assertTrue(downsampler.peekFinalized() != null);
            Assert.assertFalse(downsampler.hasPendingItems());
            Assert.assertTrue(downsampler.peekPending() == null);
        }
        else {
            Assert.assertFalse(downsampler.hasFinalizedItems() || downsampler.hasPendingItems());
            Assert.assertTrue(downsampler.peekFinalized() == null && downsampler.peekPending() == null);
        }

        List<SAMRecord> downsampledReads = downsampler.consumeFinalizedItems();
        Assert.assertFalse(downsampler.hasFinalizedItems() || downsampler.hasPendingItems());
        Assert.assertTrue(downsampler.peekFinalized() == null && downsampler.peekPending() == null);

        if ( test.numStacks == 0 ) {
            Assert.assertTrue(downsampledReads.isEmpty());
        }
        else {
            List<Integer> downsampledStackSizes = getDownsampledStackSizesAndVerifySortedness(downsampledReads);

            Assert.assertEquals(downsampledStackSizes.size(), test.numStacks);
            Assert.assertEquals(downsampledStackSizes, test.expectedStackSizes);

            int numReadsActuallyEliminated = test.totalInitialReads - downsampledReads.size();
            int numReadsReportedEliminated = downsampler.getNumberOfDiscardedItems();
            Assert.assertEquals(numReadsActuallyEliminated, numReadsReportedEliminated);
        }

        downsampler.resetStats();
        Assert.assertEquals(downsampler.getNumberOfDiscardedItems(), 0);
    }

    private List<Integer> getDownsampledStackSizesAndVerifySortedness( List<SAMRecord> downsampledReads ) {
        List<Integer> stackSizes = new ArrayList<Integer>();

        if ( downsampledReads.isEmpty() ) {
            return stackSizes;
        }

        Iterator<SAMRecord> iter = downsampledReads.iterator();
        Assert.assertTrue(iter.hasNext());

        SAMRecord previousRead = iter.next();
        int currentStackSize = 1;

        while ( iter.hasNext() ) {
            SAMRecord currentRead = iter.next();

            if ( currentRead.getReferenceIndex() > previousRead.getReferenceIndex() || currentRead.getAlignmentStart() > previousRead.getAlignmentStart() ) {
                stackSizes.add(currentStackSize);
                currentStackSize = 1;
            }
            else if ( currentRead.getReferenceIndex() < previousRead.getReferenceIndex() || currentRead.getAlignmentStart() < previousRead.getAlignmentStart() ) {
                Assert.fail(String.format("Reads are out of order: %s %s", previousRead, currentRead));
            }
            else {
                currentStackSize++;
            }

            previousRead = currentRead;
        }

        stackSizes.add(currentStackSize);
        return stackSizes;
    }

    @Test
    public void testSimplePositionalDownsamplerSignalNoMoreReadsBefore() {
        ReadsDownsampler<SAMRecord> downsampler = new SimplePositionalDownsampler<SAMRecord>(1000);

        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000000);

        Collection<SAMRecord> readStack = new ArrayList<SAMRecord>();
        readStack.addAll(ArtificialSAMUtils.createStackOfIdenticalArtificialReads(50, header, "foo", 0, 1, 100));
        downsampler.submit(readStack);

        Assert.assertFalse(downsampler.hasFinalizedItems());
        Assert.assertTrue(downsampler.peekFinalized() == null);
        Assert.assertTrue(downsampler.hasPendingItems());
        Assert.assertTrue(downsampler.peekPending() != null);

        SAMRecord laterRead = ArtificialSAMUtils.createArtificialRead(header, "foo", 0, 2, 100);
        downsampler.signalNoMoreReadsBefore(laterRead);

        Assert.assertTrue(downsampler.hasFinalizedItems());
        Assert.assertTrue(downsampler.peekFinalized() != null);
        Assert.assertFalse(downsampler.hasPendingItems());
        Assert.assertTrue(downsampler.peekPending() == null);

        List<SAMRecord> downsampledReads = downsampler.consumeFinalizedItems();

        Assert.assertEquals(downsampledReads.size(), readStack.size());
    }

    @Test
    public void testBasicUnmappedReadsSupport() {
        ReadsDownsampler<SAMRecord> downsampler = new SimplePositionalDownsampler<SAMRecord>(100);

        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000000);

        Collection<SAMRecord> readStack = new ArrayList<SAMRecord>();
        readStack.addAll(ArtificialSAMUtils.createStackOfIdenticalArtificialReads(200, header, "foo", SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX,
                                                                                  SAMRecord.NO_ALIGNMENT_START, 100));
        for ( SAMRecord read : readStack ) {
            Assert.assertTrue(read.getReadUnmappedFlag());
        }

        downsampler.submit(readStack);
        downsampler.signalEndOfInput();

        List<SAMRecord> downsampledReads = downsampler.consumeFinalizedItems();

        // Unmapped reads should not get downsampled at all by the SimplePositionalDownsampler
        Assert.assertEquals(downsampledReads.size(), readStack.size());

        for ( SAMRecord read: downsampledReads ) {
            Assert.assertTrue(read.getReadUnmappedFlag());
        }
    }

    @Test
    public void testMixedMappedAndUnmappedReadsSupport() {
        ReadsDownsampler<SAMRecord> downsampler = new SimplePositionalDownsampler<SAMRecord>(100);

        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000000);

        Collection<SAMRecord> mappedReadStack = new ArrayList<SAMRecord>();
        mappedReadStack.addAll(ArtificialSAMUtils.createStackOfIdenticalArtificialReads(200, header, "foo", 0, 1, 100));
        for ( SAMRecord read : mappedReadStack ) {
            Assert.assertFalse(read.getReadUnmappedFlag());
        }

        Collection<SAMRecord> unmappedReadStack = new ArrayList<SAMRecord>();
        unmappedReadStack.addAll(ArtificialSAMUtils.createStackOfIdenticalArtificialReads(200, header, "foo", SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX,
                                                                                          SAMRecord.NO_ALIGNMENT_START, 100));
        for ( SAMRecord read : unmappedReadStack ) {
            Assert.assertTrue(read.getReadUnmappedFlag());
        }

        downsampler.submit(mappedReadStack);
        downsampler.submit(unmappedReadStack);
        downsampler.signalEndOfInput();

        List<SAMRecord> downsampledReads = downsampler.consumeFinalizedItems();

        // Unmapped reads should not get downsampled at all by the SimplePositionalDownsampler
        Assert.assertEquals(downsampledReads.size(), 300);
        Assert.assertEquals(downsampler.getNumberOfDiscardedItems(), 100);

        int count = 1;
        for ( SAMRecord read: downsampledReads ) {
            if ( count <= 100 ) {
                Assert.assertFalse(read.getReadUnmappedFlag());
            }
            else {
                Assert.assertTrue(read.getReadUnmappedFlag());
            }

            count++;
        }
    }

    @Test
    public void testGATKSAMRecordSupport() {
        ReadsDownsampler<GATKSAMRecord> downsampler = new SimplePositionalDownsampler<GATKSAMRecord>(1000);

        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000000);

        List<GATKSAMRecord> reads = new ArrayList<GATKSAMRecord>();
        for ( int i = 0; i < 10; i++ ) {
            reads.add(ArtificialSAMUtils.createArtificialRead(header, "foo", 0, 10, 20 * i + 10));
        }

        downsampler.submit(reads);
        downsampler.signalEndOfInput();
        List<GATKSAMRecord> downsampledReads = downsampler.consumeFinalizedItems();

        Assert.assertEquals(downsampledReads.size(), 10);
    }
}
