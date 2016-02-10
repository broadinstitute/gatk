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

package org.broadinstitute.gatk.utils.locusiterator;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.downsampling.DownsampleType;
import org.broadinstitute.gatk.utils.downsampling.DownsamplingMethod;
import org.broadinstitute.gatk.utils.NGSPlatform;
import org.broadinstitute.gatk.utils.QualityUtils;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.pileup.PileupElement;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileup;
import org.broadinstitute.gatk.utils.sam.ArtificialBAMBuilder;
import org.broadinstitute.gatk.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.gatk.utils.sam.GATKSAMReadGroupRecord;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

/**
 * testing of the new (non-legacy) version of LocusIteratorByState
 */
public class LocusIteratorByStateUnitTest extends LocusIteratorByStateBaseTest {
    private static final boolean DEBUG = false;
    protected LocusIteratorByState li;

    @Test(enabled = !DEBUG)
    public void testUnmappedAndAllIReadsPassThrough() {
        final int readLength = 10;
        GATKSAMRecord mapped1 = ArtificialSAMUtils.createArtificialRead(header,"mapped1",0,1,readLength);
        GATKSAMRecord mapped2 = ArtificialSAMUtils.createArtificialRead(header,"mapped2",0,1,readLength);
        GATKSAMRecord unmapped = ArtificialSAMUtils.createArtificialRead(header,"unmapped",0,1,readLength);
        GATKSAMRecord allI = ArtificialSAMUtils.createArtificialRead(header,"allI",0,1,readLength);

        unmapped.setReadUnmappedFlag(true);
        unmapped.setCigarString("*");
        allI.setCigarString(readLength + "I");

        List<GATKSAMRecord> reads = Arrays.asList(mapped1, unmapped, allI, mapped2);

        // create the iterator by state with the fake reads and fake records
        li = makeLTBS(reads, DownsamplingMethod.NONE, true);

        Assert.assertTrue(li.hasNext());
        AlignmentContext context = li.next();
        ReadBackedPileup pileup = context.getBasePileup();
        Assert.assertEquals(pileup.depthOfCoverage(), 2, "Should see only 2 reads in pileup, even with unmapped and all I reads");

        final List<GATKSAMRecord> rawReads = li.transferReadsFromAllPreviousPileups();
        Assert.assertEquals(rawReads, reads, "Input and transferred read lists should be the same, and include the unmapped and all I reads");
    }

    @Test(enabled = true && ! DEBUG)
    public void testXandEQOperators() {
        final byte[] bases1 = new byte[] {'A','A','A','A','A','A','A','A','A','A'};
        final byte[] bases2 = new byte[] {'A','A','A','C','A','A','A','A','A','C'};

        GATKSAMRecord r1 = ArtificialSAMUtils.createArtificialRead(header,"r1",0,1,10);
        r1.setReadBases(bases1);
        r1.setBaseQualities(new byte[] {20,20,20,20,20,20,20,20,20,20});
        r1.setCigarString("10M");

        GATKSAMRecord r2 = ArtificialSAMUtils.createArtificialRead(header,"r2",0,1,10);
        r2.setReadBases(bases2);
        r2.setBaseQualities(new byte[] {20,20,20,20,20,20,20,20,20,20,20,20});
        r2.setCigarString("3=1X5=1X");

        GATKSAMRecord r3 = ArtificialSAMUtils.createArtificialRead(header,"r3",0,1,10);
        r3.setReadBases(bases2);
        r3.setBaseQualities(new byte[] {20,20,20,20,20,20,20,20,20,20,20,20});
        r3.setCigarString("3=1X5M1X");

        GATKSAMRecord r4  = ArtificialSAMUtils.createArtificialRead(header,"r4",0,1,10);
        r4.setReadBases(bases2);
        r4.setBaseQualities(new byte[] {20,20,20,20,20,20,20,20,20,20});
        r4.setCigarString("10M");

        List<GATKSAMRecord> reads = Arrays.asList(r1, r2, r3, r4);

        // create the iterator by state with the fake reads and fake records
        li = makeLTBS(reads);

        while (li.hasNext()) {
            AlignmentContext context = li.next();
            ReadBackedPileup pileup = context.getBasePileup();
            Assert.assertEquals(pileup.depthOfCoverage(), 4);
        }
    }

    @Test(enabled = true && ! DEBUG)
    public void testIndelsInRegularPileup() {
        final byte[] bases = new byte[] {'A','A','A','A','A','A','A','A','A','A'};
        final byte[] indelBases = new byte[] {'A','A','A','A','C','T','A','A','A','A','A','A'};

        GATKSAMRecord before = ArtificialSAMUtils.createArtificialRead(header,"before",0,1,10);
        before.setReadBases(bases);
        before.setBaseQualities(new byte[] {20,20,20,20,20,20,20,20,20,20});
        before.setCigarString("10M");

        GATKSAMRecord during = ArtificialSAMUtils.createArtificialRead(header,"during",0,2,10);
        during.setReadBases(indelBases);
        during.setBaseQualities(new byte[] {20,20,20,20,20,20,20,20,20,20,20,20});
        during.setCigarString("4M2I6M");

        GATKSAMRecord after  = ArtificialSAMUtils.createArtificialRead(header,"after",0,3,10);
        after.setReadBases(bases);
        after.setBaseQualities(new byte[] {20,20,20,20,20,20,20,20,20,20});
        after.setCigarString("10M");

        List<GATKSAMRecord> reads = Arrays.asList(before, during, after);

        // create the iterator by state with the fake reads and fake records
        li = makeLTBS(reads);

        boolean foundIndel = false;
        while (li.hasNext()) {
            AlignmentContext context = li.next();
            ReadBackedPileup pileup = context.getBasePileup().getBaseFilteredPileup(10);
            for (PileupElement p : pileup) {
                if (p.isBeforeInsertion()) {
                    foundIndel = true;
                    Assert.assertEquals(p.getLengthOfImmediatelyFollowingIndel(), 2, "Wrong event length");
                    Assert.assertEquals(p.getBasesOfImmediatelyFollowingInsertion(), "CT", "Inserted bases are incorrect");
                    break;
               }
            }

         }

         Assert.assertTrue(foundIndel,"Indel in pileup not found");
    }

    @Test(enabled = false && ! DEBUG)
    public void testWholeIndelReadInIsolation() {
        final int firstLocus = 44367789;

        GATKSAMRecord indelOnlyRead = ArtificialSAMUtils.createArtificialRead(header, "indelOnly", 0, firstLocus, 76);
        indelOnlyRead.setReadBases(Utils.dupBytes((byte)'A',76));
        indelOnlyRead.setBaseQualities(Utils.dupBytes((byte) '@', 76));
        indelOnlyRead.setCigarString("76I");

        List<GATKSAMRecord> reads = Arrays.asList(indelOnlyRead);

        // create the iterator by state with the fake reads and fake records
        li = makeLTBS(reads);

        // Traditionally, reads that end with indels bleed into the pileup at the following locus.  Verify that the next pileup contains this read
        // and considers it to be an indel-containing read.
        Assert.assertTrue(li.hasNext(),"Should have found a whole-indel read in the normal base pileup without extended events enabled");
        AlignmentContext alignmentContext = li.next();
        Assert.assertEquals(alignmentContext.getLocation().getStart(), firstLocus, "Base pileup is at incorrect location.");
        ReadBackedPileup basePileup = alignmentContext.getBasePileup();
        Assert.assertEquals(basePileup.getReads().size(),1,"Pileup is of incorrect size");
        Assert.assertSame(basePileup.getReads().get(0), indelOnlyRead, "Read in pileup is incorrect");
    }

    /**
     * Test to make sure that reads supporting only an indel (example cigar string: 76I) do
     * not negatively influence the ordering of the pileup.
     */
    @Test(enabled = true && ! DEBUG)
    public void testWholeIndelRead() {
        final int firstLocus = 44367788, secondLocus = firstLocus + 1;

        GATKSAMRecord leadingRead = ArtificialSAMUtils.createArtificialRead(header,"leading",0,firstLocus,76);
        leadingRead.setReadBases(Utils.dupBytes((byte)'A',76));
        leadingRead.setBaseQualities(Utils.dupBytes((byte)'@',76));
        leadingRead.setCigarString("1M75I");

        GATKSAMRecord indelOnlyRead = ArtificialSAMUtils.createArtificialRead(header,"indelOnly",0,secondLocus,76);
        indelOnlyRead.setReadBases(Utils.dupBytes((byte) 'A', 76));
        indelOnlyRead.setBaseQualities(Utils.dupBytes((byte)'@',76));
        indelOnlyRead.setCigarString("76I");

        GATKSAMRecord fullMatchAfterIndel = ArtificialSAMUtils.createArtificialRead(header,"fullMatch",0,secondLocus,76);
        fullMatchAfterIndel.setReadBases(Utils.dupBytes((byte)'A',76));
        fullMatchAfterIndel.setBaseQualities(Utils.dupBytes((byte)'@',76));
        fullMatchAfterIndel.setCigarString("75I1M");

        List<GATKSAMRecord> reads = Arrays.asList(leadingRead, indelOnlyRead, fullMatchAfterIndel);

        // create the iterator by state with the fake reads and fake records
        li = makeLTBS(reads, null, false);
        int currentLocus = firstLocus;
        int numAlignmentContextsFound = 0;

        while(li.hasNext()) {
            AlignmentContext alignmentContext = li.next();
            Assert.assertEquals(alignmentContext.getLocation().getStart(),currentLocus,"Current locus returned by alignment context is incorrect");

            if(currentLocus == firstLocus) {
                List<GATKSAMRecord> readsAtLocus = alignmentContext.getBasePileup().getReads();
                Assert.assertEquals(readsAtLocus.size(),1,"Wrong number of reads at locus " + currentLocus);
                Assert.assertSame(readsAtLocus.get(0),leadingRead,"leadingRead absent from pileup at locus " + currentLocus);
            }
            else if(currentLocus == secondLocus) {
                List<GATKSAMRecord> readsAtLocus = alignmentContext.getBasePileup().getReads();
                Assert.assertEquals(readsAtLocus.size(),1,"Wrong number of reads at locus " + currentLocus);
                Assert.assertSame(readsAtLocus.get(0),fullMatchAfterIndel,"fullMatchAfterIndel absent from pileup at locus " + currentLocus);
            }

            currentLocus++;
            numAlignmentContextsFound++;
        }

        Assert.assertEquals(numAlignmentContextsFound, 2, "Found incorrect number of alignment contexts");
    }

    /**
     * Test to make sure that reads supporting only an indel (example cigar string: 76I) are represented properly
     */
    @Test(enabled = false && ! DEBUG)
    public void testWholeIndelReadRepresentedTest() {
        final int firstLocus = 44367788, secondLocus = firstLocus + 1;

        GATKSAMRecord read1 = ArtificialSAMUtils.createArtificialRead(header,"read1",0,secondLocus,1);
        read1.setReadBases(Utils.dupBytes((byte) 'A', 1));
        read1.setBaseQualities(Utils.dupBytes((byte) '@', 1));
        read1.setCigarString("1I");

        List<GATKSAMRecord> reads = Arrays.asList(read1);

        // create the iterator by state with the fake reads and fake records
        li = makeLTBS(reads, null, false);

        while(li.hasNext()) {
            AlignmentContext alignmentContext = li.next();
            ReadBackedPileup p = alignmentContext.getBasePileup();
            Assert.assertTrue(p.getNumberOfElements() == 1);
            // TODO -- fix tests
//            PileupElement pe = p.iterator().next();
//            Assert.assertTrue(pe.isBeforeInsertion());
//            Assert.assertFalse(pe.isAfterInsertion());
//            Assert.assertEquals(pe.getBasesOfImmediatelyFollowingInsertion(), "A");
        }

        GATKSAMRecord read2 = ArtificialSAMUtils.createArtificialRead(header,"read2",0,secondLocus,10);
        read2.setReadBases(Utils.dupBytes((byte) 'A', 10));
        read2.setBaseQualities(Utils.dupBytes((byte) '@', 10));
        read2.setCigarString("10I");

        reads = Arrays.asList(read2);

        // create the iterator by state with the fake reads and fake records
        li = makeLTBS(reads, null, false);

        while(li.hasNext()) {
            AlignmentContext alignmentContext = li.next();
            ReadBackedPileup p = alignmentContext.getBasePileup();
            Assert.assertTrue(p.getNumberOfElements() == 1);
            // TODO -- fix tests
//            PileupElement pe = p.iterator().next();
//            Assert.assertTrue(pe.isBeforeInsertion());
//            Assert.assertFalse(pe.isAfterInsertion());
//            Assert.assertEquals(pe.getBasesOfImmediatelyFollowingInsertion(), "AAAAAAAAAA");
        }
    }


    /////////////////////////////////////////////
    // get event length and bases calculations //
    /////////////////////////////////////////////

    @DataProvider(name = "IndelLengthAndBasesTest")
    public Object[][] makeIndelLengthAndBasesTest() {
        final String EVENT_BASES = "ACGTACGTACGT";
        final List<Object[]> tests = new LinkedList<Object[]>();

        for ( int eventSize = 1; eventSize < 10; eventSize++ ) {
            for ( final CigarOperator indel : Arrays.asList(CigarOperator.D, CigarOperator.I) ) {
                final String cigar = String.format("2M%d%s1M", eventSize, indel.toString());
                final String eventBases = indel == CigarOperator.D ? "" : EVENT_BASES.substring(0, eventSize);
                final int readLength = 3 + eventBases.length();

                GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "read", 0, 1, readLength);
                read.setReadBases(("TT" + eventBases + "A").getBytes());
                final byte[] quals = new byte[readLength];
                for ( int i = 0; i < readLength; i++ )
                    quals[i] = (byte)(i % QualityUtils.MAX_SAM_QUAL_SCORE);
                read.setBaseQualities(quals);
                read.setCigarString(cigar);

                tests.add(new Object[]{read, indel, eventSize, eventBases.equals("") ? null : eventBases});
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = true && ! DEBUG, dataProvider = "IndelLengthAndBasesTest")
    public void testIndelLengthAndBasesTest(GATKSAMRecord read, final CigarOperator op, final int eventSize, final String eventBases) {
        // create the iterator by state with the fake reads and fake records
        li = makeLTBS(Arrays.asList((GATKSAMRecord)read), null, false);

        Assert.assertTrue(li.hasNext());

        final PileupElement firstMatch = getFirstPileupElement(li.next());

        Assert.assertEquals(firstMatch.getLengthOfImmediatelyFollowingIndel(), 0, "Length != 0 for site not adjacent to indel");
        Assert.assertEquals(firstMatch.getBasesOfImmediatelyFollowingInsertion(), null, "Getbases of following event should be null at non-adajenct event");

        Assert.assertTrue(li.hasNext());

        final PileupElement pe = getFirstPileupElement(li.next());

        if ( op == CigarOperator.D )
            Assert.assertTrue(pe.isBeforeDeletionStart());
        else
            Assert.assertTrue(pe.isBeforeInsertion());

        Assert.assertEquals(pe.getLengthOfImmediatelyFollowingIndel(), eventSize, "Length of event failed");
        Assert.assertEquals(pe.getBasesOfImmediatelyFollowingInsertion(), eventBases, "Getbases of following event failed");
    }

    private PileupElement getFirstPileupElement(final AlignmentContext context) {
        final ReadBackedPileup p = context.getBasePileup();
        Assert.assertEquals(p.getNumberOfElements(), 1);
        return p.iterator().next();
    }

    ////////////////////////////////////////////
    // comprehensive LIBS/PileupElement tests //
    ////////////////////////////////////////////

    @DataProvider(name = "MyLIBSTest")
    public Object[][] makeLIBSTest() {
        final List<Object[]> tests = new LinkedList<Object[]>();

//        tests.add(new Object[]{new LIBSTest("2=2D2=2X", 1)});
//        return tests.toArray(new Object[][]{});

        return createLIBSTests(
                Arrays.asList(1, 2),
                Arrays.asList(1, 2, 3, 4));

//        return createLIBSTests(
//                Arrays.asList(2),
//                Arrays.asList(3));
    }

    @Test(enabled = ! DEBUG, dataProvider = "MyLIBSTest")
    public void testLIBS(LIBSTest params) {
        // create the iterator by state with the fake reads and fake records
        final GATKSAMRecord read = params.makeRead();
        li = makeLTBS(Arrays.asList((GATKSAMRecord)read), null, false);
        final LIBS_position tester = new LIBS_position(read);

        int bpVisited = 0;
        int lastOffset = 0;
        while ( li.hasNext() ) {
            bpVisited++;

            AlignmentContext alignmentContext = li.next();
            ReadBackedPileup p = alignmentContext.getBasePileup();
            Assert.assertEquals(p.getNumberOfElements(), 1);
            PileupElement pe = p.iterator().next();

            Assert.assertEquals(p.getNumberOfDeletions(), pe.isDeletion() ? 1 : 0, "wrong number of deletions in the pileup");
            Assert.assertEquals(p.getNumberOfMappingQualityZeroReads(), pe.getRead().getMappingQuality() == 0 ? 1 : 0, "wront number of mapq reads in the pileup");

            tester.stepForwardOnGenome();

            if ( ! hasNeighboringPaddedOps(params.getElements(), pe.getCurrentCigarOffset()) ) {
                Assert.assertEquals(pe.isBeforeDeletionStart(), tester.isBeforeDeletionStart, "before deletion start failure");
                Assert.assertEquals(pe.isAfterDeletionEnd(), tester.isAfterDeletionEnd, "after deletion end failure");
            }

            Assert.assertEquals(pe.isBeforeInsertion(), tester.isBeforeInsertion, "before insertion failure");
            Assert.assertEquals(pe.isAfterInsertion(), tester.isAfterInsertion, "after insertion failure");
            Assert.assertEquals(pe.isNextToSoftClip(), tester.isNextToSoftClip, "next to soft clip failure");

            Assert.assertTrue(pe.getOffset() >= lastOffset, "Somehow read offsets are decreasing: lastOffset " + lastOffset + " current " + pe.getOffset());
            Assert.assertEquals(pe.getOffset(), tester.getCurrentReadOffset(), "Read offsets are wrong at " + bpVisited);

            Assert.assertEquals(pe.getCurrentCigarElement(), read.getCigar().getCigarElement(tester.currentOperatorIndex), "CigarElement index failure");
            Assert.assertEquals(pe.getOffsetInCurrentCigar(), tester.getCurrentPositionOnOperatorBase0(), "CigarElement index failure");

            Assert.assertEquals(read.getCigar().getCigarElement(pe.getCurrentCigarOffset()), pe.getCurrentCigarElement(), "Current cigar element isn't what we'd get from the read itself");

            Assert.assertTrue(pe.getOffsetInCurrentCigar() >= 0, "Offset into current cigar too small");
            Assert.assertTrue(pe.getOffsetInCurrentCigar() < pe.getCurrentCigarElement().getLength(), "Offset into current cigar too big");

            Assert.assertEquals(pe.getOffset(), tester.getCurrentReadOffset(), "Read offset failure");
            lastOffset = pe.getOffset();
        }

        final int expectedBpToVisit = read.getAlignmentEnd() - read.getAlignmentStart() + 1;
        Assert.assertEquals(bpVisited, expectedBpToVisit, "Didn't visit the expected number of bp");
    }

    // ------------------------------------------------------------
    //
    // Tests for keeping reads
    //
    // ------------------------------------------------------------

    @DataProvider(name = "LIBS_ComplexPileupTests")
    public Object[][] makeLIBS_ComplexPileupTests() {
        final List<Object[]> tests = new LinkedList<Object[]>();

        for ( final int downsampleTo : Arrays.asList(-1, 1, 2, 5, 10, 30)) {
            for ( final int nReadsPerLocus : Arrays.asList(1, 10, 60) ) {
                for ( final int nLoci : Arrays.asList(1, 10, 25) ) {
                    for ( final int nSamples : Arrays.asList(1, 2, 10) ) {
                        for ( final boolean keepReads : Arrays.asList(true, false) ) {
                            for ( final boolean grabReadsAfterEachCycle : Arrays.asList(true, false) ) {
//        for ( final int downsampleTo : Arrays.asList(1)) {
//            for ( final int nReadsPerLocus : Arrays.asList(1) ) {
//                for ( final int nLoci : Arrays.asList(1) ) {
//                    for ( final int nSamples : Arrays.asList(1) ) {
//                        for ( final boolean keepReads : Arrays.asList(true) ) {
//                            for ( final boolean grabReadsAfterEachCycle : Arrays.asList(true) ) {
                                tests.add(new Object[]{nReadsPerLocus, nLoci, nSamples,
                                        keepReads, grabReadsAfterEachCycle,
                                        downsampleTo});
                            }
                        }
                    }
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = true && ! DEBUG, dataProvider = "LIBS_ComplexPileupTests")
    public void testLIBS_ComplexPileupTests(final int nReadsPerLocus,
                                            final int nLoci,
                                            final int nSamples,
                                            final boolean keepReads,
                                            final boolean grabReadsAfterEachCycle,
                                            final int downsampleTo) {
        //logger.warn(String.format("testLIBSKeepSubmittedReads %d %d %d %b %b %b", nReadsPerLocus, nLoci, nSamples, keepReads, grabReadsAfterEachCycle, downsample));
        final int readLength = 10;

        final boolean downsample = downsampleTo != -1;
        final DownsamplingMethod downsampler = downsample
                ? new DownsamplingMethod(DownsampleType.BY_SAMPLE, downsampleTo, null)
                : new DownsamplingMethod(DownsampleType.NONE, null, null);

        final ArtificialBAMBuilder bamBuilder = new ArtificialBAMBuilder(header.getSequenceDictionary(), nReadsPerLocus, nLoci);
        bamBuilder.createAndSetHeader(nSamples).setReadLength(readLength).setAlignmentStart(1);

        final List<GATKSAMRecord> reads = bamBuilder.makeReads();
        li = new LocusIteratorByState(new FakeCloseableIterator<GATKSAMRecord>(reads.iterator()),
                downsampler, true, keepReads,
                genomeLocParser,
                bamBuilder.getSamples());

        final Set<GATKSAMRecord> seenSoFar = new HashSet<GATKSAMRecord>();
        final Set<GATKSAMRecord> keptReads = new HashSet<GATKSAMRecord>();
        int bpVisited = 0;
        while ( li.hasNext() ) {
            bpVisited++;
            final AlignmentContext alignmentContext = li.next();
            final ReadBackedPileup p = alignmentContext.getBasePileup();

            AssertWellOrderedPileup(p);

            if ( downsample ) {
                // just not a safe test
                //Assert.assertTrue(p.getNumberOfElements() <= maxDownsampledCoverage * nSamples, "Too many reads at locus after downsampling");
            } else {
                final int minPileupSize = nReadsPerLocus * nSamples;
                Assert.assertTrue(p.getNumberOfElements() >= minPileupSize);
            }

            // the number of reads starting here
            int nReadsStartingHere = 0;
            for ( final GATKSAMRecord read : p.getReads() )
                if ( read.getAlignmentStart() == alignmentContext.getPosition() )
                    nReadsStartingHere++;

            // we can have no more than maxDownsampledCoverage per sample
            final int maxCoveragePerLocus = downsample ? downsampleTo : nReadsPerLocus;
            Assert.assertTrue(nReadsStartingHere <= maxCoveragePerLocus * nSamples);

            seenSoFar.addAll(p.getReads());
            if ( keepReads && grabReadsAfterEachCycle ) {
                final List<GATKSAMRecord> locusReads = li.transferReadsFromAllPreviousPileups();


                if ( downsample ) {
                    // with downsampling we might have some reads here that were downsampled away
                    // in the pileup.  We want to ensure that no more than the max coverage per sample is added
                    Assert.assertTrue(locusReads.size() >= nReadsStartingHere);
                    Assert.assertTrue(locusReads.size() <= maxCoveragePerLocus * nSamples);
                } else {
                    Assert.assertEquals(locusReads.size(), nReadsStartingHere);
                }
                keptReads.addAll(locusReads);

                // check that all reads we've seen so far are in our keptReads
                for ( final GATKSAMRecord read : seenSoFar ) {
                    Assert.assertTrue(keptReads.contains(read), "A read that appeared in a pileup wasn't found in the kept reads: " + read);
                }
            }

            if ( ! keepReads )
                Assert.assertTrue(li.getReadsFromAllPreviousPileups().isEmpty(), "Not keeping reads but the underlying list of reads isn't empty");
        }

        if ( keepReads && ! grabReadsAfterEachCycle )
            keptReads.addAll(li.transferReadsFromAllPreviousPileups());

        if ( ! downsample ) { // downsampling may drop loci
            final int expectedBpToVisit = nLoci + readLength - 1;
            Assert.assertEquals(bpVisited, expectedBpToVisit, "Didn't visit the expected number of bp");
        }

        if ( keepReads ) {
            // check we have the right number of reads
            final int totalReads = nLoci * nReadsPerLocus * nSamples;
            if ( ! downsample ) { // downsampling may drop reads
                Assert.assertEquals(keptReads.size(), totalReads, "LIBS didn't keep the right number of reads during the traversal");

                // check that the order of reads is the same as in our read list
                for ( int i = 0; i < reads.size(); i++ ) {
                    final GATKSAMRecord inputRead = reads.get(i);
                    final GATKSAMRecord keptRead = reads.get(i);
                    Assert.assertSame(keptRead, inputRead, "Input reads and kept reads differ at position " + i);
                }
            } else {
                Assert.assertTrue(keptReads.size() <= totalReads, "LIBS didn't keep the right number of reads during the traversal");
            }

            // check uniqueness
            final Set<String> readNames = new HashSet<String>();
            for ( final GATKSAMRecord read : keptReads ) {
                Assert.assertFalse(readNames.contains(read.getReadName()), "Found duplicate reads in the kept reads");
                readNames.add(read.getReadName());
            }

            // check that all reads we've seen are in our keptReads
            for ( final GATKSAMRecord read : seenSoFar ) {
                Assert.assertTrue(keptReads.contains(read), "A read that appeared in a pileup wasn't found in the kept reads: " + read);
            }

            if ( ! downsample ) {
                // check that every read in the list of keep reads occurred at least once in one of the pileups
                for ( final GATKSAMRecord keptRead : keptReads ) {
                    Assert.assertTrue(seenSoFar.contains(keptRead), "There's a read " + keptRead + " in our keptReads list that never appeared in any pileup");
                }
            }
        }
    }

    private void AssertWellOrderedPileup(final ReadBackedPileup pileup) {
        if ( ! pileup.isEmpty() ) {
            int leftMostPos = -1;

            for ( final PileupElement pe : pileup ) {
                Assert.assertTrue(pileup.getLocation().getContig().equals(pe.getRead().getReferenceName()), "ReadBackedPileup contains an element " + pe + " that's on a different contig than the pileup itself");
                Assert.assertTrue(pe.getRead().getAlignmentStart() >= leftMostPos,
                        "ReadBackedPileup contains an element " + pe + " whose read's alignment start " + pe.getRead().getAlignmentStart()
                                + " occurs before the leftmost position we've seen previously " + leftMostPos);
            }
        }
    }

    // ---------------------------------------------------------------------------
    // make sure that downsampling isn't holding onto a bazillion reads
    //
    @DataProvider(name = "LIBS_NotHoldingTooManyReads")
    public Object[][] makeLIBS_NotHoldingTooManyReads() {
        final List<Object[]> tests = new LinkedList<Object[]>();

        for ( final int downsampleTo : Arrays.asList(1, 10)) {
            for ( final int nReadsPerLocus : Arrays.asList(100, 1000, 10000, 100000) ) {
                for ( final int payloadInBytes : Arrays.asList(0, 1024, 1024*1024) ) {
                    tests.add(new Object[]{nReadsPerLocus, downsampleTo, payloadInBytes});
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = true && ! DEBUG, dataProvider = "LIBS_NotHoldingTooManyReads")
//    @Test(enabled = true, dataProvider = "LIBS_NotHoldingTooManyReads", timeOut = 100000)
    public void testLIBS_NotHoldingTooManyReads(final int nReadsPerLocus, final int downsampleTo, final int payloadInBytes) {
        logger.warn(String.format("testLIBS_NotHoldingTooManyReads %d %d %d", nReadsPerLocus, downsampleTo, payloadInBytes));
        final int readLength = 10;

        final SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 100000);
        final int nSamples = 1;
        final List<String> samples = new ArrayList<String>(nSamples);
        for ( int i = 0; i < nSamples; i++ ) {
            final GATKSAMReadGroupRecord rg = new GATKSAMReadGroupRecord("rg" + i);
            final String sample = "sample" + i;
            samples.add(sample);
            rg.setSample(sample);
            rg.setPlatform(NGSPlatform.ILLUMINA.getDefaultPlatform());
            header.addReadGroup(rg);
        }

        final boolean downsample = downsampleTo != -1;
        final DownsamplingMethod downsampler = downsample
                ? new DownsamplingMethod(DownsampleType.BY_SAMPLE, downsampleTo, null)
                : new DownsamplingMethod(DownsampleType.NONE, null, null);

        // final List<GATKSAMRecord> reads = ArtificialSAMUtils.createReadStream(nReadsPerLocus, nLoci, header, 1, readLength);

        final WeakReadTrackingIterator iterator = new WeakReadTrackingIterator(nReadsPerLocus, readLength, payloadInBytes, header);

        li = new LocusIteratorByState(iterator,
                downsampler, true, false,
                genomeLocParser,
                samples);

        while ( li.hasNext() ) {
            final AlignmentContext next = li.next();
            Assert.assertTrue(next.getBasePileup().getNumberOfElements() <= downsampleTo, "Too many elements in pileup " + next);
            // TODO -- assert that there are <= X reads in memory after GC for some X
        }
    }

    private static class WeakReadTrackingIterator implements Iterator<GATKSAMRecord> {
        final int nReads, readLength, payloadInBytes;
        int readI = 0;
        final SAMFileHeader header;

        private WeakReadTrackingIterator(int nReads, int readLength, final int payloadInBytes, final SAMFileHeader header) {
            this.nReads = nReads;
            this.readLength = readLength;
            this.header = header;
            this.payloadInBytes = payloadInBytes;
        }

        @Override public boolean hasNext() { return readI < nReads; }
        @Override public void remove() { throw new UnsupportedOperationException("no remove"); }

        @Override
        public GATKSAMRecord next() {
            readI++;
            return makeRead();
        }

        private GATKSAMRecord makeRead() {
            final SAMReadGroupRecord rg = header.getReadGroups().get(0);
            final String readName = String.format("%s.%d.%s", "read", readI, rg.getId());
            final GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(header, readName, 0, 1, readLength);
            read.setReadGroup(new GATKSAMReadGroupRecord(rg));
            if ( payloadInBytes > 0 )
                // add a payload byte array to push memory use per read even higher
                read.setAttribute("PL", new byte[payloadInBytes]);
            return read;
        }
    }

    // ---------------------------------------------------------------------------
    //
    // make sure that adapter clipping is working properly in LIBS
    //
    // ---------------------------------------------------------------------------
    @DataProvider(name = "AdapterClippingTest")
    public Object[][] makeAdapterClippingTest() {
        final List<Object[]> tests = new LinkedList<Object[]>();

        final int start = 10;
        for ( final int goodBases : Arrays.asList(10, 20, 30) ) {
            for ( final int nClips : Arrays.asList(0, 1, 2, 10)) {
                for ( final boolean onLeft : Arrays.asList(true, false) ) {
                    final int readLength = nClips + goodBases;
                    GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "read1" , 0, start, readLength);
                    read.setProperPairFlag(true);
                    read.setReadPairedFlag(true);
                    read.setReadUnmappedFlag(false);
                    read.setMateUnmappedFlag(false);
                    read.setReadBases(Utils.dupBytes((byte) 'A', readLength));
                    read.setBaseQualities(Utils.dupBytes((byte) '@', readLength));
                    read.setCigarString(readLength + "M");

                    if ( onLeft ) {
                        read.setReadNegativeStrandFlag(true);
                        read.setMateNegativeStrandFlag(false);
                        read.setMateAlignmentStart(start + nClips);
                        read.setInferredInsertSize(readLength);
                        tests.add(new Object[]{nClips, goodBases, 0, read});
                    } else {
                        read.setReadNegativeStrandFlag(false);
                        read.setMateNegativeStrandFlag(true);
                        read.setMateAlignmentStart(start - 1);
                        read.setInferredInsertSize(goodBases - 1);
                        tests.add(new Object[]{0, goodBases, nClips, read});
                    }
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = true, dataProvider = "AdapterClippingTest")
    public void testAdapterClipping(final int nClipsOnLeft, final int nReadContainingPileups, final int nClipsOnRight, final GATKSAMRecord read) {

        li = new LocusIteratorByState(new FakeCloseableIterator<>(Collections.singletonList(read).iterator()),
                DownsamplingMethod.NONE, true, false,
                genomeLocParser,
                LocusIteratorByState.sampleListForSAMWithoutReadGroups());

        int expectedPos = read.getAlignmentStart() + nClipsOnLeft;
        int nPileups = 0;
        while ( li.hasNext() ) {
            final AlignmentContext next = li.next();
            Assert.assertEquals(next.getLocation().getStart(), expectedPos);
            nPileups++;
            expectedPos++;
        }

        final int nExpectedPileups = nReadContainingPileups;
        Assert.assertEquals(nPileups, nExpectedPileups, "Wrong number of pileups seen");
    }
}
