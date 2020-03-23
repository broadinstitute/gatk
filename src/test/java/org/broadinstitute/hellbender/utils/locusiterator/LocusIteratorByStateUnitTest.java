package org.broadinstitute.hellbender.utils.locusiterator;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.utils.NGSPlatform;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.downsampling.DownsampleType;
import org.broadinstitute.hellbender.utils.downsampling.DownsamplingMethod;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.ArtificialBAMBuilder;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

public final class LocusIteratorByStateUnitTest extends LocusIteratorByStateBaseTest {

    @Test(expectedExceptions = NoSuchElementException.class)
    public void testIteratingBeyondElements() {
        final int readLength = 3;
        final GATKRead mapped1 = ArtificialReadUtils.createArtificialRead(header, "mapped1", 0, 1, readLength);

        final List<GATKRead> reads = Arrays.asList(mapped1);

        final LocusIteratorByState li;
        li = makeLIBS(reads, DownsamplingMethod.NONE, true, header);
        for (int i = 0; i < readLength; i++) {
            Assert.assertTrue(li.hasNext());
            Assert.assertNotNull(li.next());
        }
        Assert.assertFalse(li.hasNext());
        li.next();
    }

    @Test
    public void testAdvance() {
        final int readLength = "ACTG".length();
        final GATKRead mapped1 = ArtificialReadUtils.createArtificialRead(header, "mapped1", 0, 1, readLength);
        mapped1.setBases("ACTG".getBytes());

        final List<GATKRead> reads = Arrays.asList(mapped1);

        final LocusIteratorByState li;
        li = makeLIBS(reads, DownsamplingMethod.NONE, true, header);
        final AlignmentContext alignmentContext2 = li.advanceToLocus(2, true);
        Assert.assertEquals(alignmentContext2.getPosition(), 2);

        final AlignmentContext alignmentContext3 = li.advanceToLocus(3, true);
        Assert.assertEquals(alignmentContext3.getPosition(), 3);

        final AlignmentContext alignmentContext10 = li.advanceToLocus(10, true);
        Assert.assertNull(alignmentContext10);

    }

    @Test(enabled = false)
    public void testUnmappedAndAllIReadsPassThrough() {
        final int readLength = 10;
        final GATKRead mapped1 = ArtificialReadUtils.createArtificialRead(header,"mapped1",0,1,readLength);
        final GATKRead mapped2 = ArtificialReadUtils.createArtificialRead(header,"mapped2",0,1,readLength);
        final GATKRead unmapped = ArtificialReadUtils.createArtificialRead(header,"unmapped",0,1,readLength);
        final GATKRead allI = ArtificialReadUtils.createArtificialRead(header,"allI",0,1,readLength);

        unmapped.setIsUnmapped();
        unmapped.setCigar("*");
        allI.setCigar(readLength + "I");

        final List<GATKRead> reads = Arrays.asList(mapped1, unmapped, allI, mapped2);

        // create the iterator by state with the fake reads and fake records
        final LocusIteratorByState li;
        li = makeLIBS(reads, DownsamplingMethod.NONE, true, header);

        Assert.assertTrue(li.hasNext());
        AlignmentContext context = li.next();
        ReadPileup pileup = context.getBasePileup();
        Assert.assertEquals(pileup.size(), 2, "Should see only 2 reads in pileup, even with unmapped and all I reads");

        final List<GATKRead> rawReads = li.transferReadsFromAllPreviousPileups();
        Assert.assertEquals(rawReads, reads, "Input and transferred read lists should be the same, and include the unmapped and all I reads");
    }

    @Test
    public void testAllIReadsPassThrough() {
        final int readLength = 10;
        final GATKRead mapped1 = ArtificialReadUtils.createArtificialRead(header,"mapped1",0,1,readLength);
        final GATKRead mapped2 = ArtificialReadUtils.createArtificialRead(header,"mapped2",0,1,readLength);
        final GATKRead allI = ArtificialReadUtils.createArtificialRead(header,"allI",0,1,readLength);

        allI.setCigar(readLength + "I");

        final List<GATKRead> reads = Arrays.asList(mapped1, allI, mapped2);

        // create the iterator by state with the fake reads and fake records
        final LocusIteratorByState li;
        li = makeLIBS(reads, DownsamplingMethod.NONE, true, header);

        Assert.assertTrue(li.hasNext());
        final AlignmentContext context = li.next();
        final ReadPileup pileup = context.getBasePileup();
        Assert.assertEquals(pileup.size(), 2, "Should see only 2 reads in pileup, even with all I reads");

        final List<GATKRead> rawReads = li.transferReadsFromAllPreviousPileups();
        Assert.assertEquals(rawReads, reads, "Input and transferred read lists should be the same, and include and all I reads");
    }

    @Test
    public void testXandEQOperators() {
        final byte[] bases1 = {'A','A','A','A','A','A','A','A','A','A'};
        final byte[] bases2 = {'A','A','A','C','A','A','A','A','A','C'};

        final GATKRead r1 = ArtificialReadUtils.createArtificialRead(header,"r1",0,1,10);
        r1.setBases(bases1);
        r1.setBaseQualities(new byte[] {20,20,20,20,20,20,20,20,20,20});
        r1.setCigar("10M");

        final GATKRead r2 = ArtificialReadUtils.createArtificialRead(header,"r2",0,1,10);
        r2.setBases(bases2);
        r2.setBaseQualities(new byte[] {20,20,20,20,20,20,20,20,20,20,20,20});
        r2.setCigar("3=1X5=1X");

        final GATKRead r3 = ArtificialReadUtils.createArtificialRead(header,"r3",0,1,10);
        r3.setBases(bases2);
        r3.setBaseQualities(new byte[] {20,20,20,20,20,20,20,20,20,20,20,20});
        r3.setCigar("3=1X5M1X");

        final GATKRead r4  = ArtificialReadUtils.createArtificialRead(header,"r4",0,1,10);
        r4.setBases(bases2);
        r4.setBaseQualities(new byte[] {20,20,20,20,20,20,20,20,20,20});
        r4.setCigar("10M");

        final List<GATKRead> reads = Arrays.asList(r1, r2, r3, r4);

        // create the iterator by state with the fake reads and fake records
        final LocusIteratorByState li;
        li = makeLIBS(reads, header);

        while (li.hasNext()) {
            final AlignmentContext context = li.next();
            final ReadPileup pileup = context.getBasePileup();
            Assert.assertEquals(pileup.size(), 4);
        }
    }

    @Test
    public void testIndelsInRegularPileup() {
        final byte[] bases = {'A','A','A','A','A','A','A','A','A','A'};
        final byte[] indelBases = {'A','A','A','A','C','T','A','A','A','A','A','A'};

        final GATKRead before = ArtificialReadUtils.createArtificialRead(header,"before",0,1,10);
        before.setBases(bases);
        before.setBaseQualities(new byte[] {20,20,20,20,20,20,20,20,20,20});
        before.setCigar("10M");

        final GATKRead during = ArtificialReadUtils.createArtificialRead(header,"during",0,2,10);
        during.setBases(indelBases);
        during.setBaseQualities(new byte[] {20,20,20,20,20,20,20,20,20,20,20,20});
        during.setCigar("4M2I6M");

        final GATKRead after  = ArtificialReadUtils.createArtificialRead(header,"after",0,3,10);
        after.setBases(bases);
        after.setBaseQualities(new byte[] {20,20,20,20,20,20,20,20,20,20});
        after.setCigar("10M");

        final List<GATKRead> reads = Arrays.asList(before, during, after);

        // create the iterator by state with the fake reads and fake records
        final LocusIteratorByState li;
        li = makeLIBS(reads, header);

        boolean foundIndel = false;
        while (li.hasNext()) {
            final AlignmentContext context = li.next();
            final ReadPileup pileup = context.getBasePileup().makeFilteredPileup(pe -> pe.getQual() >= 10);
            for (final PileupElement p : pileup) {
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

    /**
     * Test to make sure that reads supporting only an indel (example cigar string: 76I) do
     * not negatively influence the ordering of the pileup.
     */
    @Test
    public void testWholeIndelRead() {
        final int firstLocus = 44367788, secondLocus = firstLocus + 1;

        final GATKRead leadingRead = ArtificialReadUtils.createArtificialRead(header,"leading",0,firstLocus,76);
        leadingRead.setBases(Utils.dupBytes((byte)'A',76));
        leadingRead.setBaseQualities(Utils.dupBytes((byte)'@',76));
        leadingRead.setCigar("1M75I");

        final GATKRead indelOnlyRead = ArtificialReadUtils.createArtificialRead(header,"indelOnly",0,secondLocus,76);
        indelOnlyRead.setBases(Utils.dupBytes((byte) 'A', 76));
        indelOnlyRead.setBaseQualities(Utils.dupBytes((byte)'@',76));
        indelOnlyRead.setCigar("76I");

        final GATKRead fullMatchAfterIndel = ArtificialReadUtils.createArtificialRead(header,"fullMatch",0,secondLocus,76);
        fullMatchAfterIndel.setBases(Utils.dupBytes((byte)'A',76));
        fullMatchAfterIndel.setBaseQualities(Utils.dupBytes((byte)'@',76));
        fullMatchAfterIndel.setCigar("75I1M");

        final List<GATKRead> reads = Arrays.asList(leadingRead, indelOnlyRead, fullMatchAfterIndel);

        // create the iterator by state with the fake reads and fake records
        final LocusIteratorByState li;
        li = makeLIBS(reads, null, false, header);
        int currentLocus = firstLocus;
        int numAlignmentContextsFound = 0;

        while(li.hasNext()) {
            final AlignmentContext alignmentContext = li.next();
            Assert.assertEquals(alignmentContext.getLocation().getStart(),currentLocus,"Current locus returned by alignment context is incorrect");

            if(currentLocus == firstLocus) {
                final List<GATKRead> readsAtLocus = alignmentContext.getBasePileup().getReads();
                Assert.assertEquals(readsAtLocus.size(),1,"Wrong number of reads at locus " + currentLocus);
                Assert.assertSame(readsAtLocus.get(0),leadingRead,"leadingRead absent from pileup at locus " + currentLocus);
            }
            else if(currentLocus == secondLocus) {
                final List<GATKRead> readsAtLocus = alignmentContext.getBasePileup().getReads();
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
    @Test
    public void testWholeIndelReadRepresentedTest() {
        final int firstLocus = 44367788, secondLocus = firstLocus + 1;

        final GATKRead read1 = ArtificialReadUtils.createArtificialRead(header,"read1",0,secondLocus,1);
        read1.setBases(Utils.dupBytes((byte) 'A', 1));
        read1.setBaseQualities(Utils.dupBytes((byte) '@', 1));
        read1.setCigar("1I");

        List<GATKRead> reads = Arrays.asList(read1);

        // create the iterator by state with the fake reads and fake records
        LocusIteratorByState li;
        li = makeLIBS(reads, null, false, header);

        while(li.hasNext()) {
            final AlignmentContext alignmentContext = li.next();
            final ReadPileup p = alignmentContext.getBasePileup();
            Assert.assertTrue(p.size() == 1);
            PileupElement pe = p.iterator().next();
            Assert.assertTrue(pe.isBeforeInsertion());
            Assert.assertFalse(pe.isAfterInsertion());
            Assert.assertEquals(pe.getBasesOfImmediatelyFollowingInsertion(), "A");
        }

        final GATKRead read2 = ArtificialReadUtils.createArtificialRead(header,"read2",0,secondLocus,10);
        read2.setBases(Utils.dupBytes((byte) 'A', 10));
        read2.setBaseQualities(Utils.dupBytes((byte) '@', 10));
        read2.setCigar("10I");

        reads = Arrays.asList(read2);

        // create the iterator by state with the fake reads and fake records
        li = makeLIBS(reads, null, false, header);

        while(li.hasNext()) {
            final AlignmentContext alignmentContext = li.next();
            final ReadPileup p = alignmentContext.getBasePileup();
            Assert.assertTrue(p.size() == 1);
            PileupElement pe = p.iterator().next();
            Assert.assertTrue(pe.isBeforeInsertion());
            Assert.assertFalse(pe.isAfterInsertion());
            Assert.assertEquals(pe.getBasesOfImmediatelyFollowingInsertion(), "AAAAAAAAAA");
        }
    }

    /**
     * Test to make sure that if there are reads with Ns are keeped
     */
    @Test
    public void testKeepingNs() {
        final int firstLocus = 44367788, secondLocus = firstLocus + 1;

        final GATKRead read = ArtificialReadUtils.createArtificialRead(header,"read1",0,secondLocus,1);
        read.setBases(Utils.dupBytes((byte) 'N', 1));
        read.setBaseQualities(Utils.dupBytes((byte) '@', 1));
        read.setCigar("1I");

        // create the iterator by state with the fake reads and fake records
        LocusIteratorByState libs = makeLIBSwithNs(Collections.singletonList(read), header);

        while(libs.hasNext()) {
            final AlignmentContext alignmentContext = libs.next();
            final ReadPileup rp = alignmentContext.getBasePileup();
            Assert.assertEquals(rp.size(), 1);
            final PileupElement pe = rp.iterator().next();
            Assert.assertEquals(pe.getBase(), (byte) 'N');
        }

    }

    /////////////////////////////////////////////
    // get event length and bases calculations //
    /////////////////////////////////////////////

    @DataProvider(name = "IndelLengthAndBasesTest")
    public Object[][] makeIndelLengthAndBasesTest() {
        final String EVENT_BASES = "ACGTACGTACGT";
        final List<Object[]> tests = new LinkedList<>();

        for ( int eventSize = 1; eventSize < 10; eventSize++ ) {
            for ( final CigarOperator indel : Arrays.asList(CigarOperator.D, CigarOperator.I) ) {
                final String cigar = String.format("2M%d%s1M", eventSize, indel.toString());
                final String eventBases = indel == CigarOperator.D ? "" : EVENT_BASES.substring(0, eventSize);
                final int readLength = 3 + eventBases.length();

                final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "read", 0, 1, readLength);
                read.setBases(("TT" + eventBases + "A").getBytes());
                final byte[] quals = new byte[readLength];
                for ( int i = 0; i < readLength; i++ )
                    quals[i] = (byte)(i % QualityUtils.MAX_SAM_QUAL_SCORE);
                read.setBaseQualities(quals);
                read.setCigar(cigar);

                tests.add(new Object[]{read, indel, eventSize, eventBases.equals("") ? null : eventBases});
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = true, dataProvider = "IndelLengthAndBasesTest")
    public void testIndelLengthAndBasesTest(final GATKRead read, final CigarOperator op, final int eventSize, final String eventBases) {
        // create the iterator by state with the fake reads and fake records
        final LocusIteratorByState li;
        li = makeLIBS(Arrays.asList(read), null, false, header);

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
        final ReadPileup p = context.getBasePileup();
        Assert.assertEquals(p.size(), 1);
        return p.iterator().next();
    }

    ////////////////////////////////////////////
    // comprehensive LIBS/PileupElement tests //
    ////////////////////////////////////////////

    @DataProvider(name = "MyLIBSTest")
    public Object[][] makeLIBSTest() {
        return createLIBSTests(
                Arrays.asList(1, 2),
                Arrays.asList(1, 2, 3, 4));

    }

    @Test(enabled = true, dataProvider = "MyLIBSTest")
    public void testLIBS(final LIBSTest params) {
        // create the iterator by state with the fake reads and fake records
        final GATKRead read = params.makeRead();
        final LocusIteratorByState li;
        li = makeLIBS(Arrays.asList(read), null, false, header);
        final LIBS_position tester = new LIBS_position(read);

        int bpVisited = 0;
        int lastOffset = 0;
        while ( li.hasNext() ) {
            bpVisited++;

            final AlignmentContext alignmentContext = li.next();
            final ReadPileup p = alignmentContext.getBasePileup();
            Assert.assertEquals(p.size(), 1);
            final PileupElement pe = p.iterator().next();

            Assert.assertEquals(p.getNumberOfElements(el->el.isDeletion()), pe.isDeletion() ? 1 : 0, "wrong number of deletions in the pileup");
            Assert.assertEquals(p.getNumberOfElements(el->el.getRead().getMappingQuality() == 0), pe.getRead().getMappingQuality() == 0 ? 1 : 0, "wront number of mapq reads in the pileup");

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

        final int expectedBpToVisit = read.getEnd() - read.getStart() + 1;
        Assert.assertEquals(bpVisited, expectedBpToVisit, "Didn't visit the expected number of bp");
    }

    // ------------------------------------------------------------
    //
    // Tests for keeping reads
    //
    // ------------------------------------------------------------

    @DataProvider(name = "LIBS_ComplexPileupTests")
    public Object[][] makeLIBS_ComplexPileupTests() {
        final List<Object[]> tests = new LinkedList<>();

        for ( final int downsampleTo : Arrays.asList(-1, 1, 2, 5, 10, 30)) {
            for ( final int nReadsPerLocus : Arrays.asList(1, 10, 60) ) {
                for ( final int nLoci : Arrays.asList(1, 10, 25) ) {
                    for ( final int nSamples : Arrays.asList(1) ) {
                        for ( final boolean keepReads : Arrays.asList(true, false) ) {
                            for ( final boolean grabReadsAfterEachCycle : Arrays.asList(true, false) ) {
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

    @Test(enabled = true, dataProvider = "LIBS_ComplexPileupTests")
    public void testLIBS_ComplexPileupTests(final int nReadsPerLocus,
                                            final int nLoci,
                                            final int nSamples,
                                            final boolean keepReads,
                                            final boolean grabReadsAfterEachCycle,
                                            final int downsampleTo) {
        final int readLength = 10;

        final boolean downsample = downsampleTo != -1;
        final DownsamplingMethod downsampler = downsample
                ? new DownsamplingMethod(DownsampleType.BY_SAMPLE, downsampleTo, null)
                : new DownsamplingMethod(DownsampleType.NONE, null, null);

        final ArtificialBAMBuilder bamBuilder = new ArtificialBAMBuilder(header.getSequenceDictionary(), nReadsPerLocus, nLoci);
        bamBuilder.createAndSetHeader(nSamples).setReadLength(readLength).setAlignmentStart(1);

        final List<GATKRead> reads = bamBuilder.makeReads();
        final LocusIteratorByState li;
        li = new LocusIteratorByState(
                new FakeCloseableIterator<>(reads.iterator()),
                downsampler,
                keepReads,
                bamBuilder.getSamples(),
                bamBuilder.getHeader(),
                true
        );

        final Set<GATKRead> seenSoFar = new LinkedHashSet<>();
        final Set<GATKRead> keptReads = new LinkedHashSet<>();
        int bpVisited = 0;
        while ( li.hasNext() ) {
            bpVisited++;
            final AlignmentContext alignmentContext = li.next();
            final ReadPileup p = alignmentContext.getBasePileup();

            AssertWellOrderedPileup(p);

            if ( downsample ) {
                // just not a safe test
                //Assert.assertTrue(p.getNumberOfElements() <= maxDownsampledCoverage * nSamples, "Too many reads at locus after downsampling");
            } else {
                final int minPileupSize = nReadsPerLocus * nSamples;
                Assert.assertTrue(p.size() >= minPileupSize);
            }

            // the number of reads starting here
            int nReadsStartingHere = 0;
            for ( final GATKRead read : p.getReads() )
                if ( read.getStart() == alignmentContext.getPosition() )
                    nReadsStartingHere++;

            // we can have no more than maxDownsampledCoverage per sample
            final int maxCoveragePerLocus = downsample ? downsampleTo : nReadsPerLocus;
            Assert.assertTrue(nReadsStartingHere <= maxCoveragePerLocus * nSamples);

            seenSoFar.addAll(p.getReads());
            if ( keepReads && grabReadsAfterEachCycle ) {
                final List<GATKRead> locusReads = li.transferReadsFromAllPreviousPileups();


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
                for ( final GATKRead read : seenSoFar ) {
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
                    final GATKRead inputRead = reads.get(i);
                    final GATKRead keptRead = reads.get(i);
                    Assert.assertSame(keptRead, inputRead, "Input reads and kept reads differ at position " + i);
                }
            } else {
                Assert.assertTrue(keptReads.size() <= totalReads, "LIBS didn't keep the right number of reads during the traversal");
            }

            // check uniqueness
            final Set<String> readNames = new LinkedHashSet<>();
            for ( final GATKRead read : keptReads ) {
                Assert.assertFalse(readNames.contains(read.getName()), "Found duplicate reads in the kept reads");
                readNames.add(read.getName());
            }

            // check that all reads we've seen are in our keptReads
            for ( final GATKRead read : seenSoFar ) {
                Assert.assertTrue(keptReads.contains(read), "A read that appeared in a pileup wasn't found in the kept reads: " + read);
            }

            if ( ! downsample ) {
                // check that every read in the list of keep reads occurred at least once in one of the pileups
                for ( final GATKRead keptRead : keptReads ) {
                    Assert.assertTrue(seenSoFar.contains(keptRead), "There's a read " + keptRead + " in our keptReads list that never appeared in any pileup");
                }
            }
        }
    }

    private void AssertWellOrderedPileup(final ReadPileup pileup) {
        if ( ! pileup.isEmpty() ) {
            final int leftMostPos = -1;

            for ( final PileupElement pe : pileup ) {
                Assert.assertTrue(pileup.getLocation().getContig().equals(pe.getRead().getContig()), "ReadPileup contains an element " + pe + " that's on a different contig than the pileup itself");
                Assert.assertTrue(pe.getRead().getStart() >= leftMostPos,
                        "ReadPileup contains an element " + pe + " whose read's alignment start " + pe.getRead().getStart()
                                + " occurs before the leftmost position we've seen previously " + leftMostPos);
            }
        }
    }

    // ---------------------------------------------------------------------------
    // make sure that downsampling isn't holding onto a bazillion reads
    //
    @DataProvider(name = "LIBS_NotHoldingTooManyReads")
    public Object[][] makeLIBS_NotHoldingTooManyReads() {
        final List<Object[]> tests = new LinkedList<>();

        for ( final int downsampleTo : Arrays.asList(1, 10)) {
            for ( final int nReadsPerLocus : Arrays.asList(100, 1000, 10000, 100000) ) {
                for ( final int payloadInBytes : Arrays.asList(0, 1024, 1024*1024) ) {
                    tests.add(new Object[]{nReadsPerLocus, downsampleTo, payloadInBytes});
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = true, dataProvider = "LIBS_NotHoldingTooManyReads")
    public void testLIBS_NotHoldingTooManyReads(final int nReadsPerLocus, final int downsampleTo, final int payloadInBytes) {
        logger.warn(String.format("testLIBS_NotHoldingTooManyReads %d %d %d", nReadsPerLocus, downsampleTo, payloadInBytes));
        final int readLength = 10;

        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(1, 1, 100000);
        final int nSamples = 1;
        final List<String> samples = new ArrayList<>(nSamples);
        for ( int i = 0; i < nSamples; i++ ) {
            final SAMReadGroupRecord rg = new SAMReadGroupRecord("rg" + i);
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

        final WeakReadTrackingIterator iterator = new WeakReadTrackingIterator(nReadsPerLocus, readLength, payloadInBytes, header);

        final LocusIteratorByState li;
        li = new LocusIteratorByState(
                iterator,
                downsampler,
                false,
                samples,
                header,
                true
        );

        while ( li.hasNext() ) {
            final AlignmentContext next = li.next();
            Assert.assertTrue(next.getBasePileup().size() <= downsampleTo, "Too many elements in pileup " + next);
            // TODO -- assert that there are <= X reads in memory after GC for some X
        }
    }

    private static class WeakReadTrackingIterator implements Iterator<GATKRead> {
        final int nReads, readLength, payloadInBytes;
        int readI = 0;
        final SAMFileHeader header;

        private WeakReadTrackingIterator(final int nReads, final int readLength, final int payloadInBytes, final SAMFileHeader header) {
            this.nReads = nReads;
            this.readLength = readLength;
            this.header = header;
            this.payloadInBytes = payloadInBytes;
        }

        @Override
        public boolean hasNext() { return readI < nReads; }

        @Override
        public void remove() { throw new UnsupportedOperationException("no remove"); }

        @Override
        public GATKRead next() {
            readI++;
            return makeRead();
        }

        private GATKRead makeRead() {
            final SAMReadGroupRecord rg = header.getReadGroups().get(0);
            final String readName = String.format("%s.%d.%s", "read", readI, rg.getId());
            final GATKRead read = ArtificialReadUtils.createArtificialRead(header, readName, 0, 1, readLength);
            read.setReadGroup(rg.getId());
            if ( payloadInBytes > 0 ){
                // add a payload byte array to push memory use per read even higher
                read.setAttribute("PL", new byte[payloadInBytes]);
            }
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
        final List<Object[]> tests = new LinkedList<>();

        final int start = 10;
        for ( final int goodBases : Arrays.asList(10, 20, 30) ) {
            for ( final int nClips : Arrays.asList(0, 1, 2, 10)) {
                for ( final boolean onLeft : Arrays.asList(true, false) ) {
                    final int readLength = nClips + goodBases;
                    final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "read1" , 0, start, readLength);
                    read.setIsProperlyPaired(true);
                    read.setIsPaired(true);
                    read.setBases(Utils.dupBytes((byte) 'A', readLength));
                    read.setBaseQualities(Utils.dupBytes((byte) '@', readLength));
                    read.setCigar(readLength + "M");

                    if ( onLeft ) {
                        read.setIsReverseStrand(true);
                        read.setMateIsReverseStrand(false);
                        read.setMatePosition(read.getContig(), start + nClips);
                        read.setFragmentLength(readLength);
                        tests.add(new Object[]{nClips, goodBases, 0, read});
                    } else {
                        read.setIsReverseStrand(false);
                        read.setMateIsReverseStrand(true);
                        read.setMatePosition(read.getContig(), start - 1);
                        read.setFragmentLength(goodBases);
                        tests.add(new Object[]{0, goodBases, nClips, read});
                    }
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = true, dataProvider = "AdapterClippingTest")
    public void testAdapterClipping(final int nClipsOnLeft, final int nReadContainingPileups, final int nClipsOnRight, final GATKRead read) {

        final LocusIteratorByState li;
        li = new LocusIteratorByState(
                new FakeCloseableIterator<>(Collections.singletonList(read).iterator()),
                DownsamplingMethod.NONE,
                false,
                sampleListForSAMWithoutReadGroups(),
                header,
                true
        );

        int expectedPos = read.getStart() + nClipsOnLeft;
        int nPileups = 0;
        while ( li.hasNext() ) {
            final AlignmentContext next = li.next();
            Assert.assertEquals(next.getLocation().getStart(), expectedPos);
            nPileups++;
            expectedPos++;
        }

        final int nExpectedPileups = nReadContainingPileups;
        Assert.assertEquals(nPileups, nExpectedPileups, "\"Wrong number of pileups seen for " + read + " with " + nClipsOnLeft + " clipped bases.");
    }
}
