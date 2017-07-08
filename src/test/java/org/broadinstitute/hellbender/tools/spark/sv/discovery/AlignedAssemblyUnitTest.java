package org.broadinstitute.hellbender.tools.spark.sv.discovery;


import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.AlignedAssemblyOrExcuse;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignmentUtils;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class AlignedAssemblyUnitTest extends BaseTest{
    private static final String dummyRefName = "1";
    private static final int dummyRefId = Integer.valueOf(dummyRefName) - 1;
    private static final List<String> refNames = Collections.singletonList(dummyRefName);

    /**
     * These alignment records are supposed to be associated with the 4 possible types of evidence we could see for an inversion,
     *   where the chr1:101-200 bases are inverted, namely
     * INV55, where lower  contig coordinate is associated with a forward  strand lower reference coordinate, and
     *              higher contig coordinate is associated with a negative strand higher reference/contig coordinate
     * INV55, where lower  contig coordinate is associated with a forward  strand higher reference coordinate, and
     *              higher contig coordinate is associated with a negative strand lower reference/contig coordinate
     * INV33, where lower  contig coordinate is associated with a negative strand lower reference coordinate, and
     *              higher contig coordinate is associated with a forward  strand higher reference/contig coordinate
     * INV33, where lower  contig coordinate is associated with a forward  strand higher reference coordinate, and
     *              higher contig coordinate is associated with a negative strand lower reference/contig coordinate
     * Finally, one must be aware of the fact that BWA always outputs CIGAR with a '+'-strand representation,
     *   therefore we must use such in constructing the BwaMemAlignment's* @return objects stored in each array
     * @return an array of arrays, each composed of
     * [0] {@link BwaMemAlignment} object,
     * [1] expected reference interval,
     * [2] expected cigar,
     * [3] expected strandedness,
     * [4] expected start in assembled contig, 1-based, inclusive
     * [5] expected end in assembled contig, 1-based, inclusive
     * [6] expected contig length,
     * [7] expected {@link AlignedAssembly.AlignmentInterval} object (generated manually with all fields explicitly spell out and given to
     *                                      {@link AlignedAssembly.AlignmentInterval#AlignmentInterval(SimpleInterval, int, int, Cigar, boolean, int, int)}
     *                                      intended to be used for testing concordance between the two constructors)
     */
    @DataProvider(name = "AlignmentIntervalCtorTestForSimpleInversion")
    private Object[][] createInputsAndExpectedResults_BwaMemAlignmentConstruction() {

        final int[] alignmentStartsOnRef_0Based = {96, 196, 195, 95, 101, 201, 101, 201};
        final int[] alignmentStartsOnTig_0BasedInclusive = {0, 4, 0, 5, 0, 6, 0, 7};
        final int[] alignmentEndsOnTig_0BasedExclusive = {4, 8, 5, 10, 6, 12, 7, 14};
        final int[] seqLen = {8, 8, 10, 10, 12, 12, 14, 14};
        final int[] mapQualForBwaMemAlgn = {-1, 0, 10, 20, 30, 40, 50, 60};
        final boolean[] strandedness = {true, false, true, false, false, true, false, true};
        final String[] cigarStrings = {"4M4S", "4M4H", "5M5S", "5M5H", "6S6M", "6H6M", "7S7M", "7H7M"}; // each different number represent a different contig's pair of chimeric alignments
        final Cigar[] cigars = Arrays.stream(cigarStrings).map(TextCigarCodec::decode).toArray(Cigar[]::new);


        final Object[][] data = new Object[cigars.length][];
        for(int i=0; i<cigars.length; ++i) {
            final BwaMemAlignment bwaMemAlignment = new BwaMemAlignment(strandedness[i] ? 0 : SAMFlag.READ_REVERSE_STRAND.intValue(),
                    0, alignmentStartsOnRef_0Based[i], alignmentStartsOnRef_0Based[i]+cigars[i].getReferenceLength(),
                    strandedness[i] ? alignmentStartsOnTig_0BasedInclusive[i] : seqLen[i]-alignmentEndsOnTig_0BasedExclusive[i],
                    strandedness[i] ? alignmentEndsOnTig_0BasedExclusive[i] : seqLen[i]-alignmentStartsOnTig_0BasedInclusive[i],
                    mapQualForBwaMemAlgn[i], 0, 1, 1, cigarStrings[i],
                    null, null, 0, Integer.MIN_VALUE, Integer.MAX_VALUE);
            final SimpleInterval referenceInterval = new SimpleInterval(refNames.get(0), alignmentStartsOnRef_0Based[i]+1, bwaMemAlignment.getRefEnd());
            final AlignedAssembly.AlignmentInterval alignmentInterval = new AlignedAssembly.AlignmentInterval(referenceInterval, alignmentStartsOnTig_0BasedInclusive[i]+1, alignmentEndsOnTig_0BasedExclusive[i],
                    strandedness[i] ? cigars[i] : CigarUtils.invertCigar(cigars[i]),
                    strandedness[i], Math.max(SAMRecord.NO_MAPPING_QUALITY, bwaMemAlignment.getMapQual()), bwaMemAlignment.getNMismatches());
            data[i] = new Object[]{bwaMemAlignment, referenceInterval, strandedness[i] ? cigars[i] : CigarUtils.invertCigar(cigars[i]),
                    strandedness[i], alignmentStartsOnTig_0BasedInclusive[i]+1, alignmentEndsOnTig_0BasedExclusive[i], seqLen[i], mapQualForBwaMemAlgn[i], alignmentInterval};
        }
        return data;
    }

    @Test(dataProvider = "AlignmentIntervalCtorTestForSimpleInversion", groups = "sv")
    public void testConstructionFromBwaMemAlignment(final BwaMemAlignment bwaMemAlignment, final SimpleInterval expectedReferenceInterval, final Cigar expectedCigar,
                                                    final boolean expectedIsPositiveStrand, final int expectedStartOnContig_1BasedInclusive, final int expectedEndOnContig_1BasedInclusive,
                                                    final int expectedContigLength, final int expectedMapQualInBwaMemAlignment, final AlignedAssembly.AlignmentInterval expectedAlignmentInterval) {

        final AlignedAssembly.AlignmentInterval alignmentInterval = new AlignedAssembly.AlignmentInterval(bwaMemAlignment, refNames, expectedContigLength);
        Assert.assertEquals(alignmentInterval.referenceInterval, expectedReferenceInterval);
        Assert.assertEquals(alignmentInterval.cigarAlong5to3DirectionOfContig, expectedCigar);
        Assert.assertEquals(alignmentInterval.forwardStrand, expectedIsPositiveStrand);
        Assert.assertEquals(alignmentInterval.startInAssembledContig, expectedStartOnContig_1BasedInclusive);
        Assert.assertEquals(alignmentInterval.endInAssembledContig, expectedEndOnContig_1BasedInclusive);
        Assert.assertEquals(alignmentInterval.mapQual, Math.max(SAMRecord.NO_MAPPING_QUALITY,expectedMapQualInBwaMemAlignment));
        Assert.assertEquals(alignmentInterval, expectedAlignmentInterval);
    }



    @Test(dataProvider = "AlignmentIntervalCtorTestForSimpleInversion", groups = "sv")
    public void testConstructionFromSAMRecord(final BwaMemAlignment bwaMemAlignment, final SimpleInterval expectedReferenceInterval, final Cigar expectedCigar,
                                              final boolean expectedIsPositiveStrand, final int expectedStartOnContig_1BasedInclusive, final int expectedEndOnContig_1BasedInclusive,
                                              final int expectedContigLength, final int expectedMapQualInBwaMemAlignment, final AlignedAssembly.AlignmentInterval expectedAlignmentInterval) {

        final SAMRecord samRecord = BwaMemAlignmentUtils.applyAlignment("whatever", SVDiscoveryTestDataProvider.makeDummySequence(expectedContigLength, (byte)'A'), null, null, bwaMemAlignment, refNames, hg19Header, false, false);
        final AlignedAssembly.AlignmentInterval alignmentInterval = new AlignedAssembly.AlignmentInterval(samRecord);
        Assert.assertEquals(alignmentInterval.referenceInterval, expectedReferenceInterval);
        Assert.assertEquals(alignmentInterval.cigarAlong5to3DirectionOfContig, expectedCigar);
        Assert.assertEquals(alignmentInterval.forwardStrand, expectedIsPositiveStrand);
        Assert.assertEquals(alignmentInterval.startInAssembledContig, expectedStartOnContig_1BasedInclusive);
        Assert.assertEquals(alignmentInterval.endInAssembledContig, expectedEndOnContig_1BasedInclusive);
        Assert.assertEquals(alignmentInterval.mapQual, Math.max(SAMRecord.NO_MAPPING_QUALITY,expectedMapQualInBwaMemAlignment));
        Assert.assertEquals(alignmentInterval, expectedAlignmentInterval);
    }

    /**
     * @see #createInputsAndExpectedResults_BwaMemAlignmentConstruction()
     * 4 contigs of 3 assemblies pointing to the same inversion event.
     *  contig 1 belongs to assembly 1
     *  contig 2 and 3 belong to assembly 2
     *  contig 4 belongs to assembly 3.
     */
    @DataProvider(name = "AlignedAssemblySerializationTest")
    private Object[][] createInputsAndExpectedResults_Serialization() {

        final int[] alignmentStartsOnRef_0Based = {96, 196, 195, 95, 101, 201, 101, 201};
        final int[] alignmentStartsOnTig_0BasedInclusive = {0, 4, 0, 5, 0, 6, 0, 7};
        final int[] alignmentEndsOnTig_0BasedExclusive = {4, 8, 5, 10, 6, 12, 7, 14};
        final int[] seqLen = {8, 8, 10, 10, 12, 12, 14, 14};
        final int[] mapQual = {0, 1, 10, 20, 30, 40, 50, 60};
        final int[] mismatches = {0, 1, 1, 0, 2, 3, 3, 2};
        final boolean[] strandedness = {true, false, true, false, false, true, false, true};
        final String[] cigarStrings = {"4M4S", "4M4H", "5M5S", "5M5H", "6S6M", "6H6M", "7S7M", "7H7M"}; // each different number represent a different contig's pair of chimeric alignments
        final Cigar[] cigars = Arrays.stream(cigarStrings).map(TextCigarCodec::decode).toArray(Cigar[]::new);

        // these sequence are technically wrong the for the inversion event, but the test purpose is for serialization so it is irrelevant
        final byte[] dummySequenceForContigOne = SVDiscoveryTestDataProvider.makeDummySequence(seqLen[0], (byte)'A');

        final byte[] dummySequenceForContigTwo = SVDiscoveryTestDataProvider.makeDummySequence(seqLen[0], (byte)'T');
        final byte[] dummySequenceForContigThree = SVDiscoveryTestDataProvider.makeDummySequence(seqLen[0], (byte)'C');

        final byte[] dummySequenceForContigFour = SVDiscoveryTestDataProvider.makeDummySequence(seqLen[0], (byte)'G');


        final List<AlignedContig> allContigs = new ArrayList<>();

        for(int pair=0; pair<cigars.length/2; ++pair) {

            final List<AlignedAssembly.AlignmentInterval> alignmentIntervalsForSimpleInversion = new ArrayList<>(8);
            final SimpleInterval referenceIntervalLeft = new SimpleInterval(refNames.get(0), alignmentStartsOnRef_0Based[2*pair]+1, alignmentStartsOnRef_0Based[2*pair]+cigars[2*pair].getReferenceLength()+1);
            final AlignedAssembly.AlignmentInterval alignmentIntervalLeft = new AlignedAssembly.AlignmentInterval(referenceIntervalLeft, alignmentStartsOnTig_0BasedInclusive[2*pair]+1, alignmentEndsOnTig_0BasedExclusive[2*pair],
                    strandedness[2*pair] ? cigars[2*pair] : CigarUtils.invertCigar(cigars[2*pair]),
                    strandedness[2*pair], mapQual[2*pair], mismatches[2*pair]);
            alignmentIntervalsForSimpleInversion.add(alignmentIntervalLeft);
            final SimpleInterval referenceIntervalRight = new SimpleInterval(refNames.get(0), alignmentStartsOnRef_0Based[2*pair+1]+1, alignmentStartsOnRef_0Based[2*pair+1]+cigars[2*pair+1].getReferenceLength()+1);
            final AlignedAssembly.AlignmentInterval alignmentIntervalRight = new AlignedAssembly.AlignmentInterval(referenceIntervalRight, alignmentStartsOnTig_0BasedInclusive[2*pair+1]+1, alignmentEndsOnTig_0BasedExclusive[2*pair+1],
                    strandedness[2*pair+1] ? cigars[2*pair+1] : CigarUtils.invertCigar(cigars[2*pair+1]),
                    strandedness[2*pair+1], mapQual[2*pair+1], mismatches[2*pair+1]);
            alignmentIntervalsForSimpleInversion.add(alignmentIntervalRight);

            if (pair == 0) {
                allContigs.add( new AlignedContig(AlignedAssemblyOrExcuse.formatContigName(0, 0), dummySequenceForContigOne, alignmentIntervalsForSimpleInversion) );
            } else if (pair <3) {
                allContigs.add( new AlignedContig(AlignedAssemblyOrExcuse.formatContigName(1, pair-1), pair==1 ? dummySequenceForContigTwo : dummySequenceForContigThree, alignmentIntervalsForSimpleInversion) );
            } else {
                allContigs.add( new AlignedContig(AlignedAssemblyOrExcuse.formatContigName(2, 0), dummySequenceForContigFour, alignmentIntervalsForSimpleInversion) );
            }
        }

        final Object[][] data = new Object[3][];
        data[0] = new Object[]{0, new AlignedAssembly(0, allContigs.subList(0, 1))};
        data[1] = new Object[]{1, new AlignedAssembly(1, allContigs.subList(1, 3))};
        data[2] = new Object[]{2, new AlignedAssembly(2, allContigs.subList(3, 4))};

        return data;
    }

    @Test(dataProvider="AlignedAssemblySerializationTest", groups = "sv")
    public void testAlignedAssemblySerialization(final Integer assemblyID, final AlignedAssembly expectedAssembly) {

        final ByteArrayOutputStream bos = new ByteArrayOutputStream();
        final Output out = new Output(bos);
        final Kryo kryo = new Kryo();
        kryo.writeClassAndObject(out, expectedAssembly);
        out.flush();

        final ByteArrayInputStream bis = new ByteArrayInputStream(bos.toByteArray());
        final Input in = new Input(bis);
        @SuppressWarnings("unchecked")
        final AlignedAssembly roundTrip = (AlignedAssembly) kryo.readClassAndObject(in);
        Assert.assertEquals(roundTrip, expectedAssembly);
    }
}
