package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignmentUtils;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class AlignmentIntervalUnitTest extends BaseTest {
    private static final String dummyRefName = "1";
    private static final List<String> refNames = Collections.singletonList(dummyRefName);

    @Test(groups = "sv")
    public void testAlignmentIntervalOverlap() throws Exception {

        final AlignmentInterval ar1 = new AlignmentInterval(new SimpleInterval("1",1,5), 1,5, TextCigarCodec.decode("5M5H"),true, 60, 0, 100, false, false);
        final AlignmentInterval ar2 = new AlignmentInterval(new SimpleInterval("1",10,16), 5,10, TextCigarCodec.decode("4S6M"),true, 60, 0, 100, false, false);
        Assert.assertEquals(AlignmentInterval.overlapOnContig(ar1, ar2), 1);

        final AlignmentInterval ar3 = new AlignmentInterval(new SimpleInterval("1",1,5), 1,5, TextCigarCodec.decode("5M5H"),true, 60, 0, 100, false, false);
        final AlignmentInterval ar4 = new AlignmentInterval(new SimpleInterval("1",11,16), 6,10, TextCigarCodec.decode("5S5M"),true, 60, 0, 100, false, false);
        Assert.assertEquals(AlignmentInterval.overlapOnContig(ar3, ar4), 0);
    }

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
     * [7] expected {@link AlignmentInterval} object (generated manually with all fields explicitly spell out and given to
     *                                      {@link AlignmentInterval#AlignmentInterval(SimpleInterval, int, int, Cigar, boolean, int, int, int, boolean, boolean)}
     *                                      intended to be used for testing concordance between the two constructors)
     */
    @DataProvider(name = "AlignmentIntervalCtorTestForSimpleInversion")
    Object[][] createInputsAndExpectedResults_BwaMemAlignmentConstruction() {

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
            final AlignmentInterval alignmentInterval = new AlignmentInterval(referenceInterval, alignmentStartsOnTig_0BasedInclusive[i]+1, alignmentEndsOnTig_0BasedExclusive[i],
                    strandedness[i] ? cigars[i] : CigarUtils.invertCigar(cigars[i]),
                    strandedness[i], Math.max(SAMRecord.NO_MAPPING_QUALITY, bwaMemAlignment.getMapQual()), bwaMemAlignment.getNMismatches(), bwaMemAlignment.getAlignerScore(), false, false);
            data[i] = new Object[]{bwaMemAlignment, referenceInterval, strandedness[i] ? cigars[i] : CigarUtils.invertCigar(cigars[i]),
                    strandedness[i], alignmentStartsOnTig_0BasedInclusive[i]+1, alignmentEndsOnTig_0BasedExclusive[i], seqLen[i], mapQualForBwaMemAlgn[i], alignmentInterval};
        }
        return data;
    }

    @Test(dataProvider = "AlignmentIntervalCtorTestForSimpleInversion", groups = "sv")
    public void testConstructionFromBwaMemAlignment(final BwaMemAlignment bwaMemAlignment, final SimpleInterval expectedReferenceInterval, final Cigar expectedCigar,
                                                    final boolean expectedIsPositiveStrand, final int expectedStartOnContig_1BasedInclusive, final int expectedEndOnContig_1BasedInclusive,
                                                    final int expectedContigLength, final int expectedMapQualInBwaMemAlignment, final AlignmentInterval expectedAlignmentInterval) {

        final AlignmentInterval alignmentInterval = new AlignmentInterval(bwaMemAlignment, refNames, expectedContigLength);
        Assert.assertEquals(alignmentInterval.referenceSpan, expectedReferenceInterval);
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
                                              final int expectedContigLength, final int expectedMapQualInBwaMemAlignment, final AlignmentInterval expectedAlignmentInterval) {

        final SAMRecord samRecord = BwaMemAlignmentUtils.applyAlignment("whatever", SVDiscoveryTestDataProvider.makeDummySequence(expectedContigLength, (byte)'A'), null, null, bwaMemAlignment, refNames, hg19Header, false, false);
        final AlignmentInterval alignmentInterval = new AlignmentInterval(samRecord);
        Assert.assertEquals(alignmentInterval.referenceSpan, expectedReferenceInterval);
        Assert.assertEquals(alignmentInterval.cigarAlong5to3DirectionOfContig, expectedCigar);
        Assert.assertEquals(alignmentInterval.forwardStrand, expectedIsPositiveStrand);
        Assert.assertEquals(alignmentInterval.startInAssembledContig, expectedStartOnContig_1BasedInclusive);
        Assert.assertEquals(alignmentInterval.endInAssembledContig, expectedEndOnContig_1BasedInclusive);
        Assert.assertEquals(alignmentInterval.mapQual, Math.max(SAMRecord.NO_MAPPING_QUALITY,expectedMapQualInBwaMemAlignment));
        Assert.assertEquals(alignmentInterval, expectedAlignmentInterval);
    }
}
