package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class AlignmentRegionUnitTest {

    private static final List<String> refNames = Collections.singletonList("1");

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
     * [7] expected {@link AlignmentRegion} object (generated manually with all fields explicitly spell out and given to
     *                                      {@link AlignmentRegion#AlignmentRegion(String, String, SimpleInterval, Cigar, boolean, int, int, int, int)},
     *                                      intended to be used for testing concordance between the two constructors)
     */
    @DataProvider(name = "AlignmentRegionCtorTestForSimpleInversion")
    private static Object[][] createInputsAndExpectedResults() {

        final int[] alignmentStartsOnRef_0Based = {96, 196, 195, 95, 101, 201, 101, 201};
        final int[] alignmentStartsOnTig_0BasedInclusive = {0, 4, 0, 5, 0, 6, 0, 7};
        final int[] alignmentEndsOnTig_0BasedExclusive = {4, 8, 5, 10, 6, 12, 7, 14};
        final int[] seqLen = {8, 8, 10, 10, 12, 12, 14, 14};
        final boolean[] strandedness = {true, false, true, false, false, true, false, true};
        final String[] cigarStrings = {"4M4S", "4M4H", "5M5S", "5M5H", "6S6M", "6H6M", "7S7M", "7H7M"}; // each different number represent a different contig's pair of chimeric alignments
        final Cigar[] cigars = Arrays.stream(cigarStrings).map(TextCigarCodec::decode).toArray(Cigar[]::new);


        final Object[][] data = new Object[cigars.length][];
        for(int i=0; i<cigars.length; ++i) {
            final BwaMemAlignment bwaMemAlignment = new BwaMemAlignment(strandedness[i] ? 0 : SAMFlag.READ_REVERSE_STRAND.intValue(),
                    0, alignmentStartsOnRef_0Based[i], alignmentStartsOnRef_0Based[i]+cigars[i].getReferenceLength(),
                    strandedness[i] ? alignmentStartsOnTig_0BasedInclusive[i] : seqLen[i]-alignmentEndsOnTig_0BasedExclusive[i],
                    strandedness[i] ? alignmentEndsOnTig_0BasedExclusive[i] : seqLen[i]-alignmentStartsOnTig_0BasedInclusive[i],
                    60, 0, 1, 1, cigarStrings[i],
                    null, null, 0, Integer.MIN_VALUE, Integer.MAX_VALUE);
            final SimpleInterval referenceInterval = new SimpleInterval(refNames.get(0), alignmentStartsOnRef_0Based[i]+1, bwaMemAlignment.getRefEnd());
            data[i] = new Object[]{bwaMemAlignment, referenceInterval, strandedness[i] ? cigars[i] : CigarUtils.invertCigar(cigars[i]), strandedness[i], alignmentStartsOnTig_0BasedInclusive[i]+1, alignmentEndsOnTig_0BasedExclusive[i], seqLen[i],
                    new AlignmentRegion("1", "1", referenceInterval, strandedness[i] ? cigars[i] : CigarUtils.invertCigar(cigars[i]), strandedness[i], bwaMemAlignment.getMapQual(), bwaMemAlignment.getNMismatches(), alignmentStartsOnTig_0BasedInclusive[i]+1, alignmentEndsOnTig_0BasedExclusive[i])};
        }
        return data;
    }

    @Test(dataProvider = "AlignmentRegionCtorTestForSimpleInversion")
    public void testConstruction(final BwaMemAlignment bwaMemAlignment, final SimpleInterval expectedReferenceInterval, final Cigar expectedCigar,
                                 final boolean expectedIsPositiveStrand, final int expectedStartOnContig_1BasedInclusive, final int expectedEndOnContig_1BasedInclusive,
                                 final int expectedContigLength, final AlignmentRegion expectedAlignmentRegion) {

        final AlignmentRegion alignmentRegion = new AlignmentRegion("1", "1", expectedContigLength, bwaMemAlignment, refNames);
        Assert.assertEquals(alignmentRegion.referenceInterval, expectedReferenceInterval);
        Assert.assertEquals(alignmentRegion.cigarAlong5to3DirectionOfContig, expectedCigar);
        Assert.assertEquals(alignmentRegion.forwardStrand, expectedIsPositiveStrand);
        Assert.assertEquals(alignmentRegion.startInAssembledContig, expectedStartOnContig_1BasedInclusive);
        Assert.assertEquals(alignmentRegion.endInAssembledContig, expectedEndOnContig_1BasedInclusive);
        Assert.assertEquals(alignmentRegion, expectedAlignmentRegion);
    }

    @Test
    public void testParseAlignedAssembledContigLine() throws Exception {
        final String line = "100\t>contig-0 2498 0\t1\t7043012\t7044153\t+\t1141M1357S\t60\t1\t1141\t1";
        final AlignmentRegion region1 = AlignmentRegion.fromString(line.split(AlignmentRegion.STRING_REP_SEPARATOR, -1));
        Assert.assertEquals(region1.referenceInterval, new SimpleInterval("1", 7043012, 7044153));
        Assert.assertTrue(region1.forwardStrand);
        Assert.assertEquals(region1.cigarAlong5to3DirectionOfContig.toString(), "1141M1357S");
        Assert.assertEquals(region1.mapQual, 60);
        Assert.assertEquals(region1.startInAssembledContig, 1);
        Assert.assertEquals(region1.endInAssembledContig, 1141);
        Assert.assertEquals(region1.mismatches, 1);

        final String line2 = "100\tcontig-0\t1\t7044151\t7045306\t+\t1343S1155M\t60\t1344\t2498\t3";
        final AlignmentRegion region2 = AlignmentRegion.fromString(line2.split(AlignmentRegion.STRING_REP_SEPARATOR, -1));
        Assert.assertEquals(region2.referenceInterval, new SimpleInterval("1", 7044151, 7045306));
        Assert.assertTrue(region2.forwardStrand);
        Assert.assertEquals(region2.cigarAlong5to3DirectionOfContig.toString(), "1343S1155M");
        Assert.assertEquals(region2.mapQual, 60);
        Assert.assertEquals(region2.startInAssembledContig, 1344);
        Assert.assertEquals(region2.endInAssembledContig, 2498);
        Assert.assertEquals(region2.mismatches, 3);
    }
}