package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Collections;
import java.util.List;

public class AlignmentRegionUnitTest {

    private static final List<String> refNames = Collections.singletonList("1");
    private BwaMemAlignment createAlignment(final boolean fwdStrand, final int refStart, final int tigStart, final String cigar, final int nMismatches ) {
        return new BwaMemAlignment( fwdStrand ? 0 : SAMFlag.READ_REVERSE_STRAND.intValue(),
                                0, refStart, refStart+4, tigStart, tigStart+4, 60, nMismatches,
                                0, 0, cigar, null, null,
                                0, refStart+100, 104);
    }
    @Test
    public void testNormalizeStrandOnCreation() throws Exception {
        final AlignmentRegion ar1p1 = new AlignmentRegion("1", "1", 8,
                createAlignment(true, 96, 0, "4M4S", 0 ),
                refNames);
        Assert.assertEquals(ar1p1.referenceInterval, new SimpleInterval("1", 97, 100));
        Assert.assertEquals(ar1p1.startInAssembledContig, 1);
        Assert.assertEquals(ar1p1.endInAssembledContig, 4);
        Assert.assertTrue(ar1p1.forwardStrand);
        Assert.assertEquals(ar1p1.forwardStrandCigar, TextCigarCodec.decode("4M4S"));
        Assert.assertEquals(ar1p1.assembledContigLength, 8);
        Assert.assertEquals(ar1p1.mismatches, 0);
        Assert.assertEquals(ar1p1,
                new AlignmentRegion("1", "1", new SimpleInterval("1", 97, 100),
                                    TextCigarCodec.decode("4M4S"),
                                    true, 60, 0, 1, 4));

        final AlignmentRegion ar1p2 = new AlignmentRegion("1", "1", 8,
                createAlignment(false, 196, 4, "4H4M", 1 ),
                refNames);
        Assert.assertEquals(ar1p2.referenceInterval, new SimpleInterval("1", 197, 200));
        Assert.assertEquals(ar1p2.startInAssembledContig, 5);
        Assert.assertEquals(ar1p2.endInAssembledContig, 8);
        Assert.assertFalse(ar1p2.forwardStrand);
        Assert.assertEquals(ar1p2.forwardStrandCigar, TextCigarCodec.decode("4H4M"));
        Assert.assertEquals(ar1p2.assembledContigLength, 8);
        Assert.assertEquals(ar1p2.mismatches, 1);
        Assert.assertEquals(ar1p2,
                new AlignmentRegion("1", "1", new SimpleInterval("1", 197, 200),
                        TextCigarCodec.decode("4H4M"),
                        false, 60, 1, 5, 8));
    }

    @Test
    public void testParseAlignedAssembledContigLine() throws Exception {
        final String line = "100\t>contig-0 2498 0\t1\t7043012\t7044153\t+\t1141M1357S\t60\t1\t1141\t1";
        final AlignmentRegion region1 = AlignmentRegion.fromString(line.split(AlignmentRegion.STRING_REP_SEPARATOR, -1));
        Assert.assertEquals(region1.referenceInterval, new SimpleInterval("1", 7043012, 7044153));
        Assert.assertTrue(region1.forwardStrand);
        Assert.assertEquals(region1.forwardStrandCigar.toString(), "1141M1357S");
        Assert.assertEquals(region1.mapQual, 60);
        Assert.assertEquals(region1.startInAssembledContig, 1);
        Assert.assertEquals(region1.endInAssembledContig, 1141);
        Assert.assertEquals(region1.mismatches, 1);

        final String line2 = "100\tcontig-0\t1\t7044151\t7045306\t+\t1343S1155M\t60\t1344\t2498\t3";
        final AlignmentRegion region2 = AlignmentRegion.fromString(line2.split(AlignmentRegion.STRING_REP_SEPARATOR, -1));
        Assert.assertEquals(region2.referenceInterval, new SimpleInterval("1", 7044151, 7045306));
        Assert.assertTrue(region2.forwardStrand);
        Assert.assertEquals(region2.forwardStrandCigar.toString(), "1343S1155M");
        Assert.assertEquals(region2.mapQual, 60);
        Assert.assertEquals(region2.startInAssembledContig, 1344);
        Assert.assertEquals(region2.endInAssembledContig, 2498);
        Assert.assertEquals(region2.mismatches, 3);
    }
}