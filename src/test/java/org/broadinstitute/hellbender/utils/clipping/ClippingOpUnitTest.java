package org.broadinstitute.hellbender.utils.clipping;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public final class ClippingOpUnitTest extends GATKBaseTest {

    @Test (dataProvider = "SoftClipedReadsNewStart")
    public void testGetNewAlignmentStartOffset(final String preClip, final String postClip, final int expectedResult) {
        Cigar cpreClip = TextCigarCodec.decode(preClip);
        Cigar cpostClip = TextCigarCodec.decode(postClip);
        Assert.assertEquals(ClippingOp.getNewAlignmentStartOffset(cpostClip, cpreClip), expectedResult,
                "getNewAlignmentStartOffset returned "+ClippingOp.getNewAlignmentStartOffset(cpostClip, cpreClip)+
                        " when "+expectedResult+" was expected for "+preClip.toString()+" which was clipped to "+postClip.toString());
    }

    // provider fields: cigar string before clip, cigar string after clip, expected result for getNewAlignmentStartOffset() method
    @DataProvider(name = "SoftClipedReadsNewStart")
    public Object[][] makeRevertSoftClipsBeforeContig() {
        return new Object[][] {
                {"70M", "10S60M", 10},
                {"70M", "60M10S", 0},

                {"30M10N30M", "30S10N30M", 30},
                {"30M10N30M", "30S5N30M", 35},
                {"30M10N30M", "30M10N30S", 0},
                {"30M10N30M", "30S30M", 40},
                {"30M10N30M", "30M30S", 0},
                {"30M10N30M", "15S15M10N30M", 15},
                {"30M10N30M", "30M10N15M15S", 0},
                // Testing multiple sequential reference but not sequence consuming base types
                {"10N10D40M", "40M", 20},
                {"10N10D40M", "20S20M", 40},
                {"10N10D40M", "5N10D40M", 5},
                {"10N10D40M", "5D40M", 15},
                {"10N10D40M", "10N10D20M20S", 0},

                {"10S10N20M10N10M", "10S20M10N10M", 10},
                {"10S10N20M10N10M", "10S10N20M10S", 0},

                {"10S10I20M10I10M", "20S20M10I10M", 0},
                {"10S10I20M10I10M", "10S10I20M20S", 0},
                {"10S10I20M10I10M", "15S5I20M10I10M", 0},

                {"10H60M", "10H10S50M", 10},
                {"10H60M", "10H50M10S", 0},
                {"10H10S50M", "10H20S40M", 10},
                {"10X60M", "20S50M", 20},
                {"10I40N20M","10S20M",40}
        };
    }
}
