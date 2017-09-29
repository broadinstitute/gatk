package org.broadinstitute.hellbender.transformers;

import htsjdk.samtools.SAMUtils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class BaseQualityClipReadTransformerTest {

    @DataProvider(name = "sequenceStrings")
    public Object[][] getSequenceStrings() {
        return new Object[][] {
            {"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
            "CTCAAGTACAAGCTGATCCAGACCTACAGGGTGATGTCATTAGAGGCACTGATAACACACACACTATGGGGTGGGGGTGGACAGTTCCCCACTGCAATCC",
            "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"},

            {"####################################################################################################",
            "",
            ""},

            {"############################IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
            "GGGTGATGTCATTAGAGGCACTGATAACACACACACTATGGGGTGGGGGTGGACAGTTCCCCACTGCAATCC",
            "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"},

            {"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII######",
            "CTCAAGTACAAGCTGATCCAGACCTACAGGGTGATGTCATTAGAGGCACTGATAACACACACACTATGGGGTGGGGGTGGACAGTTCCCCACTG",
            "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"}
        };
    }

    @Test(dataProvider = "sequenceStrings")
    public void testApply(final String quals_in, final String bases_out, final String quals_out) throws Exception {
        final String bases_in = "CTCAAGTACAAGCTGATCCAGACCTACAGGGTGATGTCATTAGAGGCACTGATAACACACACACTATGGGGTGGGGGTGGACAGTTCCCCACTGCAATCC";
        final BaseQualityClipReadTransformer trans = new BaseQualityClipReadTransformer(15);
        final GATKRead read_in = ArtificialReadUtils.createArtificialRead(bases_in.getBytes(),SAMUtils.fastqToPhred(quals_in),"*");
        final GATKRead read_out = trans.apply(read_in);
        Assert.assertEquals(read_out.getBases(),bases_out.getBytes());
        Assert.assertEquals(read_out.getBaseQualities(),SAMUtils.fastqToPhred(quals_out));
    }

}