package org.broadinstitute.hellbender.transformers;

import htsjdk.samtools.SAMUtils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import java.util.Arrays;

public class BaseQualityReadTransformerTest {

    @DataProvider(name = "sequenceStrings")
    public Object[][] getSequenceStrings() {
        return new Object[][] {
                {"/IIIIIIIIIIIIIIIIIII", "NAAAAAAAAAAAAAAAAAAA"},
                {"IIIIIIIIIIIIIIIIIIII", "AAAAAAAAAAAAAAAAAAAA"},
                {"00000000000000000000", "AAAAAAAAAAAAAAAAAAAA"},
                {"0IIIIIIIIIIIIIIIIII",  "AAAAAAAAAAAAAAAAAAA"}
        };
    }

    @Test(dataProvider = "sequenceStrings")
    public void testTest(final String quals_in, final String test_out) throws Exception {
        final BaseQualityReadTransformer filter = new BaseQualityReadTransformer(15);
        final byte[] bases = quals_in.getBytes().clone();
        Arrays.fill(bases,(byte)'A');
        final GATKRead read_in = ArtificialReadUtils.createArtificialRead(bases, SAMUtils.fastqToPhred(quals_in), "*");
        final GATKRead test_i = filter.apply(read_in);
        Assert.assertEquals(test_out,test_i.getBasesString());
    }

}