package org.broadinstitute.hellbender.transformers;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
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
            "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
            "100M"},

            {"####################################################################################################",
            "",
            "",
            "*"},

            {"############################IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
            "GGGTGATGTCATTAGAGGCACTGATAACACACACACTATGGGGTGGGGGTGGACAGTTCCCCACTGCAATCC",
            "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
            "28H72M"},

            {"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII######",
            "CTCAAGTACAAGCTGATCCAGACCTACAGGGTGATGTCATTAGAGGCACTGATAACACACACACTATGGGGTGGGGGTGGACAGTTCCCCACTG",
            "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
            "94M6H"}
        };
    }

    @Test(dataProvider = "sequenceStrings")
    public void testApply(final String quals_in, final String basesOut, final String qualsOut, final String cigarOut) throws Exception {
        final String basesIn = "CTCAAGTACAAGCTGATCCAGACCTACAGGGTGATGTCATTAGAGGCACTGATAACACACACACTATGGGGTGGGGGTGGACAGTTCCCCACTGCAATCC";
        final BaseQualityClipReadTransformer trans = new BaseQualityClipReadTransformer(15);
        final Cigar cigar = new Cigar();
        cigar.add(new CigarElement(basesIn.length(), CigarOperator.MATCH_OR_MISMATCH));
        final GATKRead readIn = ArtificialReadUtils.createArtificialRead(basesIn.getBytes(),SAMUtils.fastqToPhred(quals_in), cigar.toString());
        final GATKRead readOut = trans.apply(readIn);
        Assert.assertEquals(readOut.getBases(),basesOut.getBytes());
        Assert.assertEquals(readOut.getBaseQualities(),SAMUtils.fastqToPhred(qualsOut));
        Assert.assertEquals(readOut.getCigar().toString(), cigarOut);
    }

}