package org.broadinstitute.hellbender.transformers;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class SimpleRepeatMaskTransformerTest {

    @DataProvider(name = "testData")
    public Object[][] getTestData() {
        return new Object[][]{
                {"GATCGATCGCGCGACTAGCTACGACTAGCTACGATCTAGCAGCTACG",
                        "GATCGATCGCGCGACTAGCTACGACTAGCTACGATCTAGCAGCTACG", 20, 29, 30},
                {"GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCATGCATGCATGCATGC",
                        "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTGCATGCATGCATGC", 20, 29, 30},
                {"GCGCGCGCGCGCGCGCGCGCGCGCGCGCGTATGCATGCATGCATGC",
                        "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNATGCATGCATGCATGC", 20, 29, 30},
                {"GCGCGCGCGCGCGCGCGCGCGCGCGCGCATATGCATGCATGCATGC",
                        "GCGCGCGCGCGCGCGCGCGCGCGCGCGCATATGCATGCATGCATGC", 20, 29, 30},
                {"AGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCATGCATGCATGCATGC",
                        "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTGCATGCATGCATGC", 20, 29, 30},
                {"TAGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCATGCATGCATGCATGC",
                        "TNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTGCATGCATGCATGC", 20, 29, 30},
                {"TAGCGCGCGCGCGCGCGCGTGCGCGCGCGCGCATGCATGCATGCATGC",
                 "TANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNATGCATGCATGCATGC", 20, 29, 30},
                {"TTATATATAAATATATATATATATATATTATGCGCGCGGC",
                        "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCGCGCGGC", 29, 20, 30},
                {"GCGCGCTTATATATAAATATATATATATATATATTAT",
                        "GCGCGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN", 29, 20, 30},
                {"AAAAAAAAAA", "NNNNNNNNNN", 10, 10, 11},
                {"AAAAAAAAGA", "AAAAAAAAGA", 10, 10, 11},
                {"CATATATTTGGGGGAT", "NNNNNNNNNNNNNNNN", 3, 8, 8},
                {"CATATATTTGGGGGGAT","NNNNNNNNNNNNNNGAT", 3, 8, 8},
                {"", "", 29, 20, 30}
        };
    }

    @Test(dataProvider = "testData")
    public void test(final String seqIn, final String seqOut, final int threshAT, final int threshGC, final int windowSize) {
        final SimpleRepeatMaskTransformer trans = new SimpleRepeatMaskTransformer(threshAT, threshGC, windowSize);
        final GATKRead read = new SAMRecordToGATKReadAdapter(new SAMRecord(null));
        read.setBases(seqIn.getBytes());
        read.setBaseQualities(new byte[seqIn.length()]);
        Assert.assertEquals(new String(trans.apply(read).getBases()), seqOut);
    }

}