package org.broadinstitute.hellbender.transformers;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CollectionUtil;
import org.broadinstitute.hellbender.utils.illumina.IlluminaAdapterPair;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.List;

public class IlluminaAdapterTrimTransformerTest {

    private static final List<String> ADAPTER_SEQUENCES = CollectionUtil.makeList(
            IlluminaAdapterPair.SINGLE_END.get5PrimeAdapter(),
            IlluminaAdapterPair.SINGLE_END.get3PrimeAdapter(),
            IlluminaAdapterPair.PAIRED_END.get5PrimeAdapter(),
            IlluminaAdapterPair.PAIRED_END.get3PrimeAdapter(),
            IlluminaAdapterPair.INDEXED.get5PrimeAdapter(),
            IlluminaAdapterPair.INDEXED.get3PrimeAdapter()
    );

    @DataProvider(name = "testData")
    public Object[][] getTestData() {
        return new Object[][]{
                {"AGCAGCTAGCTAGCTCGAATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTCGATCGA",
                        "AGCAGCTAGCTAGCTCG"},
                {"AGCAGCTAGCTAGCTCGAATGATACGGCGACCACCGAGATCTACACTCCTTCCCTACACGACGCTCTTCCGATCTCGATCGA",
                        "AGCAGCTAGCTAGCTCG"},
                {"AGCAGCTAGCTAGCTCGAATGATACGGCGACCACCGAGATCTACATTCCTTCCCTACACGACGCTCTTCCGATCTCGATCGA",
                        "AGCAGCTAGCTAGCTCG"},
                {"AGCAGCTAGCTAGCTCGAATGATACGGCGACCACCGAGATGTACATTCCTTCCCTACACGACGCTCTTCCGATCTCGATCGA",
                        "AGCAGCTAGCTAGCTCGAATGATACGGCGACCACCGAGATGTACATTCCTTCCCTACACGACGCTCTTCCGATCTCGATCGA"},
                {"GCATGCGATTCGAGTACTTACGGCATTCAGGTATCGAGATCGGAAGAG",
                        "GCATGCGATTCGAGTACTTACGGCATTCAGGTATCG"},
                {"GCATGCGATTCGAGTACTTACGGCATTCAGGTATCGAGATCGGAAGA",
                        "GCATGCGATTCGAGTACTTACGGCATTCAGGTATCGAGATCGGAAGA"},
        };
    }

    @Test(dataProvider = "testData")
    public void test(final String seqIn, final String seqOut) {
        final IlluminaAdapterTrimTransformer trans = new IlluminaAdapterTrimTransformer(2, 12, ADAPTER_SEQUENCES);
        final GATKRead read = new SAMRecordToGATKReadAdapter(new SAMRecord(null));
        read.setBases(seqIn.getBytes());
        read.setBaseQualities(new byte[seqIn.length()]);
        Assert.assertEquals(trans.apply(read).getBasesString(), seqOut);
    }
}