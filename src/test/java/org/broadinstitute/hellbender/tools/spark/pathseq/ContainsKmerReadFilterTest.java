package org.broadinstitute.hellbender.tools.spark.pathseq;

import org.broadinstitute.hellbender.tools.spark.sv.SVKmerShort;
import org.broadinstitute.hellbender.tools.spark.sv.SVKmerizer;
import org.broadinstitute.hellbender.tools.spark.utils.LargeLongHopscotchSet;
import org.broadinstitute.hellbender.tools.spark.utils.LongHopscotchSet;
import org.broadinstitute.hellbender.tools.spark.utils.LongIterator;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;

public class ContainsKmerReadFilterTest extends BaseTest {

    private final int kSize = 11;
    private LongHopscotchSet kmerSet;
    private SVKmerShort kmerMask;
    private File kmerSetFile;

    @BeforeMethod
    public void before() {
        final String kmerRef = "ATCGAGCGCTAGCGATGGCGCGCGATCGCGCTAGCGCGCTAGC";

        final byte[] maskBytes = new byte[]{3};
        kmerMask = SVKmerShort.getMask(maskBytes, kSize);

        final SVKmerizer kmerizer = new SVKmerizer(kmerRef.getBytes(), kSize, 1, new SVKmerShort(kSize));
        kmerSet = new LongHopscotchSet(kmerRef.length());
        while (kmerizer.hasNext()) {
            kmerSet.add(((SVKmerShort)kmerizer.next()).canonical(kSize).mask(kmerMask).getLong());
        }
        final LargeLongHopscotchSet largeKmerSet = new LargeLongHopscotchSet(kmerSet.size());
        LongIterator itr = kmerSet.iterator();
        while (itr.hasNext()) {
            largeKmerSet.add(itr.next());
        }
        kmerSetFile = createTempFile("kmerset",".hss");
        if (!kmerSetFile.delete()) {
            Assert.fail();
        }
        PSKmerUtils.writeKmerSet(kmerSetFile.getAbsolutePath(), new PSKmerSet(largeKmerSet, kSize, kmerMask));
    }

    @DataProvider(name = "sequenceStrings")
    public Object[][] getSequenceStrings() {
        return new Object[][]{
                {"GCGCTAGCGAT", 1, Boolean.FALSE},
                {"ATCGCTAGCGC", 1, Boolean.FALSE},
                {"GCGCTAGCGATG", 1, Boolean.FALSE},
                {"GCGCTAGCGATG", 2, Boolean.FALSE},
                {"GCGCTAGCGATGG", 3, Boolean.FALSE},
                {"GCGCTAGCGATT", 1, Boolean.FALSE},
                {"GAGCGCTAGCGTAGCGCGCTAG", 2, Boolean.FALSE},
                {"GCGCTAGCGAT", 2, Boolean.TRUE},
                {"GCGCTAGCGATA", 2, Boolean.TRUE},
                {"GCGCTAGCGAAT", 1, Boolean.TRUE},
                {"GCGCTAGCGANT", 1, Boolean.TRUE},
                {"", 1, Boolean.TRUE}
        };
    }

    @Test(dataProvider = "sequenceStrings")
    public void testTest(final String bases_in, final int kmerCountThreshold, final Boolean test_out) {
        final ContainsKmerReadFilterSpark filter = new ContainsKmerReadFilterSpark(kmerSetFile.getAbsolutePath(), kmerCountThreshold);
        final byte[] quals = new byte[bases_in.length()];
        Arrays.fill(quals, (byte) 30);
        GATKRead read_in = ArtificialReadUtils.createArtificialRead(bases_in.getBytes(), quals, "*");
        final Boolean test_i = filter.call(read_in);
        Assert.assertEquals(test_out, test_i);
        ContainsKmerReadFilter.closeKmerLib();
    }

}