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
    private File kmerSetFile;

    @BeforeMethod
    public void before() {
        final String kmerRef = "ATCGAGCGCTAGCGATGGCGCGCGATCGCGCTAGCGCGCTAGC";
        final SVKmerizer kmerizer = new SVKmerizer(kmerRef.getBytes(), kSize, 1, new SVKmerShort(kSize));
        kmerSet = new LongHopscotchSet(kmerRef.length());
        while (kmerizer.hasNext()) {
            kmerSet.add(kmerizer.next().canonical(kSize).getLong());
        }
        final LargeLongHopscotchSet largeKmerSet = new LargeLongHopscotchSet(1024 * 1024, kmerSet.size());
        LongIterator itr = kmerSet.iterator();
        while (itr.hasNext()) {
            largeKmerSet.add(itr.next());
        }
        kmerSetFile = createTempFile("kmerset",".hss");
        if (!kmerSetFile.delete()) {
            Assert.fail();
        }
        PSKmerLibUtils.writeLargeLongHopscotchSet(largeKmerSet, kmerSetFile.getAbsolutePath(), null);
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
        final byte[] mask = new byte[0];
        final ContainsKmerReadFilterSpark filter = new ContainsKmerReadFilterSpark(kmerSetFile.getAbsolutePath(), kSize,
                mask, kmerCountThreshold, null);
        final byte[] quals = bases_in.getBytes().clone();
        Arrays.fill(quals, (byte) 30);
        GATKRead read_in = ArtificialReadUtils.createArtificialRead(bases_in.getBytes(), quals, "*");
        final Boolean test_i = filter.call(read_in);
        Assert.assertEquals(test_out, test_i);
        ContainsKmerReadFilter.closeKmerLib();
    }

}