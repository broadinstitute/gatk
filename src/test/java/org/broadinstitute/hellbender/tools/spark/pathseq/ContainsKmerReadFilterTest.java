package org.broadinstitute.hellbender.tools.spark.pathseq;

import htsjdk.samtools.SAMUtils;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.ContainsKmerReadFilterSpark;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVKmer;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVKmerShort;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVKmerizer;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchSet;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;

public class ContainsKmerReadFilterTest extends BaseTest {

    private final int kSize = 11;
    private HopscotchSet<SVKmer> kmerSet;

    @BeforeMethod
    public void before() {
        final String kmerRef = "ATCGAGCGCTAGCGATGGCGCGCGATCGCGCTAGCGCGCTAGC";
        final SVKmerizer kmerizer = new SVKmerizer(kmerRef.getBytes(),kSize,new SVKmerShort(kSize));
        final ArrayList<SVKmer> kmerList = new ArrayList<>();
        while (kmerizer.hasNext()) {
            kmerList.add(kmerizer.next());
        }
        kmerSet = new HopscotchSet<>(kmerList);
    }

    @DataProvider(name = "sequenceStrings")
    public Object[][] getSequenceStrings() {
        return new Object[][] {
                {"GCGCGCGATCG", Boolean.FALSE},
                {"TAGCGATGGCGCGCGATCACGCTAG", Boolean.FALSE},
                {"TAGCGATGGCACGCGATCGAGCTAG", Boolean.TRUE},
                {"TAGCGA", Boolean.TRUE},
                {"", Boolean.TRUE}
        };
    }

    @Test(dataProvider = "sequenceStrings")
    public void testTest(final String bases_in, final Boolean test_out) throws Exception {
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        final ContainsKmerReadFilterSpark filter = new ContainsKmerReadFilterSpark(ctx.broadcast(kmerSet),kSize);
        final byte[] quals = bases_in.getBytes().clone();
        Arrays.fill(quals,(byte)'I');
        SAMUtils.fastqToPhred(quals);
        GATKRead read_in = ArtificialReadUtils.createArtificialRead(bases_in.getBytes(), quals, "*");
        final Boolean test_i = filter.call(read_in);
        Assert.assertEquals(test_out,test_i);
    }

}