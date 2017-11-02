package org.broadinstitute.hellbender.tools.walkers.validation.basicshortmutpileup;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.nio.file.Path;
import java.util.Collections;
import java.util.List;

public class AllelePileupCounterUnitTest extends GATKBaseTest {

    private static final String TEST_DREAM_BAM_DIR = largeFileTestDir + "mutect/dream_synthetic_bams/";
    private static final String TEST_DREAM_TUMOR_BAM = TEST_DREAM_BAM_DIR + "tumor.bam";
    private static final String TEST_DREAM_TUMOR_BAI = TEST_DREAM_TUMOR_BAM + ".bai";

    @Test(dataProvider = "basicCounts")
    public void testBasic(final SimpleInterval interval, final Allele ref, final Allele alt, int minBaseQuality, int gtRefCount, int gtAltCount) {

        final List<Path> bamFiles = Collections.singletonList(IOUtils.getPath(TEST_DREAM_TUMOR_BAM));
        final List<Path> baiFiles = Collections.singletonList(IOUtils.getPath(TEST_DREAM_TUMOR_BAI));
        final ReadsDataSource readsDataSource = new ReadsDataSource(bamFiles, baiFiles);
        final ReadPileup readPileup = GATKProtectedVariantContextUtils.getPileup(interval, () -> readsDataSource.query(interval));

        final AllelePileupCounter allelePileupCounter = new AllelePileupCounter(ref,
                Collections.singletonList(alt), minBaseQuality);
        allelePileupCounter.addPileup(readPileup);

        Assert.assertEquals(allelePileupCounter.getCountMap().entrySet().size(), 2);
        Assert.assertEquals(allelePileupCounter.getCountMap().get(ref).intValue(), gtRefCount);
        Assert.assertEquals(allelePileupCounter.getCountMap().get(alt).intValue(), gtAltCount);
    }

    // Chose a place in IGV and made sure the numbers made sense.
    // 20:10,022,820
    @DataProvider(name="basicCounts")
    public Object[][] createBasicTest() {
        return new Object[][]{
                {new SimpleInterval("20", 10022820, 10022820),
                        Allele.create("A", true), Allele.create("T"), 0, 19, 14},
                {new SimpleInterval("20", 10022820, 10022820),
                        Allele.create("A", true), Allele.create("T"), 20, 18, 14},
                {new SimpleInterval("20", 10022820, 10022820),
                        Allele.create("A", true), Allele.create("C"), 20, 18, 0},
                {new SimpleInterval("20", 10022820, 10022820),
                        Allele.create("A", true), Allele.create("T"), 100, 0, 0},

                // At this position, we have one low quality SNP G>A and three G>GA and one GA > G.  23 total reads.  3 support the indel,
                // and one supports the SNP (again, with low base quality).  18 total support the ref
                // two ref are BQ<20
                {new SimpleInterval("20", 10080550, 10080550),
                        Allele.create("G", true), Allele.create("A"), 100, 0, 0},
                {new SimpleInterval("20", 10080550, 10080550),
                        Allele.create("G", true), Allele.create("A"), 0, 18, 1},
                {new SimpleInterval("20", 10080550, 10080550),
                        Allele.create("G", true), Allele.create("A"), 20, 16, 0},
                {new SimpleInterval("20", 10080550, 10080550),
                        Allele.create("G", true), Allele.create("GA"), 0, 18, 3},
                {new SimpleInterval("20", 10080550, 10080550),
                        Allele.create("G", true), Allele.create("GA"), 20, 16, 3},
                {new SimpleInterval("20", 10080550, 10080550),
                        Allele.create("GA", true), Allele.create("G"), 0, 18, 1},
                {new SimpleInterval("20", 10080550, 10080550),
                        Allele.create("GA", true), Allele.create("G"), 20, 16, 1},
        };
    }
}
