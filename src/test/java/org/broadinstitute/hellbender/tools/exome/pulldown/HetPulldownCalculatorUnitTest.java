package org.broadinstitute.hellbender.tools.exome.pulldown;

import htsjdk.samtools.*;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.SamLocusIterator;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.TestResources;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.stream.IntStream;

/**
 * Unit tests for {@link HetPulldownCalculator}.  Uses BAM and SNP files generated from hg19mini using wgsim.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class HetPulldownCalculatorUnitTest extends BaseTest {
    private static final String TEST_SUB_DIR = TestResources.publicTestDir + "org/broadinstitute/hellbender/tools/exome";

    private static final File NORMAL_BAM_FILE = new File(TEST_SUB_DIR, "normal.sorted.bam");
    private static final File NORMAL_UNSORTED_BAM_FILE = new File(TEST_SUB_DIR, "normal.unsorted.bam");
    private static final File TUMOR_BAM_FILE = new File(TEST_SUB_DIR, "tumor.sorted.bam");
    private static final File SNP_FILE = new File(TEST_SUB_DIR, "common_SNP.interval_list");
    private static final File REF_FILE = new File(TestResources.hg19MiniReference);

    private static final int MINIMUM_MAPPING_QUALITY = 30;
    private static final int MINIMUM_BASE_QUALITY = 20;

    private static final HetPulldownCalculator calculator = new HetPulldownCalculator(REF_FILE, SNP_FILE,
            MINIMUM_MAPPING_QUALITY, MINIMUM_BASE_QUALITY, ValidationStringency.STRICT);

    private static SAMFileHeader normalHeader;
    private static SAMFileHeader tumorHeader;

    @BeforeClass
    public void initHeaders() throws IOException {
        try (final SamReader normalBamReader = SamReaderFactory.makeDefault().open(NORMAL_BAM_FILE);
             final SamReader tumorBamReader = SamReaderFactory.makeDefault().open(TUMOR_BAM_FILE)) {
            normalHeader = normalBamReader.getFileHeader();
            tumorHeader = tumorBamReader.getFileHeader();
        }
    }

    private static Nucleotide.Counter makeBaseCounts(final int aCount, final int cCount,
                                                     final int gCount, final int tCount) {
        final Nucleotide.Counter baseCounts = new Nucleotide.Counter();
        IntStream.range(0, aCount).forEach(i -> baseCounts.add(Nucleotide.A));
        IntStream.range(0, cCount).forEach(i -> baseCounts.add(Nucleotide.C));
        IntStream.range(0, gCount).forEach(i -> baseCounts.add(Nucleotide.G));
        IntStream.range(0, tCount).forEach(i -> baseCounts.add(Nucleotide.T));
        return baseCounts;
    }

    private static void assertEqualBaseCounts(final Nucleotide.Counter actual, final Nucleotide.Counter expected) {
        for (final Nucleotide base : HetPulldownCalculator.BASES) {
            Assert.assertEquals(actual.get(base), expected.get(base));
        }
    }


    @DataProvider(name = "inputGetPileupBaseCount")
    public Object[][] inputGetPileupBaseCount() throws IOException {
        try (final SamReader bamReader = SamReaderFactory.makeDefault().open(NORMAL_BAM_FILE)) {
            final IntervalList intervals = new IntervalList(bamReader.getFileHeader());
            intervals.add(new Interval("1", 100, 100));
            intervals.add(new Interval("1", 11000, 11000));
            intervals.add(new Interval("1", 14000, 14000));
            intervals.add(new Interval("1", 14630, 14630));

            final SamLocusIterator locusIterator = new SamLocusIterator(bamReader, intervals);

            final Nucleotide.Counter baseCounts1 = makeBaseCounts(0, 0, 0, 0);
            final Nucleotide.Counter baseCounts2 = makeBaseCounts(0, 9, 0, 0);
            final Nucleotide.Counter baseCounts3 = makeBaseCounts(12, 0, 0, 0);
            final Nucleotide.Counter baseCounts4 = makeBaseCounts(0, 0, 8, 9);

            if (!locusIterator.hasNext()) {
                throw new SAMException("Can't get locus to start iteration. Check that " + NORMAL_BAM_FILE.toString()
                        + " contains 1:0-16000.");
            }
            final SamLocusIterator.LocusInfo locus1 = locusIterator.next();
            final SamLocusIterator.LocusInfo locus2 = locusIterator.next();
            final SamLocusIterator.LocusInfo locus3 = locusIterator.next();
            final SamLocusIterator.LocusInfo locus4 = locusIterator.next();
            locusIterator.close();

            return new Object[][]{
                    {locus1, baseCounts1},
                    {locus2, baseCounts2},
                    {locus3, baseCounts3},
                    {locus4, baseCounts4}
            };
        }
    }

    @Test(dataProvider = "inputGetPileupBaseCount")
    public void testGetPileupBaseCount(final SamLocusIterator.LocusInfo locus,
                                       final Nucleotide.Counter expected) {
        final Nucleotide.Counter result = HetPulldownCalculator.getPileupBaseCounts(locus);
        assertEqualBaseCounts(result, expected);
    }

    @DataProvider(name = "inputIsPileupHetCompatible")
    public Object[][] inputIsPileupHetCompatible() {
        final Nucleotide.Counter baseCountsUsualHet = makeBaseCounts(50, 50, 0, 0);
        final Nucleotide.Counter baseCountsUsualHom = makeBaseCounts(50, 1, 0, 0);
        final Nucleotide.Counter baseCountsEdgeHom = makeBaseCounts(21, 1, 0, 8);
        final Nucleotide.Counter baseCountsEmpty = makeBaseCounts(0, 0, 0, 0);

        //if pval < pvalThreshold, expected = false
        return new Object[][]{
                {baseCountsUsualHet, 100, 0.05, true}, //pval = 1.0
                {baseCountsUsualHom, 51, 0.05, false}, //pval = 4.6185277824406525e-14
                {baseCountsEdgeHom, 30, 0.05, false},  //pval = 0.04277394525706768
                {baseCountsEdgeHom, 30, 0.04, true},   //pval = 0.04277394525706768
                {baseCountsEmpty, 0, 0.05, false},     //pval = 1.0
        };
    }

    @Test(dataProvider = "inputIsPileupHetCompatible")
    public void testIsPileupHetCompatible(final Nucleotide.Counter baseCounts, final int totalBaseCounts,
                                          final double pvalThreshold,
                                          final boolean expected) {
        final boolean result = HetPulldownCalculator.isPileupHetCompatible(baseCounts, totalBaseCounts,
                pvalThreshold);
        Assert.assertEquals(result, expected);
    }

    @DataProvider(name = "inputGetNormalHetPulldown")
    public Object[][] inputGetNormalHetPulldown() {
        final Pulldown normalHetPulldown1 = new Pulldown(normalHeader);
        normalHetPulldown1.add(new AllelicCount(new SimpleInterval("1", 11522, 11522), 7, 4));
        normalHetPulldown1.add(new AllelicCount(new SimpleInterval("1", 12098, 12098), 8, 6));
        normalHetPulldown1.add(new AllelicCount(new SimpleInterval("1", 14630, 14630), 9, 8));
        normalHetPulldown1.add(new AllelicCount(new SimpleInterval("2", 14689, 14689), 6, 9));
        normalHetPulldown1.add(new AllelicCount(new SimpleInterval("2", 14982, 14982), 6, 5));

        //changing pValThreshold from 0.05 -> 0.95 only keeps hets close to balanced
        final Pulldown normalHetPulldown2 = new Pulldown(normalHeader);
        normalHetPulldown2.add(new AllelicCount(new SimpleInterval("1", 14630, 14630), 9, 8));
        normalHetPulldown2.add(new AllelicCount(new SimpleInterval("2", 14982, 14982), 6, 5));

        return new Object[][]{
                {0.05, normalHetPulldown1},
                {0.95, normalHetPulldown2}
        };
    }

    @Test(dataProvider = "inputGetNormalHetPulldown")
    public void testGetNormalHetPulldown(final double pvalThreshold, final Pulldown expected) {
        final Pulldown result = calculator.getNormal(NORMAL_BAM_FILE, pvalThreshold, 10);
        Assert.assertEquals(result, expected);
    }

    @Test(expectedExceptions = UserException.class)
    public void testGetHetPulldownWithUnsortedBAMFile() {
        final Pulldown result = calculator.getNormal(NORMAL_UNSORTED_BAM_FILE, 0.05, 10);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testBadpValue() {
        final Pulldown result = calculator.getNormal(NORMAL_BAM_FILE, -1, 10);
    }

    @DataProvider(name = "inputGetTumorHetPulldown")
    public Object[][] inputGetTumorHetPulldown() {
        final Pulldown tumorHetPulldown = new Pulldown(normalHeader);
        tumorHetPulldown.add(new AllelicCount(new SimpleInterval("1", 11522, 11522), 7, 4));
        tumorHetPulldown.add(new AllelicCount(new SimpleInterval("1", 12098, 12098), 8, 6));
        tumorHetPulldown.add(new AllelicCount(new SimpleInterval("1", 14630, 14630), 9, 8));
        tumorHetPulldown.add(new AllelicCount(new SimpleInterval("2", 14689, 14689), 6, 9));
        tumorHetPulldown.add(new AllelicCount(new SimpleInterval("2", 14982, 14982), 6, 5));

        final IntervalList normalHetIntervals = new IntervalList(tumorHeader);
        normalHetIntervals.add(new Interval("1", 11522, 11522));
        normalHetIntervals.add(new Interval("1", 12098, 12098));
        normalHetIntervals.add(new Interval("1", 14630, 14630));
        normalHetIntervals.add(new Interval("2", 14689, 14689));
        normalHetIntervals.add(new Interval("2", 14982, 14982));

        return new Object[][]{
                {normalHetIntervals, tumorHetPulldown}
        };
    }

    @DataProvider(name = "inputGetTumorHetPulldownMin15")
    public Object[][] inputGetTumorHetPulldown15() {
        final Pulldown tumorHetPulldown = new Pulldown(normalHeader);
        tumorHetPulldown.add(new AllelicCount(new SimpleInterval("1", 14630, 14630), 9, 8));
        tumorHetPulldown.add(new AllelicCount(new SimpleInterval("2", 14689, 14689), 6, 9));

        final IntervalList normalHetIntervals = new IntervalList(tumorHeader);
        normalHetIntervals.add(new Interval("1", 14630, 14630));
        normalHetIntervals.add(new Interval("2", 14689, 14689));

        return new Object[][]{
                {normalHetIntervals, tumorHetPulldown}
        };
    }


    @Test(dataProvider = "inputGetTumorHetPulldown")
    public void testGetTumorHetPulldown(final IntervalList normalHetIntervals,
                                        final Pulldown expected) {
        final Pulldown result = calculator.getTumor(TUMOR_BAM_FILE, normalHetIntervals, 10);
        Assert.assertEquals(result, expected);
    }

    @Test(dataProvider = "inputGetTumorHetPulldownMin15")
    public void testGetTumorHetPulldown15(final IntervalList normalHetIntervals,
                                        final Pulldown expected) {
        final Pulldown result = calculator.getTumor(TUMOR_BAM_FILE, normalHetIntervals, 15);
        Assert.assertEquals(result, expected);
    }
}
