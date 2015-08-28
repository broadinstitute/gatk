package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Unit tests for {@link SNPSegmenter}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class SNPSegmenterUnitTest extends BaseTest {
    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/tools/exome/";

    //Tests that segments are correctly determined using allelic counts from SNP sites.
    //Segment-mean and target-number columns from expected segment file are not checked.
    @Test
    public void testAllelicFractionBasedSegmentation() throws IOException {
        final float minLogValue = -10.f;
        final String sampleName = "TCGA-02-0001-01C-01D-0182-01";

        final File snpFile = new File(TEST_SUB_DIR + "snps-simplified-for-allelic-fraction-segmentation.tsv");
        final List<AllelicCount> snpCounts = new AllelicCountCollection(snpFile).getCounts();

        final File resultFile = createTempFile("snp-segmenter-test-result", ".seg");
        SNPSegmenter.writeSegmentFile(snpCounts, sampleName, resultFile, minLogValue);

        final File expectedFile = new File(TEST_SUB_DIR + "snp-segmenter-test-expected.seg");

        Assert.assertTrue(resultFile.exists(), "SNPSegmenterTest output was not written to temp file: " + resultFile);

        final List<SimpleInterval> result = SegmentUtils.readIntervalsFromSegfile(resultFile);
        final List<SimpleInterval> expected = SegmentUtils.readIntervalsFromSegfile(expectedFile);

        Assert.assertEquals(result, expected);
    }

    //Tests that AllelicCounts are correctly transformed.
    @Test
    public void testTransformAllelicCounts() {

        final File snpFile = new File(TEST_SUB_DIR + "snps-simplified-for-allelic-fraction-transformation.tsv");
        final List<AllelicCount> snpCounts = new AllelicCountCollection(snpFile).getCounts();

        final List<TargetCoverage> resultTargets = snpCounts.stream()
                .map(count -> count.toTargetCoverage("snp-target")).collect(Collectors.toList());

        final List<TargetCoverage> expectedTargets = Arrays.asList(
                new TargetCoverage("snp-target", new SimpleInterval("1", 212360, 212360), -3.01056924171),
                new TargetCoverage("snp-target", new SimpleInterval("1", 241501, 241501), -1.55551872283),
                new TargetCoverage("snp-target", new SimpleInterval("1", 242173, 242173), -1.05062607307),
                new TargetCoverage("snp-target", new SimpleInterval("1", 256641, 256641), -2.45943161864),
                new TargetCoverage("snp-target", new SimpleInterval("1", 261164, 261164), -4.91288933623),
                new TargetCoverage("snp-target", new SimpleInterval("1", 267204, 267204), -1.70626879694),
                new TargetCoverage("snp-target", new SimpleInterval("1", 282282, 282282), -2.65207669658),
                new TargetCoverage("snp-target", new SimpleInterval("1", 291649, 291649), -1.34042443850),
                new TargetCoverage("snp-target", new SimpleInterval("1", 376402, 376402), -2.45943161864),
                new TargetCoverage("snp-target", new SimpleInterval("1", 408347, 408347), -1.85345133665),
                new TargetCoverage("snp-target", new SimpleInterval("1", 415813, 415813), -1.52356195606),
                new TargetCoverage("snp-target", new SimpleInterval("1", 426517, 426517), -3.58496250072),
                new TargetCoverage("snp-target", new SimpleInterval("1", 429357, 429357), -1.50250034053),
                new TargetCoverage("snp-target", new SimpleInterval("1", 455201, 455201), -1.09310940439),
                new TargetCoverage("snp-target", new SimpleInterval("1", 466369, 466369), -4.28540221886),
                new TargetCoverage("snp-target", new SimpleInterval("1", 545461, 545461), -1.30932805811),
                new TargetCoverage("snp-target", new SimpleInterval("1", 665716, 665716), -1.93859945534),
                new TargetCoverage("snp-target", new SimpleInterval("1", 679370, 679370), -1.04212547567),
                new TargetCoverage("snp-target", new SimpleInterval("1", 704935, 704935), -5.62935662008));

        Assert.assertEquals(resultTargets.size(), expectedTargets.size());
        for (int index = 0; index < expectedTargets.size(); index++) {
            Assert.assertEquals(resultTargets.get(index).getInterval(), expectedTargets.get(index).getInterval());
            Assert.assertEquals(resultTargets.get(index).getName(), expectedTargets.get(index).getName());
            Assert.assertEquals(resultTargets.get(index).getCoverage(), expectedTargets.get(index).getCoverage(),
                    0.0000000001);
        }
    }
}
