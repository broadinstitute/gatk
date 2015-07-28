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
        final double allelicFractionSkew = 1.;
        final float minLogValue = -10.f;
        final String sampleName = "TCGA-02-0001-01C-01D-0182-01";

        final File snpFile = new File(TEST_SUB_DIR + "snps-simplified-for-allelic-fraction-segmentation.tsv");
        final List<AllelicCount> snpCounts = new AllelicCountCollection(snpFile).getCounts();

        final File resultFile = createTempFile("snp-segmenter-test-result", ".seg");
        SNPSegmenter.writeSegmentFile(snpCounts, allelicFractionSkew, sampleName, resultFile, minLogValue);

        final File expectedFile = new File(TEST_SUB_DIR + "snp-segmenter-test-expected.seg");

        Assert.assertTrue(resultFile.exists(), "SNPSegmenterTest output was not written to temp file: " + resultFile);

        final List<Segment> result = SegmentUtils.readUncalledSegments(resultFile);
        final List<Segment> expected = SegmentUtils.readUncalledSegments(expectedFile);

        Assert.assertEquals(result, expected);
    }

    //Tests that AllelicCounts are correctly transformed.
    @Test
    public void testTransformAllelicCounts() {
        final double allelicFractionSkew = 0.96;

        final File snpFile = new File(TEST_SUB_DIR + "snps-simplified-for-allelic-fraction-transformation.tsv");
        final List<AllelicCount> snpCounts = new AllelicCountCollection(snpFile).getCounts();

        final List<TargetCoverage> resultTargets = snpCounts.stream()
                .map(count -> count.toTargetCoverage("snp-target", allelicFractionSkew)).collect(Collectors.toList());

        final List<TargetCoverage> expectedTargets = Arrays.asList(
                new TargetCoverage("snp-target", new SimpleInterval("1", 212360, 212360), 0.3559124087591240),
                new TargetCoverage("snp-target", new SimpleInterval("1", 241501, 241501), 0.1397938144329896),
                new TargetCoverage("snp-target", new SimpleInterval("1", 242173, 242173), 0.0372413793103448),
                new TargetCoverage("snp-target", new SimpleInterval("1", 256641, 256641), 0.3381818181818182),
                new TargetCoverage("snp-target", new SimpleInterval("1", 261164, 261164), 0.4868049792531120),
                new TargetCoverage("snp-target", new SimpleInterval("1", 267204, 267204), 0.1735483870967741),
                new TargetCoverage("snp-target", new SimpleInterval("1", 282282, 282282), 0.3209090909090909),
                new TargetCoverage("snp-target", new SimpleInterval("1", 291649, 291649), 0.1250955414012738),
                new TargetCoverage("snp-target", new SimpleInterval("1", 376402, 376402), 0.2981818181818181),
                new TargetCoverage("snp-target", new SimpleInterval("1", 408347, 408347), 0.2032704402515723),
                new TargetCoverage("snp-target", new SimpleInterval("1", 415813, 415813), 0.1721739130434782),
                new TargetCoverage("snp-target", new SimpleInterval("1", 426517, 426517), 0.4366666666666666),
                new TargetCoverage("snp-target", new SimpleInterval("1", 429357, 429357), 0.1670588235294118),
                new TargetCoverage("snp-target", new SimpleInterval("1", 455201, 455201), 0.0112499999999999),
                new TargetCoverage("snp-target", new SimpleInterval("1", 466369, 466369), 0.4287179487179487),
                new TargetCoverage("snp-target", new SimpleInterval("1", 545461, 545461), 0.1164912280701754),
                new TargetCoverage("snp-target", new SimpleInterval("1", 665716, 665716), 0.2191304347826086),
                new TargetCoverage("snp-target", new SimpleInterval("1", 679370, 679370), 0.0056115107913669),
                new TargetCoverage("snp-target", new SimpleInterval("1", 704935, 704935), 0.4597979797979797));

        Assert.assertEquals(resultTargets.size(), expectedTargets.size());
        for (int index = 0; index < expectedTargets.size(); index++) {
            Assert.assertEquals(resultTargets.get(index).getInterval(), expectedTargets.get(index).getInterval());
            Assert.assertEquals(resultTargets.get(index).getName(), expectedTargets.get(index).getName());
            Assert.assertEquals(resultTargets.get(index).getCoverage(), expectedTargets.get(index).getCoverage(),
                    0.0000000001);
        }
    }
}
