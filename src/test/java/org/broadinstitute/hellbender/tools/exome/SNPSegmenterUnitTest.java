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
    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/exome/";

    /**
     * Tests that segments are correctly determined using allelic counts from SNP sites.
     //Segment-mean and target-number columns from expected segment file are not checked.
     * @throws IOException  if temporary target file cannot be created by {@link SNPSegmenter#writeSegmentFile}
     */
    @Test
    public void testAllelicFractionBasedSegmentation() throws IOException {
        final String sampleName = "test";

        final File snpFile = new File(TEST_SUB_DIR + "snps-simplified-for-allelic-fraction-segmentation.tsv");
        final List<AllelicCount> snpCounts = new AllelicCountCollection(snpFile).getCounts();

        final File resultFile = createTempFile("snp-segmenter-test-result", ".seg");
        SNPSegmenter.writeSegmentFile(snpCounts, sampleName, resultFile);

        final File expectedFile = new File(TEST_SUB_DIR + "snp-segmenter-test-expected.seg");

        Assert.assertTrue(resultFile.exists(), "SNPSegmenterTest output was not written to temp file: " + resultFile);

        final List<SimpleInterval> result = SegmentUtils.readIntervalsFromSegfile(resultFile);
        final List<SimpleInterval> expected = SegmentUtils.readIntervalsFromSegfile(expectedFile);

        Assert.assertEquals(result, expected);
    }

    /**
     * Tests that AllelicCounts are correctly transformed to alternate-allele fractions.
     */
    @Test
    public void testTransformAllelicCountsToAltAlleleFractions() {

        final File snpFile = new File(TEST_SUB_DIR + "snps-simplified-for-allelic-fraction-transformation.tsv");
        final List<AllelicCount> snpCounts = new AllelicCountCollection(snpFile).getCounts();

        final List<Double> resultFractions = snpCounts.stream().map(count -> count.toAltAlleleFraction())
                .collect(Collectors.toList());

        final List<Double> expectedFractions = Arrays.asList(
                1.0000000000,
                0.0000000000,
                0.5172413793,
                0.8181818182,
                0.9668049793,
                0.3064516129,
                0.1590909091,
                0.6050955414,
                0.1818181818,
                0.2767295597,
                0.6521739130,
                0.9166666667,
                0.6470588235,
                0.4687500000,
                0.0512820513,
                0.5964912281,
                0.2608695652,
                0.4856115108,
                0.0202020202
        );

        Assert.assertEquals(resultFractions.size(), expectedFractions.size());
        for (int index = 0; index < expectedFractions.size(); index++) {
            Assert.assertEquals(resultFractions.get(index), expectedFractions.get(index), 0.0000000001);
        }
    }

    /**
     * Tests that AllelicCounts are correctly transformed to target coverages with minor-allele fractions.
     */
    @Test
    public void testTransformAllelicCountsToLog2MinorAlleleFractionTargetCoverages() {

        final File snpFile = new File(TEST_SUB_DIR + "snps-simplified-for-allelic-fraction-transformation.tsv");
        final List<AllelicCount> snpCounts = new AllelicCountCollection(snpFile).getCounts();

        final List<TargetCoverage> resultTargets = snpCounts.stream()
                .map(count -> count.toMinorAlleleFractionTargetCoverage("snp-target")).collect(Collectors.toList());

        final List<TargetCoverage> expectedTargets = Arrays.asList(
                new TargetCoverage("snp-target", new SimpleInterval("1", 212360, 212360), 0.00000000),
                new TargetCoverage("snp-target", new SimpleInterval("1", 241501, 241501), 0.00000000),
                new TargetCoverage("snp-target", new SimpleInterval("1", 242173, 242173), 0.48275862),
                new TargetCoverage("snp-target", new SimpleInterval("1", 256641, 256641), 0.18181818),
                new TargetCoverage("snp-target", new SimpleInterval("1", 261164, 261164), 0.03319502),
                new TargetCoverage("snp-target", new SimpleInterval("1", 267204, 267204), 0.30645161),
                new TargetCoverage("snp-target", new SimpleInterval("1", 282282, 282282), 0.15909091),
                new TargetCoverage("snp-target", new SimpleInterval("1", 291649, 291649), 0.39490446),
                new TargetCoverage("snp-target", new SimpleInterval("1", 376402, 376402), 0.18181818),
                new TargetCoverage("snp-target", new SimpleInterval("1", 408347, 408347), 0.27672956),
                new TargetCoverage("snp-target", new SimpleInterval("1", 415813, 415813), 0.34782609),
                new TargetCoverage("snp-target", new SimpleInterval("1", 426517, 426517), 0.08333333),
                new TargetCoverage("snp-target", new SimpleInterval("1", 429357, 429357), 0.35294118),
                new TargetCoverage("snp-target", new SimpleInterval("1", 455201, 455201), 0.46875000),
                new TargetCoverage("snp-target", new SimpleInterval("1", 466369, 466369), 0.05128205),
                new TargetCoverage("snp-target", new SimpleInterval("1", 545461, 545461), 0.40350877),
                new TargetCoverage("snp-target", new SimpleInterval("1", 665716, 665716), 0.26086957),
                new TargetCoverage("snp-target", new SimpleInterval("1", 679370, 679370), 0.48561151),
                new TargetCoverage("snp-target", new SimpleInterval("1", 704935, 704935), 0.02020202));

        Assert.assertEquals(resultTargets.size(), expectedTargets.size());
        for (int index = 0; index < expectedTargets.size(); index++) {
            Assert.assertEquals(resultTargets.get(index).getInterval(), expectedTargets.get(index).getInterval());
            Assert.assertEquals(resultTargets.get(index).getName(), expectedTargets.get(index).getName());
            Assert.assertEquals(resultTargets.get(index).getCoverage(), expectedTargets.get(index).getCoverage(),
                    0.00000001);
        }
    }
}
