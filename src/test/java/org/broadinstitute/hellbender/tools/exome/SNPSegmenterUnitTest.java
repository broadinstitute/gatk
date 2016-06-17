package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.List;

/**
 * Unit tests for {@link SNPSegmenter}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class SNPSegmenterUnitTest extends BaseTest {
    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/exome/";

    /**
     * Tests that segments are correctly determined using allelic counts from SNP sites.
     * Segment-mean and target-number columns from expected segment file are not checked.
     */
    @Test
    public void testAllelicFractionBasedSegmentation() {
        final String sampleName = "test";

        final File snpFile = new File(TEST_SUB_DIR, "snps-simplified-for-allelic-fraction-segmentation.tsv");
        final List<AllelicCount> snpCounts = new AllelicCountCollection(snpFile).getCounts();
        final TargetCollection<AllelicCount> snps = new HashedListTargetCollection<>(snpCounts);

        final File resultFile = createTempFile("snp-segmenter-test-result", ".seg");
        SNPSegmenter.writeSegmentFile(snps, sampleName, resultFile);

        final File expectedFile = new File(TEST_SUB_DIR, "snp-segmenter-test-expected.seg");

        Assert.assertTrue(resultFile.exists(), "SNPSegmenterTest output was not written to temp file: " + resultFile);

        final List<SimpleInterval> result = SegmentUtils.readIntervalsFromSegmentFile(resultFile);
        final List<SimpleInterval> expected = SegmentUtils.readIntervalsFromSegmentFile(expectedFile);

        Assert.assertEquals(result, expected);
    }
}
