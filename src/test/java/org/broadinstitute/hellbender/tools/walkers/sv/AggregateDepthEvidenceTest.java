package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.engine.GatkToolIntegrationTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

public class AggregateDepthEvidenceTest extends GatkToolIntegrationTest {

    // If true, update the expected outputs in tests that assert an exact or approximate match vs. prior output,
    // instead of actually running the tests.
    public static final boolean UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS = false;

    public static final double FLOAT_TOLERANCE = 1e-4;

    public static final String TOOL_TEST_DIR = toolsTestDir + "walkers/sv/AggregateDepthEvidence";
    public static final String LARGE_FILE_TEST_DIR = largeFileTestDir + "sv";

    public static final long LARGE_VARIANT_SIZE = 2500000L;
    public static final int LARGE_VARIANT_POINTS = 500;
    public static final int LARGE_VARIANT_WINDOW = 2000;
    public static final int NUM_BINS = 10;
    public static final int MAX_QUAL = 999;
    public static final double POWER_THRESHOLD = 0.8;

    public static final List<String> attributesWithJitter = Arrays.asList(
            GATKSVVCFConstants.READ_DEPTH_QUALITY_ATTRIBUTE,
            GATKSVVCFConstants.READ_DEPTH_SECOND_MAX_QUALITY_ATTRIBUTE,
            GATKSVVCFConstants.READ_DEPTH_MEDIAN_SEPARATION_ATTRIBUTE
    );

    /**
     * Make sure that someone didn't leave the {@value UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS} toggle turned on
     */
    @Test
    public void assertThatExpectedOutputUpdateToggleIsDisabled() {
        Assert.assertFalse(UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS, "The toggle to update expected outputs should not be left enabled");
    }

    @Test
    public void test() {
        final File vcfFile = new File(TOOL_TEST_DIR, "aggregate_depth_test.vcf.gz");
        final File rdFile = new File(LARGE_FILE_TEST_DIR, "aggregate_depth_test.rd.txt.gz");
        final File medianCountsFile = new File(TOOL_TEST_DIR, "aggregate_depth_test.medianCov.tsv");
        final String outputFile = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? TOOL_TEST_DIR + "/aggregate_depth_test.output.vcf.gz" : createTempFile("aggregated", ".vcf.gz").getAbsolutePath();

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addVCF(vcfFile)
                .addOutput(outputFile)
                .addReference(hg38Reference)
                .add(AggregateDepthEvidence.DEPTH_EVIDENCE_FILE_PATH_LONG_NAME, rdFile)
                .add(AggregateDepthEvidence.MEDIAN_COUNTS_FILE_PATH_LONG_NAME, medianCountsFile)
                .add(AggregateDepthEvidence.LARGE_VARIANT_SIZE_LONG_NAME, LARGE_VARIANT_SIZE)
                .add(AggregateDepthEvidence.LARGE_VARIANT_POINTS_LONG_NAME, LARGE_VARIANT_POINTS)
                .add(AggregateDepthEvidence.LARGE_VARIANT_WINDOW_LONG_NAME, LARGE_VARIANT_WINDOW)
                .add(AggregateDepthEvidence.NUM_BINS_LONG_NAME, NUM_BINS)
                .add(AggregateDepthEvidence.MAX_QUALITY_LONG_NAME, MAX_QUAL)
                .add(AggregateDepthEvidence.POWER_THRESHOLD_LONG_NAME, POWER_THRESHOLD);

        runCommandLine(args, AggregateDepthEvidence.class.getSimpleName());

        final File expectedFile = new File(TOOL_TEST_DIR, "aggregate_depth_test.output.vcf.gz");
        final Pair<VCFHeader, List<VariantContext>> expected = VariantContextTestUtils.readEntireVCFIntoMemory(expectedFile.getPath());
        final Pair<VCFHeader, List<VariantContext>> output = VariantContextTestUtils.readEntireVCFIntoMemory(outputFile);

        Assert.assertEquals(expected.getValue().size(), output.getValue().size());
        final Iterator<VariantContext> expectedIterator = expected.getValue().iterator();
        final Iterator<VariantContext> outputIterator = output.getValue().iterator();
        while (expectedIterator.hasNext()) {
            final VariantContext expectedVariant = expectedIterator.next();
            final VariantContext variant = outputIterator.next();
            VariantContextTestUtils.assertVariantContextsAreEqual(variant, expectedVariant, Collections.emptyList(), attributesWithJitter);
            for (final String attribute : attributesWithJitter) {
                final double expectedValue = expectedVariant.getAttributeAsDouble(attribute, 0);
                if (Double.isNaN(expectedValue)) {
                    Assert.assertTrue(Double.isNaN(variant.getAttributeAsDouble(attribute, 0)));
                } else {
                    final double actualValue = variant.getAttributeAsDouble(attribute, 0);
                    Assert.assertTrue(Math.abs(expectedValue - actualValue) < FLOAT_TOLERANCE);
                }
            }
        }
    }
}