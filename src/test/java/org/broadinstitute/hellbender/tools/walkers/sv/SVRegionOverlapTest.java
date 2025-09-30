package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.engine.GatkToolIntegrationTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.IntervalMergingRule;
import org.broadinstitute.hellbender.utils.IntervalSetRule;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

public class SVRegionOverlapTest extends GatkToolIntegrationTest {

    // If true, update the expected outputs in tests that assert an exact or approximate match vs. prior output,
    // instead of actually running the tests.
    public static final boolean UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS = false;

    public static final double FLOAT_TOLERANCE = 1e-1;

    // borrow from SVStratify
    public static final String TOOL_TEST_DIR = toolsTestDir + "walkers/sv/SVStratify/";
    public static final String INPUT_VCF_FILE = TOOL_TEST_DIR + "bwa_melt.chr22.vcf.gz";
    public static final String SEGDUP_REGION_FILE = TOOL_TEST_DIR + "hg38.SegDup.chr22.bed";
    public static final String REPEATMASKER_REGION_FILE = TOOL_TEST_DIR + "hg38.RM.chr22_subsampled.bed";

    public static final String SEGDUP_REGION_NAME = "SD";
    public static final String REPEATMASKER_REGION_NAME = "RM";

    public static final IntervalSetRule REGIONS_SET_RULE = IntervalSetRule.UNION;
    public static final IntervalMergingRule REGIONS_MERGING_RULE = IntervalMergingRule.OVERLAPPING_ONLY;
    public static final int REGION_PADDING = 0;
    public static final boolean SUPPRESS_ENDPOINT_COUNTS = false;
    public static final boolean SUPPRESS_OVERLAP_FRACTION = false;

    public static final List<String> attributesWithJitter = Arrays.asList(
            GATKSVVCFConstants.NUM_END_OVERLAPS_INFO_BASE + SEGDUP_REGION_NAME,
            GATKSVVCFConstants.OVERLAP_FRACTION_INFO_BASE + SEGDUP_REGION_NAME,
            GATKSVVCFConstants.NUM_END_OVERLAPS_INFO_BASE + REPEATMASKER_REGION_NAME,
            GATKSVVCFConstants.OVERLAP_FRACTION_INFO_BASE + REPEATMASKER_REGION_NAME
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
        final String outputFile = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? getToolTestDataDir() + "/sv_region_overlap_test.output.vcf.gz" : createTempFile("overlap", ".vcf.gz").getAbsolutePath();

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addVCF(INPUT_VCF_FILE)
                .addOutput(outputFile)
                .addReference(hg38Reference)
                .add(SVRegionOverlap.REGIONS_NAME_LONG_NAME, SEGDUP_REGION_NAME)
                .add(SVRegionOverlap.REGIONS_FILE_LONG_NAME, SEGDUP_REGION_FILE)
                .add(SVRegionOverlap.REGIONS_NAME_LONG_NAME, REPEATMASKER_REGION_NAME)
                .add(SVRegionOverlap.REGIONS_FILE_LONG_NAME, REPEATMASKER_REGION_FILE)
                .add(SVRegionOverlap.REGIONS_SET_RULE_LONG_NAME, REGIONS_SET_RULE)
                .add(SVRegionOverlap.REGIONS_MERGING_RULE_LONG_NAME, REGIONS_MERGING_RULE)
                .add(SVRegionOverlap.REGION_PADDING_LONG_NAME, REGION_PADDING)
                .add(SVRegionOverlap.SUPPRESS_ENDPOINT_COUNTS_LONG_NAME, SUPPRESS_ENDPOINT_COUNTS)
                .add(SVRegionOverlap.SUPPRESS_OVERLAP_FRACTION_LONG_NAME, SUPPRESS_OVERLAP_FRACTION);

        runCommandLine(args, SVRegionOverlap.class.getSimpleName());

        final File expectedFile = new File(getToolTestDataDir(), "sv_region_overlap_test.output.vcf.gz");
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