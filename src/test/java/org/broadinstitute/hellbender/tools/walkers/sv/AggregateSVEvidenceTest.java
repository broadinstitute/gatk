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

public class AggregateSVEvidenceTest extends GatkToolIntegrationTest {

    // If true, update the expected outputs in tests that assert an exact or approximate match vs. prior output,
    // instead of actually running the tests.
    public static final boolean UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS = false;

    public static final String TOOL_TEST_DIR = toolsTestDir + "walkers/sv/AggregateSVEvidence";
    public static final String LARGE_FILE_TEST_DIR = largeFileTestDir + "sv";

    public static final double FLOAT_TOLERANCE = 1e-4;

    public static final int PE_INNER_WINDOW = 50;
    public static final int PE_OUTER_WINDOW = 500;
    public static final int SR_WINDOW = 200;
    public static final int SR_CROSSOVER = 20;

    public static final List<String> attributesWithJitter = Arrays.asList(
            GATKSVVCFConstants.FIRST_SPLIT_QUALITY_ATTRIBUTE,
            GATKSVVCFConstants.FIRST_SPLIT_CARRIER_SIGNAL_ATTRIBUTE,
            GATKSVVCFConstants.SECOND_SPLIT_QUALITY_ATTRIBUTE,
            GATKSVVCFConstants.SECOND_SPLIT_CARRIER_SIGNAL_ATTRIBUTE,
            GATKSVVCFConstants.TOTAL_SPLIT_QUALITY_ATTRIBUTE,
            GATKSVVCFConstants.TOTAL_SPLIT_CARRIER_SIGNAL_ATTRIBUTE,
            GATKSVVCFConstants.DISCORDANT_PAIR_QUALITY_ATTRIBUTE,
            GATKSVVCFConstants.DISCORDANT_PAIR_CARRIER_SIGNAL_ATTRIBUTE,
            GATKSVVCFConstants.PESR_QUALITY_ATTRIBUTE,
            GATKSVVCFConstants.PESR_CARRIER_SIGNAL_ATTRIBUTE,
            GATKSVVCFConstants.BAF_HET_RATIO_ATTRIBUTE,
            GATKSVVCFConstants.BAF_DEL_LOGLIK_ATTRIBUTE,
            GATKSVVCFConstants.BAF_KS_STAT_ATTRIBUTE,
            GATKSVVCFConstants.BAF_KS_Q_ATTRIBUTE
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
        final File vcfFile = new File(TOOL_TEST_DIR, "aggregate_sv_test.vcf.gz");
        final File peFile = new File(LARGE_FILE_TEST_DIR, "aggregate_sv_test.pe.txt.gz");
        final File srFile = new File(LARGE_FILE_TEST_DIR, "aggregate_sv_test.sr.txt.gz");
        final File bafFile = new File(LARGE_FILE_TEST_DIR, "aggregate_sv_test.baf.txt.gz");
        final File coverageFile = new File(TOOL_TEST_DIR, "aggregate_sv_test.medianCov.tsv");
        final File ploidyTable = new File(TOOL_TEST_DIR, "1kg_ref_panel_v1.ploidy_table.tsv");
        final String outputFile = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? TOOL_TEST_DIR + "/aggregate_sv_test.expected.vcf.gz" : createTempFile("aggregated", ".vcf.gz").getAbsolutePath();

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addVCF(vcfFile)
                .addOutput(outputFile)
                .addReference(hg38Reference)
                .add(AggregateSVEvidence.DISCORDANT_PAIRS_LONG_NAME, peFile)
                .add(AggregateSVEvidence.SPLIT_READ_LONG_NAME, srFile)
                .add(AggregateSVEvidence.BAF_LONG_NAME, bafFile)
                .add(AggregateSVEvidence.MEDIAN_COVERAGE_LONG_NAME, coverageFile)
                .add(SVCluster.PLOIDY_TABLE_LONG_NAME, ploidyTable)
                .add(AggregateSVEvidence.PE_INNER_WINDOW_LONG_NAME, PE_INNER_WINDOW)
                .add(AggregateSVEvidence.PE_OUTER_WINDOW_LONG_NAME, PE_OUTER_WINDOW)
                .add(AggregateSVEvidence.SR_WINDOW_LONG_NAME, SR_WINDOW)
                .add(AggregateSVEvidence.SR_CROSSOVER_LONG_NAME, SR_CROSSOVER);

        runCommandLine(args, AggregateSVEvidence.class.getSimpleName());

        final File expectedFile = new File(TOOL_TEST_DIR, "aggregate_sv_test.expected.vcf.gz");
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
                    Assert.assertTrue(Math.abs(expectedValue - variant.getAttributeAsDouble(attribute, 0)) < FLOAT_TOLERANCE);
                }
            }
        }
    }
}