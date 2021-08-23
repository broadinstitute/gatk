package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.engine.GatkToolIntegrationTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

public class AggregateSVEvidenceTest extends GatkToolIntegrationTest {

    public static final String TOOL_TEST_DIR = toolsTestDir + "walkers/sv/AggregatePairedEndAndSplitReadEvidence";

    public static final int PE_INNER_WINDOW = 50;
    public static final int PE_OUTER_WINDOW = 500;
    public static final int SR_WINDOW = 200;
    public static final int SR_CROSSOVER = 20;

    @Test
    public void test() {
        final File vcfFile = new File(TOOL_TEST_DIR, "test_hg38.vcf.gz");
        final File peFile = new File(TOOL_TEST_DIR, "test_hg38.pe.txt.gz");
        final File srFile = new File(TOOL_TEST_DIR, "test_hg38.sr.txt.gz");
        final File coverageFile = new File(TOOL_TEST_DIR, "test_hg38.sample_coverage.tsv");
        final File outputFile = createTempFile("aggregated", ".vcf.gz");

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addVCF(vcfFile)
                .addOutput(outputFile)
                .addReference(hg38Reference)
                .add(AggregateSVEvidence.DISCORDANT_PAIRS_LONG_NAME, peFile)
                .add(AggregateSVEvidence.SPLIT_READ_LONG_NAME, srFile)
                .add(AggregateSVEvidence.SAMPLE_COVERAGE_LONG_NAME, coverageFile)
                .add(AggregateSVEvidence.PE_INNER_WINDOW_LONG_NAME, PE_INNER_WINDOW)
                .add(AggregateSVEvidence.PE_OUTER_WINDOW_LONG_NAME, PE_OUTER_WINDOW)
                .add(AggregateSVEvidence.SR_WINDOW_LONG_NAME, SR_WINDOW)
                .add(AggregateSVEvidence.SR_CROSSOVER_LONG_NAME, SR_CROSSOVER);

        runCommandLine(args, AggregateSVEvidence.class.getSimpleName());

        final File expectedFile = new File(TOOL_TEST_DIR, "test_hg38.expected.vcf.gz");
        final Pair<VCFHeader, List<VariantContext>> expected = VariantContextTestUtils.readEntireVCFIntoMemory(expectedFile.getPath());
        final Pair<VCFHeader, List<VariantContext>> output = VariantContextTestUtils.readEntireVCFIntoMemory(outputFile.getPath());

        Assert.assertEquals(expected.getValue().size(), output.getValue().size());
        final Iterator<VariantContext> expectedIterator = expected.getValue().iterator();
        final Iterator<VariantContext> outputIterator = output.getValue().iterator();
        while (expectedIterator.hasNext()) {
            final VariantContext expectedVariant = expectedIterator.next();
            final VariantContext variant = outputIterator.next();
            VariantContextTestUtils.assertVariantContextsAreEqual(variant, expectedVariant, Collections.emptyList(), Collections.emptyList());
        }
    }
}