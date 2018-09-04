package org.broadinstitute.hellbender.tools.funcotator.metadata;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Collections;
import java.util.List;

public class SamplePairExtractorUnitTest extends GATKBaseTest {
    private static final String EXAMPLE_TUMOR_SAMPLE_NAME = toolsTestDir + "funcotator/NA12878.ob_filtered.vcf";
    private static final String EXAMPLE_M2_TUMOR_SAMPLE_NAME = toolsTestDir + "funcotator/M2_output_vcf_sample_name_test.vcf";
    private static final String EXAMPLE_M2_2_TUMOR_SAMPLE_NAME = toolsTestDir + "funcotator/SM-74P4M-filtered_sample_name_test.vcf";
    private static final String EXAMPLE_GERMLINE_NAME = toolsTestDir + "funcotator/0816201804HC0_R01C01.pik3ca.vcf";

    @DataProvider
    public Object[][] provideHeaders() {
        return new Object[][] {
                {EXAMPLE_TUMOR_SAMPLE_NAME, Collections.singletonList(new TumorNormalPair("TUMOR", "NORMAL"))},
                {EXAMPLE_M2_TUMOR_SAMPLE_NAME,  Collections.singletonList(new TumorNormalPair("S1234T", "S1234N"))},
                {EXAMPLE_M2_2_TUMOR_SAMPLE_NAME,  Collections.singletonList(new TumorNormalPair("SM-74P4M", "SM-74NEG"))},
                {EXAMPLE_GERMLINE_NAME,  Collections.singletonList(new TumorNormalPair("0816201804HC0_R01C01", ""))}
        };
    }

    @Test(dataProvider = "provideHeaders")
    public void testExtractTumorNormalPairs(final String vcfFilename, final List<Pair<String, String>> gtPairs) {
        final Pair<VCFHeader, List<VariantContext>> samples = VariantContextTestUtils.readEntireVCFIntoMemory(vcfFilename);
        final List<TumorNormalPair> pairs = SamplePairExtractor.extractPossibleTumorNormalPairs(samples.getLeft());
        Assert.assertEquals(pairs, gtPairs);
    }
}
