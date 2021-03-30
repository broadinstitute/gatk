package org.broadinstitute.hellbender.tools.variantdb.nextgen;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.List;

public class ExtractCohortRemoveAnnotationsUnitTest extends GATKBaseTest{
    private static final String ALLELE_SPECIFIC_DIRECTORY = toolsTestDir + "walkers/annotator/allelespecific";
    private static final File ORIGINAL_TEST_FILE =  new File(ALLELE_SPECIFIC_DIRECTORY + "/GenotypeGVCFs.output.vcf");

    @Test
    public void testRemoveAnnotations() {
        VCFHeader testVCFHeader = VariantContextTestUtils.getVCFHeader(ORIGINAL_TEST_FILE.getPath());
        ExtractCohortEngine engine = new ExtractCohortEngine(
                null,
                null,
                testVCFHeader,
                null,
                null,
                null,
                null,
                "spec-ops-aou.kc_high_cov_ccdg.exported_cohort_100_test",
                null,
                null,
                null,
                0,
                false,
                0.0,
                0.0,
                null,
                null,
                null,
                false);
        List<VariantContext> variantContexts = VariantContextTestUtils.getVariantContexts(ORIGINAL_TEST_FILE); // list variantContexts from VCF file
        VariantContext expectedVC = engine.removeAnnotations(variantContexts.get(0)); // single variantContext -- with annotations removed
        Assert.assertFalse(expectedVC.hasAttribute(GATKVCFConstants.FISHER_STRAND_KEY), "removeAnnotations did not remove FS annotation.");
        Assert.assertFalse(expectedVC.hasAttribute(GATKVCFConstants.AS_QUAL_BY_DEPTH_KEY), "removeAnnotations did not remove AS_QD annotation.");
        Assert.assertFalse(expectedVC.hasAttribute(GATKVCFConstants.STRAND_ODDS_RATIO_KEY), "removeAnnotations did not remove SOR annotation.");
    }
}
