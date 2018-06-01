package org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific;

import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.tools.walkers.annotator.Annotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Created by emeryj on 8/11/17.
 */

public class AS_RMSMappingQualityUnitTest extends ReducibleAnnotationBaseTest {

    @Override
    protected List<Annotation> getAnnotationsToUse() {
        return Collections.singletonList(new AS_RMSMappingQuality());
    }

    @Override
    protected String getRawKey() {
        return GATKVCFConstants.AS_RAW_RMS_MAPPING_QUALITY_KEY;
    }

    @Override
    protected String getKey() {
        return GATKVCFConstants.AS_RMS_MAPPING_QUALITY_KEY;
    }

    @Test
    public void testFinalizeAnnotations() throws Exception {
        final FeatureInput<VariantContext> dbSNPBinding = null;
        final List<FeatureInput<VariantContext>> features = Collections.emptyList();
        final VariantAnnotatorEngine vae = new VariantAnnotatorEngine(Collections.singletonList(new AS_RMSMappingQuality()), dbSNPBinding, features, false);
        final Allele refAllele = Allele.create("A", true);
        final Allele altAllele1 = Allele.create("T");
        final Allele altAllele2 = Allele.create("C");
        final Genotype genotype = new GenotypeBuilder("sample2", Arrays.asList(refAllele, altAllele1, altAllele2))
                .AD(new int[]{2,80,9}).make();

        final VariantContext vc = new VariantContextBuilder(new VariantContextBuilder())
                .alleles(Arrays.asList(refAllele, altAllele1, altAllele2))
                .chr("1").start(15L).stop(15L)
                .attribute(GATKVCFConstants.AS_RAW_RMS_MAPPING_QUALITY_KEY, "400.00|285.00|385.00")
                .genotypes(genotype)
                .make();
        final VariantContext result = vae.finalizeAnnotations(vc, vc);
        Assert.assertNull(result.getAttribute(GATKVCFConstants.AS_RAW_RMS_MAPPING_QUALITY_KEY));
        Assert.assertNotNull(result.getAttribute(GATKVCFConstants.AS_RMS_MAPPING_QUALITY_KEY));
        Assert.assertEquals(result.getAttribute(GATKVCFConstants.AS_RMS_MAPPING_QUALITY_KEY),String.format("%.2f",Math.sqrt(285.00/80)) + "," + String.format("%.2f",Math.sqrt(385.00/9)));
    }

}