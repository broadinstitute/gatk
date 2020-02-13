package org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.tools.walkers.annotator.Annotation;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Created by emeryj on 8/11/17.
 */
public class AS_QualByDepthUnitTest extends ReducibleAnnotationBaseTest {
    @Override
    protected List<Annotation> getAnnotationsToUse() {
        return Collections.singletonList(new AS_QualByDepth());
    }

    @Override
    protected String getRawKey() {
        return GATKVCFConstants.AS_QUAL_KEY;
    }

    @Override
    protected String getKey() {
        return GATKVCFConstants.AS_QUAL_KEY;
    }

    @Test
    public void testParseQualList() {
        final int NON_REF_INDEX = 1; //finalized QD values won't have an entry for ref
        final String goodQualList = "|234|0"; //combined VCs from GenomicsDB have zero value for NON-REF
        final String trickyQualList = "|234|"; //older single-sample GVCFs don't have a value for NON-REF -- GenomicsDB assigns an empty value, which is fair
        final VariantContext withNonRefValue = new VariantContextBuilder(GATKVariantContextUtils.makeFromAlleles("good", "chr1", 10001, Arrays.asList("A","T", Allele.NON_REF_STRING))).
                attribute(GATKVCFConstants.AS_RAW_QUAL_APPROX_KEY, goodQualList).make();
        final VariantContext withNoNonRefValue = new VariantContextBuilder(GATKVariantContextUtils.makeFromAlleles("good", "chr1", 10001, Arrays.asList("A","T", Allele.NON_REF_STRING))).
                attribute(GATKVCFConstants.AS_RAW_QUAL_APPROX_KEY, trickyQualList).make();
        Assert.assertEquals(AS_QualByDepth.parseQualList(withNonRefValue).size(), withNonRefValue.getAlternateAlleles().size());
        final List<Integer> qualsFromNoNonRefList = AS_QualByDepth.parseQualList(withNoNonRefValue);
        Assert.assertEquals(qualsFromNoNonRefList.size(), withNonRefValue.getAlternateAlleles().size());
        Assert.assertEquals((int)qualsFromNoNonRefList.get(NON_REF_INDEX), 0);  //int cast because Integer results in ambiguous method call for assert
    }

}