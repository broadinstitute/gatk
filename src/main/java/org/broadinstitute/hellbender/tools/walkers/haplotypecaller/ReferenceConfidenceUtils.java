package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.apache.commons.collections4.ListUtils;

import java.util.Arrays;
import java.util.List;

public class ReferenceConfidenceUtils {

    public static VariantContext addNonRefSymbolicAllele(final VariantContext mergedVC) {
        final List<Allele> alleleList = ListUtils.union(mergedVC.getAlleles(), Arrays.asList(Allele.NON_REF_ALLELE));
        return new VariantContextBuilder(mergedVC).alleles(alleleList).make();
    }
}
