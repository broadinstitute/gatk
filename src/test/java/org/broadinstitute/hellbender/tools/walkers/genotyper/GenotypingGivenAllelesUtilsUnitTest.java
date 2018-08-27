package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class GenotypingGivenAllelesUtilsUnitTest {

    @Test
    public void testComposeGivenAllelesVariantContextFromVariantList() {

        final SimpleInterval loc = new SimpleInterval("20", 68, 68);

        final List<Allele> sameLocDelAlleles1 = Arrays.asList(Allele.create("GT", true), Allele.create("G"));
        final VariantContext sameLocDelVc1 = new VariantContextBuilder("a", "20", 68, 69, sameLocDelAlleles1).make();

        final List<Allele> sameLocSnpAlleles1 = Arrays.asList(Allele.create("G", true), Allele.create("C"));
        final VariantContext sameLocSnpVc1 = new VariantContextBuilder("a", "20", 68, 68, sameLocSnpAlleles1).make();

        final List<Allele> spanningDelAlleles1 = Arrays.asList(Allele.create("ATGTA", true), Allele.create("A"));
        final VariantContext spanningDelVc1 = new VariantContextBuilder("a", "20", 66, 70, spanningDelAlleles1).make();

        final List<Allele> expectedAlleles = Arrays.asList(Allele.create("GT", true), Allele.create("G"), Allele.create("CT"));
        final VariantContext expectedVC = new VariantContextBuilder("a", "20", 68, 69, expectedAlleles).make();

        final VariantContext variantContext =
                GenotypingGivenAllelesUtils.composeGivenAllelesVariantContextFromVariantList(Arrays.asList(spanningDelVc1, sameLocDelVc1, sameLocSnpVc1),
                        loc,
                        true);

        VariantContextTestUtils.assertVariantContextsAreEqual(variantContext, expectedVC, Collections.emptyList());
    }
}