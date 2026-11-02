package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class PermutectDatasetEngineUnitTest {

    // a variant on contig 20 with the given start, ref allele and alt alleles
    private static VariantContext variant(final int start, final String ref, final String... alts) {
        final List<Allele> alleles = new ArrayList<>();
        alleles.add(Allele.create(ref, true));
        Arrays.stream(alts).map(alt -> Allele.create(alt, false)).forEach(alleles::add);
        return new VariantContextBuilder("test", "20", start, start + ref.length() - 1, alleles).make();
    }

    private static VariantContext filtered(final VariantContext vc) {
        return new VariantContextBuilder(vc).filter("FILTERED").make();
    }

    @DataProvider
    public Object[][] truthMatchingData() {
        final VariantContext mnv = variant(100, "AT", "GC");
        final VariantContext firstSnv = variant(100, "A", "G");
        final VariantContext secondSnv = variant(101, "T", "C");

        // {call, truth records overlapping the call, whether each alt allele of the call matches the truth}
        return new Object[][] {
                // alleles match the same allele in the truth, after padding to a common ref allele
                {firstSnv, List.of(firstSnv), List.of(true)},
                {firstSnv, List.of(variant(100, "A", "C")), List.of(false)},
                {firstSnv, List.of(variant(100, "AT", "GT")), List.of(true)},
                {mnv, List.of(mnv), List.of(true)},
                {variant(100, "ATC", "A"), List.of(variant(100, "ATC", "A")), List.of(true)},

                // an MNV matches truth that splits it into SNVs, if the truth has each of its substituted bases, unfiltered
                {mnv, List.of(firstSnv, secondSnv), List.of(true)},
                {mnv, List.of(firstSnv), List.of(false)},
                {mnv, List.of(secondSnv), List.of(false)},
                {mnv, List.of(firstSnv, filtered(secondSnv)), List.of(false)},
                {mnv, List.of(firstSnv, variant(101, "T", "A")), List.of(false)},
                {variant(100, "ACG", "TCA"), List.of(variant(100, "A", "T"), variant(102, "G", "A")), List.of(true)},

                // an SNV matches only the same SNV, not an MNV that contains it or truth starting before it
                {firstSnv, List.of(mnv), List.of(false)},
                {secondSnv, List.of(mnv), List.of(false)},
                {firstSnv, List.of(variant(99, "CA", "C")), List.of(false)},

                // at a multiallelic site, an MNV and an SNV padded to the MNV's ref allele each match truth SNVs
                {variant(100, "AT", "GC", "GT"), List.of(firstSnv, secondSnv), List.of(true, true)},
                {variant(100, "AT", "GC", "GT"), List.of(firstSnv), List.of(false, true)}
        };
    }

    @Test(dataProvider = "truthMatchingData")
    public void testAltAllelesMatchingTruth(final VariantContext call, final List<VariantContext> truthVCs, final List<Boolean> expected) {
        Assert.assertEquals(PermutectDatasetEngine.altAllelesMatchingTruth(call, truthVCs), expected);
    }
}
