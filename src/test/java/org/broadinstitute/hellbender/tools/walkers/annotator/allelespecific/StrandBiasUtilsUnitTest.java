package org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific;

import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.*;

public class StrandBiasUtilsUnitTest extends GATKBaseTest {

    private final String sample1 = "NA1";
    private static final String CONTIG = "1";

    private static final Allele REF = Allele.create("T", true);
    private static final Allele ALT = Allele.create("A", false);

    private VariantContext makeVC(final long position) {
        final GenotypesContext testGC = GenotypesContext.create(2);
        // sample1 -> A/T with GQ 30
        testGC.add(new GenotypeBuilder(sample1).alleles(Arrays.asList(REF, ALT)).GQ(30).make());

        return (new VariantContextBuilder())
                .alleles(Arrays.asList(REF, ALT)).chr(CONTIG).start(position).stop(position).genotypes(testGC).make();
    }

    /**
     * Test for issue #6766
     */
    @Test
    public void testMakeRawAnnotationStringsNullPointerExceptionIssue6766() {
        VariantContext vc = makeVC(1);
        Map<Allele, List<Integer>> perAlleleValues = new LinkedHashMap<>();
        perAlleleValues.put(REF, Arrays.asList(2, 5));
        perAlleleValues.put(ALT, null);
        try {
            StrandBiasUtils.makeRawAnnotationString(vc.getAlleles(), perAlleleValues);
        } catch (NullPointerException npe) {
            Assert.fail("This code should not be throwing a NullPointerException");
        }
    }
}
