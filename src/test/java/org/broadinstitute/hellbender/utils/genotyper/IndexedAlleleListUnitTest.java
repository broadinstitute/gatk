package org.broadinstitute.hellbender.utils.genotyper;

import htsjdk.variant.variantcontext.Allele;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.List;

import static org.broadinstitute.hellbender.utils.genotyper.AlleleListUnitTester.assertAlleleList;

/**
 * Tests {@link org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleListUnitTest}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class IndexedAlleleListUnitTest {

    @Test
    public void testEmptyConstructor() {
        final AlleleList<Allele> subject = new IndexedAlleleList<>();
        assertAlleleList(subject, Collections.<Allele>emptyList());
    }

    @Test(dataProvider= "alleleCountMaxAlleleLengthData")
    public void testArrayConstructor(final int alleleCount, final int maxAlleleLength) {
        final Allele[] alleles = AlleleListUnitTester.generateRandomAlleles(alleleCount, maxAlleleLength);

        final LinkedHashSet<Allele> nonRepeatedAlleles = new LinkedHashSet<>(Arrays.asList(alleles));
        final IndexedAlleleList<Allele> subject = new IndexedAlleleList<>(alleles);
        assertAlleleList(subject, Arrays.asList(nonRepeatedAlleles.toArray(new Allele[nonRepeatedAlleles.size()])));
    }

    @Test(dataProvider= "alleleCountMaxAlleleLengthData")
    public void testCollectionConstructor(final int alleleCount, final int maxAlleleLength) {
        final Allele[] alleles = AlleleListUnitTester.generateRandomAlleles(alleleCount, maxAlleleLength);

        final List<Allele> alleleList = Arrays.asList(alleles);
        final LinkedHashSet<Allele> nonRepeatedAlleles = new LinkedHashSet<>(Arrays.asList(alleles));
        final IndexedAlleleList<Allele> subject = new IndexedAlleleList<>(alleleList);
        assertAlleleList(subject, Arrays.asList(nonRepeatedAlleles.toArray(new Allele[nonRepeatedAlleles.size()])));
    }

    private static final int[] SAMPLE_COUNT = { 0, 1, 5, 10, 20};

    private static final int[] MAX_ALLELE_LENGTH = { 1, 2, 3, 10 };

    @DataProvider(name="alleleCountMaxAlleleLengthData")
    public Object[][] alleleCountMaxAlleleLengthData() {
        final Object[][] result = new Object[SAMPLE_COUNT.length * MAX_ALLELE_LENGTH.length][];
        int nextIndex = 0;
        for (int i = 0; i < SAMPLE_COUNT.length; i++)
            for (int j = 0; j < MAX_ALLELE_LENGTH.length; j++)
                result[nextIndex++] = new Object[] { SAMPLE_COUNT[i], MAX_ALLELE_LENGTH[j]};
        return result;
    }
}
