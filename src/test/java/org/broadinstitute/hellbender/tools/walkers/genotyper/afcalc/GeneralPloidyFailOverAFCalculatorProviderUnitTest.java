package org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc;

import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeCalculationArgumentCollection;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

/**
 * Tests {@link GeneralPloidyFailOverAFCalculatorProvider}
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class GeneralPloidyFailOverAFCalculatorProviderUnitTest extends GATKBaseTest {

    private final static int[] PLOIDIES = new int[] { AFCalculatorImplementation.UNBOUND_PLOIDY,1,2,3,4,10 };
    private final static int[] MAX_ALT_ALLELES = new int[] { AFCalculatorImplementation.UNBOUND_ALTERNATIVE_ALLELE_COUNT,1,2,3,4,10};

    @Test(dataProvider= "getMatrixOfPlodiesAndMaxAltAlleles")
    public void testAFCalculatorProvider(final int ploidy, final int maxAltAlleles) {
        final GenotypeCalculationArgumentCollection args = new GenotypeCalculationArgumentCollection();
        args.MAX_ALTERNATE_ALLELES = maxAltAlleles;
        args.samplePloidy = ploidy;

        final GeneralPloidyFailOverAFCalculatorProvider provider = new GeneralPloidyFailOverAFCalculatorProvider(args);

        final AFCalculator calculator = provider.getInstance(ploidy,maxAltAlleles);
        Assert.assertNotNull(calculator);
        final AFCalculatorImplementation implementation = AFCalculatorImplementation.fromCalculatorClass(calculator.getClass());
        Assert.assertTrue(implementation.usableForParams(ploidy,maxAltAlleles));
        for (final int PLOIDY : PLOIDIES) {
            for (final int MAX_ALT_ALLELE : MAX_ALT_ALLELES) {
                if (implementation.usableForParams(PLOIDY, MAX_ALT_ALLELE)) {
                    Assert.assertSame(provider.getInstance(PLOIDY, MAX_ALT_ALLELE), calculator);
                } else {
                    final AFCalculator failOver = provider.getInstance(PLOIDY, MAX_ALT_ALLELE);
                    Assert.assertNotNull(failOver);
                    final AFCalculatorImplementation failOverImplementation = AFCalculatorImplementation.fromCalculatorClass(failOver.getClass());
                    Assert.assertTrue(failOverImplementation.usableForParams(PLOIDY, MAX_ALT_ALLELE));
                    Assert.assertEquals(failOverImplementation, AFCalculatorImplementation.EXACT_GENERAL_PLOIDY);
                }
            }
        }
    }

    @DataProvider(name="getMatrixOfPlodiesAndMaxAltAlleles")
    public Object[][] getMatrixOfPlodiesAndMaxAltAlleles() {
        final Object[][] result = new Object[PLOIDIES.length * MAX_ALT_ALLELES.length][];
        int idx = 0;
        for (final int PLOIDY : PLOIDIES) {
            for (final int MAX_ALT_ALLELE : MAX_ALT_ALLELES) {
                result[idx++] = new Object[]{PLOIDY, MAX_ALT_ALLELE};
            }
        }
        return result;
    }

}
