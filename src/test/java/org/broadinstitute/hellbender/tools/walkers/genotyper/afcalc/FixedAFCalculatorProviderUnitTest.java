package org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc;

import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeCalculationArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.genotyper.StandardCallerArgumentCollection;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;

/**
 * Tests {@link org.broadinstitute.gatk.tools.walkers.genotyper.afcalc.FixedAFCalculatorProvider}
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class FixedAFCalculatorProviderUnitTest {

    @Test(dataProvider="nonThreadSafeConstructorsData")
    public void testNonThreadSafeConstructors(final int ploidy, final int maxAltAlleles, final AFCalculatorImplementation preferred) {
        final GenotypeCalculationArgumentCollection args = new GenotypeCalculationArgumentCollection();
        args.MAX_ALTERNATE_ALLELES = maxAltAlleles;
        args.samplePloidy = ploidy;
        final StandardCallerArgumentCollection callerArgs = new StandardCallerArgumentCollection();
        if (preferred != null ) {
            callerArgs.requestedAlleleFrequencyCalculationModel = preferred;
        }
        callerArgs.genotypeArgs = args;
        final FixedAFCalculatorProvider providerCallerArgs = new FixedAFCalculatorProvider(callerArgs, true);
        final FixedAFCalculatorProvider providerCallerArgsNoVerify = new FixedAFCalculatorProvider(callerArgs, false);
        final FixedAFCalculatorProvider providerGenotypingArgs = new FixedAFCalculatorProvider(args, true);

        Assert.assertNotNull(providerCallerArgs.getInstance(ploidy, maxAltAlleles));
        Assert.assertNotNull(providerCallerArgsNoVerify.getInstance(ploidy, maxAltAlleles));
        Assert.assertTrue(AFCalculatorImplementation.fromCalculatorClass(providerCallerArgs.getInstance(ploidy, maxAltAlleles).getClass()).usableForParams(ploidy, maxAltAlleles));
        Assert.assertTrue(AFCalculatorImplementation.fromCalculatorClass(providerCallerArgsNoVerify.getInstance(ploidy, maxAltAlleles).getClass()).usableForParams(ploidy, maxAltAlleles));
        Assert.assertNotNull(providerGenotypingArgs.getInstance(ploidy, maxAltAlleles));
        Assert.assertTrue(AFCalculatorImplementation.fromCalculatorClass(providerGenotypingArgs.getInstance(ploidy, maxAltAlleles).getClass()).usableForParams(ploidy, maxAltAlleles));

        final VariantContext vc= new VariantContextBuilder().chr("chr1").alleles("A", "T").make();
        Assert.assertEquals(providerCallerArgs.getInstance(vc, ploidy, maxAltAlleles), providerCallerArgs.getInstance(ploidy, maxAltAlleles));//equal because there's no samples in vc
        Assert.assertEquals(providerCallerArgsNoVerify.getInstance(vc, ploidy, maxAltAlleles), providerCallerArgsNoVerify.getInstance(ploidy, maxAltAlleles));//equal because there's no samples in vc

        if (preferred != null && preferred.usableForParams(ploidy,maxAltAlleles)) {
            Assert.assertEquals(AFCalculatorImplementation.fromCalculatorClass(providerCallerArgs.getInstance(ploidy, maxAltAlleles).getClass()), preferred);
            Assert.assertEquals(AFCalculatorImplementation.fromCalculatorClass(providerCallerArgsNoVerify.getInstance(ploidy, maxAltAlleles).getClass()), preferred);
        }
    }


    private static final int[] PLOIDIES = { 1,2,3,4,10 };
    private static final int[] MAX_ALT_ALLELES = { 1,2,3,4,10};

    @DataProvider(name="nonThreadSafeConstructorsData")
    public Object[][] nonThreadSafeConstructorsData() {
        final Object[][] result = new Object[PLOIDIES.length * MAX_ALT_ALLELES.length * (AFCalculatorImplementation.values().length + 1)][];
        int idx = 0;
        for (int i = 0; i < PLOIDIES.length; i++) {
            for (int j = 0; j < MAX_ALT_ALLELES.length; j++) {
                result[idx++] = new Object[] { PLOIDIES[i], MAX_ALT_ALLELES[j], null };
                for (final AFCalculatorImplementation impl : AFCalculatorImplementation.values()) {
                    result[idx++] = new Object[]{PLOIDIES[i], MAX_ALT_ALLELES[j], impl};
                }
            }
        }
        return result;
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testPloidyError() throws Exception {
        final GenotypeCalculationArgumentCollection args = new GenotypeCalculationArgumentCollection();
        args.samplePloidy = -2;
        final StandardCallerArgumentCollection callerArgs = new StandardCallerArgumentCollection();
        callerArgs.genotypeArgs = args;

        new FixedAFCalculatorProvider(callerArgs, false);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testMaxAltAllelesError() throws Exception {
        final GenotypeCalculationArgumentCollection args = new GenotypeCalculationArgumentCollection();
        args.samplePloidy = 2;
        args.MAX_ALTERNATE_ALLELES = -2;
        final StandardCallerArgumentCollection callerArgs = new StandardCallerArgumentCollection();
        callerArgs.genotypeArgs = args;

        new FixedAFCalculatorProvider(callerArgs, false);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testTesInstanceInvalidPloidyError() throws Exception {
        final GenotypeCalculationArgumentCollection args = new GenotypeCalculationArgumentCollection();
        args.samplePloidy = 2;
        args.MAX_ALTERNATE_ALLELES = 2;
        final StandardCallerArgumentCollection callerArgs = new StandardCallerArgumentCollection();
        callerArgs.genotypeArgs = args;

        final FixedAFCalculatorProvider p = new FixedAFCalculatorProvider(callerArgs, true);
        p.getInstance(5, 2);
    }

    @Test
    public void testTesInstanceInvalidPloidyError_noVerify() throws Exception {
        final GenotypeCalculationArgumentCollection args = new GenotypeCalculationArgumentCollection();
        args.samplePloidy = 2;
        args.MAX_ALTERNATE_ALLELES = 2;
        final StandardCallerArgumentCollection callerArgs = new StandardCallerArgumentCollection();
        callerArgs.genotypeArgs = args;

        final FixedAFCalculatorProvider p = new FixedAFCalculatorProvider(callerArgs, false);
        p.getInstance(5, 2); //this passes - no validation
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testTesInstanceInvalidAlleleNumberError() throws Exception {
        final GenotypeCalculationArgumentCollection args = new GenotypeCalculationArgumentCollection();
        args.samplePloidy = 2;
        args.MAX_ALTERNATE_ALLELES = 2;
        final StandardCallerArgumentCollection callerArgs = new StandardCallerArgumentCollection();
        callerArgs.genotypeArgs = args;

        final FixedAFCalculatorProvider p = new FixedAFCalculatorProvider(callerArgs, true);
        p.getInstance(2, 18);
    }

    @Test
    public void testTesInstanceInvalidAlleleNumberError_noVerify() throws Exception {
        final GenotypeCalculationArgumentCollection args = new GenotypeCalculationArgumentCollection();
        args.samplePloidy = 2;
        args.MAX_ALTERNATE_ALLELES = 2;
        final StandardCallerArgumentCollection callerArgs = new StandardCallerArgumentCollection();
        callerArgs.genotypeArgs = args;

        final FixedAFCalculatorProvider p = new FixedAFCalculatorProvider(callerArgs, false);
        p.getInstance(2, 18); //this passes - no validation
    }
}
