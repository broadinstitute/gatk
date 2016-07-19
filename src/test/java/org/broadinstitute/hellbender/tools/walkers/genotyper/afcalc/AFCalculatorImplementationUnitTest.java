package org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import static org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc.AFCalculatorImplementation.*;

public final class AFCalculatorImplementationUnitTest {

    @DataProvider(name = "AFCalculatorImplementation")
    public Iterator<Object[]> AFCalculatorImplementation() {
        final List<Object[]> list = new ArrayList<>();
        list.add(new Object[]{2, 2,   EXACT_ORIGINAL,  EXACT_ORIGINAL});
        list.add(new Object[]{2, 2,   EXACT_REFERENCE, EXACT_REFERENCE});
        list.add(new Object[]{2, 6,   EXACT_REFERENCE, EXACT_REFERENCE});
        list.add(new Object[]{2, 10,  EXACT_REFERENCE, EXACT_REFERENCE});
        list.add(new Object[]{2, 100, EXACT_REFERENCE, EXACT_REFERENCE});
        list.add(new Object[]{1, 6,   EXACT_REFERENCE, EXACT_GENERAL_PLOIDY});
        list.add(new Object[]{1, 10,  EXACT_REFERENCE, EXACT_GENERAL_PLOIDY});
        list.add(new Object[]{1, 100, EXACT_REFERENCE, EXACT_GENERAL_PLOIDY});

        list.add(new Object[]{2, 3, EXACT_ORIGINAL, EXACT_INDEPENDENT});
        return list.iterator();
    }

    @Test(dataProvider = "AFCalculatorImplementation")
    public void testPickBestOne(final int ploidy, final int ac, final AFCalculatorImplementation preferred, final AFCalculatorImplementation expected) {
        Assert.assertEquals(expected, AFCalculatorImplementation.bestValue(ploidy, ac, preferred));
    }

    @DataProvider(name = "impls")
    public Iterator<Object[]> impls() {
        final List<Object[]> list = new ArrayList<>();
        list.add(new Object[]{EXACT_ORIGINAL, OriginalDiploidExactAFCalculator.class});
        list.add(new Object[]{EXACT_GENERAL_PLOIDY, GeneralPloidyExactAFCalculator.class});
        list.add(new Object[]{EXACT_INDEPENDENT, IndependentAllelesDiploidExactAFCalculator.class});
        list.add(new Object[]{EXACT_REFERENCE, ReferenceDiploidExactAFCalculator.class});
        return list.iterator();
    }

    @Test(dataProvider = "impls")
    public void instance(final AFCalculatorImplementation impl, final Class<? extends AFCalculatorImplementation> clazz) {
        Assert.assertEquals(impl.newInstance().getClass(), clazz);
    }

    @Test
    public void testFromCalcClass() throws Exception {
        Assert.assertEquals(EXACT_INDEPENDENT, fromCalculatorClass(IndependentAllelesDiploidExactAFCalculator.class));
        Assert.assertEquals(EXACT_REFERENCE, fromCalculatorClass(ReferenceDiploidExactAFCalculator.class));
        Assert.assertEquals(EXACT_ORIGINAL, fromCalculatorClass(OriginalDiploidExactAFCalculator.class));
        Assert.assertEquals(EXACT_GENERAL_PLOIDY, fromCalculatorClass(GeneralPloidyExactAFCalculator.class));
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testFromCalcClassNull() throws Exception {
        fromCalculatorClass(null);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testFromCalcClassAbstract() throws Exception {
        fromCalculatorClass(ExactAFCalculator.class);
    }
}
