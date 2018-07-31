package org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc;

import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.variant.HomoSapiensConstants;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public final class GetGLsTest extends GATKBaseTest.TestDataProvider {
    GenotypesContext GLs;
    int numAltAlleles;
    final AFCalculator calc;
    final int[] expectedACs;
    final double[] priors;
    final String priorName;

    GetGLsTest(final AFCalculator calc, int numAltAlleles, List<Genotype> arg, final double[] priors, final String priorName) {
        super(GetGLsTest.class);
        GLs = GenotypesContext.create(new ArrayList<>(arg));
        this.numAltAlleles = numAltAlleles;
        this.calc = calc;
        this.priors = priors;
        this.priorName = priorName;

        expectedACs = new int[numAltAlleles+1];
        for ( int alleleI = 0; alleleI < expectedACs.length; alleleI++ ) {
            expectedACs[alleleI] = 0;
            final Allele allele = getAlleles().get(alleleI);
            for ( Genotype g : arg ) {
                expectedACs[alleleI] += Collections.frequency(g.getAlleles(), allele);
            }
        }
    }

    public AFCalculationResult execute() {
        return getCalc().getLog10PNonRef(getVC(), HomoSapiensConstants.DEFAULT_PLOIDY, numAltAlleles, getPriors());
    }

    public AFCalculationResult executeRef() {
        final AFCalculator ref = AFCalculatorImplementation.EXACT_REFERENCE.newInstance();
        return ref.getLog10PNonRef(getVC(), HomoSapiensConstants.DEFAULT_PLOIDY, numAltAlleles, getPriors());
    }

    public double[] getPriors() {
        return priors;
    }

    public AFCalculator getCalc() {
        return calc;
    }

    public VariantContext getVC() {
        VariantContextBuilder builder = new VariantContextBuilder("test", "1", 1, 1, getAlleles());
        builder.genotypes(GLs);
        return builder.make();
    }

    public List<Allele> getAlleles() {
        return Arrays.asList(Allele.create("A", true),
                Allele.create("C"),
                Allele.create("G"),
                Allele.create("T")).subList(0, numAltAlleles+1);
    }

    public int getExpectedAltAC(final int alleleI) {
        return expectedACs[alleleI+1];
    }

    public String toString() {
        return String.format("%s model=%s prior=%s input=%s", super.toString(), calc.getClass().getSimpleName(),
                priorName, GLs.size() > 5 ? String.format("%d samples", GLs.size()) : GLs);
    }
}
