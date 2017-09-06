package org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc;

import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;

public final class GeneralPloidyAFCalculationModelUnitTest extends BaseTest {

    static double[] AA1, AB1, BB1;
    static double[] AA2, AB2, AC2, BB2, BC2, CC2;
    static double[] A4_1, B4_1, C4_1, D4_1, E4_1,F4_1;
    static double[] A4_400, B4_310, C4_220, D4_130, E4_121, F4_013;
    static final int numSamples = 4;
    static final int samplePloidy = 4;   // = 2*samplesPerPool

    @BeforeClass
    public void before() {
        // legacy diploid cases
        AA1 = new double[]{-5.0, -20.0, -20.0};
        AB1 = new double[]{-20.0, 0.0, -20.0};
        BB1 = new double[]{-20.0, -20.0, 0.0};

        // diploid, nAlleles = 3. Ordering is [2 0 0] [1 1 0] [0 2 0] [1 0 1] [0 1 1] [0 0 2], ie AA AB BB AC BC CC
        AA2 = new double[]{0.0, -20.0, -20.0, -20.0, -20.0, -20.0};
        AB2 = new double[]{-20.0, 0.0, -20.0, -20.0, -20.0, -20.0};
        AC2 = new double[]{-20.0, -20.0, -20.0, 0.0, -20.0, -20.0};
        BB2 = new double[]{-20.0, -20.0, 0.0, -20.0, -20.0, -20.0};
        BC2 = new double[]{-20.0, -20.0, -20.0, -20.0, 0.0, -20.0};
        CC2 = new double[]{-20.0, -20.0, -20.0, -20.0, -20.0, 0.0};
        
        // pool (i.e. polyploid cases)
        // NAlleles = 2, ploidy=4
        // ordering is [4 0] [3 1] [2 2 ] [1 3] [0 4]

        A4_1 = new double[]{-3.0, -20.0, -20.0, -20.0, -20.0};
        B4_1 = new double[]{-20.0, 0.0, -20.0, -20.0, -20.0};
        C4_1 = new double[]{-20.0, -20.0, 0.0, -20.0, -20.0};
        D4_1 = new double[]{-20.0, -20.0, 0.0,   0.0, -20.0};
        E4_1 = new double[]{-20.0, -20.0, 0.0,   0.0, -20.0};
        F4_1 = new double[]{-20.0, -20.0, -20.0,   -20.0, 0.0};

        // NAlleles = 3, ploidy = 4
        // ordering is [4 0 0] [3 1 0] [2 2 0] [1 3 0] [0 4 0] [3 0 1] [2 1 1] [1 2 1] [0 3 1] [2 0 2] [1 1 2] [0 2 2] [1 0 3] [0 1 3] [0 0 4]
        A4_400 = new double[]{0.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0};
        B4_310 = new double[]{-20.0, 0.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0};
        C4_220 = new double[]{-20.0, -20.0, 0.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0};
        D4_130 = new double[]{-20.0, -20.0, -20.0,   0.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0};
        E4_121 = new double[]{-20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0,   0.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0};
        F4_013 = new double[]{-20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, 0.0, -20.0};

    }

    private class GetGLsTest extends TestDataProvider {
        GenotypesContext GLs;
        int numAltAlleles;
        String name;
        int ploidy;
        private GetGLsTest(String name, int numAltAlleles, int ploidy, Genotype... arg) {
            super(GetGLsTest.class, name);
            GLs = GenotypesContext.create(arg);
            this.name = name;
            this.numAltAlleles = numAltAlleles;
            this.ploidy = ploidy;
        }

        public String toString() {
            return String.format("%s input=%s", super.toString(), GLs);
        }
    }

    private static Genotype createGenotype(String name, double[] gls, int ploidy) {
        Allele[] alleles = new Allele[ploidy];
        
        for (int i=0; i < ploidy; i++)
            alleles[i] = Allele.NO_CALL;
        
        return new GenotypeBuilder(name, Arrays.asList(alleles)).PL(gls).make();
    }                              

    @DataProvider(name = "getGLs")
    public Object[][] createGLsData() {

        // bi-allelic diploid case
        new GetGLsTest("B0", 1, 2, createGenotype("AA1", AA1,2), createGenotype("AA2", AA1,2), createGenotype("AA3", AA1,2));
        new GetGLsTest("B1", 1, 2, createGenotype("AA1", AA1,2), createGenotype("AA2", AA1,2), createGenotype("AB", AB1,2));
        new GetGLsTest("B2", 1, 2, createGenotype("AA1", AA1,2), createGenotype("BB", BB1,2), createGenotype("AA2", AA1,2));
        new GetGLsTest("B3a", 1, 2, createGenotype("AB", AB1,2), createGenotype("AA", AA1,2), createGenotype("BB", BB1,2));
        new GetGLsTest("B3b", 1, 2, createGenotype("AB1", AB1,2), createGenotype("AB2", AB1,2), createGenotype("AB3", AB1,2));
        new GetGLsTest("B4", 1, 2, createGenotype("BB1", BB1,2), createGenotype("BB2", BB1,2), createGenotype("AA", AA1,2));
        new GetGLsTest("B5", 1, 2, createGenotype("BB1", BB1,2), createGenotype("AB", AB1,2), createGenotype("BB2", BB1,2));
        new GetGLsTest("B6", 1, 2, createGenotype("BB1", BB1,2), createGenotype("BB2", BB1,2), createGenotype("BB3", BB1,2));

        // tri-allelic diploid case
        new GetGLsTest("B1C0", 2, 2, createGenotype("AA1", AA2,2), createGenotype("AA2", AA2,2), createGenotype("AB", AB2,2));
        new GetGLsTest("B0C1", 2, 2, createGenotype("AA1", AA2,2), createGenotype("AA2", AA2,2), createGenotype("AC", AC2,2));
        new GetGLsTest("B1C1a", 2,2,  createGenotype("AA", AA2,2), createGenotype("AB", AB2,2), createGenotype("AC", AC2,2));
        new GetGLsTest("B1C1b", 2,2,  createGenotype("AA1", AA2,2), createGenotype("AA2", AA2,2), createGenotype("BC", BC2,2));
        new GetGLsTest("B2C1", 2, 2, createGenotype("AB1", AB2,2), createGenotype("AB2", AB2,2), createGenotype("AC", AC2,2));
        new GetGLsTest("B3C2a", 2, 2, createGenotype("AB", AB2,2), createGenotype("BC1", BC2,2), createGenotype("BC2", BC2,2));
        new GetGLsTest("B3C2b", 2, 2, createGenotype("AB", AB2,2), createGenotype("BB", BB2,2), createGenotype("CC", CC2,2));

        // bi-allelic pool case
        new GetGLsTest("P0", 1, samplePloidy, createGenotype("A4_1", A4_1,samplePloidy), createGenotype("A4_1", A4_1,samplePloidy), createGenotype("A4_1", A4_1,samplePloidy));
        new GetGLsTest("P1", 1, samplePloidy,createGenotype("A4_1", A4_1,samplePloidy), createGenotype("B4_1", B4_1,samplePloidy), createGenotype("A4_1", A4_1,samplePloidy));
        new GetGLsTest("P2a", 1,samplePloidy, createGenotype("A4_1", A4_1,samplePloidy), createGenotype("C4_1", C4_1,samplePloidy), createGenotype("A4_1", A4_1,samplePloidy));
        new GetGLsTest("P2b", 1, samplePloidy,createGenotype("B4_1", B4_1,samplePloidy), createGenotype("B4_1", B4_1,samplePloidy), createGenotype("A4_1", A4_1,samplePloidy));
        new GetGLsTest("P4", 1, samplePloidy,createGenotype("A4_1", A4_1,samplePloidy), createGenotype("C4_1", C4_1,samplePloidy), createGenotype("C4_1", C4_1,samplePloidy));
        new GetGLsTest("P6", 1, samplePloidy,createGenotype("A4_1", A4_1,samplePloidy), createGenotype("F4_1", F4_1,samplePloidy), createGenotype("C4_1", C4_1,samplePloidy));
        new GetGLsTest("P8", 1, samplePloidy,createGenotype("A4_1", A4_1,samplePloidy), createGenotype("F4_1", F4_1,samplePloidy), createGenotype("F4_1", F4_1,samplePloidy));

        // multi-allelic pool case
        new GetGLsTest("B1C3", 2, samplePloidy,createGenotype("A4_400", A4_400,samplePloidy), createGenotype("A4_400", A4_400,samplePloidy), createGenotype("F4_013", F4_013,samplePloidy));
        new GetGLsTest("B3C9", 2, samplePloidy,createGenotype("F4_013", F4_013,samplePloidy), createGenotype("F4_013", F4_013,samplePloidy), createGenotype("F4_013", F4_013,samplePloidy));
        new GetGLsTest("B6C0", 2, samplePloidy,createGenotype("B4_310", B4_310,samplePloidy), createGenotype("C4_220", C4_220,samplePloidy), createGenotype("D4_130", D4_130,samplePloidy));
        new GetGLsTest("B6C4", 2, samplePloidy,createGenotype("D4_130", D4_130,samplePloidy), createGenotype("E4_121", E4_121,samplePloidy), createGenotype("F4_013", F4_013,samplePloidy));
        new GetGLsTest("B4C7", 2, samplePloidy,createGenotype("F4_013", F4_013,samplePloidy), createGenotype("E4_121", E4_121,samplePloidy), createGenotype("F4_013", F4_013,samplePloidy));
        new GetGLsTest("B2C3", 2, samplePloidy,createGenotype("A4_400", A4_400,samplePloidy), createGenotype("F4_013", F4_013,samplePloidy), createGenotype("B4_310", B4_310,samplePloidy));

        return GetGLsTest.getTests(GetGLsTest.class);
    }

    @Test(dataProvider = "getGLs")
    public void testGLs(GetGLsTest cfg) {
        final int len = GenotypeLikelihoods.numLikelihoods(1 + cfg.numAltAlleles, cfg.ploidy * cfg.GLs.size());
        double[] priors = new double[len];  // flat priors

        final GeneralPloidyExactAFCalculator calc = new GeneralPloidyExactAFCalculator();
        calc.combineSinglePools(cfg.GLs, cfg.ploidy,cfg.numAltAlleles + 1, priors);
        int nameIndex = 1;

        for ( int allele = 0; allele < cfg.numAltAlleles; allele++, nameIndex+=2 ) {
            int expectedAlleleCount = Integer.valueOf(cfg.name.substring(nameIndex, nameIndex + 1));
            int calculatedAlleleCount = calc.getAltAlleleCountOfMAP(allele);
            Assert.assertEquals(calculatedAlleleCount, expectedAlleleCount);
        }
    }

}
