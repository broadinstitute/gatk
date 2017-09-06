package org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc;

import htsjdk.variant.variantcontext.*;
import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.lang.math.NumberUtils;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypingEngine;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;


public final class AFCalculationUnitTest extends BaseTest {
    static Allele A = Allele.create("A", true);
    static Allele C = Allele.create("C");
    static Allele G = Allele.create("G");

    static int sampleNameCounter = 0;
    static Genotype AA1, AB1, BB1, NON_INFORMATIVE1;
    static Genotype AA2, AB2, AC2, BB2, BC2, CC2, NON_INFORMATIVE2;
    final double[] FLAT_3SAMPLE_PRIORS = MathUtils.normalizeLog10(new double[2 * 3 + 1]);  // flat priors

    final private static boolean INCLUDE_BIALLELIC = true;
    final private static boolean INCLUDE_TRIALLELIC = true;
    final private static boolean Guillermo_FIXME = false; // TODO -- can only be enabled when GdA fixes bug
    final private static boolean DEBUG_ONLY = false;

    protected static List<AFCalculator> createAFCalculators(final List<AFCalculatorImplementation> calcs, final int maxAltAlleles, final int ploidy) {
        final List<AFCalculator> AFCalculators = new LinkedList<>();

        for ( final AFCalculatorImplementation calc : calcs ) {
            if (calc.usableForParams(ploidy,maxAltAlleles))
                AFCalculators.add(calc.newInstance());
            else
                throw new IllegalStateException("cannot use " + calc + " calculator instance with combination " + maxAltAlleles + " " + ploidy);
        }

        return AFCalculators;
    }

    @BeforeClass
    public void before() {
        AA1 = makePL(Arrays.asList(A, A), 0, 20, 20);
        AB1 = makePL(Arrays.asList(A, C), 20, 0, 20);
        BB1 = makePL(Arrays.asList(C, C), 20, 20, 0);
        NON_INFORMATIVE1 = makePL(Arrays.asList(Allele.NO_CALL, Allele.NO_CALL), 0, 0, 0);

        AA2 = makePL(Arrays.asList(A, A), 0, 20, 20, 20, 20, 20);
        AB2 = makePL(Arrays.asList(A, C), 20, 0, 20, 20, 20, 20);
        BB2 = makePL(Arrays.asList(C, C), 20, 20, 0, 20, 20, 20);
        AC2 = makePL(Arrays.asList(A, G), 20, 20, 20, 0, 20, 20);
        BC2 = makePL(Arrays.asList(C, G), 20, 20, 20, 20, 0, 20);
        CC2 = makePL(Arrays.asList(G, G), 20, 20, 20, 20, 20, 0);
        NON_INFORMATIVE2 = makePL(Arrays.asList(Allele.NO_CALL, Allele.NO_CALL), 0, 0, 0, 0, 0, 0);
    }

    protected static Genotype makePL(final List<Allele> expectedGT, int ... pls) {
        GenotypeBuilder gb = new GenotypeBuilder("sample" + sampleNameCounter++);
        gb.alleles(expectedGT);
        gb.PL(pls);
        return gb.make();
    }


    private static final int MAX_ALT_ALLELES = 2;
    private static final int PLOIDY = 2;


    @DataProvider(name = "wellFormedGLs")
    public Object[][] createSimpleGLsData() {
        final List<Genotype> biAllelicSamples = Arrays.asList(AA1, AB1, BB1);
        final List<Genotype> triAllelicSamples = Arrays.asList(AA2, AB2, BB2, AC2, BC2, CC2);

        for ( final int nSamples : Arrays.asList(1, 2, 3, 4) ) {
            List<AFCalculator> calcs = createAFCalculators(Arrays.asList(AFCalculatorImplementation.values()), MAX_ALT_ALLELES, PLOIDY);

            final int nPriorValues = 2*nSamples+1;
            final double[] flatPriors = MathUtils.normalizeLog10(new double[nPriorValues]);  // flat priors
            final double[] humanPriors = new double[nPriorValues];
            GenotypingEngine.computeAlleleFrequencyPriors(nPriorValues - 1, humanPriors, 0.001, new ArrayList<>());

            for ( final double[] priors : Arrays.asList(flatPriors, humanPriors) ) { // , humanPriors) ) {
                for ( AFCalculator model : calcs ) {
                    final String priorName = priors == humanPriors ? "human" : "flat";

                    // bi-allelic
                    if ( INCLUDE_BIALLELIC && nSamples <= biAllelicSamples.size() )
                        for ( List<Genotype> genotypes : Utils.makePermutations(biAllelicSamples, nSamples, true) )
                            new GetGLsTest(model, 1, genotypes, priors, priorName);

                    // tri-allelic
                    if ( INCLUDE_TRIALLELIC && ( ! priorName.equals("human") || Guillermo_FIXME ) && ! ( model instanceof OriginalDiploidExactAFCalculator) ) // || model != generalCalc ) )
                        for ( List<Genotype> genotypes : Utils.makePermutations(triAllelicSamples, nSamples, true) )
                            new GetGLsTest(model, 2, genotypes, priors, priorName);
                }
            }
        }

        return GetGLsTest.getTests(GetGLsTest.class);
    }

    @Test(enabled = true && ! DEBUG_ONLY, dataProvider = "wellFormedGLs")
    public void testBiallelicGLs(GetGLsTest cfg) {
        if ( cfg.getAlleles().size() == 2 )
            testResultSimple(cfg);
    }

    @Test(enabled = true && ! DEBUG_ONLY, dataProvider = "wellFormedGLs")
    public void testTriallelicGLs(GetGLsTest cfg) {
        if ( cfg.getAlleles().size() > 2 )
            testResultSimple(cfg);
    }

    private static class NonInformativeData {
        final Genotype nonInformative;
        final List<Genotype> called;
        final int nAltAlleles;

        private NonInformativeData(List<Genotype> called, Genotype nonInformative, int nAltAlleles) {
            this.called = called;
            this.nonInformative = nonInformative;
            this.nAltAlleles = nAltAlleles;
        }
    }

    @DataProvider(name = "GLsWithNonInformative")
    public Object[][] makeGLsWithNonInformative() {
        List<Object[]> tests = new ArrayList<>();

        final List<NonInformativeData> nonInformativeTests = new LinkedList<>();
        nonInformativeTests.add(new NonInformativeData(Arrays.asList(AB1), NON_INFORMATIVE1, 1));
        nonInformativeTests.add(new NonInformativeData(Arrays.asList(AB2), NON_INFORMATIVE2, 2));
        nonInformativeTests.add(new NonInformativeData(Arrays.asList(AB2, BC2), NON_INFORMATIVE2, 2));

        for ( final int nNonInformative : Arrays.asList(1, 10, 100) ) {
            for ( final NonInformativeData testData : nonInformativeTests ) {
                final List<Genotype> samples = new ArrayList<>();
                samples.addAll(testData.called);
                samples.addAll(Collections.nCopies(nNonInformative, testData.nonInformative));

                final int nSamples = samples.size();
                List<AFCalculator> calcs = createAFCalculators(Arrays.asList(AFCalculatorImplementation.values()), MAX_ALT_ALLELES, PLOIDY);

                final double[] priors = MathUtils.normalizeLog10(new double[2*nSamples+1]);  // flat priors

                for ( AFCalculator model : calcs ) {
                    if ( testData.nAltAlleles > 1 && model instanceof OriginalDiploidExactAFCalculator)
                        continue;

                    final GetGLsTest onlyInformative = new GetGLsTest(model, testData.nAltAlleles, testData.called, priors, "flat");

                    for ( int rotation = 0; rotation < nSamples; rotation++ ) {
                        Collections.rotate(samples, 1);
                        final GetGLsTest withNonInformative = new GetGLsTest(model, testData.nAltAlleles, samples, priors, "flat");
                        tests.add(new Object[]{onlyInformative, withNonInformative});
                    }
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = true && ! DEBUG_ONLY, dataProvider = "GLsWithNonInformative", dependsOnMethods = {"testBiallelicGLs", "testTriallelicGLs"})
    public void testGLsWithNonInformative(GetGLsTest onlyInformative, GetGLsTest withNonInformative) {
        final AFCalculationResult expected = onlyInformative.execute();
        final AFCalculationResult actual = withNonInformative.execute();

        testResultSimple(withNonInformative);
        compareAFCalcResults(actual, expected, onlyInformative.getCalc(), onlyInformative.numAltAlleles, true);
    }

    private void testResultSimple(final GetGLsTest cfg) {
        final AFCalculationResult refResultTracker = cfg.executeRef();
        final AFCalculationResult resultTracker = cfg.execute();
        try {
            compareAFCalcResults(resultTracker, refResultTracker, cfg.getCalc(), cfg.numAltAlleles, true);
        } catch (Throwable t) {
            cfg.execute();
            throw new RuntimeException(t);
        }
        Assert.assertNotNull(resultTracker.getAllelesUsedInGenotyping());
        Assert.assertTrue(cfg.getAlleles().containsAll(resultTracker.getAllelesUsedInGenotyping()), "Result object has alleles not in our initial allele list");

        for ( int altAlleleI = 0; altAlleleI < cfg.numAltAlleles; altAlleleI++ ) {
            int expectedAlleleCount = cfg.getExpectedAltAC(altAlleleI);
            int calcAC_MLE = resultTracker.getAlleleCountsOfMLE()[altAlleleI];

            final Allele allele = cfg.getAlleles().get(altAlleleI+1);
            Assert.assertEquals(calcAC_MLE, expectedAlleleCount, "MLE AC not equal to expected AC for allele " + allele);
        }
    }

    private void compareAFCalcResults(final AFCalculationResult actual, final AFCalculationResult expected, final AFCalculator calc, final int maxAltAlleles, final boolean onlyPosteriorsShouldBeEqual) {
        // note we cannot really test the multi-allelic case because we actually meaningfully differ among the models here
        final double TOLERANCE = maxAltAlleles > 1 ? 1000 : 0.1; // much tighter constraints on bi-allelic results

        if ( ! onlyPosteriorsShouldBeEqual ) {
            Assert.assertEquals(actual.getLog10PriorOfAFEq0(), expected.getLog10PriorOfAFEq0(), TOLERANCE, "Priors AF == 0");
            Assert.assertEquals(actual.getLog10PriorOfAFGT0(), expected.getLog10PriorOfAFGT0(), TOLERANCE, "Priors AF > 0");
            Assert.assertEquals(actual.getLog10LikelihoodOfAFEq0(), expected.getLog10LikelihoodOfAFEq0(), TOLERANCE, "Likelihoods AF == 0");
            Assert.assertEquals(actual.getLog10LikelihoodOfAFGT0(), expected.getLog10LikelihoodOfAFGT0(), TOLERANCE, "Likelihoods AF > 0");
        }
        Assert.assertEquals(actual.getLog10PosteriorOfAFEq0(), expected.getLog10PosteriorOfAFEq0(), TOLERANCE, "Posteriors AF == 0");
        Assert.assertEquals(actual.getLog10PosteriorOfAFGT0(), expected.getLog10PosteriorOfAFGT0(), TOLERANCE, "Posteriors AF > 0");
        Assert.assertTrue(Arrays.equals(actual.getAlleleCountsOfMLE(), expected.getAlleleCountsOfMLE()), "MLE ACs ");
        Assert.assertEquals(actual.getAllelesUsedInGenotyping(), expected.getAllelesUsedInGenotyping(), "Alleles used in genotyping");

        for ( final Allele a : expected.getAllelesUsedInGenotyping() ) {
            if ( ! a.isReference() ) {
                Assert.assertEquals(actual.getAlleleCountAtMLE(a), expected.getAlleleCountAtMLE(a), "MLE AC for allele " + a);
                // TODO -- enable me when IndependentAllelesDiploidExactAFCalc works properly
//                if ( ! ( calc instanceof GeneralPloidyExactAFCalc ) )
//                    // TODO -- delete when general ploidy works properly with multi-allelics
//                    Assert.assertEquals(actual.isPolymorphic(a, 0.0), expected.isPolymorphic(a, 0.0), "isPolymorphic with thread 0.0 for allele " + a);
            }
        }
    }

    @Test(enabled = true && ! DEBUG_ONLY, dataProvider = "Models")
    public void testLargeGLs(final ExactAFCalculator calc) {
        final Genotype BB = makePL(Arrays.asList(C, C), 20000000, 20000000, 0);
        GetGLsTest cfg = new GetGLsTest(calc, 1, Arrays.asList(BB, BB, BB), FLAT_3SAMPLE_PRIORS, "flat");

        final AFCalculationResult resultTracker = cfg.execute();

        int calculatedAlleleCount = resultTracker.getAlleleCountsOfMLE()[0];
        Assert.assertEquals(calculatedAlleleCount, 6);
    }

    @Test(enabled = true && ! DEBUG_ONLY, dataProvider = "Models")
    public void testMismatchedGLs(final ExactAFCalculator calc) {
        final Genotype AB = makePL(Arrays.asList(A, C), 2000, 0, 2000, 2000, 2000, 2000);
        final Genotype AC = makePL(Arrays.asList(A, G), 100, 100, 100, 0, 100, 100);
        GetGLsTest cfg = new GetGLsTest(calc, 2, Arrays.asList(AB, AC), FLAT_3SAMPLE_PRIORS, "flat");

        final AFCalculationResult resultTracker = cfg.execute();

        Assert.assertEquals(resultTracker.getAlleleCountsOfMLE()[0], 1);
        Assert.assertEquals(resultTracker.getAlleleCountsOfMLE()[1], 1);
    }

    // --------------------------------------------------------------------------------
    //
    // Code to test that the pNonRef value is meaningful
    //
    // --------------------------------------------------------------------------------

    private static class PNonRefData {
        final Genotype g;
        final double pNonRef, tolerance;
        final boolean canScale;
        final List<AFCalculatorImplementation> badModels;
        final VariantContext vc;

        private PNonRefData(final VariantContext vc, Genotype g, double pNonRef, double tolerance, final boolean canScale) {
            this(vc, g, pNonRef, tolerance, canScale, Collections.emptyList());
        }

        private PNonRefData(final VariantContext vc, Genotype g, double pNonRef, double tolerance, final boolean canScale, final List<AFCalculatorImplementation> badModels) {
            this.g = g;
            this.pNonRef = pNonRef;
            this.tolerance = tolerance;
            this.canScale = canScale;
            this.badModels = badModels;
            this.vc = vc;
        }

        public PNonRefData scale(final int scaleFactor) {
            if ( canScale ) {
                final int[] PLs = new int[g.getPL().length];
                for ( int i = 0; i < PLs.length; i++ ) PLs[i] = g.getPL()[i] * ((int) Math.log10(scaleFactor)+1);
                final Genotype scaledG = new GenotypeBuilder(g).PL(PLs).make();
                final double scaledPNonRef = pNonRef < 0.5 ? pNonRef / scaleFactor : 1 - ((1-pNonRef) / scaleFactor);
                return new PNonRefData(vc, scaledG, scaledPNonRef, tolerance, true);
            } else {
                return this;
            }
        }
    }

    @DataProvider(name = "PNonRef")
    public Object[][] makePNonRefTest() {
        List<Object[]> tests = new ArrayList<>();

        final List<Allele> AA = Arrays.asList(A, A);
        final List<Allele> AC = Arrays.asList(A, C);
        final List<Allele> CC = Arrays.asList(C, C);
        final List<Allele> AG = Arrays.asList(A, G);
        final List<Allele> GG = Arrays.asList(G, G);
        final List<Allele> CG = Arrays.asList(C, G);

        final VariantContext vc2 = new VariantContextBuilder("x","1", 1, 1, Arrays.asList(A, C)).make();
        final VariantContext vc3 = new VariantContextBuilder("x","1", 1, 1, Arrays.asList(A, C, G)).make();
        final AFCalculatorTestBuilder.PriorType priorType = AFCalculatorTestBuilder.PriorType.flat;

        final double TOLERANCE = 0.5;

        final List<PNonRefData> initialPNonRefData = Arrays.asList(
                // bi-allelic sites
                new PNonRefData(vc2, makePL(AA, 0, 10, 10), 0.1666667, TOLERANCE, true),
                new PNonRefData(vc2, makePL(AA, 0, 1, 10), 0.4721084, TOLERANCE, false),
                new PNonRefData(vc2, makePL(AA, 0, 1, 1), 0.6136992, TOLERANCE, false),
                new PNonRefData(vc2, makePL(AA, 0, 5, 5), 0.3874259, TOLERANCE, false),
                new PNonRefData(vc2, makePL(AC, 10, 0, 10), 0.9166667, TOLERANCE, true),
                new PNonRefData(vc2, makePL(CC, 10, 10, 0), 0.9166667, TOLERANCE, true),

                // tri-allelic sites -- cannot scale because of the naivety of our scaling estimator
                new PNonRefData(vc3, makePL(AA, 0, 10, 10, 10, 10, 10), 0.3023255813953489, TOLERANCE * 2, false), // more tolerance because constrained model is a bit inaccurate
                new PNonRefData(vc3, makePL(AC, 10, 0, 10, 10, 10, 10), 0.9166667, TOLERANCE, false),
                new PNonRefData(vc3, makePL(CC, 10, 10, 0, 10, 10, 10), 0.9166667, TOLERANCE, false),
                new PNonRefData(vc3, makePL(AG, 10, 10, 10, 0, 10, 10), 0.9166667, TOLERANCE, false),
                new PNonRefData(vc3, makePL(CG, 10, 10, 10, 10, 0, 10), 0.80, TOLERANCE, false),
                new PNonRefData(vc3, makePL(GG, 10, 10, 10, 10, 10, 0), 0.9166667, TOLERANCE, false)
        );

        for ( AFCalculatorImplementation modelType : Arrays.asList(AFCalculatorImplementation.EXACT_REFERENCE, AFCalculatorImplementation.EXACT_INDEPENDENT) ) {
            for ( int nNonInformative = 0; nNonInformative < 3; nNonInformative++ ) {
                for ( final PNonRefData rootData : initialPNonRefData ) {
                    for ( int plScale = 1; plScale <= 100000; plScale *= 10 ) {
                        if ( ! rootData.badModels.contains(modelType) && (plScale == 1 || rootData.canScale) ) {
                            final PNonRefData data = rootData.scale(plScale);
                            tests.add(new Object[]{data.vc, modelType, priorType, Arrays.asList(data.g), data.pNonRef, data.tolerance, nNonInformative});
                        }
                    }
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = true && ! DEBUG_ONLY, dataProvider = "PNonRef")
    public void testPNonRef(final VariantContext vcRoot,
                             AFCalculatorImplementation modelType,
                             AFCalculatorTestBuilder.PriorType priorType,
                             final List<Genotype> genotypes,
                             final double expectedPNonRef,
                             final double tolerance,
                             final int nNonInformative) {
        final AFCalculatorTestBuilder testBuilder
                = new AFCalculatorTestBuilder(1, vcRoot.getNAlleles()-1, modelType, priorType);

        final VariantContextBuilder vcb = new VariantContextBuilder(vcRoot);
        vcb.genotypes(genotypes);

        final AFCalculationResult resultTracker = testBuilder.makeModel().getLog10PNonRef(vcb.make(), PLOIDY, MAX_ALT_ALLELES, testBuilder.makePriors());

        Assert.assertEquals(resultTracker.getLog10PosteriorOfAFGT0(), Math.log10(expectedPNonRef), tolerance,
                "Actual pNonRef not within tolerance " + tolerance + " of expected");
    }

    @DataProvider(name = "PNonRefBiallelicSystematic")
    public Object[][] makePNonRefBiallelicSystematic() {
        List<Object[]> tests = new ArrayList<>();

        final List<Integer> bigNonRefPLs = Arrays.asList(0, 1, 2, 3, 4, 5, 10, 15, 20, 25, 50, 100, 1000);
        final List<List<Integer>> bigDiploidPLs = removeBadPLs(Utils.makePermutations(bigNonRefPLs, 3, true));

        for ( AFCalculatorImplementation modelType : AFCalculatorImplementation.values() ) {

            if ( false ) { // for testing only
                tests.add(new Object[]{modelType, toGenotypes(Arrays.asList(Arrays.asList(0, 100, 0)))});
            } else {
                if ( modelType == AFCalculatorImplementation.EXACT_GENERAL_PLOIDY ) continue; // TODO -- GENERAL_PLOIDY DOESN'T WORK

                // test all combinations of PLs for 1 sample
                for ( final List<List<Integer>> PLsPerSample : Utils.makePermutations(bigDiploidPLs, 1, true) ) {
                    tests.add(new Object[]{modelType, toGenotypes(PLsPerSample)});
                }


                final List<List<Integer>> smallDiploidPLs = new LinkedList<>();
                for ( final int nonRefPL : Arrays.asList(5, 10, 20, 30) ) {
                    for ( int i = 0; i < 2; i++ ) {
                        List<Integer> pls = new ArrayList<>(Collections.nCopies(3, nonRefPL));
                        pls.set(i, 0);
                        smallDiploidPLs.add(pls);
                    }
                }

                for ( final List<List<Integer>> PLsPerSample : Utils.makePermutations(smallDiploidPLs, 5, false) ) {
                    tests.add(new Object[]{modelType, toGenotypes(PLsPerSample)});
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    final List<List<Integer>> removeBadPLs(List<List<Integer>> listOfPLs) {
        List<List<Integer>> clean = new LinkedList<>();

        for ( final List<Integer> PLs : listOfPLs ) {
            int x = PLs.get(0);
            boolean bad = false;
            for ( int pl1 : PLs )
                if ( pl1 > x )
                    bad = true;
                else
                    x = pl1;
            if ( ! bad ) clean.add(PLs);
        }

        return clean;
    }

    private List<Genotype> toGenotypes(final List<List<Integer>> PLsPerSample) {
        final List<Allele> nocall = Arrays.asList(Allele.NO_CALL, Allele.NO_CALL);
        final List<Genotype> genotypes = new ArrayList<>(PLsPerSample.size());

        for ( final List<Integer> PLs : PLsPerSample ) {
            final int[] pls = ArrayUtils.toPrimitive(PLs.toArray(new Integer[3]));
            final int min = NumberUtils.min(pls);
            for ( int i = 0; i < pls.length; i++ ) pls[i] -= min;
            genotypes.add(makePL(nocall, pls));
        }

        return genotypes;
    }

    @Test(enabled = true && ! DEBUG_ONLY, dataProvider = "PNonRefBiallelicSystematic")
    public void PNonRefBiallelicSystematic(AFCalculatorImplementation modelType, final List<Genotype> genotypes) {
        //logger.warn("Running " + modelType + " with " + genotypes);
        final AFCalculatorTestBuilder refBuilder = new AFCalculatorTestBuilder(genotypes.size(), 1, AFCalculatorImplementation.EXACT_REFERENCE, AFCalculatorTestBuilder.PriorType.human);
        final AFCalculatorTestBuilder testBuilder = new AFCalculatorTestBuilder(genotypes.size(), 1, modelType, AFCalculatorTestBuilder.PriorType.human);

        final VariantContextBuilder vcb = new VariantContextBuilder("x", "1", 1, 1, Arrays.asList(A, C));
        vcb.genotypes(genotypes);

        final AFCalculationResult refResult = refBuilder.makeModel().getLog10PNonRef(vcb.make(), PLOIDY, MAX_ALT_ALLELES, testBuilder.makePriors());
        final AFCalculationResult testResult = testBuilder.makeModel().getLog10PNonRef(vcb.make(), PLOIDY, MAX_ALT_ALLELES, testBuilder.makePriors());

        final double tolerance = 1e-3;
        Assert.assertEquals(testResult.getLog10PosteriorOfAFGT0(), refResult.getLog10PosteriorOfAFGT0(), tolerance,
                "Actual pNonRef not within tolerance " + tolerance + " of expected");
        Assert.assertEquals(testResult.getAlleleCountsOfMLE(), refResult.getAlleleCountsOfMLE(),
                "Actual MLE " + Utils.join(",", testResult.getAlleleCountsOfMLE()) + " not equal to expected " + Utils.join(",", refResult.getAlleleCountsOfMLE()));
    }

    // --------------------------------------------------------------------------------
    //
    // Test priors
    //
    // --------------------------------------------------------------------------------

    @DataProvider(name = "Models")
    public Object[][] makeModels() {
        List<Object[]> tests = new ArrayList<>();

        for ( final AFCalculatorImplementation calc : AFCalculatorImplementation.values() ) {
            if ( calc.usableForParams(2, 4) )
                tests.add(new Object[]{AFCalculatorImplementation.bestValue(2, AFCalculatorImplementation.UNBOUND_ALTERNATIVE_ALLELE_COUNT, null).newInstance()});
        }

        return tests.toArray(new Object[][]{});
    }


    @Test(enabled = true, dataProvider =  "Models")
    public void testNoPrior(final AFCalculator model) {
        for ( int REF_PL = 10; REF_PL <= 20; REF_PL += 10 ) {
            final Genotype AB = makePL(Arrays.asList(A, C), REF_PL, 0, 10000);

            final double[] flatPriors = new double[]{0.0,0.0,0.0};
            final double[] noPriors = new double[3];
            // test that function computeAlleleFrequency correctly operates when the flat prior option is set
            // computeAlleleFrequencyPriors takes linear priors
            final ArrayList<Double> inputPrior = new ArrayList<>();
            inputPrior.add(1.0/3);
            inputPrior.add(1.0/3);
            GenotypingEngine.computeAlleleFrequencyPriors(2, noPriors, 0.0, inputPrior);

            GetGLsTest cfgFlatPrior = new GetGLsTest(model, 1, Arrays.asList(AB), flatPriors, "flatPrior");
            GetGLsTest cfgNoPrior = new GetGLsTest(model, 1, Arrays.asList(AB), flatPriors, "noPrior");
            final AFCalculationResult resultTrackerFlat = cfgFlatPrior.execute();
            final AFCalculationResult resultTrackerNoPrior = cfgNoPrior.execute();

            final double pRefWithNoPrior = AB.getLikelihoods().getAsVector()[0];
            final double pHetWithNoPrior = AB.getLikelihoods().getAsVector()[1]  - Math.log10(0.5);
            final double nonRefPost = Math.pow(10, pHetWithNoPrior) / (Math.pow(10, pRefWithNoPrior) + Math.pow(10, pHetWithNoPrior));
            final double log10NonRefPost = Math.log10(nonRefPost);

            if ( ! Double.isInfinite(log10NonRefPost) ) {
                // check that the no-prior and flat-prior constructions yield same result
                Assert.assertEquals(resultTrackerFlat.getLog10PosteriorOfAFGT0(), resultTrackerNoPrior.getLog10PosteriorOfAFGT0());
            }

        }
    }

    @Test(enabled = true && !DEBUG_ONLY, dataProvider = "Models")
    public void testBiallelicPriors(final AFCalculator model) {

        for ( int REF_PL = 10; REF_PL <= 20; REF_PL += 10 ) {
            final Genotype AB = makePL(Arrays.asList(A, C), REF_PL, 0, 10000);

            for ( int log10NonRefPrior = 1; log10NonRefPrior < 10*REF_PL; log10NonRefPrior += 1 ) {
                final double refPrior = 1 - QualityUtils.qualToErrorProb(log10NonRefPrior);
                final double nonRefPrior = (1-refPrior) / 2;
                final double[] priors = MathUtils.normalizeLog10(MathUtils.toLog10(new double[]{refPrior, nonRefPrior, nonRefPrior}));
                if ( ! Double.isInfinite(priors[1]) ) {
                    GetGLsTest cfg = new GetGLsTest(model, 1, Arrays.asList(AB), priors, "pNonRef" + log10NonRefPrior);
                    final AFCalculationResult resultTracker = cfg.execute();
                    final int actualAC = resultTracker.getAlleleCountsOfMLE()[0];

                    final double pRefWithPrior = AB.getLikelihoods().getAsVector()[0] + priors[0];
                    final double pHetWithPrior = AB.getLikelihoods().getAsVector()[1] + priors[1] - Math.log10(0.5);
                    final double nonRefPost = Math.pow(10, pHetWithPrior) / (Math.pow(10, pRefWithPrior) + Math.pow(10, pHetWithPrior));
                    final double log10NonRefPost = Math.log10(nonRefPost);

                    if ( ! Double.isInfinite(log10NonRefPost) )
                        Assert.assertEquals(resultTracker.getLog10PosteriorOfAFGT0(), log10NonRefPost, 1e-2);

                    if ( nonRefPost >= 0.9 )
                        Assert.assertTrue(resultTracker.isPolymorphic(C, -1));

                    final int expectedMLEAC = 1; // the MLE is independent of the prior
                    Assert.assertEquals(actualAC, expectedMLEAC,
                            "actual AC with priors " + log10NonRefPrior + " not expected "
                                    + expectedMLEAC + " priors " + Utils.join(",", priors));
                }
            }
        }
    }

    @Test(enabled = true && ! DEBUG_ONLY, dataProvider = "Models")

    // --------------------------------------------------------------------------------
    //
    // Test that polymorphic sites (bi and tri) are properly called
    //
    // --------------------------------------------------------------------------------

    @DataProvider(name = "polyTestProvider")
    public Object[][] makePolyTestProvider() {
        List<Object[]> tests = new ArrayList<>();

        // list of all high-quality models in the system
        final List<AFCalculatorImplementation> models = Arrays.asList(
                AFCalculatorImplementation.DEFAULT,
                AFCalculatorImplementation.EXACT_REFERENCE,
                AFCalculatorImplementation.EXACT_INDEPENDENT);

        // note that we cannot use small PLs here or the thresholds are hard to set
        for ( final int nonTypePLs : Arrays.asList(100, 1000) ) {
            for ( final AFCalculatorImplementation model : models ) {
                for ( final int allele1AC : Arrays.asList(0, 1, 2, 10, 100, 1000, 10000) ) {
                    for ( final int nSamples : Arrays.asList(1, 10, 100, 1000, 10000) ) {
//        for ( final int nonTypePLs : Arrays.asList(10) ) {
//            for ( final AFCalcFactory.Calculation model : models ) {
//                for ( final int allele1AC : Arrays.asList(100) ) {
//                    for ( final int nSamples : Arrays.asList(1000) ) {
                        if ( nSamples < allele1AC ) continue;

                        final double pPerSample = Math.pow(10, nonTypePLs / -10.0);
                        final double errorFreq = pPerSample * nSamples;
                        final boolean poly1 = allele1AC > errorFreq && (nonTypePLs * allele1AC) > 30;

                        // bi-allelic tests
                        {
                            final AFCalculatorTestBuilder testBuilder
                                    = new AFCalculatorTestBuilder(nSamples, 1, model, AFCalculatorTestBuilder.PriorType.human);
                            final List<Integer> ACs = Arrays.asList(allele1AC);
                            tests.add(new Object[]{testBuilder, ACs, nonTypePLs, Arrays.asList(poly1)});
                        }

                        // multi-allelic tests
                        for ( final int allele2AC : Arrays.asList(0, 1, 2, 10, 20, 50) ) {
                            if ( nSamples < allele2AC || allele1AC + allele2AC > nSamples || nSamples > 100 || nSamples == 1)
                                continue;

                            final AFCalculatorTestBuilder testBuilder
                                    = new AFCalculatorTestBuilder(nSamples, 2, model, AFCalculatorTestBuilder.PriorType.human);
                            final List<Integer> ACs = Arrays.asList(allele1AC, allele2AC);
                            final boolean poly2 = allele2AC > errorFreq && (nonTypePLs * allele2AC) > 90;
                            tests.add(new Object[]{testBuilder, ACs, nonTypePLs, Arrays.asList(poly1, poly2)});
                        }
                    }
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = true && ! DEBUG_ONLY, dataProvider = "polyTestProvider")
    public void testCallingGeneral(final AFCalculatorTestBuilder testBuilder, final List<Integer> ACs, final int nonTypePL, final List<Boolean> expectedPoly ) {
        testCalling(testBuilder, ACs, nonTypePL, expectedPoly);
    }

    @DataProvider(name = "polyTestProviderLotsOfAlleles")
    public Object[][] makepolyTestProviderLotsOfAlleles() {
        List<Object[]> tests = new ArrayList<>();

        // list of all high-quality models in the system
        final List<AFCalculatorImplementation> models = Arrays.asList(AFCalculatorImplementation.EXACT_INDEPENDENT);

        final List<Integer> alleleCounts = Arrays.asList(0, 1, 2, 3, 4, 5, 10, 20);

        final int nonTypePLs = 1000;
        final int nAlleles = 4;
        for ( final AFCalculatorImplementation model : models ) {
            for ( final List<Integer> ACs : Utils.makePermutations(alleleCounts, nAlleles, true) ) {
                final List<Boolean> isPoly = new ArrayList<>(ACs.size());
                for ( final int ac : ACs ) {
                    isPoly.add(ac > 0);
                }

                final double acSum = ACs.stream().mapToDouble(i->i).sum();
                for ( final int nSamples : Arrays.asList(1, 10, 100) ) {
                    if ( nSamples < acSum ) {
                        continue;
                    }
                    final AFCalculatorTestBuilder testBuilder
                            = new AFCalculatorTestBuilder(nSamples, nAlleles, model, AFCalculatorTestBuilder.PriorType.human);
                    tests.add(new Object[]{testBuilder, ACs, nonTypePLs, isPoly});
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = true && ! DEBUG_ONLY, dataProvider = "polyTestProviderLotsOfAlleles")
    public void testCallingLotsOfAlleles(final AFCalculatorTestBuilder testBuilder, final List<Integer> ACs, final int nonTypePL, final List<Boolean> expectedPoly ) {
        testCalling(testBuilder, ACs, nonTypePL, expectedPoly);
    }

    private void testCalling(final AFCalculatorTestBuilder testBuilder, final List<Integer> ACs, final int nonTypePL, final List<Boolean> expectedPoly) {
        final AFCalculator calc = testBuilder.makeModel();
        final double[] priors = testBuilder.makePriors();
        final VariantContext vc = testBuilder.makeACTest(ACs, 0, nonTypePL);
        final AFCalculationResult result = calc.getLog10PNonRef(vc, PLOIDY, testBuilder.numAltAlleles, priors);

        boolean anyPoly = false;
        for ( final boolean onePoly : expectedPoly ) anyPoly = anyPoly || onePoly;

        if ( anyPoly )
            Assert.assertTrue(result.getLog10PosteriorOfAFGT0() > -1);

        for ( int altI = 1; altI < result.getAllelesUsedInGenotyping().size(); altI++ ) {
            final int i = altI - 1;
            final Allele alt = result.getAllelesUsedInGenotyping().get(altI);

            // must be getCalledChrCount because we cannot ensure that the VC made has our desired ACs
            Assert.assertEquals(result.getAlleleCountAtMLE(alt), vc.getCalledChrCount(alt));
            Assert.assertEquals(result.isPolymorphic(alt, -1), (boolean) expectedPoly.get(i), "isPolymorphic for allele " + alt + " " + result.getLog10PosteriorOfAFEq0ForAllele(alt));
        }
    }
}