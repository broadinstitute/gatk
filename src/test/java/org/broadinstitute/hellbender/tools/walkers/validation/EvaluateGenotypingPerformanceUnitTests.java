package org.broadinstitute.hellbender.tools.walkers.validation;

import com.google.common.collect.ImmutableMap;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Random;

public class EvaluateGenotypingPerformanceUnitTests extends GATKBaseTest {

    @Test
    void testPearsonCorrelationAggregator() {
        final Random rand = new Random(234);
        final int nEntries = 1000;
        final double[] x = new double[nEntries];
        final double[] y = new double[nEntries];

        final EvaluateGenotypingPerformance.PearsonCorrelationAggregator aggregator = new EvaluateGenotypingPerformance.PearsonCorrelationAggregator();

        for (int i = 0; i < nEntries; i++) {
            x[i] = rand.nextDouble();
            y[i] = rand.nextDouble();

            aggregator.addEntry(x[i], y[i]);
        }

        final PearsonsCorrelation pearsonsCorrelation = new PearsonsCorrelation();
        Assert.assertEquals(aggregator.getCorrelation(), pearsonsCorrelation.correlation(x, y), 1e-10);

    }

    @DataProvider(name = "getConcordanceStateDataProvider")
    Object[][] getConcordanceStateDataProvider() {
        return new Object[][] {
                {Arrays.asList(Allele.ALT_A, Allele.REF_G), Arrays.asList(Allele.REF_G, Allele.ALT_A), false, ConcordanceState.TRUE_POSITIVE},
                {Arrays.asList(Allele.ALT_A, Allele.REF_G), Arrays.asList(Allele.REF_G, Allele.ALT_A), true, ConcordanceState.FILTERED_FALSE_NEGATIVE},
                {Arrays.asList(Allele.ALT_A, Allele.ALT_A), Arrays.asList(Allele.REF_G, Allele.ALT_A), false, ConcordanceState.FALSE_POSITIVE},
                {Arrays.asList(Allele.ALT_A, Allele.ALT_A), Arrays.asList(Allele.REF_G, Allele.ALT_A), true, ConcordanceState.FILTERED_FALSE_NEGATIVE},
                {Arrays.asList(Allele.ALT_A, Allele.ALT_A), Arrays.asList(Allele.ALT_G, Allele.ALT_A), false, ConcordanceState.FALSE_POSITIVE},
                {Arrays.asList(Allele.ALT_A, Allele.ALT_A), Arrays.asList(Allele.ALT_A, Allele.ALT_A), false, ConcordanceState.TRUE_POSITIVE},
                {Arrays.asList(Allele.REF_A, Allele.REF_A), Arrays.asList(Allele.REF_A, Allele.REF_A), false, ConcordanceState.TRUE_NEGATIVE},
                {Arrays.asList(Allele.REF_A, Allele.REF_A), Arrays.asList(Allele.REF_A, Allele.REF_A), true, ConcordanceState.TRUE_NEGATIVE},
                {Arrays.asList(Allele.REF_A, Allele.REF_A), Arrays.asList(Allele.REF_A, Allele.ALT_C), false, ConcordanceState.FALSE_POSITIVE},
                {Arrays.asList(Allele.REF_A, Allele.REF_A), Arrays.asList(Allele.REF_A, Allele.ALT_C), true, ConcordanceState.FILTERED_TRUE_NEGATIVE},

                //differing ploidy, generally from different treatment of chromosome X in males

                {Arrays.asList(Allele.REF_A, Allele.REF_A), Arrays.asList(Allele.REF_A), false, ConcordanceState.TRUE_NEGATIVE},
                {Arrays.asList(Allele.ALT_A, Allele.ALT_A), Arrays.asList(Allele.REF_A), false, ConcordanceState.FALSE_NEGATIVE},
                {Arrays.asList(Allele.REF_A, Allele.REF_A), Arrays.asList(Allele.ALT_A), false, ConcordanceState.FALSE_POSITIVE},
                {Arrays.asList(Allele.ALT_A, Allele.ALT_A), Arrays.asList(Allele.ALT_A), false, ConcordanceState.TRUE_POSITIVE},
                {Arrays.asList(Allele.REF_A), Arrays.asList(Allele.REF_A, Allele.REF_A), false, ConcordanceState.TRUE_NEGATIVE},
                {Arrays.asList(Allele.REF_A), Arrays.asList(Allele.ALT_A, Allele.REF_A), false, ConcordanceState.FALSE_POSITIVE},
                {Arrays.asList(Allele.REF_A), Arrays.asList(Allele.ALT_A, Allele.REF_A), true, ConcordanceState.FILTERED_TRUE_NEGATIVE},
                {Arrays.asList(Allele.ALT_A), Arrays.asList(Allele.REF_A, Allele.REF_A), false, ConcordanceState.FALSE_NEGATIVE},
                {Arrays.asList(Allele.ALT_A), Arrays.asList(Allele.REF_A, Allele.ALT_A), false, ConcordanceState.FALSE_POSITIVE},
                {Arrays.asList(Allele.ALT_A), Arrays.asList(Allele.ALT_A, Allele.ALT_A), false, ConcordanceState.TRUE_POSITIVE},
                {Arrays.asList(Allele.ALT_A), Arrays.asList(Allele.ALT_A, Allele.ALT_A), true, ConcordanceState.FILTERED_FALSE_NEGATIVE}
        };
    }

    @Test(dataProvider = "getConcordanceStateDataProvider")
    void testGetConcordanceState(final List<Allele> truthAlleles, final List<Allele> evalAlleles, final boolean evalWasFiltered, final ConcordanceState expectedConcordanceState ) {
        final GenotypeBuilder genotypeBuilder = new GenotypeBuilder();
        final Genotype truthGenotype = genotypeBuilder.alleles(truthAlleles).make();
        final Genotype evalGenotype = genotypeBuilder.alleles(evalAlleles).make();
        final EvaluateGenotypingPerformance evaluateGenotypingPerformance = new EvaluateGenotypingPerformance();
        evaluateGenotypingPerformance.allowDifferingPloidy = true;
        final ConcordanceState observedConcordanceState = evaluateGenotypingPerformance.getConcordanceState(truthGenotype, evalGenotype, evalWasFiltered);

        Assert.assertEquals(observedConcordanceState, expectedConcordanceState);
    }

    @DataProvider(name = "getConcordanceStateDifferingPloidyExceptionDataProvider")
    Object[][] getConcordanceStateDifferingPloidyExceptionDataProvider() {
        return new Object[][] {
                {Arrays.asList(Allele.REF_A, Allele.REF_A), Arrays.asList(Allele.REF_A)},
                {Arrays.asList(Allele.ALT_A, Allele.ALT_A), Arrays.asList(Allele.REF_A)},
                {Arrays.asList(Allele.REF_A, Allele.REF_A), Arrays.asList(Allele.ALT_A)},
                {Arrays.asList(Allele.ALT_A, Allele.ALT_A), Arrays.asList(Allele.ALT_A)},
                {Arrays.asList(Allele.REF_A), Arrays.asList(Allele.REF_A, Allele.REF_A)},
                {Arrays.asList(Allele.REF_A), Arrays.asList(Allele.ALT_A, Allele.REF_A)},
                {Arrays.asList(Allele.REF_A), Arrays.asList(Allele.ALT_A, Allele.REF_A)},
                {Arrays.asList(Allele.ALT_A), Arrays.asList(Allele.REF_A, Allele.REF_A)},
                {Arrays.asList(Allele.ALT_A), Arrays.asList(Allele.REF_A, Allele.ALT_A)},
                {Arrays.asList(Allele.ALT_A), Arrays.asList(Allele.ALT_A, Allele.ALT_A)},
                {Arrays.asList(Allele.ALT_A), Arrays.asList(Allele.ALT_A, Allele.ALT_A)}
        };
    }

    @Test(dataProvider = "getConcordanceStateDifferingPloidyExceptionDataProvider", expectedExceptions = GATKException.class)
    void testGetConcordanceStateDifferingPloidyException(final List<Allele> truthAlleles, final List<Allele> evalAlleles) {
        final GenotypeBuilder genotypeBuilder = new GenotypeBuilder();
        final Genotype truthGenotype = genotypeBuilder.alleles(truthAlleles).make();
        final Genotype evalGenotype = genotypeBuilder.alleles(evalAlleles).make();

        final EvaluateGenotypingPerformance evaluateGenotypingPerformance = new EvaluateGenotypingPerformance();
        evaluateGenotypingPerformance.getConcordanceState(truthGenotype, evalGenotype, false);
    }

    @DataProvider(name = "buildAFMapDataProvider")
    Object[][] buildAFMapDataProvider() {
        return new Object[][]{
                {
                        ImmutableMap.builder().
                                put("annotation1", Arrays.asList(0.3, 0.4)).
                                put("annotation2", Arrays.asList(0.1, 0.05)).
                                build(),
                        0,
                        ImmutableMap.builder().
                                put("annotation1", 0.3).
                                put("annotation2", 0.1).
                                build()
                },
                {
                        ImmutableMap.builder().
                                put("annotation1", Arrays.asList(0.3, 0.4)).
                                put("annotation2", Arrays.asList(0.1, 0.05)).
                                build(),
                        1,
                        ImmutableMap.builder().
                                put("annotation1", 0.4).
                                put("annotation2", 0.05).
                                build()
                }
        };
    }

    @Test(dataProvider = "buildAFMapDataProvider")
    void testBuildAFMap(final Map<String, List<Double>> annotationsMap, final int altAlleleIndex, final Map<String, Double> expectedMap) {
        final VariantContextBuilder vb = new VariantContextBuilder("snp1", "chr1", 100, 100, Arrays.asList(Allele.REF_A, Allele.ALT_G, Allele.ALT_C));
        for(final Map.Entry<String, List<Double>> entry : annotationsMap.entrySet()) {
            vb.attribute(entry.getKey(), entry.getValue());
        }

        final VariantContext vc = vb.make();

        final Map<String, Double> actualMap = EvaluateGenotypingPerformance.buildAFMapForIndex(vc, annotationsMap.keySet(), altAlleleIndex);
        Assert.assertEquals(actualMap, expectedMap);
    }

    @Test
    void testBuildAFMapSkipsMissingAnnotations() {
        final VariantContext vc = new VariantContextBuilder("snp1", "chr1", 100, 100, Arrays.asList(Allele.REF_A, Allele.ALT_G, Allele.ALT_C))
                .attribute("AF", Arrays.asList(".", "0.05"))
                .make();

        Assert.assertTrue(EvaluateGenotypingPerformance.buildAFMapForIndex(vc, Arrays.asList("AF"), 0).isEmpty());
        Assert.assertEquals(EvaluateGenotypingPerformance.buildAFMapForIndex(vc, Arrays.asList("AF"), 1), ImmutableMap.of("AF", 0.05));
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    void testRejectsMultiallelicEvaluation() {
        final VariantContext eval = new VariantContextBuilder("eval", "chr1", 100, 100, Arrays.asList(Allele.REF_A, Allele.ALT_G, Allele.ALT_C)).make();
        EvaluateGenotypingPerformance.validateBiallelicEval(eval);
    }

    @Test
    void testAcceptsBiallelicEvaluation() {
        final VariantContext eval = new VariantContextBuilder("eval", "chr1", 100, 100, Arrays.asList(Allele.REF_A, Allele.ALT_G)).make();
        EvaluateGenotypingPerformance.validateBiallelicEval(eval);
    }

    @DataProvider(name = "projectGenotypeDataProvider")
    Object[][] projectGenotypeDataProvider() {
        //truth is A>C,G, scoring the C allele
        return new Object[][] {
                {Arrays.asList(Allele.REF_A, Allele.REF_A), Arrays.asList(Allele.REF_A, Allele.REF_A)},
                {Arrays.asList(Allele.REF_A, Allele.ALT_C), Arrays.asList(Allele.REF_A, Allele.ALT_C)},
                {Arrays.asList(Allele.ALT_C, Allele.ALT_C), Arrays.asList(Allele.ALT_C, Allele.ALT_C)},
                {Arrays.asList(Allele.REF_A, Allele.ALT_G), Arrays.asList(Allele.REF_A, Allele.REF_A)},
                {Arrays.asList(Allele.ALT_G, Allele.ALT_G), Arrays.asList(Allele.REF_A, Allele.REF_A)},
                {Arrays.asList(Allele.ALT_C, Allele.ALT_G), Arrays.asList(Allele.ALT_C, Allele.REF_A)},
                {Arrays.asList(Allele.ALT_G), Arrays.asList(Allele.REF_A)},
                {Arrays.asList(Allele.ALT_C), Arrays.asList(Allele.ALT_C)}
        };
    }

    @Test(dataProvider = "projectGenotypeDataProvider")
    void testProjectGenotypeOntoAltAllele(final List<Allele> alleles, final List<Allele> expectedAlleles) {
        final Genotype genotype = new GenotypeBuilder("sample", alleles).make();
        final Genotype projected = EvaluateGenotypingPerformance.projectGenotypeOntoAltAllele(genotype, Allele.REF_A, Allele.ALT_C);
        Assert.assertEquals(projected.getAlleles(), expectedAlleles);
        Assert.assertEquals(projected.getSampleName(), "sample");
    }

    @Test
    void testProjectNullGenotype() {
        Assert.assertNull(EvaluateGenotypingPerformance.projectGenotypeOntoAltAllele(null, Allele.REF_A, Allele.ALT_C));
    }

    @DataProvider(name = "getDosageFracDataProvider")
    Object[][] getDosageFracDataProvider() {
        final List<Allele> multiallelic = Arrays.asList(Allele.REF_A, Allele.ALT_C, Allele.ALT_G);
        final List<Allele> biallelic = Arrays.asList(Allele.REF_A, Allele.ALT_C);
        return new Object[][] {
                //from genotype call
                {multiallelic, Arrays.asList(Allele.ALT_C, Allele.ALT_G), null, Allele.ALT_C, null, 0.5},
                {multiallelic, Arrays.asList(Allele.ALT_C, Allele.ALT_G), null, Allele.ALT_G, null, 0.5},
                {multiallelic, Arrays.asList(Allele.ALT_G, Allele.ALT_G), null, Allele.ALT_C, null, 0.0},
                {biallelic, Arrays.asList(Allele.ALT_C, Allele.ALT_C), null, Allele.ALT_C, null, 1.0},
                {biallelic, Arrays.asList(Allele.ALT_C), null, Allele.ALT_C, null, 1.0},
                //from dosage field, with one value per alt allele
                {multiallelic, Arrays.asList(Allele.ALT_C, Allele.ALT_G), "0.9,1.1", Allele.ALT_C, "DS", 0.45},
                {multiallelic, Arrays.asList(Allele.ALT_C, Allele.ALT_G), "0.9,1.1", Allele.ALT_G, "DS", 0.55},
                {biallelic, Arrays.asList(Allele.REF_A, Allele.ALT_C), "0.8", Allele.ALT_C, "DS", 0.4},
                //dosage field requested but missing from genotype, so fall back to genotype call
                {biallelic, Arrays.asList(Allele.REF_A, Allele.ALT_C), null, Allele.ALT_C, "DS", 0.5}
        };
    }

    @Test(dataProvider = "getDosageFracDataProvider")
    void testGetDosageFrac(final List<Allele> vcAlleles, final List<Allele> genotypeAlleles, final String dosage, final Allele altAllele,
                           final String dosageField, final double expectedDosageFrac) {
        final VariantContext vc = new VariantContextBuilder("vc", "chr1", 100, 100, vcAlleles).make();
        final GenotypeBuilder genotypeBuilder = new GenotypeBuilder("sample", genotypeAlleles);
        if (dosage != null) {
            genotypeBuilder.attribute("DS", dosage);
        }
        Assert.assertEquals(EvaluateGenotypingPerformance.getDosageFrac(genotypeBuilder.make(), vc, altAllele, dosageField), expectedDosageFrac, 1e-10);
    }

    @Test
    void testGetDosageFracNullGenotype() {
        final VariantContext vc = new VariantContextBuilder("vc", "chr1", 100, 100, Arrays.asList(Allele.REF_A, Allele.ALT_C)).make();
        Assert.assertEquals(EvaluateGenotypingPerformance.getDosageFrac(null, vc, Allele.ALT_C, "DS"), 0.0);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    void testGetDosageFracMissingAlleleDosage() {
        final VariantContext vc = new VariantContextBuilder("vc", "chr1", 100, 100, Arrays.asList(Allele.REF_A, Allele.ALT_C, Allele.ALT_G)).make();
        final Genotype genotype = new GenotypeBuilder("sample", Arrays.asList(Allele.ALT_C, Allele.ALT_G)).attribute("DS", "0.9").make();
        EvaluateGenotypingPerformance.getDosageFrac(genotype, vc, Allele.ALT_G, "DS");
    }

    @Test
    void testMatchesSplitEvalsToMultiallelicTruth() {
        final VariantContext truth = new VariantContextBuilder("truth", "chr1", 100, 100, Arrays.asList(Allele.REF_A, Allele.ALT_C, Allele.ALT_G)).make();
        final VariantContext evalC = new VariantContextBuilder("eval", "chr1", 100, 100, Arrays.asList(Allele.REF_A, Allele.ALT_C)).make();
        final VariantContext evalG = new VariantContextBuilder("eval", "chr1", 100, 100, Arrays.asList(Allele.REF_A, Allele.ALT_G)).make();
        final VariantContext evalT = new VariantContextBuilder("eval", "chr1", 100, 100, Arrays.asList(Allele.REF_A, Allele.ALT_T)).make();

        final EvaluateGenotypingPerformance evaluateGenotypingPerformance = new EvaluateGenotypingPerformance();
        //both split eval records should be compared to the same multiallelic truth record
        Assert.assertTrue(evaluateGenotypingPerformance.allowTruthToMatchMultipleEvals());
        Assert.assertTrue(evaluateGenotypingPerformance.shouldVariantsBeMatched(truth, evalC));
        Assert.assertTrue(evaluateGenotypingPerformance.shouldVariantsBeMatched(truth, evalG));
        Assert.assertFalse(evaluateGenotypingPerformance.shouldVariantsBeMatched(truth, evalT));
    }

    @Test
    void testGetACBin() {
        Assert.assertNull(EvaluateGenotypingPerformance.getACBin(0, 20));
        Assert.assertEquals(EvaluateGenotypingPerformance.getACBin(1, 20), Integer.valueOf(0));
        Assert.assertEquals(EvaluateGenotypingPerformance.getACBin(20, 20), Integer.valueOf(19));
        Assert.assertEquals(EvaluateGenotypingPerformance.getACBin(21, 20), Integer.valueOf(19));
    }

    @DataProvider(name = "getBinDataProvider")
    Object[][] getBinDataProvider() {
        return new Object[][] {
                {12, 0.01, 0.075, 5},
                {12, 0.01, 0.002, 0},
                {12, 0.01, 0, 0},
                {12, 0.01, 1, 11},
                {12, 0.01,0.86, 11},
                {12, 0.01,0.5, 10},

                {6,0.3,0.7, 4},
                {21,0.001,0.024, 10}
        };
    }

    @Test(dataProvider = "getBinDataProvider")
    void testGetBinDataProvider(final int nBins, final double firstBinRightEdge, final double val, final int expectedBin) {
        final EvaluateGenotypingPerformance evaluateGenotypingPerformance = new EvaluateGenotypingPerformance();
        evaluateGenotypingPerformance.nBins = nBins;
        evaluateGenotypingPerformance.firstBinRightEdge = firstBinRightEdge;
        final int actualBin = evaluateGenotypingPerformance.getBin(val);
        Assert.assertEquals(actualBin, expectedBin);
    }

}