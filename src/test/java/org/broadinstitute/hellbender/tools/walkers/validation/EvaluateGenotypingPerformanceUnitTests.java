package org.broadinstitute.hellbender.tools.walkers.validation;

import com.google.common.collect.ImmutableMap;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
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

    @Test
    void testGetACForAltAllele() {
        final List<Allele> alleles = Arrays.asList(Allele.REF_A, Allele.ALT_C, Allele.ALT_G);
        final VariantContext withAC = new VariantContextBuilder("vc", "chr1", 100, 100, alleles).attribute("AC", Arrays.asList(3, 7)).make();
        Assert.assertEquals(EvaluateGenotypingPerformance.getACForAltAllele(withAC, 0), Integer.valueOf(3));
        Assert.assertEquals(EvaluateGenotypingPerformance.getACForAltAllele(withAC, 1), Integer.valueOf(7));

        //missing AC
        final VariantContext noAC = new VariantContextBuilder("vc", "chr1", 100, 100, alleles).make();
        Assert.assertNull(EvaluateGenotypingPerformance.getACForAltAllele(noAC, 0));

        //AC with fewer values than alt alleles
        final VariantContext shortAC = new VariantContextBuilder("vc", "chr1", 100, 100, alleles).attribute("AC", "3").make();
        Assert.assertEquals(EvaluateGenotypingPerformance.getACForAltAllele(shortAC, 0), Integer.valueOf(3));
        Assert.assertNull(EvaluateGenotypingPerformance.getACForAltAllele(shortAC, 1));
    }

    @DataProvider(name = "getGPBinDataProvider")
    Object[][] getGPBinDataProvider() {
        return new Object[][] {
                {0.0, 10, 0},
                {0.05, 10, 0},
                {0.1, 10, 1},
                {0.95, 10, 9},
                {1.0, 10, 9},
                {0.5, 1, 0},
                {1.0, 1, 0}
        };
    }

    @Test(dataProvider = "getGPBinDataProvider")
    void testGetGPBin(final double gp, final int nBins, final int expectedBin) {
        Assert.assertEquals(EvaluateGenotypingPerformance.getGPBin(gp, nBins), expectedBin);
    }

    @DataProvider(name = "getGenotypeProbabilitiesDataProvider")
    Object[][] getGenotypeProbabilitiesDataProvider() {
        final List<Allele> diploid = Arrays.asList(Allele.REF_A, Allele.ALT_C);
        final List<Allele> haploid = Collections.singletonList(Allele.ALT_C);
        return new Object[][] {
                //VCF-style comma-separated string
                {new GenotypeBuilder("s", diploid).attribute("GP", "0.1,0.7,0.2").make(), new double[] {0.1, 0.7, 0.2}},
                //BCF-style typed lists
                {new GenotypeBuilder("s", diploid).attribute("GP", Arrays.asList(0.1, 0.7, 0.2)).make(), new double[] {0.1, 0.7, 0.2}},
                {new GenotypeBuilder("s", diploid).attribute("GP", Arrays.asList(0.1f, 0.7f, 0.2f)).make(), new double[] {0.1f, 0.7f, 0.2f}},
                //haploid
                {new GenotypeBuilder("s", haploid).attribute("GP", "0.3,0.7").make(), new double[] {0.3, 0.7}},
                {new GenotypeBuilder("s", haploid).attribute("GP", Arrays.asList(0.3, 0.7)).make(), new double[] {0.3, 0.7}},
                //no GP
                {new GenotypeBuilder("s", diploid).make(), null},
                {null, null},
                //missing values
                {new GenotypeBuilder("s", diploid).attribute("GP", ".,.,.").make(), null},
                {new GenotypeBuilder("s", diploid).attribute("GP", "0.1,.,0.2").make(), null},
                //number of values does not match ploidy
                {new GenotypeBuilder("s", diploid).attribute("GP", "0.3,0.7").make(), null},
                {new GenotypeBuilder("s", haploid).attribute("GP", "0.1,0.7,0.2").make(), null},
                {new GenotypeBuilder("s", diploid).attribute("GP", "0.1,0.2,0.3,0.4").make(), null}
        };
    }

    @Test(dataProvider = "getGenotypeProbabilitiesDataProvider")
    void testGetGenotypeProbabilities(final Genotype genotype, final double[] expected) {
        final VariantContext vc = new VariantContextBuilder("vc", "chr1", 100, 100, Arrays.asList(Allele.REF_A, Allele.ALT_C)).make();
        final double[] actual = EvaluateGenotypingPerformance.getGenotypeProbabilities(genotype, vc);
        if (expected == null) {
            Assert.assertNull(actual);
        } else {
            Assert.assertEquals(actual, expected, 1e-9);
        }
    }

    @DataProvider(name = "getGenotypeProbabilitiesBadInputDataProvider")
    Object[][] getGenotypeProbabilitiesBadInputDataProvider() {
        final List<Allele> diploid = Arrays.asList(Allele.REF_A, Allele.ALT_C);
        return new Object[][] {
                {new GenotypeBuilder("s", diploid).attribute("GP", "0,30,60").make()}, //phred-scaled
                {new GenotypeBuilder("s", diploid).attribute("GP", Arrays.asList(0, 30, 60)).make()},
                {new GenotypeBuilder("s", diploid).attribute("GP", "-0.1,0.6,0.5").make()},
                {new GenotypeBuilder("s", diploid).attribute("GP", "a,b,c").make()}
        };
    }

    @Test(dataProvider = "getGenotypeProbabilitiesBadInputDataProvider", expectedExceptions = UserException.BadInput.class)
    void testGetGenotypeProbabilitiesBadInput(final Genotype genotype) {
        final VariantContext vc = new VariantContextBuilder("vc", "chr1", 100, 100, Arrays.asList(Allele.REF_A, Allele.ALT_C)).make();
        EvaluateGenotypingPerformance.getGenotypeProbabilities(genotype, vc);
    }

    private static List<List<EvaluateGenotypingPerformance.GPCalibrationMetrics>> makeCalibrations(final int nBins) {
        final List<List<EvaluateGenotypingPerformance.GPCalibrationMetrics>> calibrations = new ArrayList<>();
        for (int i = 0; i < nBins; i++) {
            final List<EvaluateGenotypingPerformance.GPCalibrationMetrics> thisBin = new ArrayList<>();
            for (final GenotypeType genotypeType : EvaluateGenotypingPerformance.GP_GENOTYPE_TYPES) {
                thisBin.add(new EvaluateGenotypingPerformance.GPCalibrationMetrics("s", (i + 0.5) / nBins, VariantContext.Type.SNP, genotypeType));
            }
            calibrations.add(thisBin);
        }
        return calibrations;
    }

    @Test
    void testUpdateGPCalibrationDiploid() {
        final VariantContext vc = new VariantContextBuilder("vc", "chr1", 100, 100, Arrays.asList(Allele.REF_A, Allele.ALT_C)).make();
        final List<List<EvaluateGenotypingPerformance.GPCalibrationMetrics>> calibrations = makeCalibrations(2);
        final Genotype eval = new GenotypeBuilder("s", Arrays.asList(Allele.REF_A, Allele.ALT_C)).attribute("GP", "0.1,0.8,0.1").make();
        final Genotype truth = new GenotypeBuilder("s", Arrays.asList(Allele.REF_A, Allele.ALT_C)).make();
        EvaluateGenotypingPerformance.updateGPCalibration(calibrations, eval, truth, vc, 2);

        //HOM_REF and HOM_VAR (gp 0.1) in bin 0, incorrect; HET (gp 0.8) in bin 1, correct
        assertCalibration(calibrations.get(0).get(0), 0, 1);
        assertCalibration(calibrations.get(0).get(1), 0, 0);
        assertCalibration(calibrations.get(0).get(2), 0, 1);
        assertCalibration(calibrations.get(1).get(0), 0, 0);
        assertCalibration(calibrations.get(1).get(1), 1, 1);
        assertCalibration(calibrations.get(1).get(2), 0, 0);
    }

    @Test
    void testUpdateGPCalibrationHaploid() {
        final VariantContext vc = new VariantContextBuilder("vc", "chr1", 100, 100, Arrays.asList(Allele.REF_A, Allele.ALT_C)).make();
        final List<List<EvaluateGenotypingPerformance.GPCalibrationMetrics>> calibrations = makeCalibrations(2);
        final Genotype eval = new GenotypeBuilder("s", Collections.singletonList(Allele.ALT_C)).attribute("GP", "0.2,0.8").make();
        final Genotype truth = new GenotypeBuilder("s", Collections.singletonList(Allele.ALT_C)).make();
        EvaluateGenotypingPerformance.updateGPCalibration(calibrations, eval, truth, vc, 2);

        //REF (gp 0.2) maps to HOM_REF in bin 0, incorrect; ALT (gp 0.8) maps to HOM_VAR in bin 1, correct; HET is never incremented
        assertCalibration(calibrations.get(0).get(0), 0, 1);
        assertCalibration(calibrations.get(1).get(2), 1, 1);
        for (final List<EvaluateGenotypingPerformance.GPCalibrationMetrics> bin : calibrations) {
            assertCalibration(bin.get(1), 0, 0);
        }
        assertCalibration(calibrations.get(0).get(2), 0, 0);
        assertCalibration(calibrations.get(1).get(0), 0, 0);
    }

    @Test
    void testUpdateGPCalibrationNullTruthIsHomRef() {
        final VariantContext vc = new VariantContextBuilder("vc", "chr1", 100, 100, Arrays.asList(Allele.REF_A, Allele.ALT_C)).make();
        final List<List<EvaluateGenotypingPerformance.GPCalibrationMetrics>> calibrations = makeCalibrations(1);
        final Genotype eval = new GenotypeBuilder("s", Arrays.asList(Allele.REF_A, Allele.REF_A)).attribute("GP", "0.9,0.1,0.0").make();
        EvaluateGenotypingPerformance.updateGPCalibration(calibrations, eval, null, vc, 1);
        assertCalibration(calibrations.get(0).get(0), 1, 1);
        assertCalibration(calibrations.get(0).get(1), 0, 1);
        assertCalibration(calibrations.get(0).get(2), 0, 1);
    }

    private static void assertCalibration(final EvaluateGenotypingPerformance.GPCalibrationMetrics metrics, final int expectedCorrect, final int expectedTotal) {
        Assert.assertEquals(metrics.CORRECT, expectedCorrect, metrics.GTYPE + " at " + metrics.GP_BIN_CENTER);
        Assert.assertEquals(metrics.TOTAL, expectedTotal, metrics.GTYPE + " at " + metrics.GP_BIN_CENTER);
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