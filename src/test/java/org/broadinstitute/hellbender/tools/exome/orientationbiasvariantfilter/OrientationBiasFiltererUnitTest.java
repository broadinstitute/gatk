package org.broadinstitute.hellbender.tools.exome.orientationbiasvariantfilter;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.tools.picard.analysis.artifacts.Transition;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

public class OrientationBiasFiltererUnitTest extends BaseTest {
    private static final File TEST_RESOURCE_DIR = new File("src/test/resources/org/broadinstitute/hellbender/tools/exome/orientationbiasvariantfilter/");
    public static final String smallM2Vcf = TEST_RESOURCE_DIR.getAbsolutePath() + "/small_m2.vcf";
    public static final String smallM2VcfMore = TEST_RESOURCE_DIR.getAbsolutePath() + "/small_m2_more_variants.vcf";

    @Test
    public void testAnnotateVariantContextWithPreprocessingValues() {
        final FeatureDataSource<VariantContext> featureDataSource = new FeatureDataSource<>(new File(smallM2Vcf));
        SortedSet<Transition> relevantTransitions = new TreeSet<>();
        relevantTransitions.add(Transition.transitionOf('G', 'T'));

        final Map<Transition, Double> preAdapterQFakeScoreMap = new HashMap<>();
        final double amGTPreAdapterQ = 20.0;
        preAdapterQFakeScoreMap.put(Transition.transitionOf('G', 'T'), amGTPreAdapterQ);  // preAdapterQ suppression will do nothing.

        for (final VariantContext vc : featureDataSource) {
            final VariantContext updatedVariantContext = OrientationBiasFilterer.annotateVariantContextWithPreprocessingValues(vc, relevantTransitions, preAdapterQFakeScoreMap);

            final Genotype genotypeTumor = updatedVariantContext.getGenotype("TUMOR");
            final Genotype genotypeNormal = updatedVariantContext.getGenotype("NORMAL");

            Assert.assertNotEquals(genotypeTumor.getExtendedAttribute(OrientationBiasFilterConstants.FOB, VCFConstants.EMPTY_ALLELE), VCFConstants.EMPTY_ALLELE);
            Assert.assertNotEquals(genotypeTumor.getExtendedAttribute(OrientationBiasFilterConstants.P_ARTIFACT_FIELD_NAME, VCFConstants.EMPTY_ALLELE), VCFConstants.EMPTY_ALLELE);

            assertArtifact(amGTPreAdapterQ, genotypeTumor, Transition.transitionOf('G', 'T'));

            // The NORMAL is always ref/ref in the example file.
            assertNormal(genotypeNormal);
        }
    }

    private boolean assertArtifact(double amPreAdapterQ, final Genotype genotypeTumor, final Transition transition) {
        final Transition transitionComplement = transition.complement();

        boolean result = false;

        // Check whether this genotype is reverse complement or actual artifact mode
        if (genotypeTumor.getAllele(0).basesMatch(String.valueOf(transition.ref())) && genotypeTumor.getAllele(1).basesMatch(String.valueOf(transition.call()))) {
            // not complement (i.e. artifact mode)
            Assert.assertTrue(genotypeTumor.getExtendedAttribute(OrientationBiasFilterConstants.PRE_ADAPTER_METRIC_RC_FIELD_NAME).equals(OrientationBiasFilterer.PRE_ADAPTER_METRIC_NOT_ARTIFACT_SCORE));
            Assert.assertTrue(genotypeTumor.getExtendedAttribute(OrientationBiasFilterConstants.PRE_ADAPTER_METRIC_FIELD_NAME).equals(amPreAdapterQ));
            Assert.assertEquals(genotypeTumor.getExtendedAttribute(OrientationBiasFilterConstants.IS_ORIENTATION_BIAS_ARTIFACT_MODE), String.valueOf(true));
            Assert.assertEquals(genotypeTumor.getExtendedAttribute(OrientationBiasFilterConstants.IS_ORIENTATION_BIAS_RC_ARTIFACT_MODE), String.valueOf(false));
            result = true;

        } else if (genotypeTumor.getAllele(0).basesMatch(String.valueOf(transitionComplement.ref())) && genotypeTumor.getAllele(1).basesMatch(String.valueOf(transitionComplement.call()))) {
            //complement
            Assert.assertTrue(genotypeTumor.getExtendedAttribute(OrientationBiasFilterConstants.PRE_ADAPTER_METRIC_RC_FIELD_NAME).equals(amPreAdapterQ));
            Assert.assertTrue(genotypeTumor.getExtendedAttribute(OrientationBiasFilterConstants.PRE_ADAPTER_METRIC_FIELD_NAME).equals(OrientationBiasFilterer.PRE_ADAPTER_METRIC_NOT_ARTIFACT_SCORE));
            Assert.assertEquals(genotypeTumor.getExtendedAttribute(OrientationBiasFilterConstants.IS_ORIENTATION_BIAS_ARTIFACT_MODE), String.valueOf(false));
            Assert.assertEquals(genotypeTumor.getExtendedAttribute(OrientationBiasFilterConstants.IS_ORIENTATION_BIAS_RC_ARTIFACT_MODE), String.valueOf(true));
            result = true;
        }

        return result;
    }

    @Test
    public void testAnnotateVariantContextWithPreprocessingValuesMultiArtifact() {
        final FeatureDataSource<VariantContext> featureDataSource = new FeatureDataSource<>(new File(smallM2VcfMore));
        SortedSet<Transition> relevantTransitions = new TreeSet<>();
        relevantTransitions.add(Transition.transitionOf('G', 'T'));
        relevantTransitions.add(Transition.transitionOf('C', 'T'));

        final Map<Transition, Double> preAdapterQFakeScoreMap = new HashMap<>();
        final double amGTPreAdapterQ = 20.0;
        final double amCTPreAdapterQ = 25.0;
        preAdapterQFakeScoreMap.put(relevantTransitions.first(), amGTPreAdapterQ);  // preAdapterQ suppression will do nothing.
        preAdapterQFakeScoreMap.put(relevantTransitions.last(), amCTPreAdapterQ);  // preAdapterQ suppression will do nothing.

        for (final VariantContext vc : featureDataSource) {
            final VariantContext updatedVariantContext = OrientationBiasFilterer.annotateVariantContextWithPreprocessingValues(vc, relevantTransitions, preAdapterQFakeScoreMap);

            final Genotype genotypeTumor = updatedVariantContext.getGenotype("TUMOR");
            final Genotype genotypeNormal = updatedVariantContext.getGenotype("NORMAL");

            // This is mostly just to make sure that nobody breaks the test itself.  I.e. that this test will test all tumor genotype paths be artifact or non-artifact.
            boolean wasGenotypeTumorTested = false;

            // Check whether this genotype is reverse complement or actual artifact mode
            wasGenotypeTumorTested |= assertArtifact(amGTPreAdapterQ, genotypeTumor, relevantTransitions.first());
            wasGenotypeTumorTested |= assertArtifact(amCTPreAdapterQ, genotypeTumor, relevantTransitions.last());

            // Check any variants that are not an artifact mode but are SNP
            if (!OrientationBiasUtils.isGenotypeInTransitionsWithComplement(genotypeTumor, relevantTransitions)) {
                assertNotTransition(genotypeTumor);
                wasGenotypeTumorTested = true;
            } else {

                // Check attributes common to all variants in artifact mode
                Assert.assertNotEquals(genotypeTumor.getExtendedAttribute(OrientationBiasFilterConstants.FOB, VCFConstants.EMPTY_ALLELE), VCFConstants.EMPTY_ALLELE);
                Assert.assertNotEquals(genotypeTumor.getExtendedAttribute(OrientationBiasFilterConstants.P_ARTIFACT_FIELD_NAME, VCFConstants.EMPTY_ALLELE), VCFConstants.EMPTY_ALLELE);
            }

            // The NORMAL is always ref/ref in the example file.
            assertNormal(genotypeNormal);

            Assert.assertTrue(wasGenotypeTumorTested, "The test seems to be broken...  A variant context was tested, but it had no tumor genotype.");
        }
    }

    private void assertNotTransition(final Genotype genotype) {
        Assert.assertEquals(genotype.getExtendedAttribute(OrientationBiasFilterConstants.IS_ORIENTATION_BIAS_ARTIFACT_MODE), String.valueOf(false));
        Assert.assertEquals(genotype.getExtendedAttribute(OrientationBiasFilterConstants.IS_ORIENTATION_BIAS_RC_ARTIFACT_MODE), String.valueOf(false));
        Assert.assertEquals(genotype.getExtendedAttribute(OrientationBiasFilterConstants.PRE_ADAPTER_METRIC_RC_FIELD_NAME), OrientationBiasFilterer.PRE_ADAPTER_METRIC_NOT_ARTIFACT_SCORE);
        Assert.assertEquals(genotype.getExtendedAttribute(OrientationBiasFilterConstants.PRE_ADAPTER_METRIC_FIELD_NAME), OrientationBiasFilterer.PRE_ADAPTER_METRIC_NOT_ARTIFACT_SCORE);
        Assert.assertEquals(genotype.getExtendedAttribute(OrientationBiasFilterConstants.FOB), VCFConstants.EMPTY_ALLELE);
        Assert.assertEquals(genotype.getExtendedAttribute(OrientationBiasFilterConstants.P_ARTIFACT_FIELD_NAME), VCFConstants.EMPTY_ALLELE);
    }

    private void assertNormal(final Genotype genotypeNormal) {
        Assert.assertEquals(genotypeNormal.getExtendedAttribute(OrientationBiasFilterConstants.IS_ORIENTATION_BIAS_ARTIFACT_MODE), String.valueOf(false));
        Assert.assertEquals(genotypeNormal.getExtendedAttribute(OrientationBiasFilterConstants.IS_ORIENTATION_BIAS_RC_ARTIFACT_MODE), String.valueOf(false));
        Assert.assertNull(genotypeNormal.getExtendedAttribute(OrientationBiasFilterConstants.PRE_ADAPTER_METRIC_RC_FIELD_NAME));
        Assert.assertNull(genotypeNormal.getExtendedAttribute(OrientationBiasFilterConstants.PRE_ADAPTER_METRIC_FIELD_NAME));
        Assert.assertNull(genotypeNormal.getExtendedAttribute(OrientationBiasFilterConstants.FOB));
        Assert.assertNull(genotypeNormal.getExtendedAttribute(OrientationBiasFilterConstants.P_ARTIFACT_FIELD_NAME));
    }

    @Test
    public void testAnnotateVariantContextWithFilterValuesMultiArtifact() {
        final FeatureDataSource<VariantContext> featureDataSource = new FeatureDataSource<>(new File(smallM2VcfMore));
        SortedSet<Transition> relevantTransitions = new TreeSet<>();
        relevantTransitions.add(Transition.transitionOf('G', 'T'));
        relevantTransitions.add(Transition.transitionOf('C', 'T'));

        final Map<Transition, Double> preAdapterQFakeScoreMap = new HashMap<>();
        final double amGTPreAdapterQ = 20.0;
        final double amCTPreAdapterQ = 25.0;
        preAdapterQFakeScoreMap.put(relevantTransitions.first(), amGTPreAdapterQ);  // preAdapterQ suppression will do nothing.
        preAdapterQFakeScoreMap.put(relevantTransitions.last(), amCTPreAdapterQ);  // preAdapterQ suppression will do nothing.
        final List<VariantContext> updatedVariants = new ArrayList<>();

        for (final VariantContext vc : featureDataSource) {
            final VariantContext updatedVariantContext = OrientationBiasFilterer.annotateVariantContextWithPreprocessingValues(vc, relevantTransitions, preAdapterQFakeScoreMap);
            updatedVariants.add(updatedVariantContext);
        }
        final List<String> sampleNames = updatedVariants.get(0).getSampleNamesOrderedByName();

        // Create a mapping from sample name to a genotype->variant context map
        final Map<String, SortedMap<Genotype, VariantContext>> sampleNameToVariants = OrientationBiasFilterer.createSampleToGenotypeVariantContextSortedMap(sampleNames, updatedVariants);
    }



    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testTotalNumToCutCalculationFDRThreshNegative() {
        OrientationBiasFilterer.calculateTotalNumToCut(-0.5, 200, Collections.nCopies(20, 1.0) );
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testTotalNumToCutCalculationFDRThreshZero() {
        OrientationBiasFilterer.calculateTotalNumToCut(0.0, 200, Collections.nCopies(20, 1.0) );
    }

    @Test
    public void testCreateSampleToGenotypeVCMap() {

        // Setup the test
        final FeatureDataSource<VariantContext> featureDataSource = new FeatureDataSource<>(new File(smallM2VcfMore));
        SortedSet<Transition> relevantTransitions = new TreeSet<>();
        relevantTransitions.add(Transition.transitionOf('G', 'T'));
        relevantTransitions.add(Transition.transitionOf('C', 'T'));

        final Map<Transition, Double> preAdapterQFakeScoreMap = new HashMap<>();
        final double amGTPreAdapterQ = 20.0;
        final double amCTPreAdapterQ = 25.0;
        preAdapterQFakeScoreMap.put(relevantTransitions.first(), amGTPreAdapterQ);  // preAdapterQ suppression will do nothing.
        preAdapterQFakeScoreMap.put(relevantTransitions.last(), amCTPreAdapterQ);  // preAdapterQ suppression will do nothing.
        final List<VariantContext> updatedVariants = new ArrayList<>();

        for (final VariantContext vc : featureDataSource) {
            final VariantContext updatedVariantContext = OrientationBiasFilterer.annotateVariantContextWithPreprocessingValues(vc, relevantTransitions, preAdapterQFakeScoreMap);
            updatedVariants.add(updatedVariantContext);
        }
        final List<String> sampleNames = updatedVariants.get(0).getSampleNamesOrderedByName();

        // Do the test
        // Create a mapping from sample name to a genotype->variant context map with the second map sorted by p_artifact (i.e. OrientationBiasFilterConstants.P_ARTIFACT_FIELD_NAME)
        final Map<String, SortedMap<Genotype, VariantContext>> sampleNameToVariants = OrientationBiasFilterer.createSampleToGenotypeVariantContextSortedMap(sampleNames, updatedVariants);
        Assert.assertEquals(sampleNameToVariants.keySet().size(), 2);
        Assert.assertTrue(sampleNameToVariants.keySet().contains("TUMOR"));
        Assert.assertTrue(sampleNameToVariants.keySet().contains("NORMAL"));

        Assert.assertEquals(sampleNameToVariants.get("TUMOR").keySet().size(), 8);

        // None of the normal genotypes should have a pvalue, so cannot/shouldn't be added to the sorted map
        Assert.assertEquals(sampleNameToVariants.get("NORMAL").keySet().size(), 0);

        // Check that the sorted map is getting smaller (or same) values of p_artifact and not staying put.
        double previousPArtifact = Double.POSITIVE_INFINITY;
        for (final Genotype genotypeTumor : sampleNameToVariants.get("TUMOR").keySet()) {
            final Double pArtifact = OrientationBiasUtils.getGenotypeDouble(genotypeTumor, OrientationBiasFilterConstants.P_ARTIFACT_FIELD_NAME, Double.POSITIVE_INFINITY);
            Assert.assertNotNull(pArtifact);
            Assert.assertTrue(pArtifact <= previousPArtifact);
            Assert.assertNotEquals(pArtifact, Double.POSITIVE_INFINITY);
        }
    }

    @Test(dataProvider = "testTotalNumToCutCalculation")
    public void testTotalNumToCutCalculation(Pair<Double, List<Double>> fdrThreshAndPArtifact, int gt) {
        final double fdrThresh = fdrThreshAndPArtifact.getLeft();
        final List<Double> pArtifactScores = fdrThreshAndPArtifact.getRight();
        final int unfilteredGenotypeCount = pArtifactScores.size();

        final int guess = OrientationBiasFilterer.calculateTotalNumToCut(fdrThresh, unfilteredGenotypeCount, pArtifactScores);
        Assert.assertEquals(guess, gt);
    }

    @DataProvider(name = "testTotalNumToCutCalculation")
    public Object[][] testTotalNumToCutCalculation() {

        // This test does not include the suppression factor (based on preAdapterQ)
        final List<Double> fakePArtifactValues = new ArrayList<>();
        fakePArtifactValues.addAll(Collections.nCopies(99, 1.0));
        fakePArtifactValues.addAll(Collections.nCopies(1, 0.001));
        fakePArtifactValues.addAll(Collections.nCopies(100, 0.0));

        final List<Double> fakePArtifactValuesBigger = new ArrayList<>();
        fakePArtifactValuesBigger.addAll(fakePArtifactValues);
        fakePArtifactValuesBigger.addAll(Collections.nCopies(1000, 0.0));

        final Object [] [] result = new Object[][]{
                // fdrThresh, Array of doubles for p_Artifact -- including non-artifact
                // Results double checked against matlab implementation, manually
                {Pair.of(0.01, fakePArtifactValues), 99},
                {Pair.of(0.01, fakePArtifactValuesBigger), 100}
        };

        return result;
    }
}
