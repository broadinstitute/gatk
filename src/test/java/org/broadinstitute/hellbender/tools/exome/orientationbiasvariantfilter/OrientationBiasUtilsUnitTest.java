package org.broadinstitute.hellbender.tools.exome.orientationbiasvariantfilter;


import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.tools.picard.analysis.artifacts.Transition;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;
import java.util.stream.IntStream;


public class OrientationBiasUtilsUnitTest extends BaseTest {
    private static final File TEST_RESOURCE_DIR = new File("src/test/resources/org/broadinstitute/hellbender/tools/exome/orientationbiasvariantfilter/");
    public static final String smallM2VcfMore = TEST_RESOURCE_DIR.getAbsolutePath() + "/small_m2_more_variants.vcf";

    @Test(dataProvider = "BasicTransitionsWithRC")
    public void testFindReverseComplementTransitions(List<Transition> transition, List<Transition> reverseComplement) {

        final List<Transition> result = OrientationBiasUtils.createReverseComplementTransitions(transition);
        IntStream.range(0, result.size()).forEachOrdered(i -> Assert.assertEquals(result.get(i), reverseComplement.get(i)));
    }



    @DataProvider(name = "BasicTransitionsWithRC")
    public Object[][] basicTransitionsWithRC() {
        return new Object[][]{
                {Collections.singletonList(Transition.transitionOf('C', 'A')), Collections.singletonList(Transition.transitionOf('G', 'T'))},
                {Collections.singletonList(Transition.transitionOf('A', 'T')), Collections.singletonList(Transition.transitionOf('T', 'A'))},
                {Arrays.asList(Transition.transitionOf('A', 'T'), Transition.transitionOf('C', 'T')), Arrays.asList(Transition.transitionOf('T', 'A'), Transition.transitionOf('G', 'A'))}
        };
    }

    @Test
    public void testCountNumTransition() {

        // Setup the test
        final FeatureDataSource<VariantContext> featureDataSource = new FeatureDataSource<>(new File(smallM2VcfMore));
        SortedSet<Transition> relevantTransitions = new TreeSet<>();
        relevantTransitions.add(Transition.transitionOf('T', 'A'));
        final List<VariantContext> variantContexts = new ArrayList<>();
        for (final VariantContext vc : featureDataSource) {
            variantContexts.add(vc);
        }

        // Should be one, since one of the variants was filtered.
        Assert.assertEquals(OrientationBiasUtils.calculateNumTransition("TUMOR", variantContexts, relevantTransitions.first()), 1);
        Assert.assertEquals(OrientationBiasUtils.calculateNumTransition("NORMAL", variantContexts, relevantTransitions.first()), 0);
    }

    @Test
    public void testCountNumUnfiltered() {

        // Setup the test
        final FeatureDataSource<VariantContext> featureDataSource = new FeatureDataSource<>(new File(smallM2VcfMore));
        final List<VariantContext> variantContexts = new ArrayList<>();
        for (final VariantContext vc : featureDataSource) {
            variantContexts.add(vc);
        }
        Assert.assertEquals(OrientationBiasUtils.calculateUnfilteredNonRefGenotypeCount(variantContexts, "TUMOR"), 9);
        Assert.assertEquals(OrientationBiasUtils.calculateUnfilteredNonRefGenotypeCount(variantContexts, "NORMAL"), 0);
    }

    @Test(dataProvider = "GenotypeFiltering")
    public void testUpdatingGenotypeFilter(String inFilter, String addFilter, String updatedFilter) {
        Assert.assertEquals(OrientationBiasUtils.addFilterToGenotype(inFilter, addFilter), updatedFilter);
    }

    @DataProvider(name = "GenotypeFiltering")
    public Object[][] genotypeFiltering() {
        return new Object[][]{
                {VCFConstants.PASSES_FILTERS_v4, "Foo", "Foo"},
                {"Foo", "Bar", "Foo;Bar"},
                {"Foo", "Bar;Baz", "Foo;Bar;Baz"},
                {"Foo;Baz", "Bar", "Foo;Baz;Bar"},
                {null, "Foo", "Foo"},
                {VCFConstants.UNFILTERED, "Foo", "Foo"},
        };
    }
}
