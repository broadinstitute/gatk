package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.VariantContextTestUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.function.Function;


public class GenotypeGVCFsUnitTest extends BaseTest {

    public static final int MIN_DP = 15;
    public static final int DP = 3;
    public static final Function<GenotypeBuilder, GenotypeBuilder> ADD_DP = b -> b.DP(DP);
    public static final Allele REF = Allele.create("A", true);
    public static final Allele ALT = Allele.create("C", false);

    @DataProvider(name = "variantContexts")
    public Object[][] variantContexts() {
        return new Object[][]{
                {
                        getHetWithGenotype(Collections.singletonList(new GenotypeBuilder("SAMPLE", Arrays.asList(ALT, ALT))
                                .DP(DP)
                                .attribute(GATKVCFConstants.MIN_DP_FORMAT_KEY, MIN_DP)
                                .attribute(GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY, new int[]{1, 1, 1, 1})
                                .attribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY, "bad")
                                .PL(new int[]{50, 20, 0}).make())),
                        false,
                        Collections.singletonList(new GenotypeBuilder("SAMPLE", Arrays.asList(ALT, ALT))
                                .DP(MIN_DP)
                                .attribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY, GenotypeGVCFs.PHASED_HOM_VAR_STRING)
                                .AD(new int[]{MIN_DP, 0})
                                .PL(new int[]{50, 20, 0}).make())
                },
                {
                        getHetWithGenotype(Collections.singletonList(new GenotypeBuilder("SAMPLE", Arrays.asList(REF,REF))
                            .GQ(20)
                            .DP(MIN_DP)
                            .PL(new int[]{0,15,50}).make())),
                        true,
                        Collections.singletonList(new GenotypeBuilder("SAMPLE", Arrays.asList(REF, REF))
                                .DP(MIN_DP)
                                .AD(new int[]{MIN_DP, 0})
                                .attribute(GATKVCFConstants.REFERENCE_GENOTYPE_QUALITY, 20).make())
                }
        };
    }


    @Test(dataProvider= "variantContexts")
    public void testCleanupGenotypeAnnotations(VariantContext vc, boolean createRefGTs,  List<Genotype> expected){
        final List<Genotype> genotypes = GenotypeGVCFs.cleanupGenotypeAnnotations(vc, createRefGTs);
        VariantContextTestUtils.assertGenotypesAreEqual(genotypes.get(0), expected.get(0));
    }

    private VariantContext getHetWithGenotype(List<Genotype> genotypes){
        return new VariantContextBuilder().alleles(Arrays.asList(REF, ALT))
                .loc("1",100,100)
                .genotypes(genotypes).make();
    }

    @SafeVarargs
    private static List<Genotype> generateGenotypes(Function<GenotypeBuilder, GenotypeBuilder>... customizations){
        final List<Genotype> genotypes = new ArrayList<>();
        for(int i = 0; i< customizations.length; i++) {
            final Function<GenotypeBuilder, GenotypeBuilder> customization = customizations[i];
            final GenotypeBuilder builder = new GenotypeBuilder("Sample_" + i);
            final GenotypeBuilder applied = customization.apply(builder);
            genotypes.add(applied.make());
        }
        return genotypes;
    }

    @DataProvider
    public Object[][] getMinDPData(){
        final Function<GenotypeBuilder, GenotypeBuilder> addMinDP = b -> b.attribute(GATKVCFConstants.MIN_DP_FORMAT_KEY, MIN_DP);
        return new Object[][]{
                {getHetWithGenotype( generateGenotypes(addMinDP)), MIN_DP}, //MIN_DP only
                {getHetWithGenotype( generateGenotypes(ADD_DP)), DP}, // DP only
                {getHetWithGenotype( generateGenotypes(addMinDP.andThen(ADD_DP))), MIN_DP} //MIN_DP replaces DP
        };
    }

    @Test(dataProvider = "getMinDPData")
    public void testMinDPReplacedWithDP(VariantContext vc, int expectedDepth){
        Assert.assertEquals(GenotypeGVCFs.cleanupGenotypeAnnotations(vc, false).get(0).getDP(), expectedDepth);
        Assert.assertNull(GenotypeGVCFs.cleanupGenotypeAnnotations(vc, false).get(0).getExtendedAttribute(GATKVCFConstants.MIN_DP_FORMAT_KEY));
    }

    @Test
    public void testSBRemoved(){
        final VariantContext vcWithSB = getHetWithGenotype(generateGenotypes(b -> b.attribute(GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY, new int[]{6, 11, 11, 10})));
        Assert.assertNotNull(vcWithSB.getGenotype("Sample_0").getAnyAttribute(GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY));
        final Genotype afterCleanup = GenotypeGVCFs.cleanupGenotypeAnnotations(vcWithSB, true).get(0);
        Assert.assertNull(afterCleanup.getExtendedAttribute(GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY));
    }

    @Test
    public void testADCreated(){
        final VariantContext noAD = getHetWithGenotype(generateGenotypes(ADD_DP));
        Assert.assertNull(noAD.getGenotype("Sample_0").getAD());
        final Genotype afterCleanup = GenotypeGVCFs.cleanupGenotypeAnnotations(noAD, true).get(0);
        Assert.assertEquals(afterCleanup.getAD(), new int[]{DP, 0});
    }

    @Test
    public void testPhasingUpdate(){
        final VariantContext withPhasing = getHetWithGenotype(generateGenotypes(b ->  b.attribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY, "bad")
                .alleles(Arrays.asList(ALT,ALT)),
                b -> b.attribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY, "something").alleles(Arrays.asList(REF,REF))));
        final Genotype homVarAfterCleanup = GenotypeGVCFs.cleanupGenotypeAnnotations(withPhasing, true).get(0);
        Assert.assertEquals(homVarAfterCleanup.getAnyAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY), GenotypeGVCFs.PHASED_HOM_VAR_STRING);
        final Genotype homRefAfterCleaning = GenotypeGVCFs.cleanupGenotypeAnnotations(withPhasing, true).get(1);
        Assert.assertEquals(homRefAfterCleaning.getAnyAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY), "something");
    }

    @Test
    public void testRGQ0IsNoCall(){
        final List<Allele> noCall = GATKVariantContextUtils.noCallAlleles(2);
        final VariantContext gq0 = getHetWithGenotype(generateGenotypes(b -> b.GQ(0).DP(10).alleles(noCall), // GQ = 0
                                                                        (b -> b.DP(10).alleles(noCall))));  //no GQ
        final List<Genotype> genotypes = GenotypeGVCFs.cleanupGenotypeAnnotations(gq0, true);
        for( Genotype genotype : genotypes ){
            Assert.assertEquals(genotype.getAlleles(), noCall);
        }
    }

    @DataProvider
    public Object[][] getVariantsForIsProperlyPolymorphic(){
        return new Object[][]{
                {new VariantContextBuilder("test", "1", 1, 1, Collections.singleton(REF)).make(), false},
                {new VariantContextBuilder("test", "1", 1, 1, Arrays.asList(REF, ALT)).make(), true},
                {new VariantContextBuilder("test", "1", 1, 1, Arrays.asList(REF, GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE)).make(), false},
                {new VariantContextBuilder("test", "1", 1, 1, Arrays.asList(REF, Allele.SPAN_DEL)).make(), false},
                {new VariantContextBuilder("test", "1", 1, 2, Arrays.asList(REF, Allele.create("<SOME_SYMBOLIC_ALLELE>"))).make(), false},
                {new VariantContextBuilder("test", "1", 1, 1, Arrays.asList(REF, ALT, Allele.create("G"))).make(), true},
                {null, false}
        };
    }

    @Test(dataProvider = "getVariantsForIsProperlyPolymorphic")
    public void testIsProperlyPolymorphic(VariantContext vc, boolean expected){
        Assert.assertEquals(GenotypeGVCFs.isProperlyPolymorphic(vc), expected);
    }

    @DataProvider
    public Object[][] getSpanningAndNonSpanningAlleles(){
        return new Object[][]{
                {Allele.SPAN_DEL, true},
                {GATKVCFConstants.SPANNING_DELETION_SYMBOLIC_ALLELE_DEPRECATED, true},
                {GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE, false},
                {Allele.NO_CALL, false},
                {Allele.create("A", true), false}
        };
    }

    @Test(dataProvider = "getSpanningAndNonSpanningAlleles")
    public void testIsSpanningDeletion(Allele allele, boolean expected){
        Assert.assertEquals(GenotypeGVCFs.isSpanningDeletion(allele), expected);
    }
}
