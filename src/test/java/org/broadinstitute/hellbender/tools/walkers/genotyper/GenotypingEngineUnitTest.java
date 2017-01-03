package org.broadinstitute.hellbender.tools.walkers.genotyper;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc.GeneralPloidyFailOverAFCalculatorProvider;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import org.testng.Assert;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class GenotypingEngineUnitTest extends BaseTest {

    private GenotypingEngine<UnifiedArgumentCollection> genotypingEngine;
    private final Allele refAllele = Allele.create("TCCTTCCTTCCCTCCCTCCCTC", true);
    private final Allele altT = Allele.create("T");
    private List<Allele> allelesDel;
    private  final List<Allele> gtAlleles = GATKVariantContextUtils.noCallAlleles(2);

    @BeforeTest
    public void init() {
        final GenotypeCalculationArgumentCollection genotypeArgs = new GenotypeCalculationArgumentCollection();
        final UnifiedArgumentCollection uac = new UnifiedArgumentCollection();
        uac.genotypeArgs = new GenotypeCalculationArgumentCollection(genotypeArgs);
        final SampleList samples = new IndexedSampleList("test");
        genotypingEngine = new MinimalGenotypingEngine(uac, samples, new GeneralPloidyFailOverAFCalculatorProvider(genotypeArgs));

        allelesDel = new ArrayList<>(Arrays.asList(refAllele,
                Allele.create("TCTTTCCTTCCCTCCCTCCCTCCCTCCCTTCCTTCCCTCCCTCCCTC"),
                Allele.create("TCCCTCCCTCCCTTCCTTCCCTCCCTCCCTC"),
                altT,
                GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE));

        final int deletionSize = refAllele.length() - altT.length();

        final int start = 1;
        final VariantContext deletionVC = new VariantContextBuilder("testDeletion", "1", start, start + deletionSize, allelesDel).make();
        genotypingEngine.recordDeletion(deletionSize, deletionVC);
    }

    @DataProvider(name="testCoveredByDeletionData")
    public Object[][] testCoveredByDeletionData() {
        return new Object[][] {
                {"test1", new SimpleInterval("1", 1, 1), false},
                {"test2", new SimpleInterval("1", 2, 2), true},
                {"test3", new SimpleInterval("1", 40, 40), false},
                {"test4", new SimpleInterval("2", 9, 9), false}
        };
    }

    @Test(dataProvider = "testCoveredByDeletionData")
    public void testCoveredByDeletion(final String test, final SimpleInterval location, final boolean isCovered) {
        final VariantContext vc = new VariantContextBuilder("test", location.getContig(), location.getStart(), location.getEnd(), allelesDel).make();
        Assert.assertEquals(isCovered, genotypingEngine.isVcCoveredByDeletion(vc), test + " failed");
    }

    @Test(dependsOnMethods={"testCoveredByDeletion"})  // Want to run testCoveredByDeletion first to ensure clean list or recorded deletions
    public void testCalculateGenotypes() {
        genotypingEngine.clearUpstreamDeletionsLoc();

         // Remove deletion
        final List<Genotype> genotypes = Arrays.asList(
                new GenotypeBuilder("sample1").alleles(gtAlleles).PL(new double[]{0, 0, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}).make(),  // first alt
                new GenotypeBuilder("sample2").alleles(gtAlleles).PL(new double[]{0, 0, 0, 0, 0, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0}).make()); // second alt
        final VariantContext vc = new VariantContextBuilder("test", "1",1, refAllele.length(), allelesDel).genotypes(genotypes).make();
        final VariantContext vcOut = genotypingEngine.calculateGenotypes(vc, GenotypeLikelihoodsCalculationModel.INDEL, null);
        Assert.assertFalse(vcOut.getAlleles().contains(altT));

        // Make sure the spanning deletion is removed since the deletion was removed
        final Allele refAlleleSpanDel = Allele.create("C", true);
        final List<Allele> vcAllelesSpanDel = new ArrayList<>(Arrays.asList(refAlleleSpanDel, Allele.SPAN_DEL, GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE));
        final List<Genotype> genotypesSpanDel = Arrays.asList(
                new GenotypeBuilder("sample1").alleles(gtAlleles).PL(new double[]{0, 0, 100, 0, 0, 0}).make(), // first alt
                new GenotypeBuilder("sample2").alleles(gtAlleles).PL(new double[]{0, 0, 0, 0, 0, 0}).make());
        final VariantContext vcSpanDel = new VariantContextBuilder("test1", "1",2, 2 + refAlleleSpanDel.length() - 1, vcAllelesSpanDel).
                genotypes(genotypesSpanDel).make();
        final VariantContext vcOut1 = genotypingEngine.calculateGenotypes(vcSpanDel, GenotypeLikelihoodsCalculationModel.INDEL, null);
        Assert.assertFalse(vcOut1.getAlleles().contains(Allele.SPAN_DEL));
    }
}
