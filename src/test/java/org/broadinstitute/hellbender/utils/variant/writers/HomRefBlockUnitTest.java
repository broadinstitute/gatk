package org.broadinstitute.hellbender.utils.variant.writers;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.variant.HomoSapiensConstants;
import org.testng.Assert;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class HomRefBlockUnitTest extends BaseTest {
    private static final String SAMPLE_NAME = "foo";
    private static final Allele REF = Allele.create("A", true);
    VariantContext vc;

    @BeforeMethod
    public void setUp() throws Exception {
        vc = new VariantContextBuilder(SAMPLE_NAME, "20", 1, 1, getAlleles()).make();
    }

    private List<Allele> getAlleles() {
        return Arrays.asList(REF, Allele.create("C"));
    }

    @Test
    public void testBasicConstruction() {
        final HomRefBlock band = new HomRefBlock(vc, 10, 20, HomoSapiensConstants.DEFAULT_PLOIDY);
        Assert.assertSame(band.getStartingVC(), vc);
        Assert.assertEquals(band.getRef(), vc.getReference());
        Assert.assertEquals(band.getGQLowerBound(), 10);
        Assert.assertEquals(band.getGQUpperBound(), 20);
        Assert.assertEquals(band.withinBounds(1), false);
        Assert.assertEquals(band.withinBounds(10), true);
        Assert.assertEquals(band.withinBounds(11), true);
        Assert.assertEquals(band.withinBounds(20), false);
        Assert.assertEquals(band.withinBounds(21), false);
    }

    @Test
    public void testMinMedian() {
        //TODO - might be better to make this test use a data provider?
        final HomRefBlock band = new HomRefBlock(vc, 10, 20, HomoSapiensConstants.DEFAULT_PLOIDY);
        final GenotypeBuilder gb = new GenotypeBuilder("NA12878");
        gb.alleles(vc.getAlleles());

        int pos = vc.getStart();
        band.add(pos++, gb.DP(10).GQ(11).PL(new int[]{0,11,100}).make());
        Assert.assertEquals(band.getEnd(), pos - 1);
        assertValues(band, 10, 10);

        band.add(pos++, gb.DP(11).GQ(10).PL(new int[]{0,10,100}).make());
        Assert.assertEquals(band.getEnd(), pos - 1);
        assertValues(band, 10, 11);

        band.add(pos++, gb.DP(12).GQ(12).PL(new int[]{0,12,100}).make());
        Assert.assertEquals(band.getEnd(), pos - 1);
        assertValues(band, 10, 11);

        band.add(pos++, gb.DP(13).GQ(15).PL(new int[]{0,15,100}).make());
        Assert.assertEquals(band.getEnd(), pos - 1);
        band.add(pos++, gb.DP(14).GQ(16).PL(new int[]{0,16,100}).make());
        Assert.assertEquals(band.getEnd(), pos - 1);
        band.add(pos++, gb.DP(15).GQ(17).PL(new int[]{0,17,100}).make());
        Assert.assertEquals(band.getEnd(), pos - 1);
        band.add(pos++, gb.DP(16).GQ(18).PL(new int[]{0,18,100}).make());
        Assert.assertEquals(band.getEnd(), pos - 1);
        assertValues(band, 10, 13);
        Assert.assertEquals(band.getSize(), pos - vc.getStart());
        Assert.assertTrue(Arrays.equals(band.getMinPLs(), new int[]{0,10,100}));
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testBadAdd() {
        final HomRefBlock band = new HomRefBlock(vc, 10, 20, HomoSapiensConstants.DEFAULT_PLOIDY);
        final GenotypeBuilder gb = new GenotypeBuilder("NA12878");

        band.add(vc.getStart() + 10, gb.DP(10).GQ(11).PL(new int[]{0,10,100}).make());
    }

    private void assertValues(final HomRefBlock band, final int minDP, final int medianDP) {
        Assert.assertEquals(band.getMinDP(), minDP);
        Assert.assertEquals(band.getMedianDP(), medianDP);
    }


    @DataProvider(name = "ContiguousData")
    public Object[][] makeContiguousData() {
        List<Object[]> tests = new ArrayList<>();

        for ( final String chrMod : Arrays.asList("", ".mismatch") ) {
            for ( final int offset : Arrays.asList(-10, -1, 0, 1, 10) ) {
                final boolean equals = chrMod.equals("") && offset == 0;
                tests.add(new Object[]{vc.getContig() + chrMod, vc.getStart() + offset, equals});
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "ContiguousData")
    public void testIsContiguous(final String contig, final int pos, final boolean expected) {
        final HomRefBlock band = new HomRefBlock(vc, 10, 20, HomoSapiensConstants.DEFAULT_PLOIDY);
        final VariantContext testVC = new VariantContextBuilder(vc).chr(contig).start(pos).stop(pos).make();
        Assert.assertEquals(band.isContiguous(testVC), expected);
    }

    @DataProvider(name = "minPLs")
    public Object[][] getMinPLs(){
        return new Object[][] {
                {new int[] {0, 1, 5}, 1},
                {new int[] {10, 50, 30}, 20},
                {new int[] {0, 10, 20, 10, 5, 50}, 5}
        };
    }

    @DataProvider(name = "badMinPls")
    public Object[][] getBadMinPLs(){
        return new Object[][] {
                {new int[] {0}}, // too few PLs
                {new int[] {20, 1, 0}}, //first should be lowest since this is a homRefBlock
                {new int[] {}},
        };
    }

    @Test(dataProvider = "minPLs")
    public void testGenotypeQualityFromPls(int[] minPLs, int expected){
        Assert.assertEquals(HomRefBlock.genotypeQualityFromPLs(minPLs), expected);
    }

    @Test(dataProvider = "badMinPls", expectedExceptions = GATKException.class)
    public void testGenotypeQualityFromPLsBadPLs(int[] minPLs){
        HomRefBlock.genotypeQualityFromPLs(minPLs); //this should explode
    }

    @Test
    public void testToVariantContext(){
        final GenotypeBuilder gb = new GenotypeBuilder(SAMPLE_NAME, vc.getAlleles());
        final Genotype genotype1 = gb.GQ(15).DP(6).PL(new int[]{0, 10, 100}).make();
        final Genotype genotype2 = gb.GQ(17).DP(10).PL(new int[]{0, 5, 80}).make();

        HomRefBlock block = new HomRefBlock(vc, 10, 20, HomoSapiensConstants.DEFAULT_PLOIDY);
        block.add(vc.getEnd(), genotype1);
        block.add(vc.getEnd() + 1, genotype2);

        final VariantContext newVc = block.toVariantContext(SAMPLE_NAME);
        Assert.assertEquals(newVc.getGenotypes().size(), 1);
        final Genotype genotype = newVc.getGenotypes().get(0);
        Assert.assertEquals(genotype.getDP(),8); //dp should be median of the added DPs
        Assert.assertEquals(genotype.getGQ(), 5); //GQ should have been recalculated with the minPls
        Assert.assertTrue(genotype.getAlleles().stream().allMatch(a -> a.equals(REF)));
    }
}
