package org.broadinstitute.hellbender.utils.variant.writers;

import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.variant.HomoSapiensConstants;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class HomRefBlockUnitTest extends GATKBaseTest {
    private static final String SAMPLE_NAME = "foo";
    private static final Allele REF = Allele.create("A", true);
    private static final int START = 1001;
    private static final int END = START + 1000;

    private static List<Allele> getAlleles() {
        return Arrays.asList(REF, Allele.create("C"));
    }

    private static Genotype getGenotype() {
        return new GenotypeBuilder(SAMPLE_NAME, Arrays.asList(REF, REF)).DP(10).make();
    }

    private static VariantContext getVariantContext() {
        return new VariantContextBuilder(SAMPLE_NAME, "20", START, START, getAlleles()).genotypes(getGenotype()).make();
    }

    @Test
    public void testBasicConstruction() {
        final VariantContext vc = getVariantContext();
        final GVCFBlock band = getHomRefBlock(vc);
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
        final VariantContext vc = getVariantContext();
        final HomRefBlock band = getHomRefBlock(vc);
        assertValues(band, 10, 10);
        Assert.assertEquals(band.getEnd(), START);

        int pos = band.getStart() + 1;

        final GenotypeBuilder gb = new GenotypeBuilder(SAMPLE_NAME);
        gb.alleles(vc.getAlleles());

        band.add(pos++, gb.DP(11).GQ(10).PL(new int[]{0, 10, 100}).make());
        Assert.assertEquals(band.getEnd(), pos - 1);  //-1 because we already did a postfix add above
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
        Assert.assertEquals(band.getMinPLs(), new int[]{0, 10, 100});
    }

    @DataProvider
    public static Object[][] badAdditions() {
        final VariantContext vc = getVariantContext();
        return new Object[][]{
                {vc.getStart(), getValidGenotypeBuilder().PL((int[])null).make()}, //no PLs
                {END + 2, getValidGenotypeBuilder().make()}, //bad start
                {vc.getStart() - 1, getValidGenotypeBuilder().make()}, //bad start
                {vc.getStart(), getValidGenotypeBuilder().GQ(1).make()}, // GQ out of bounds
                {vc.getStart(), getValidGenotypeBuilder().GQ(100).make()}, // GQ out of bounds
                {vc.getStart(), getValidGenotypeBuilder().alleles(Arrays.asList(REF, REF, REF)).make()}, //wrong ploidy
                {vc.getStart(), getValidGenotypeBuilder().noPL().noGQ().make()}, //no PL and no GQ
                {vc.getStart(), null}, //null genotype
        };
    }

    private static GenotypeBuilder getValidGenotypeBuilder() {
        return new GenotypeBuilder(SAMPLE_NAME)
                .PL(getPLArray())
                .GQ(getPLArray()[1])
                .alleles(getVariantContext().getAlleles());
    }

    @Test(dataProvider = "badAdditions", expectedExceptions = IllegalArgumentException.class)
    public void testBadAdd(int start, Genotype gb) {
        final HomRefBlock block = getHomRefBlock(getVariantContext());
        block.add(START, END, getValidGenotypeBuilder().make());
        block.add(start, gb);
    }

    @Test(expectedExceptions = GATKException.class)
    public void testCantAddDifferentNumbersOfPls(){
        final VariantContext vc = getVariantContext();
        final GVCFBlock band = getHomRefBlock(getVariantContext());
        band.add(vc.getStart()+1, getValidGenotypeBuilder().make() );
        band.add(vc.getStart() + 2, getValidGenotypeBuilder().PL(new int[] {1,2,4,5,6}).make() );
    }

    @Test
    public void testNoPLs() {
        //add VC with no PLs to block with PLs
        final VariantContext vc = getVariantContext();
        final HomRefBlock band = getHomRefBlock(getVariantContext());
        Assert.assertNull(band.getMinPLs());
        Assert.assertEquals(band.getMinGQ(), -1);
        band.add(vc.getStart() + 1, getValidGenotypeBuilder().make() );
        Assert.assertNotNull(band.getMinPLs());
        Assert.assertEquals(band.getMinPLs()[1], getPLArray()[1]);
        Assert.assertEquals(band.getMinGQ(), getPLArray()[1]);
        band.add(vc.getStart() + 2, getValidGenotypeBuilder().noPL().make());
        Assert.assertEquals(band.getSize(), 3);
        Assert.assertEquals(band.getMinPLs()[1], getPLArray()[1]);
        Assert.assertEquals(band.getMinGQ(), getPLArray()[1]);

        //add VC with no PLs to block with no PLs
        final HomRefBlock band2 = getHomRefBlock(getVariantContext());
        Assert.assertNull(band2.getMinPLs());
        Assert.assertEquals(band2.getMinGQ(), -1);
        band2.add(vc.getStart() + 1, getValidGenotypeBuilder().noPL().make() );
        Assert.assertEquals(band2.getMinGQ(), getPLArray()[1]);
        Assert.assertNull(band2.getMinPLs());
        band2.add(vc.getStart() + 2, getValidGenotypeBuilder().noPL().make());
        Assert.assertNull(band2.getMinPLs());
        Assert.assertEquals(band2.getMinGQ(), getPLArray()[1]);
        Assert.assertTrue(band2.getSize() == 3);
    }

    private static int[] getPLArray() {
        return new int[]{0,15,100};
    }

    private static void assertValues(final GVCFBlock band, final int minDP, final int medianDP) {
        Assert.assertEquals(band.getMinDP(), minDP);
        Assert.assertEquals(band.getMedianDP(), medianDP);
    }


    @DataProvider(name = "ContiguousData")
    public Object[][] makeContiguousData() {
        final List<Object[]> tests = new ArrayList<>();
        final VariantContext vc = getVariantContext();

        for ( final String chrMod : Arrays.asList("", ".mismatch") ) {
            for ( final int offset : Arrays.asList(-10, -1, 0, 1, 10) ) {
                final boolean equals = chrMod.isEmpty() && (Math.abs(offset) <= 1); //allow adding of VCs with overlap or adjacent start
                tests.add(new Object[]{vc.getContig() + chrMod, vc.getStart() + offset, equals});
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "ContiguousData")
    public void testIsContiguous(final String contig, final int pos, final boolean expected) {
        final VariantContext vc = getVariantContext();
        final GVCFBlock band = getHomRefBlock(vc);
        final VariantContext testVC = new VariantContextBuilder(vc).chr(contig).start(pos).stop(pos).make();
        Assert.assertEquals(band.isContiguous(testVC), expected);
    }

    @Test
    public void testToVariantContext(){
        final VariantContext vc = getVariantContext();
        final GenotypeBuilder gb = new GenotypeBuilder(SAMPLE_NAME, vc.getAlleles());
        final Genotype genotype1 = gb.GQ(15).DP(6).PL(new int[]{0, 10, 100}).make();
        final Genotype genotype2 = gb.GQ(17).DP(10).PL(new int[]{0, 5, 80}).make();

        final HomRefBlock block = getHomRefBlock(vc);
        block.add(vc.getEnd() + 1, genotype1);
        block.add(vc.getEnd() + 2, genotype2);

        final VariantContext newVc = block.toVariantContext(SAMPLE_NAME, false);
        Assert.assertEquals(newVc.getGenotypes().size(), 1);
        final Genotype genotype = newVc.getGenotypes().get(0);
        Assert.assertEquals(genotype.getDP(),10); //dp should be median of the added DPs and the initial DP (from getVariantContext())
        Assert.assertEquals(genotype.getGQ(), 5); //GQ should have been recalculated with the minPls
        Assert.assertTrue(genotype.getAlleles().stream().allMatch(a -> a.equals(REF)));
    }

    public static HomRefBlock getHomRefBlock(VariantContext vc) {
        return new HomRefBlock(vc, 10, 20, HomoSapiensConstants.DEFAULT_PLOIDY);
    }

    @Test
    public void testAddGQGreaterThanMaxGenotypeQual() {
        final VariantContext vc = getVariantContext();
        final GVCFBlock block90_100 = new HomRefBlock(vc, 90, 100, HomoSapiensConstants.DEFAULT_PLOIDY);

        // Test that adding a Genotype with GQ > 99 succeeds (ie., doesn't throw).
        // Internally, HomRefBlock should treat this GQ as 99.
        block90_100.add(vc.getStart()+1, new GenotypeBuilder(SAMPLE_NAME, vc.getAlleles()).GQ(150).DP(10).PL(new int[]{0, 10, 100}).make());
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testConstructorThrowsOnUpperGQBoundTooLarge() {
        final GVCFBlock block = new HomRefBlock(getVariantContext(), 90, 101, HomoSapiensConstants.DEFAULT_PLOIDY);
    }
}
