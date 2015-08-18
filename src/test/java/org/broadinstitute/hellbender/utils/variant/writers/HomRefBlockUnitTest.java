/*package org.broadinstitute.hellbender.engine.writers;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
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
    VariantContext vc;

    @BeforeMethod
    public void setUp() throws Exception {
        vc = new VariantContextBuilder("foo", "20", 1, 1, Arrays.asList(Allele.create("A", true), Allele.create("C"))).make();
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
        Assert.assertEquals(band.getStop(), pos - 1);
        assertValues(band, 10, 10, 11, 11);

        band.add(pos++, gb.DP(11).GQ(10).PL(new int[]{0,10,100}).make());
        Assert.assertEquals(band.getStop(), pos - 1);
        assertValues(band, 10, 11, 10, 11);

        band.add(pos++, gb.DP(12).GQ(12).PL(new int[]{0,12,100}).make());
        Assert.assertEquals(band.getStop(), pos - 1);
        assertValues(band, 10, 11, 10, 11);

        band.add(pos++, gb.DP(13).GQ(15).PL(new int[]{0,15,100}).make());
        Assert.assertEquals(band.getStop(), pos - 1);
        band.add(pos++, gb.DP(14).GQ(16).PL(new int[]{0,16,100}).make());
        Assert.assertEquals(band.getStop(), pos - 1);
        band.add(pos++, gb.DP(15).GQ(17).PL(new int[]{0,17,100}).make());
        Assert.assertEquals(band.getStop(), pos - 1);
        band.add(pos++, gb.DP(16).GQ(18).PL(new int[]{0,18,100}).make());
        Assert.assertEquals(band.getStop(), pos - 1);
        assertValues(band, 10, 13, 10, 15);
        Assert.assertEquals(band.getSize(), pos - vc.getStart());
        Assert.assertTrue(Arrays.equals(band.getMinPLs(), new int[]{0,10,100}));
    }

    @Test
    public void testBigGQIsCapped() {
        final HomRefBlock band = new HomRefBlock(vc, 10, 20, HomoSapiensConstants.DEFAULT_PLOIDY);
        final GenotypeBuilder gb = new GenotypeBuilder("NA12878");
        gb.alleles(vc.getAlleles());

        band.add(vc.getStart(), gb.DP(1000).GQ(1000).PL(new int[]{0,10,100}).make());
        assertValues(band, 1000, 1000, 99, 99);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testBadAdd() {
        final HomRefBlock band = new HomRefBlock(vc, 10, 20, HomoSapiensConstants.DEFAULT_PLOIDY);
        final GenotypeBuilder gb = new GenotypeBuilder("NA12878");

        band.add(vc.getStart() + 10, gb.DP(10).GQ(11).PL(new int[]{0,10,100}).make());
    }

    private void assertValues(final HomRefBlock band, final int minDP, final int medianDP, final int minGQ, final int medianGQ) {
        Assert.assertEquals(band.getMinDP(), minDP);
        Assert.assertEquals(band.getMedianDP(), medianDP);
        Assert.assertEquals(band.getMinGQ(), minGQ);
        Assert.assertEquals(band.getMedianGQ(), medianGQ);
    }


    @DataProvider(name = "ContiguousData")
    public Object[][] makeContiguousData() {
        List<Object[]> tests = new ArrayList<>();

        for ( final String chrMod : Arrays.asList("", ".mismatch") ) {
            for ( final int offset : Arrays.asList(-10, -1, 0, 1, 10) ) {
                final boolean equals = chrMod.equals("") && offset == 0;
                tests.add(new Object[]{vc.getChr() + chrMod, vc.getStart() + offset, equals});
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

    @Test
    public void testToVCFHeaderLine() {
        final HomRefBlock band = new HomRefBlock(vc, 10, 20, HomoSapiensConstants.DEFAULT_PLOIDY);
        Assert.assertEquals(band.toVCFHeaderLine().getKey(), "GVCFBlock10-20", "Wrong key for HomRefBlock " + band);
    }
}
*/
