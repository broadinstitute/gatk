package org.broadinstitute.hellbender.utils;

import org.apache.commons.lang3.tuple.Pair;
import org.testng.Assert;
import org.testng.annotations.Test;

import static org.broadinstitute.hellbender.utils.MathUtils.promote;

public final class MannWhitneyUUnitTest {

    @Test
    public void testMWU() {
//        final int[] i0  = new int[]{0,6,7,8,9,10};
//        final int[] i00 = new int[]{1,2,3,4,5,11};
//
//        Assert.assertEquals(z(todouble(i0), todouble(i00), false), 1.072, 0.001);
//        Assert.assertEquals(z(todouble(i0), todouble(i00), true), 1.072, 0.001);
//        Assert.assertEquals(z(todouble(i00), todouble(i0), false), -1.072, 0.001);
//        Assert.assertEquals(z(todouble(i00), todouble(i0), true),  -1.072, 0.001);


        final int[] i1 = {2,4,5,6,8};
        final int[] i2 = {1,3,7,9,10,11,12,13};

        Assert.assertEquals(z(promote(i1), promote(i2), false), -1.3805, 0.001);

        Assert.assertEquals(z(promote(i1), promote(i2), true), -1.3805, 0.001);
        Assert.assertEquals(z(promote(i2), promote(i1), false),  1.3805, 0.001);
        Assert.assertEquals(z(promote(i2), promote(i1), true),  1.3805, 0.001);

        final int[] i3 = {0,2,4};
        final int[] i4 = {1,5,6,7,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34};

        Assert.assertEquals(z(promote(i3), promote(i4), false), -2.7240756662204, 0.001);
        Assert.assertEquals(z(promote(i3), promote(i4), true),  -2.7240756662204, 0.001);
    }

    private static double z(final double[] d1, final double[] d2, final boolean dither){
        MannWhitneyU mwu = new MannWhitneyU(dither);
        for ( double dp : d1 ) {
            mwu.add(dp,MannWhitneyU.USet.SET1);
        }
        for ( double dp : d2 ) {
            mwu.add(dp,MannWhitneyU.USet.SET2);
        }
        final Pair<Double, Double> pair = mwu.runOneSidedTest(MannWhitneyU.USet.SET1);
        return pair.getLeft();
    }

}
