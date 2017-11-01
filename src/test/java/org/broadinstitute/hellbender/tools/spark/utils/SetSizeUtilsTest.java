package org.broadinstitute.hellbender.tools.spark.utils;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

public class SetSizeUtilsTest extends GATKBaseTest {

    private static boolean isPrime(final long iii) {
        if (iii % 2 == 0 || iii % 3 == 0) return false;
        long iFact = 5;
        while (iFact * iFact <= iii) {
            if (iii % iFact == 0 || iii % (iFact + 2) == 0) return false;
            iFact += 6;
        }
        return true;
    }

    @Test
    void legalCapacitiesTest() {
        final int[] caps = SetSizeUtils.legalSizes;
        final int nCaps = caps.length;
        // test that they're spaced properly -- each is supposed to be about sqrt(2) bigger than the previous one
        for (int idx = 1; idx < nCaps; ++idx) {
            final double err = Math.abs(1. - Math.sqrt(2.) * caps[idx - 1] / caps[idx]);
            Assert.assertTrue(err < .015, "testing capacity " + caps[idx] + " at index " + idx);
        }
        // test that they're all primes
        for (int idx = 0; idx < nCaps; ++idx) {
            Assert.assertTrue(isPrime(caps[idx]), "testing capacity " + caps[idx] + " at index " + idx);
        }
    }

    @Test(expectedExceptions = Exception.class)
    void getLegalSizeBelowTest() {
        Assert.assertEquals(SetSizeUtils.getLegalSizeBelow(1000),719);
        Assert.assertEquals(SetSizeUtils.getLegalSizeBelow(2147483630),2147483629);
        Assert.assertTrue(SetSizeUtils.getLegalSizeBelow(2039) < 2039);

        //Illegal max size
        SetSizeUtils.getLegalSizeBelow(1);
    }

    @Test(expectedExceptions = Exception.class)
    void getLegalSizeAboveTest() {
        Assert.assertEquals(SetSizeUtils.getLegalSizeAbove(1000),1021);
        Assert.assertEquals(SetSizeUtils.getLegalSizeAbove(2147483628),2147483629);
        Assert.assertTrue(SetSizeUtils.getLegalSizeAbove(2039) > 2039);

        //Illegal min size
        SetSizeUtils.getLegalSizeAbove(2147483630);
    }

}