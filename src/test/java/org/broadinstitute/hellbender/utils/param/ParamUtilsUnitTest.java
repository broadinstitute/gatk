package org.broadinstitute.hellbender.utils.param;

import org.testng.Assert;
import org.testng.annotations.Test;

/**
 * Created by lichtens on 8/25/15.
 */
public class ParamUtilsUnitTest {

    @Test
    public void testInRangeSuccess(){
        Assert.assertTrue(4 == ParamUtils.inRange(4, 3, 6, "Range calculation did not work properly"), "Did not return proper value");
        Assert.assertTrue(4.1 == ParamUtils.inRange(4.1, 3, 6, "Range calculation did not work properly"), "Did not return proper value");
        Assert.assertTrue(4.1 == ParamUtils.inRange(4.1, -3, 6, "Range calculation did not work properly"), "Did not return proper value");
        Assert.assertTrue(0 == ParamUtils.inRange(0, -3, 6, "Range calculation did not work properly"), "Did not return proper value");
        Assert.assertTrue(0.0 == ParamUtils.inRange(0.0, -3, 6, "Range calculation did not work properly"), "Did not return proper value");
        Assert.assertTrue(0 == ParamUtils.inRange(0, -3, 6, "Range calculation did not work properly"), "Did not return proper value");
        Assert.assertTrue(0.0 == ParamUtils.inRange(0.0, -3, 6, "Range calculation did not work properly"), "Did not return proper value");
        Assert.assertTrue(-1 == ParamUtils.inRange(-1, -3, 6, "Range calculation did not work properly"), "Did not return proper value");
        Assert.assertTrue(-1.5 == ParamUtils.inRange(-1.5, -3, 6, "Range calculation did not work properly"), "Did not return proper value");
        Assert.assertTrue(-1 == ParamUtils.inRange(-1, -3, 6, "Range calculation did not work properly"), "Did not return proper value");
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInRangeFailureAllPositiveInt(){
        ParamUtils.inRange(4, 7, 10, "Range calculation did not work properly");
    }
    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInRangeFailureAllPositiveDouble(){
        ParamUtils.inRange(4.1, 7.2, 10, "Range calculation did not work properly");
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInRangeFailureMinNegativeDouble(){
        ParamUtils.inRange(45.2, -7.2, 10, "Range calculation did not work properly");
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInRangeFailureValNegativeDouble(){
        ParamUtils.inRange(-10.3, -7.2, 10, "Range calculation did not work properly");
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInRangeFailureValNegativeInt(){
        ParamUtils.inRange(-10, -7, 10, "Range calculation did not work properly");
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInRangeFailureValNegativeInt2(){
        ParamUtils.inRange(-10, -7.2, 10, "Range calculation did not work properly");
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInRangeFailureValMinMaxUnreasonable(){
        ParamUtils.inRange(5, 10, -7, "Range calculation did not work properly.  Min was greater than max, so will always throw exception.");
    }

    @Test
    public void testInRangeLongDoubleDouble(){
        ParamUtils.inRange(5, -1.2, 7.1, "Range calculation did not work properly!");
    }

    @Test(expectedExceptions = IllegalArgumentException.class, expectedExceptionsMessageRegExp = "That is infinite!")
    public void testInfinite(){
        ParamUtils.isFinite(Double.POSITIVE_INFINITY, "That is infinite!");
    }

    @Test(expectedExceptions = IllegalArgumentException.class, expectedExceptionsMessageRegExp = "That is infinite!")
    public void testNegInfinite(){
        ParamUtils.isFinite(Double.NEGATIVE_INFINITY, "That is infinite!");
    }

    @Test(expectedExceptions = IllegalArgumentException.class, expectedExceptionsMessageRegExp = "That is negative!")
    public void testNegativeIsNotPositive(){
        ParamUtils.isPositiveOrZero(Double.NEGATIVE_INFINITY, "That is negative!");
    }

    @Test(expectedExceptions = IllegalArgumentException.class, expectedExceptionsMessageRegExp = "That is negative!")
    public void testNegativeIsNotPositive2(){
        ParamUtils.isPositiveOrZero(-4.2, "That is negative!");
    }

    @Test
    public void testLogB(){
        Assert.assertEquals(ParamUtils.logb(4, 2), 2.0);
        Assert.assertEquals(ParamUtils.logb(4.1, 2), 2.0356, 1e-3);
        Assert.assertEquals(ParamUtils.logb(4.1, 2.1), 1.9018, 1e-3);
        Assert.assertEquals(ParamUtils.logb(4, 2.1), 1.8685, 1e-3);

    }

    @Test
    public void testLog2(){
        Assert.assertEquals(2.0, ParamUtils.log2(4));
        Assert.assertEquals(10.0, ParamUtils.log2(1024));
        Assert.assertEquals(-13.2877, ParamUtils.log2(.0001), 1e-3);
    }

    @Test
    public void testLogBOfBaseZero(){
        Assert.assertEquals(-0.0, ParamUtils.logb(337, -0));
        Assert.assertEquals(-0.0, ParamUtils.logb(337, 0));
        Assert.assertTrue(Double.isNaN(ParamUtils.logb(-337, 0)));
        Assert.assertTrue(Double.isNaN(ParamUtils.logb(-337, -0)));
        Assert.assertTrue(Double.isNaN(ParamUtils.logb(0, 0)));
        Assert.assertTrue(Double.isNaN(ParamUtils.logb(-0, -0)));
        Assert.assertTrue(Double.isNaN(ParamUtils.logb(-0, 0)));
        Assert.assertTrue(Double.isNaN(ParamUtils.logb(0, -0)));
    }

    @Test
    public void testLogBOfZero(){
        Assert.assertEquals(ParamUtils.logb(0, 337), Double.NEGATIVE_INFINITY);
        Assert.assertEquals(ParamUtils.logb(-0, 337), Double.NEGATIVE_INFINITY);
    }

    @Test
    public void testLogBOfNegativeNumber() {
        // Note you cannot do Double.NaN == val  see: http://stackoverflow.com/questions/1456566/how-do-you-test-to-see-if-a-double-is-equal-to-nan
        Assert.assertTrue(Double.isNaN(ParamUtils.logb(-4, 5)));
        Assert.assertTrue(Double.isNaN(ParamUtils.logb(4, -5.1)));
    }

    @Test
    public void testLog2OfNegativeNumber() {
        Assert.assertTrue(Double.isNaN(ParamUtils.log2(-4)));
        Assert.assertTrue(Double.isNaN(ParamUtils.log2(-5.1)));
    }
}
