package org.broadinstitute.hellbender.utils;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public final class UnsignedTypeUtilTest {

    @DataProvider(name="uByteData")
    public Object[][] uByteData() {
        return new Object[][]{
            {(byte)0xFF, 255},
            {(byte)0xFE, 254},
            {(byte)0x80, 128},
            {(byte)0x7F, 127},
            {(byte)0x70, 112},
            {(byte)0x02, 2},
            {(byte)0x01, 1},
            {(byte)0x00, 0},
        };
    }

    @DataProvider(name="uShortData")
    public Object[][] uShortData() {
        return new Object[][]{
                {(short)0xFFFF, 65535},
                {(short)0xFFFE, 65534},
                {(short)0x8021, 32801},
                {(short)0x7FFF, 32767},
                {(short)0x5545, 21829},
                {(short)0x0002, 2},
                {(short)0x0001, 1},
                {(short)0x0000, 0}
        };
    }

    @DataProvider(name="uIntData")
    public Object[][] uIntData() {
        return new Object[][]{
                {0xFFFFFFFF,  4294967295L},
                {0xFFFFFFFE,  4294967294L},
                {0x81014000,  2164342784L},
                {0x7FFFFFFF,  2147483647L},
                {0x10502100,  273686784L},
                {0x00000002, 2L},
                {0x00000001, 1L},
                {0x00000000, 0L}
        };
    }

    /** Convert an unsigned byte to a signed int */
    @Test(dataProvider="uByteData")
    public void uByteToIntTest(final byte unsignedByte, final int expectedInt) {
        Assert.assertEquals(UnsignedTypeUtil.uByteToInt(unsignedByte), expectedInt);
    }

    /** Convert an unsigned byte to a signed short */
    @Test(dataProvider="uByteData")
    public void uByteToShortTest(final byte unsignedByte, final int expectedInt) {
        final short expectedShort = (short) expectedInt;
        Assert.assertEquals(UnsignedTypeUtil.uByteToShort(unsignedByte), expectedShort);
    }

    /** Convert an unsigned short to an Int */
    @Test(dataProvider="uShortData")
    public void uShortToIntTest(final short unsignedShort, final int expectedInt) {
        Assert.assertEquals(UnsignedTypeUtil.uShortToInt(unsignedShort), expectedInt);
    }

    /** Convert an unsigned int to a long */
    @Test(dataProvider="uIntData")
    public void uIntToLongTest(final int unsignedInt, final long expectedLong) {
        Assert.assertEquals(UnsignedTypeUtil.uIntToLong(unsignedInt), expectedLong);
    }
}
