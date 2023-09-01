package org.broadinstitute.hellbender.utils.read;

import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;

public class FlowBasedKeyCodecUnitTest extends GATKBaseTest {

    @DataProvider(name = "testData")
    public Object[][] getTestData() {

        final Object[][]        testData = {

                // trivial cases
                { "T", "TGCA", "1", false },
                { "TT", "TGCA", "2", false },
                { "TGCA", "TGCA", "1,1,1,1", false },
                { "TA", "TGCA", "1,0,0,1", false },
                { "TTAATG", "TGCA", "2,0,0,2,1,1", false },

                // clipping
                { "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
                          , "TGCA", "130", false },

                // N processing
                { "TNTA", "TGCA", "3,0,0,1", false},
                { "TTNA", "TGCA", "3,0,0,1", false},
                { "TTAN", "TGCA", "2,0,0,2", false},
                { "TTAN", "TGCA", "2,0,0,2", false},
                { "NTTA", "TGCA", "3,0,0,1", false},
                { "NGGA", "TGCA", "1,2,0,1", false},
                { "NTGGA", "TGCA", "2,2,0,1", false},

                { "NT*GGA", "TGCA", "2,2,0,1", true}
        };

        return testData;
    }

    @Test(dataProvider = "testData")
    public void testBase2Key(final String basesAsString, final String flowOrder,
                                    final String expectedKeyAsString, final boolean expectException) {

        // int version
        try {
            final int[] intKey = FlowBasedKeyCodec.baseArrayToKey(basesAsString.getBytes(), flowOrder);
            Assert.assertNotNull(intKey);
            final String        intKeyAsString = StringUtils.join(intKey, ',');
            Assert.assertEquals(intKeyAsString, expectedKeyAsString);
        } catch (Exception e) {
            if ( !expectException )
                throw e;
        }
    }

    @DataProvider(name="makeReadArrayTests")
    public Object[][] makeByteArrayTests(){
        List<Object[]> tests = new ArrayList<>();
        tests.add(new Object[]{new byte[]{'T','T','T','A','T','G','C'}, new byte[]{10,10,10,10,10,10,10}, "ACTG", (byte)0, new byte[]{0,0,10,10,10,10,10,10,10,10}});
        tests.add(new Object[]{new byte[]{'T','T','T','A','T','G','C'}, new byte[]{10,10,10,10,10,10,10}, "ACTG", (byte)10, new byte[]{10,10,10,10,10,10,10,10,10,10}});
        tests.add(new Object[]{new byte[]{'T','T','T','A','T','G','C'}, new byte[]{10,5,10,10,10,10,10}, "ACTG", (byte)0, new byte[]{0,0,5,5,10,10,10,10,10,10}});
        tests.add(new Object[]{new byte[]{'T','T','T','A','T','G','C'}, new byte[]{10,25,10,10,10,10,10}, "ACTG", (byte)0, new byte[]{0,0,10,10,10,10,10,10,10,10}});
        tests.add(new Object[]{new byte[]{'T','T','T','A','T','G','C'}, new byte[]{1,2,3,4,5,6,7}, "ACTG", (byte)0, new byte[]{0,0,1,1,4,4,5,6,6,7}});

        return tests.toArray(new Object[][]{});
    }

    @Test (dataProvider = "makeReadArrayTests")
    public void testBaseArray2KeySpace(final byte[] readBases, final byte[] qualArray, final String flowOrder, final byte defualtQual, final byte[] expected) {
        final int[] flowBases = FlowBasedKeyCodec.baseArrayToKey(readBases, flowOrder);

        final byte[] result = FlowBasedKeyCodec.baseArrayToKeySpace(readBases, flowBases.length, qualArray, defualtQual, flowOrder);

        Assert.assertEquals(flowBases.length, result.length, "Read bases in flow space and baseArrayToKeySpace do not match in length");
        Assert.assertEquals(result, expected);
    }

}
