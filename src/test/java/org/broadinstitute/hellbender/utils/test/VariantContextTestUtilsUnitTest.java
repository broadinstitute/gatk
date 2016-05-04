package org.broadinstitute.hellbender.utils.test;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class VariantContextTestUtilsUnitTest {

    @DataProvider(name="valuesToNormalize")
    public Object[][] getValuesToNormalize(){
        final Object aSpecificObject = new Object();
        return new Object[][] {
                {"-2.172e+00", -2.172},
                {"-2.172", -2.172},
                {"-2.172e+01", -21.72},
                {-21.72, -21.72},
                {10, 10},
                {"SomeValue", "SomeValue"},
                {aSpecificObject,  aSpecificObject}

        };
    }

    @Test(dataProvider = "valuesToNormalize")
    public void testNormalizeScientificNotation(Object toNormalize, Object expected){
        Assert.assertEquals(VariantContextTestUtils.normalizeScientificNotation(toNormalize), expected);
    }
}
