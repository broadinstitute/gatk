package org.broadinstitute.hellbender.utils.test;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class VariantContextTestUtilsUnitTest {

    @DataProvider(name="scientificNotationValuesToNormalize")
    public Object[][] getScientificNotationValuesToNormalize(){
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

    @Test(dataProvider = "scientificNotationValuesToNormalize")
    public void testNormalizeScientificNotation(Object toNormalize, Object expected){
        Assert.assertEquals(VariantContextTestUtils.normalizeScientificNotation(toNormalize), expected);
    }

    @DataProvider(name="integerValuesToNormalize")
    public Object[][] getIntegerValuesToNormalize(){
        final Object aSpecificObject = new Object();
        return new Object[][] {
                {"27", new Integer(27)},
                {"-27", new Integer(-27)},
                {"0", new Integer(0)},
                {"-27014", new Integer(-27014)},
                {1, 1},
                {1, 1},
                {-1, -1},
                {1, 1},
                {"-2.172e+00", "-2.172e+00"},
                {"-2.172", "-2.172"},
                {-21.72, -21.72},
                {10, 10},
                {"SomeValue", "SomeValue"},
                {aSpecificObject,  aSpecificObject}
        };
    }

    @Test(dataProvider = "integerValuesToNormalize")
    public void testNormalizeInteger(Object toNormalize, Object expected){
        Assert.assertEquals(VariantContextTestUtils.normalizeToInteger(toNormalize), expected);
    }

}
