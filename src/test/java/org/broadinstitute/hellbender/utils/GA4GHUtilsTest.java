package org.broadinstitute.hellbender.utils;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import static org.testng.Assert.*;

public class GA4GHUtilsTest extends GATKBaseTest {
    
    @DataProvider(name = "testSha")
    public Object[][] testSha() {
        return new Object[][]{
                {"", "z4PhNX7vuL3xVChQ1m2AB9Yg5AULVxXc"},
                {"ACGT", "aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2"}
        };
    }

    @Test(dataProvider= "testSha")
    public void testEncodeAsSha512t24u(String input, String expected) {
        String result = GA4GHUtils.encodeAsSha512t24u(input);
        assertEquals(result, expected);
    }
}