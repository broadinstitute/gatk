package org.broadinstitute.hellbender.utils;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import static org.testng.Assert.*;

public class GA4GHUtilsTest extends GATKBaseTest {

    @DataProvider(name = "testSha")
    public Object[][] testSha() {
        return new Object[][]{
                //values from https://github.com/ga4gh/vrs-python/blob/main/src/ga4gh/core/digests.py
                {"", "z4PhNX7vuL3xVChQ1m2AB9Yg5AULVxXc"},
                {"ACGT", "aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2"}
        };
    }

    @Test(dataProvider= "testSha")
    public void testEncodeAsSha512t24u(String input, String expected) {
        assertEquals(GA4GHUtils.encodeAsSha512t24u(input), expected);
    }

    @Test
    public void testComposeGA4GHIdentifier() {
        final String prefix = "VA";
        final String digest = "z4PhNX7vuL3xVChQ1m2AB9Yg5AULVxXc";
        final String expected = GA4GHUtils.NAMESPACE + GA4GHUtils.NAMESPACE_SEP + prefix + GA4GHUtils.GA4GH_PREFIX_SEP + digest;
        assertEquals(GA4GHUtils.composeGA4GHIdentifier(prefix, digest), expected);
    }
}