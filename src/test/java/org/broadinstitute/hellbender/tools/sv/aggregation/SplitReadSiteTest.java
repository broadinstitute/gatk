package org.broadinstitute.hellbender.tools.sv.aggregation;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Collections;

import static org.testng.Assert.*;

public class SplitReadSiteTest {

    @Test
    public void testGetters() {
        final SplitReadSite site = new SplitReadSite("chr1", 100, true,
                Collections.singletonMap("TEST_SAMPLE", 5),
                new EvidenceStatUtils.PoissonTestResult(0.1, 0.2, 0.3));
        Assert.assertEquals(site.getContig(), "chr1");
        Assert.assertEquals(site.getPosition(), 100);
        Assert.assertEquals(site.getStrand(), true);
        Assert.assertEquals(site.getCount("TEST_SAMPLE"), 5);
        Assert.assertEquals(site.getCount("NOT_TEST_SAMPLE"), 0);
        Assert.assertEquals(site.getP(), 0.1);
        Assert.assertEquals(site.getCarrierSignal(), 0.2);
        Assert.assertEquals(site.getBackgroundSignal(), 0.3);
    }

    @DataProvider(name = "testEqualsData")
    public Object[][] testEqualsData() {
        return new Object[][]{
                {
                        new SplitReadSite("chr1", 100, true, Collections.singletonMap("TEST_SAMPLE", 5), new EvidenceStatUtils.PoissonTestResult(0.1, 0.2, 0.3)),
                        new SplitReadSite("chr1", 100, true, Collections.singletonMap("TEST_SAMPLE", 5), new EvidenceStatUtils.PoissonTestResult(0.1, 0.2, 0.3)),
                        true
                },
                {
                        new SplitReadSite("chr1", 100, true, Collections.singletonMap("TEST_SAMPLE", 5), new EvidenceStatUtils.PoissonTestResult(0.1, 0.2, 0.3)),
                        new SplitReadSite("chr1", 100, true, Collections.singletonMap("TEST_SAMPLE", 5), new EvidenceStatUtils.PoissonTestResult(0.2, 0.2, 0.3)),
                        false
                },
                {
                        new SplitReadSite("chr1", 100, true, Collections.singletonMap("TEST_SAMPLE", 5), new EvidenceStatUtils.PoissonTestResult(0.1, 0.2, 0.3)),
                        new SplitReadSite("chr2", 100, true, Collections.singletonMap("TEST_SAMPLE", 5), new EvidenceStatUtils.PoissonTestResult(0.1, 0.2, 0.3)),
                        false
                },
                {
                        new SplitReadSite("chr1", 100, true, Collections.singletonMap("TEST_SAMPLE", 5), new EvidenceStatUtils.PoissonTestResult(0.1, 0.2, 0.3)),
                        new SplitReadSite("chr1", 200, true, Collections.singletonMap("TEST_SAMPLE", 5), new EvidenceStatUtils.PoissonTestResult(0.1, 0.2, 0.3)),
                        false
                },
                {
                        new SplitReadSite("chr1", 100, true, Collections.singletonMap("TEST_SAMPLE", 5), new EvidenceStatUtils.PoissonTestResult(0.1, 0.2, 0.3)),
                        new SplitReadSite("chr1", 100, false, Collections.singletonMap("TEST_SAMPLE", 5), new EvidenceStatUtils.PoissonTestResult(0.1, 0.2, 0.3)),
                        false
                },
                {
                        new SplitReadSite("chr1", 100, true, Collections.singletonMap("TEST_SAMPLE", 5), new EvidenceStatUtils.PoissonTestResult(0.1, 0.2, 0.3)),
                        new SplitReadSite("chr1", 100, true, Collections.singletonMap("TEST_SAMPLE", 3), new EvidenceStatUtils.PoissonTestResult(0.1, 0.2, 0.3)),
                        false
                },
                {
                        new SplitReadSite("chr1", 100, true, Collections.singletonMap("TEST_SAMPLE", 5), new EvidenceStatUtils.PoissonTestResult(0.1, 0.2, 0.3)),
                        new SplitReadSite("chr1", 100, true, Collections.singletonMap("TEST_SAMPLE_2", 5), new EvidenceStatUtils.PoissonTestResult(0.1, 0.2, 0.3)),
                        false
                },

        };
    }

    @Test(dataProvider= "testEqualsData")
    public void testEquals(SplitReadSite a, SplitReadSite b, boolean expected) {
        Assert.assertEquals(a.equals(b), expected);
        if (expected) {
            Assert.assertEquals(a.hashCode(), b.hashCode());
        } else {
            Assert.assertNotEquals(a.hashCode(), b.hashCode());
        }
    }
}