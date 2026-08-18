package org.broadinstitute.hellbender.tools.copynumber;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.HashSet;
import java.util.Set;

/**
 * Unit tests for {@link GermlineCNVCaller}.
 */
public final class GermlineCNVCallerUnitTest extends GATKBaseTest {

    @Test
    public void testGenerateRetrySeedIsNonNegative() {
        for (int i = 0; i < 1000; i++) {
            final int seed = GermlineCNVCaller.generateRetrySeed();
            Assert.assertTrue(seed >= 0, "Retry seed must be non-negative (valid for numpy RandomState), but got: " + seed);
        }
    }

    @Test
    public void testGenerateRetrySeedIsLessThanMaxInt() {
        for (int i = 0; i < 1000; i++) {
            final int seed = GermlineCNVCaller.generateRetrySeed();
            Assert.assertTrue(seed < Integer.MAX_VALUE, "Retry seed must be less than Integer.MAX_VALUE, but got: " + seed);
        }
    }

    @Test
    public void testGenerateRetrySeedIsNonDeterministic() {
        final Set<Integer> seeds = new HashSet<>();
        for (int i = 0; i < 100; i++) {
            seeds.add(GermlineCNVCaller.generateRetrySeed());
        }
        Assert.assertTrue(seeds.size() > 1, "Retry seeds should not all be identical across calls");
    }
}
