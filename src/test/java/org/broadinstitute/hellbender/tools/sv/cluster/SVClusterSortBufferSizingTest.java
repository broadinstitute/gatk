package org.broadinstitute.hellbender.tools.sv.cluster;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

/**
 * Unit tests for {@link SVClusterWalker#genotypeHeavyMaxRecordsInRam(int, int)}, the pure arithmetic that
 * sample-scales the in-RAM window of genotype-heavy sort buffers in low-mem mode.
 */
public final class SVClusterSortBufferSizingTest {

    @DataProvider(name = "sizingData")
    public Object[][] sizingData() {
        final long budget = SVClusterWalker.LOW_MEM_GENOTYPE_BUFFER_BUDGET; // 40_000_000
        final int floor = SVClusterWalker.LOW_MEM_SORT_BUFFER_FLOOR;        // 100
        return new Object[][]{
                // {maxRecordsInRam, numSamples, expected}
                // No samples (e.g. sites-only): unchanged.
                {10000, 0, 10000},
                // Small cohort: budget/n far exceeds maxRecordsInRam -> clamps to maxRecordsInRam (unchanged).
                {10000, 100, 10000},
                {10000, (int) (budget / 10000) - 1, 10000}, // just below the scale-down threshold
                // Scale-down branch: budget/n is between floor and maxRecordsInRam.
                {10000, 5_000, (int) (budget / 5_000)},     // 12M/5k = 2400
                {10000, 50_000, (int) (budget / 50_000)},   // 12M/50k = 240
                // Floor clamp: budget/n drops below the floor.
                {10000, 250_000, floor},                    // 12M/250k = 48 -> floored to 50
                {10000, 1_000_000, floor},                  // 12M/1M = 12 -> floored to 50
                {10000, (int) (budget / floor) + 1, floor}, // just past where budget/n < floor
                // maxRecordsInRam itself below the scaled value caps the result (the cap=2 forced-spill case).
                {2, 250_000, 2},
        };
    }

    @Test(dataProvider = "sizingData")
    public void testGenotypeHeavyMaxRecordsInRam(final int maxRecordsInRam, final int numSamples, final int expected) {
        Assert.assertEquals(
                SVClusterWalker.genotypeHeavyMaxRecordsInRam(maxRecordsInRam, numSamples),
                expected,
                "maxRecordsInRam=" + maxRecordsInRam + ", numSamples=" + numSamples);
    }

    @Test
    public void testNeverExceedsMaxOrDropsBelowFloorForPositiveSamples() {
        final int max = 10000;
        for (int n = 1; n <= 2_000_000; n += 7919) { // sweep a prime stride across the range
            final int result = SVClusterWalker.genotypeHeavyMaxRecordsInRam(max, n);
            Assert.assertTrue(result <= max, "exceeded maxRecordsInRam at n=" + n + ": " + result);
            Assert.assertTrue(result >= Math.min(max, SVClusterWalker.LOW_MEM_SORT_BUFFER_FLOOR),
                    "dropped below floor at n=" + n + ": " + result);
        }
    }
}
