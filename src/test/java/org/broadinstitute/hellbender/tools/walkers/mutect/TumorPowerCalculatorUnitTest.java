package org.broadinstitute.hellbender.tools.walkers.mutect;

import org.testng.Assert;
import org.testng.annotations.Test;

import static org.testng.Assert.*;

/**
 * Created by tsato on 6/19/16.
 */
public class TumorPowerCalculatorUnitTest {
    @Test
    public void testCachedPowerCalculation() throws Exception {
        TumorPowerCalculator tpc = new TumorPowerCalculator(0.001, 2.0, 0.0);
        final double epsilon = 0.01;
        Assert.assertEquals(tpc.cachedPowerCalculation(100, 0.2), 1.0, epsilon);
        Assert.assertEquals(tpc.cachedPowerCalculation(30,0.1), 0.8864, epsilon);
        Assert.assertEquals(tpc.cachedPowerCalculation(0,0.02), 0.0, epsilon);
        Assert.assertEquals(tpc.cachedPowerCalculation(5, 0.01), 0.0520, epsilon);
    }
}