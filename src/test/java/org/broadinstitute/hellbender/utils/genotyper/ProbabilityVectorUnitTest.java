package org.broadinstitute.hellbender.utils.genotyper;

import org.testng.Assert;
import org.testng.annotations.Test;

public class ProbabilityVectorUnitTest {
    @Test
    public void basicPVCreationTest() {

        // basic sanity checks for compression
        ProbabilityVector pv = new ProbabilityVector(new double[]{0.0, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY});
        Assert.assertEquals(pv.getMaxVal(), 0);

        pv = new ProbabilityVector(new double[]{Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY,0.0, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY});
        Assert.assertEquals(pv.getMinVal(), 2);
        Assert.assertEquals(pv.getMaxVal(), 2);

        pv = new ProbabilityVector(new double[]{Double.NEGATIVE_INFINITY,0.0,0.0, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY});
        Assert.assertEquals(pv.getMinVal(), 1);
        Assert.assertEquals(pv.getMaxVal(), 2);

    
        pv = new ProbabilityVector(new double[]{-9.0,-10.0,-11.0,-20.0,0.0,-100.0,-3.0});
        Assert.assertEquals(pv.getMinVal(), 0);
        Assert.assertEquals(pv.getMaxVal(), 6);
        Assert.assertEquals(pv.getProbabilityVector().length, 7);

        pv = new ProbabilityVector(new double[]{-30.1,-20.0,-100.0});
        Assert.assertEquals(pv.getMinVal(), 1);
        Assert.assertEquals(pv.getMaxVal(), 1);
        Assert.assertEquals(pv.getProbabilityVector().length, 1);

        pv = new ProbabilityVector(new double[]{-3.0,-4.0,-5.0});
        Assert.assertEquals(pv.getMinVal(), 0);
        Assert.assertEquals(pv.getMaxVal(), 2);
        Assert.assertEquals(pv.getProbabilityVector().length, 3);

    }
    
    @Test
    public void testDotProduct() {
        ProbabilityVector v1 = new ProbabilityVector(new double[]{-3.0,-4.0,-5.0});
        ProbabilityVector v2 = new ProbabilityVector(new double[]{Double.NEGATIVE_INFINITY,0.0, Double.NEGATIVE_INFINITY,});
        Assert.assertEquals(v1.logDotProduct(v2), -4.0, 1e-3);

        v1 = new ProbabilityVector(new double[]{-3.0,-3.0,});
        v2 = new ProbabilityVector(new double[]{-3.0,-3.0,});
        // log10(10^-6+10^-6) = -5.69897
        Assert.assertEquals(v1.logDotProduct(v2), -5.69897, 1e-3);

        // test non-overlapping range case
        v1 = new ProbabilityVector(new double[]{-3.0,-3.0, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY});
        v2 = new ProbabilityVector(new double[]{Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY,-3.0,-3.0,});
        Assert.assertEquals(v1.logDotProduct(v2), Double.NEGATIVE_INFINITY, 1e-3);

    }
    
}
