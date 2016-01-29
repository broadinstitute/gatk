package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.Test;

/**
 * Created by davidben on 1/25/16.
 */
public class AllelicCountUnitTest {

    @Test
    public void testEstimateMinorAlleleFraction() {
        final SimpleInterval interval = new SimpleInterval("contig", 1, 2);
        final double tolerance = 0.0001;

        //exact MLE values obtained from Sage Math
        Assert.assertEquals(new AllelicCount(interval, 10, 10).estimateMinorAlleleFraction(), 0.5, tolerance);
        Assert.assertEquals(new AllelicCount(interval, 12, 19).estimateMinorAlleleFraction(), 0.39936538237877861, tolerance);
        Assert.assertEquals(new AllelicCount(interval, 40, 12).estimateMinorAlleleFraction(), 0.23076923128457277, tolerance);
        Assert.assertEquals(new AllelicCount(interval, 53, 67).estimateMinorAlleleFraction(), 0.44744832677710078, tolerance);
        Assert.assertEquals(new AllelicCount(interval, 55, 45).estimateMinorAlleleFraction(), 0.5, tolerance);
    }

    @Test
    public void testEstimateAltAlleleFraction() {
        final SimpleInterval interval = new SimpleInterval("contig", 1, 2);
        final double tolerance = 0.0001;

        Assert.assertEquals(new AllelicCount(interval, 10, 10).estimateAltAlleleFraction(), 0.5, tolerance);
        Assert.assertEquals(new AllelicCount(interval, 20, 10).estimateAltAlleleFraction(), 0.3333333, tolerance);
        Assert.assertEquals(new AllelicCount(interval, 10, 20).estimateAltAlleleFraction(), 0.6666666, tolerance);
        Assert.assertEquals(new AllelicCount(interval, 10, 90).estimateAltAlleleFraction(), 0.9, tolerance);
    }

    @Test
    public void testToMinorAlleleFractionTargetCoverage() {
        final SimpleInterval interval = new SimpleInterval("contig", 1, 2);
        final String sample = "sample";
        final TargetCoverage result = new AllelicCount(interval, 10, 10).toMinorAlleleFractionTargetCoverage(sample);

        final TargetCoverage expected = new TargetCoverage(sample, interval, 0.5);
        Assert.assertEquals(result.getInterval(), expected.getInterval());
        Assert.assertEquals(result.getName(), expected.getName());
        Assert.assertEquals(result.getCoverage(), expected.getCoverage());
    }
}