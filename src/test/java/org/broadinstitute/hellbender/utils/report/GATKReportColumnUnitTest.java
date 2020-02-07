package org.broadinstitute.hellbender.utils.report;

import org.testng.Assert;
import org.testng.annotations.Test;

import static org.testng.Assert.*;

public class GATKReportColumnUnitTest {

    @Test
    public void testFormatValue() {
        //String.format("%.8f") produced strings that didn't parse for non-finite Doubles
        final GATKReportColumn testCol = new GATKReportColumn("values", "%.8f");
        Assert.assertTrue(Double.isNaN(Double.parseDouble(testCol.formatValue(Double.NaN))));
        Assert.assertTrue(Double.isInfinite(Double.parseDouble(testCol.formatValue(Double.POSITIVE_INFINITY))));
        Assert.assertTrue(Double.isInfinite(Double.parseDouble(testCol.formatValue(Double.NEGATIVE_INFINITY))));
    }
}