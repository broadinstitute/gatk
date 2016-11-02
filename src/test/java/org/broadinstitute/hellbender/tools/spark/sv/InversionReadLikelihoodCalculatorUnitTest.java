package org.broadinstitute.hellbender.tools.spark.sv;

import org.apache.commons.lang3.StringUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;

/**
 * Created by shuang on 9/16/16.
 */
public class InversionReadLikelihoodCalculatorUnitTest {

    @Test
    public void testFormatDebugLikelihoodArray(){
        final List<SVDummyAllele> alleleList = Arrays.asList(new SVDummyAllele("A".getBytes(), true), new SVDummyAllele("T".getBytes(), false), new SVDummyAllele("C".getBytes(), false), new SVDummyAllele("G".getBytes(), false));
        final int rc = 4;
        final double[] rll = new double[alleleList.size()*rc];
        Arrays.fill(rll, -1.0);
        final String message = InversionReadLikelihoodCalculator.formatDebugLikelihoodArray(rll, alleleList, rc);
        Assert.assertEquals(StringUtils.countMatches(message, "A*"), 1);
        Assert.assertEquals(StringUtils.countMatches(message, "T"), 1);
        Assert.assertEquals(StringUtils.countMatches(message, "C"), 1);
        Assert.assertEquals(StringUtils.countMatches(message, "G"), 1);

        Assert.assertEquals(StringUtils.countMatches(message, "\n"), alleleList.size()*2);
        Assert.assertEquals(StringUtils.countMatches(message, ","), alleleList.size()*(rc-1));
    }
}
