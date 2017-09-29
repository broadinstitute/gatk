package org.broadinstitute.hellbender.tools.spark.pathseq;

import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

public class PSUtilsTest extends BaseTest {

    @Override
    public String getTestedClassName() {
        return PSUtils.class.getSimpleName();
    }

    @DataProvider(name = "maskData")
    public Object[][] getMaskData() {
        return new Object[][]{
                {"", 31, new byte[]{}},
                {"0,1", 31, new byte[]{0,1}},
                {"0,1,", 31, new byte[]{0,1}},
                {"0,1,10,8", 31, new byte[]{0,1,10,8}},
                {"30", 31, new byte[]{30}}
        };
    }

    @Test(dataProvider = "maskData")
    public void testParseMask(final String maskArg, final int kSize, final byte[] expected) {
        final byte[] result = PSUtils.parseMask(maskArg, kSize);
        Assert.assertNotNull(result);
        Assert.assertEquals(result, expected);
    }

    @DataProvider(name = "badMaskData")
    public Object[][] getBadMaskData() {
        return new Object[][]{
                {"0,-1", 31},
                {"0,31", 31}
        };
    }

    @Test(dataProvider = "badMaskData", expectedExceptions = IllegalArgumentException.class)
    public void testParseMaskException(final String maskArg, final int kSize) {
        PSUtils.parseMask(maskArg, kSize);
    }

    @Test(groups = "spark")
    public void testPrimaryReads() {
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        final List<GATKRead> readList = new ArrayList<>(2);

        final GATKRead primaryRead = ArtificialReadUtils.createRandomRead(101);
        readList.add(primaryRead);

        final GATKRead secondaryRead = ArtificialReadUtils.createRandomRead(101);
        secondaryRead.setIsSecondaryAlignment(true);
        readList.add(secondaryRead);

        final GATKRead supplementaryRead = ArtificialReadUtils.createRandomRead(101);
        supplementaryRead.setIsSupplementaryAlignment(true);
        readList.add(supplementaryRead);

        final JavaRDD<GATKRead> readRDD = ctx.parallelize(readList);
        final List<GATKRead> result = PSUtils.primaryReads(readRDD).collect();

        Assert.assertTrue(result.contains(primaryRead));
        Assert.assertTrue(result.size() == 1);
    }

    @Test
    public void testLogItemizedWarning() {
        final Collection<String> items = new ArrayList<>(3);
        PSUtils.logItemizedWarning(logger, items, "Test warning statement");
        items.add("x");
        items.add("y");
        items.add("z");
        PSUtils.logItemizedWarning(logger, items, "Test warning statement");
    }

}