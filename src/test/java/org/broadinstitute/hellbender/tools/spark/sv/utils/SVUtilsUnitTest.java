package org.broadinstitute.hellbender.tools.spark.sv.utils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.spark.sv.TestUtilsForSV;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;

/**
 * Unit tests for SVUtils.
 */
public class SVUtilsUnitTest extends GATKBaseTest {

    @Test(groups = "sv")
    void hashMapCapacityTest() {
        Assert.assertEquals(SVUtils.hashMapCapacity(150),201);
    }

    @Test(groups = "sv")
    public void testConversionBetweenLocatableAndSvInterval() {
        SimpleInterval simpleInterval = new SimpleInterval("20", 10000, 10000);
        SVInterval svInterval = new SVInterval(0, 10000, 10001);
        Assert.assertEquals(SVUtils.convertLocatable(simpleInterval, TestUtilsForSV.b37_seqDict), svInterval);
        Assert.assertEquals(SVUtils.convertSVInterval(svInterval, TestUtilsForSV.b37_seqDict), simpleInterval);
    }

    @DataProvider
    private Object[][] forTestConvertFromLocatableExpectingException() {
        final List<Object[]> data = new ArrayList<>(3);
        data.add(new Object[]{new SimpleInterval("2", 10000, 10000), TestUtilsForSV.b37_seqDict});
        data.add(new Object[]{new SimpleInterval("chr20", 10000, 10000), TestUtilsForSV.b37_seqDict});
        data.add(new Object[]{new SimpleInterval("20", 10000, 10000), TestUtilsForSV.b38_seqDict_chr20_chr21});
        return data.toArray(new Object[data.size()][]);
    }
    @Test(groups = "sv", dataProvider = "forTestConvertFromLocatableExpectingException", expectedExceptions = IllegalArgumentException.class)
    public void testConvertFromLocatableExpectingException(final Locatable locatable, final SAMSequenceDictionary dictionary) {
        SVUtils.convertLocatable(locatable, dictionary);
    }

    @DataProvider
    private Object[][] forTestConvertFromSVIntervalExpectingException() {
        final List<Object[]> data = new ArrayList<>(2);
        data.add(new Object[]{new SVInterval(3, 10000, 10000), TestUtilsForSV.b37_seqDict}); // invalid chromosome
        data.add(new Object[]{new SVInterval(0, 10000, 10000), TestUtilsForSV.b37_seqDict}); // length zero
        return data.toArray(new Object[data.size()][]);
    }
    @Test(groups = "sv", dataProvider = "forTestConvertFromSVIntervalExpectingException", expectedExceptions = IllegalArgumentException.class)
    public void testConvertFromSVIntervalExpectingException(final SVInterval svInterval, final SAMSequenceDictionary dictionary) {
        SVUtils.convertSVInterval(svInterval, dictionary);
    }
}
