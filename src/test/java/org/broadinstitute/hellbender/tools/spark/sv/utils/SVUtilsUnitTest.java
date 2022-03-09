package org.broadinstitute.hellbender.tools.spark.sv.utils;

import htsjdk.samtools.*;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SVInterval;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.broadinstitute.hellbender.tools.walkers.sv.SVAnnotateUnitTest.createSequenceDictionary;

/**
 * Unit tests for SVUtils.
 */
public class SVUtilsUnitTest extends GATKBaseTest {

    @Test(groups = "sv")
    void hashMapCapacityTest() {
        Assert.assertEquals(SVUtils.hashMapCapacity(150),201);
    }

    @DataProvider
    private Object[][] forGetSamRecordComparatorExpectingException() {
        final List<Object[]> data = new ArrayList<>(3);
        data.add(new Object[]{SAMFileHeader.SortOrder.unsorted});
        data.add(new Object[]{SAMFileHeader.SortOrder.duplicate});
        data.add(new Object[]{SAMFileHeader.SortOrder.unknown});
        return data.toArray(new Object[data.size()][]);
    }
    @Test(groups = "sv", expectedExceptions = UserException.class, dataProvider = "forGetSamRecordComparatorExpectingException")
    void testGetSamRecordComparatorExpectingException(final SAMFileHeader.SortOrder sortOrder) {
        SVUtils.getSamRecordComparator(sortOrder);
    }

    @DataProvider
    private Object[][] forGetSamRecordComparator() {
        final List<Object[]> data = new ArrayList<>(2);
        data.add(new Object[]{SAMFileHeader.SortOrder.coordinate, SAMRecordCoordinateComparator.class});
        data.add(new Object[]{SAMFileHeader.SortOrder.queryname, SAMRecordQueryNameComparator.class});
        return data.toArray(new Object[data.size()][]);
    }
    @Test(groups = "sv", dataProvider = "forGetSamRecordComparator")
    void testGetSamRecordComparator(final SAMFileHeader.SortOrder sortOrder, final Class<SAMRecordComparator> expectedClass) {
        Assert.assertEquals(SVUtils.getSamRecordComparator(sortOrder).getClass(), expectedClass);
    }

    @DataProvider(name="locatableToSVInterval")
    public Object[][] getLocatableToSVIntervalData() {
        return new Object[][] {
                { new SVInterval(0, 1, 101), new SimpleInterval("chr1", 1, 100)},
                { new SVInterval(1, 201, 202), new SimpleInterval("chr2", 201, 201)},
                { null, new SimpleInterval("chr3", 2, 3)}
        };
    }


    @Test(dataProvider = "locatableToSVInterval")
    public void testLocatableToSVInterval(
            final SVInterval expectedSVInterval,
            final SimpleInterval locatable)
    {
        final SAMSequenceDictionary sequenceDictionary = createSequenceDictionary(Arrays.asList("chr1", "chr2"));
        try {
            Assert.assertEquals(SVUtils.locatableToSVInterval(locatable, sequenceDictionary), expectedSVInterval);
        } catch (IllegalArgumentException e) {
            Assert.assertNull(expectedSVInterval);
        }
    }
}
