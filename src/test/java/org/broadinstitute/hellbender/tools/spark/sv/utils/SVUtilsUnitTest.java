package org.broadinstitute.hellbender.tools.spark.sv.utils;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecordComparator;
import htsjdk.samtools.SAMRecordCoordinateComparator;
import htsjdk.samtools.SAMRecordQueryNameComparator;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.UserException;
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
}
