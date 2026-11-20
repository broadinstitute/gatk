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

    @DataProvider(name = "fnvByteArray64")
    public Object[][] getFnvByteArray64Data() {
        return new Object[][] {
                {new byte[0], new long[0]},
                {new byte[] {0}, new long[] {0}},
                {new byte[] {0x7f}, new long[] {0x7fL}},
                {new byte[] {(byte) 0x80}, new long[] {0x80L}},
                {new byte[] {(byte) 0xff}, new long[] {0xffL}},
                {new byte[] {(byte) 0xf0, (byte) 0xf1, (byte) 0xf2, (byte) 0xf3}, new long[] {0xf0f1f2f3L}},
                {new byte[] {(byte) 0x80, 1, 2, 3, 4, 5, 6, 7}, new long[] {0x8001020304050607L}},
                {new byte[] {1, (byte) 0xff, 3, 4, 5, 6, 7, 8}, new long[] {0x01ff030405060708L}},
                {new byte[] {1, 2, 3, 4, 5, 6, 7, (byte) 0x80}, new long[] {0x0102030405060780L}},
                {new byte[] {1, 2, 3, 4, 5, 6, 7, 8, (byte) 0x80}, new long[] {0x0102030405060708L, 0x80L}},
                {new byte[] {1, 2, 3, 4, 5, 6, 7, 8, (byte) 0xff, 1, 2, 3, 4, 5, 6, (byte) 0x80, (byte) 0xff, 1},
                        new long[] {0x0102030405060708L, 0xff01020304050680L, 0xff01L}},
                {new byte[] {'A', 'C', 'G', 'T'}, new long[] {0x41434754L}},
                {new byte[] {'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'N'}, new long[] {0x4143475441434754L, 0x4eL}}
        };
    }

    @Test(groups = "sv", dataProvider = "fnvByteArray64")
    public void testFnvByteArray64(final byte[] bytes, final long[] packedChunks) {
        final long seed = 0x0123456789abcdefL;
        long expected = seed;
        // Preserve the existing byte-array default seed and zero-padding of partial chunks.
        long expectedDefault = 1099511628211L;
        for (final long chunk : packedChunks) {
            expected = SVUtils.fnvLong64(expected, chunk);
            expectedDefault = SVUtils.fnvLong64(expectedDefault, chunk);
        }
        Assert.assertEquals(SVUtils.fnvByteArray64(seed, bytes), expected);
        Assert.assertEquals(SVUtils.fnvByteArray64(bytes), expectedDefault);
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
