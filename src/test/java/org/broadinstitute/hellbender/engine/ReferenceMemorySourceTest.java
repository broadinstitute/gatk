package org.broadinstitute.hellbender.engine;

import com.google.api.client.util.Lists;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.testng.Assert;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * Tests for ReferenceMemorySource
 */
public class ReferenceMemorySourceTest {

    private final byte[] ref1Bytes = "CATCAT".getBytes();
    private final SimpleInterval interval1 = new SimpleInterval("chr1", 100, 105);
    private final ReferenceBases ref1 = new ReferenceBases(ref1Bytes, interval1);
    private ReferenceMemorySource memorySource;

    @DataProvider(name = "data")
    public Object[][] createReferenceMemorySourceTestData() {
        final SimpleInterval subInterval = new SimpleInterval("chr1", 101, 103);
        final byte[] siBytes = "ATC".getBytes();
        final SimpleInterval subInterval2 = new SimpleInterval("chr1", 101, 105);
        final byte[] si2Bytes = "ATCAT".getBytes();
        return new Object[][] {
            { interval1, ref1Bytes },
            { subInterval, siBytes },
            { subInterval2, si2Bytes }
        };
    }

    @DataProvider(name = "badIntervals")
    public Object[][] createBadIntervals() {
        return new Object[][] {
            // before
            { new SimpleInterval("chr1",99,103) },
            // after
            { new SimpleInterval("chr1",103,106) },
            // both
            { new SimpleInterval("chr1",99,106) }
        };
    }

    @BeforeTest
    public void init() {
        List<SAMSequenceRecord> l = new ArrayList<>();
        l.add(new SAMSequenceRecord("chr1",1000000));
        SAMSequenceDictionary seqDir = new SAMSequenceDictionary(l);
        memorySource = new ReferenceMemorySource(ref1, seqDir);
    }

    @Test(dataProvider="data")
    public void testQuery(SimpleInterval interval, byte[] bytes) throws Exception {
        checkEquals(memorySource.query(interval), bytes);
    }

    @Test(dataProvider="data")
    public void testQueryAndPrefetch(SimpleInterval interval, byte[] bytes) throws Exception {
        Assert.assertEquals(memorySource.queryAndPrefetch(interval).getBases(), bytes);
    }

    @Test(dataProvider="badIntervals", expectedExceptions = java.lang.IllegalArgumentException.class)
    public void testQueryOutOfBounds(SimpleInterval interval) {
        // we want to explode right away, not after going through the iterator for a while.
        memorySource.query(interval);
    }

    @Test(dataProvider="badIntervals", expectedExceptions = java.lang.IllegalArgumentException.class)
    public void testQueryAndPrefetchOutOfBounds(SimpleInterval interval) {
        memorySource.queryAndPrefetch(interval);
    }

    private void checkEquals(Iterator<Byte> actual, byte[] expected) {
        for (int i=0; i<expected.length; i++) {
            Assert.assertTrue(actual.hasNext());
            Assert.assertEquals(actual.next().byteValue(), expected[i]);
        }
        Assert.assertFalse(actual.hasNext());
    }
}