package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class CalledSVGraphEventTest extends GATKBaseTest {

    private static final List<SAMSequenceRecord> TEST_RECORDS = Arrays.asList(new SAMSequenceRecord("0", 10000),
            new SAMSequenceRecord("1", 10000), new SAMSequenceRecord("2", 10000),
            new SAMSequenceRecord("3", 10000));
    private static final SAMSequenceDictionary TEST_DICT = new SAMSequenceDictionary(TEST_RECORDS);

    @DataProvider(name="eventData")
    private Object[][] eventDataGenerator() {
        final List<Object[]> data = new ArrayList<>(20);

        data.add(new Object[]{CalledSVGraphEvent.Type.DEL, 0, 1000, 2000, 0, 0, false});
        data.add(new Object[]{CalledSVGraphEvent.Type.DUP, 3, 0, 100, 3, 1, false});
        data.add(new Object[]{CalledSVGraphEvent.Type.UR, 0, 1000, 5000, 8, 1, true});

        return data.toArray(new Object[data.size()][]);
    }

    @Test(groups = "sv", dataProvider = "eventData")
    public void testConstructorAndGetters(final CalledSVGraphEvent.Type type, final int contig, final int start,
                                          final int end, final int groupId, final int pathId, final boolean resolved) {
        final SVInterval interval = new SVInterval(contig, start, end);
        /*final CalledSVGraphEvent event = new CalledSVGraphEvent(type, interval, groupId, pathId, resolved);
        Assert.assertEquals(type, event.getType());
        Assert.assertEquals(interval, event.getInterval());
        Assert.assertEquals(groupId, event.getGroupId());
        Assert.assertEquals(pathId, event.getPathId());
        Assert.assertEquals(resolved, event.isResolved());
        Assert.assertNotNull(event.bedString(TEST_DICT));*/
    }

    @Test(groups = "sv")
    public void testGetHeader() {
        final String header = CalledSVGraphEvent.bedHeader();
        Assert.assertNotNull(header);
    }

}