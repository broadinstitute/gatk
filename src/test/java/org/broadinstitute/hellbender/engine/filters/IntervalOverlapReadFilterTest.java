package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public final class IntervalOverlapReadFilterTest {

    // header for all the tests
    private static final SAMFileHeader HEADER = ArtificialReadUtils.createArtificialSamHeader();

    // creates the read filter with the test header without the interval
    private final IntervalOverlapReadFilter createReadFilter(final List<String> intervals) {
        final IntervalOverlapReadFilter filter = new IntervalOverlapReadFilter(intervals);
        filter.setHeader(HEADER);
        return filter;
    }

    private final GATKRead createRead(final String intervalString) {
        final SimpleInterval interval = new SimpleInterval(intervalString);
        final GATKRead read = ArtificialReadUtils.createArtificialRead(interval.getLengthOnReference() + "M");
        read.setPosition(interval);
        return read;
    }

    @DataProvider
    public Object[][] keepReadsData() {
        final IntervalOverlapReadFilter filter = createReadFilter(Arrays.asList("1:100-10000", "2:1-300", "3:1-10", "UNMAPPED"));
        return new Object[][] {
                {filter, ArtificialReadUtils.createArtificialUnmappedRead(HEADER, new byte[10], new byte[10])},
                {filter, createRead("1:1-150")},
                {filter, createRead("1:200-400")},
                {filter, createRead("2:100-300")},
                {filter, createRead("3:1-10")},
                {filter, createRead("3:10-100")}
        };
    }

    @Test(dataProvider = "keepReadsData")
    public void testKeepRead(final IntervalOverlapReadFilter filter, final GATKRead read) {
        Assert.assertTrue(filter.test(read));
    }

    @DataProvider
    public Object[][] filteredReadsData() {
        final IntervalOverlapReadFilter filter = createReadFilter(Arrays.asList("1:100-10000", "2:1-300", "3:1-10"));
        return new Object[][] {
                {filter, ArtificialReadUtils.createArtificialUnmappedRead(HEADER, new byte[10], new byte[10])},
                {filter, createRead("1:1-10")},
                {filter, createRead("1:10001-10010")},
                {filter, createRead("2:500-600")},
                {filter, createRead("10:1-100000")},
        };
    }

    @Test(dataProvider = "filteredReadsData")
    public void testFilterOutRead(final IntervalOverlapReadFilter filter, final GATKRead read) {
        Assert.assertFalse(filter.test(read));
    }

    @Test(dataProvider = "filteredReadsData", expectedExceptions = UserException.class)
    public void testFilterErrorsWhenProvidedBadIntervals(final IntervalOverlapReadFilter f, final GATKRead read) {
        final IntervalOverlapReadFilter filter = createReadFilter(Arrays.asList("1:-100-1000000000000", "2:1-300", "3:1-10"));
        Assert.assertFalse(filter.test(read));
    }
}
