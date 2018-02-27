package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.utils.IntHistogramTest;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.Test;
import scala.Tuple2;

import java.util.*;

public class IntervalCoverageFinderUnitTest extends GATKBaseTest {
    private final static LibraryStatistics LIBRARY_STATISTICS =
            new LibraryStatistics(IntHistogramTest.genLogNormalSample(400, 175, 10000).getCDF(),
                    60000000000L, 600000000L, 1200000000000L, 3000000000L);


    @Test
    public void testIntervalCoverageFinder() throws Exception {
        final StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection params = new StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection();
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeaderWithGroups(3, 1, 10000000, 1);
        final Set<Integer> crossContigIgnoreSet = new HashSet<>();
        final ReadMetadata readMetadata = new ReadMetadata(crossContigIgnoreSet, header, LIBRARY_STATISTICS, null, 2L, 2L, 1);
        final ArrayList<SVInterval> intervals = new ArrayList<>(2);
        final SVInterval interval1 = new SVInterval(0, 1000, 1100);
        intervals.add(interval1);
        final GATKRead read1 = ArtificialReadUtils.createArtificialRead(header, "read1", 0, 1050, 100);
        final GATKRead read2 = ArtificialReadUtils.createArtificialRead(header, "read1", 0, 1090, 100);
        final ArrayList<GATKRead> reads = new ArrayList<>();
        reads.add(read1);
        reads.add(read2);

        IntervalCoverageFinder intervalCoverageFinder = new IntervalCoverageFinder(readMetadata, intervals, reads.iterator(), new SVReadFilter(params));

        Iterator<Tuple2<Integer, int[]>> iterator = intervalCoverageFinder.iterator();
        Assert.assertTrue(iterator.hasNext());
        Tuple2<Integer, int[]> next = iterator.next();
        Assert.assertEquals(next._1.intValue(), 0);
        int[] expected = new int[100];
        for (int i = 50; i < 90; i++) {
            expected[i] = 1;
        }
        for (int i = 90; i < 100; i++) {
            expected[i] = 2;
        }

        Assert.assertTrue(Arrays.equals(next._2, expected));

        final SVInterval interval2 = new SVInterval(0, 1100, 1200);
        intervals.add(interval2);

        Assert.assertTrue(! interval1.overlaps(interval2));

        intervalCoverageFinder = new IntervalCoverageFinder(readMetadata, intervals, reads.iterator(), new SVReadFilter(params));

        iterator = intervalCoverageFinder.iterator();
        Assert.assertTrue(iterator.hasNext());
        next = iterator.next();
        Assert.assertEquals(next._1.intValue(), 0);
        expected = new int[100];
        for (int i = 50; i < 90; i++) {
            expected[i] = 1;
        }
        for (int i = 90; i < 100; i++) {
            expected[i] = 2;
        }

        Assert.assertTrue(iterator.hasNext());
        next = iterator.next();
        Assert.assertEquals(next._1.intValue(), 1);
        expected = new int[100];
        for (int i = 0; i < 50; i++) {
            expected[i] = 2;
        }
        for (int i = 50; i < 90; i++) {
            expected[i] = 1;
        }

        Assert.assertTrue(Arrays.equals(next._2, expected));
    }


    @Test(groups="sv")
    public void testGetHighCoverageIntervals() throws Exception {
        final SVInterval interval = new SVInterval(0, 100, 200);
        final int[] coverageArray = new int[interval.getLength()];
        for (int i = 10; i < 20; i++) {
            coverageArray[i] = 101;
        }
        coverageArray[15] = 401;

        for (int i = 60; i < 80; i++) {
            coverageArray[i] = 101;
        }

        for (int i = 90; i < 100; i++) {
            coverageArray[i] = 101;
        }

        final Iterator<IntervalCoverageFinder.CandidateCoverageInterval> highCoverageIntervals =
                IntervalCoverageFinder.getHighCoverageIntervalsInWindow(100, 400, interval, coverageArray);
        Assert.assertTrue(highCoverageIntervals.hasNext());
        IntervalCoverageFinder.CandidateCoverageInterval next = highCoverageIntervals.next();
        Assert.assertEquals(next.getInterval(), new SVInterval(0, 110, 120));

        Assert.assertTrue(highCoverageIntervals.hasNext());
        next = highCoverageIntervals.next();
        Assert.assertEquals(next.getInterval(), new SVInterval(0, 190, 200));
        Assert.assertFalse(next.containsMaxCoveragePeak());

        Assert.assertFalse(highCoverageIntervals.hasNext());
    }

    @Test(groups="sv")
    public void testGetHighCoverageIntervalsEntireInterval() throws Exception {
        final SVInterval interval = new SVInterval(0, 100, 200);
        final int[] coverageArray = new int[interval.getLength()];
        for (int i = 0; i < coverageArray.length; i++) {
            coverageArray[i] = 401;
        }
        final Iterator<IntervalCoverageFinder.CandidateCoverageInterval> highCoverageIntervals =
                IntervalCoverageFinder.getHighCoverageIntervalsInWindow(100, 400, interval, coverageArray);
        Assert.assertTrue(highCoverageIntervals.hasNext());
        final SVInterval next = highCoverageIntervals.next().getInterval();
        Assert.assertEquals(next, interval);
    }

}