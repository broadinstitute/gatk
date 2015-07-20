package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;

/**
 * Created by davidben on 8/3/15.
 */
public class SegmentUtilsUnitTest extends BaseTest {
    //a common set of CalledIntervals for tests
    private static final CalledInterval ci1 = new CalledInterval(new SimpleInterval("chr1", 1, 4), "call");
    private static final CalledInterval ci2 = new CalledInterval(new SimpleInterval("chr1", 5, 12), "call");
    private static final CalledInterval ci3 = new CalledInterval(new SimpleInterval("chr2", 1, 10), "call");

    @Test
    public void testReadUncalledSegments() {
        final File TEST_DIR = new File("src/test/resources/org/broadinstitute/tools/exome/caller");
        final File TEST_SEGMENTS = new File(TEST_DIR,"segments.tsv");


        List<SimpleInterval> segments = SegmentUtils.readIntervalsFromSegfile(TEST_SEGMENTS);
        Assert.assertEquals(segments.size(), 4);
        Assert.assertEquals(segments.get(0).getContig(), "chr");
        Assert.assertEquals(segments.get(1).getStart(), 300);
        Assert.assertEquals(segments.get(2).getEnd(), 550);
    }

    @Test
    public void testReadAndWriteCalledIntervals() {
        final File file = createTempFile("test",".txt");
        List<CalledInterval> calledIntervals = Arrays.asList(ci1, ci2, ci3);

        SegmentUtils.writeCalledIntervalsToSegfile(file, calledIntervals, "sample");
        List<CalledInterval> sameIntervals = SegmentUtils.readCalledIntervalsFromSegfile(file);
        Assert.assertEquals(calledIntervals, sameIntervals);
    }

    @Test
    public void testMeanTargetCoverage() {
        //a common set of targets for tests
        final TargetCoverage target1 = new TargetCoverage("myTarget", new SimpleInterval("chr1", 1, 2), 1.0);
        final TargetCoverage target2 = new TargetCoverage("myTarget", new SimpleInterval("chr1", 3, 4), 2.0);
        final TargetCoverage target3 = new TargetCoverage("myTarget", new SimpleInterval("chr1", 5, 6), 3.0);
        final TargetCoverage target4 = new TargetCoverage("myTarget", new SimpleInterval("chr2", 1, 2), 4.0);
        final TargetCoverage target5 = new TargetCoverage("myTarget", new SimpleInterval("chr2", 3, 4), 5.0);
        final TargetCoverage target6 = new TargetCoverage("myTarget", new SimpleInterval("chr2", 5, 6), 6.0);
        final HashedListTargetCollection<TargetCoverage> targets =
                new HashedListTargetCollection<>(Arrays.asList(target1, target2, target3, target4, target5, target6));

        Assert.assertEquals(SegmentUtils.meanTargetCoverage(ci1, targets), (1.0 + 2.0) / 2, 0.00001);
        Assert.assertEquals(SegmentUtils.meanTargetCoverage(ci2, targets), (3.0) / 1, 0.00001);
        Assert.assertEquals(SegmentUtils.meanTargetCoverage(ci3, targets), (4.0 + 5.0 + 6.0) / 3, 0.00001);
    }
}
