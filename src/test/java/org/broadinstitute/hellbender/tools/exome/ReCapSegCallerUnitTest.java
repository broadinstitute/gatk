package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;

public final class ReCapSegCallerUnitTest extends BaseTest{
    @Test
    public void testMakeCalls() {
        List<TargetCoverage> targetList = new ArrayList<TargetCoverage>();
        //add amplification targets
        for (int i = 0; i < 10; i++) {
            targetList.add(new TargetCoverage("arbitrary_name", new SimpleInterval("chr", 100 + 2 * i, 101 + 2 * i), 1.0));
        }
        //add deletion targets
        for (int i = 0; i < 10; i++) {
            targetList.add(new TargetCoverage("arbitrary_name", new SimpleInterval("chr", 300 + 2 * i, 301 + 2 * i), -1.0));
        }
        //add targets that don't belong to a segment
        for (int i = 1; i < 10; i++) {
            targetList.add(new TargetCoverage("arbitrary_name", new SimpleInterval("chr", 400 + 2 * i, 401 + 2 * i), 0.0));
        }
        //add obviously neutral targets with some small spread
        for (int i = -5; i < 6; i++) {
            targetList.add(new TargetCoverage("arbitrary_name", new SimpleInterval("chr", 500 + 2 * i, 501 + 2 * i), 0.01 * i));
        }
        //add spread-out targets to a neutral segment (mean near zero)
        for (int i = -5; i < 6; i++) {
            targetList.add(new TargetCoverage("arbitrary_name", new SimpleInterval("chr", 700 + 2 * i, 701 + 2 * i), 0.1 * i));
        }

        HashedListExonCollection<TargetCoverage> targets = new HashedListExonCollection<TargetCoverage>(targetList);

        List<Segment> segments = new ArrayList<Segment>();
        segments.add(new Segment("sample", new SimpleInterval("chr", 100, 200), targets)); //amplification
        segments.add(new Segment("sample", new SimpleInterval("chr", 300, 400), targets)); //deletion
        segments.add(new Segment("sample", new SimpleInterval("chr", 450, 550), targets)); //neutral
        segments.add(new Segment("sample", new SimpleInterval("chr", 650, 750), targets)); //neutral

        ReCapSegCaller.makeCalls(targets, segments);

        Assert.assertEquals(segments.get(0).getCall(), ReCapSegCaller.AMPLIFICATION_CALL);
        Assert.assertEquals(segments.get(1).getCall(), ReCapSegCaller.DELETION_CALL);
        Assert.assertEquals(segments.get(2).getCall(), ReCapSegCaller.NEUTRAL_CALL);
        Assert.assertEquals(segments.get(3).getCall(), ReCapSegCaller.NEUTRAL_CALL);
    }
}
