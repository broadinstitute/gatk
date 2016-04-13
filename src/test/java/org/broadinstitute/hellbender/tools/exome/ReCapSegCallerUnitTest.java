package org.broadinstitute.hellbender.tools.exome;

import org.apache.commons.collections4.list.SetUniqueList;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public final class ReCapSegCallerUnitTest extends BaseTest{
    @Test
    public void testMakeCalls() {

        final List<Target> targets = new ArrayList<>();
        final List<String> columnNames = Arrays.asList("Sample");
        final List<Double> coverage = new ArrayList<>();

        //add amplification targets
        for (int i = 0; i < 10; i++) {
            final SimpleInterval interval = new SimpleInterval("chr", 100 + 2 * i, 101 + 2 * i);
            targets.add(new Target(TargetUtils.createDummyTargetName(interval), interval));
            coverage.add(ParamUtils.log2(2.0));
        }
        //add deletion targets
        for (int i = 0; i < 10; i++) {
            final SimpleInterval interval = new SimpleInterval("chr", 300 + 2 * i, 301 + 2 * i);
            targets.add(new Target(TargetUtils.createDummyTargetName(interval), interval));
            coverage.add(ParamUtils.log2(0.5));
        }
        //add targets that don't belong to a segment
        for (int i = 1; i < 10; i++) {
            final SimpleInterval interval = new SimpleInterval("chr", 400 + 2 * i, 401 + 2 * i);
            targets.add(new Target(TargetUtils.createDummyTargetName(interval), interval));
            coverage.add(ParamUtils.log2(1.0));
        }
        //add obviously neutral targets with some small spread
        for (int i = -5; i < 6; i++) {
            final SimpleInterval interval = new SimpleInterval("chr", 500 + 2 * i, 501 + 2 * i);
            targets.add(new Target(TargetUtils.createDummyTargetName(interval), interval));
            coverage.add(ParamUtils.log2(0.01 * i + 1));
        }
        //add spread-out targets to a neutral segment (mean near zero)
        for (int i = -5; i < 6; i++) {
            final SimpleInterval interval = new SimpleInterval("chr", 700 + 2 * i, 701 + 2 * i);
            targets.add(new Target(TargetUtils.createDummyTargetName(interval), interval));
            coverage.add(ParamUtils.log2(0.1 * i + 1));
        }

        final RealMatrix coverageMatrix = new Array2DRowRealMatrix(targets.size(), 1);
        coverageMatrix.setColumn(0, coverage.stream().mapToDouble(x->x).toArray());
        final int n = targets.size();
        final int m = coverageMatrix.getRowDimension();
        final ReadCountCollection counts = new ReadCountCollection(targets, columnNames, coverageMatrix);

        List<ModeledSegment> segments = new ArrayList<>();
        segments.add(new ModeledSegment(new SimpleInterval("chr", 100, 200), 100, ParamUtils.log2(2.0))); //amplification
        segments.add(new ModeledSegment(new SimpleInterval("chr", 300, 400), 100, ParamUtils.log2(0.5))); //deletion
        segments.add(new ModeledSegment(new SimpleInterval("chr", 450, 550), 100, ParamUtils.log2(1))); //neutral
        segments.add(new ModeledSegment(new SimpleInterval("chr", 650, 750), 100, ParamUtils.log2(1))); //neutral

        List<ModeledSegment> calls = ReCapSegCaller.makeCalls(counts, segments);

        Assert.assertEquals(calls.get(0).getCall(), ReCapSegCaller.AMPLIFICATION_CALL);
        Assert.assertEquals(calls.get(1).getCall(), ReCapSegCaller.DELETION_CALL);
        Assert.assertEquals(calls.get(2).getCall(), ReCapSegCaller.NEUTRAL_CALL);
        Assert.assertEquals(calls.get(3).getCall(), ReCapSegCaller.NEUTRAL_CALL);
    }
}
