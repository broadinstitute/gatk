package org.broadinstitute.hellbender.tools.copynumber.coverage.caller;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.copynumber.coverage.copyratio.CopyRatio;
import org.broadinstitute.hellbender.tools.copynumber.coverage.copyratio.CopyRatioCollection;
import org.broadinstitute.hellbender.tools.copynumber.coverage.segmentation.CopyRatioSegment;
import org.broadinstitute.hellbender.tools.copynumber.coverage.segmentation.CopyRatioSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleSampleMetadata;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static org.broadinstitute.hellbender.tools.copynumber.coverage.caller.CalledCopyRatioSegment.Call.*;

public final class ReCapSegCallerUnitTest extends GATKBaseTest {
    @Test
    public void testMakeCalls() {

        final SampleMetadata sampleMetadata = new SimpleSampleMetadata("Sample");
        final List<SimpleInterval> intervals = new ArrayList<>();
        final List<Double> testData = new ArrayList<>();

        //add amplification intervals
        for (int i = 0; i < 10; i++) {
            final SimpleInterval interval = new SimpleInterval("chr", 101 + i, 101 + i);
            intervals.add(interval);
            testData.add(ParamUtils.log2(2.0));
        }
        //add deletion intervals
        for (int i = 0; i < 10; i++) {
            final SimpleInterval interval = new SimpleInterval("chr", 201 + i, 201 + i);
            intervals.add(interval);
            testData.add(ParamUtils.log2(0.5));
        }
        //add obviously neutral intervals with some small spread
        for (int i = 0; i < 10; i++) {
            final SimpleInterval interval = new SimpleInterval("chr", 301 + i, 301 + i);
            intervals.add(interval);
            testData.add(ParamUtils.log2(0.01 * (i - 5) + 1));
        }
        //add spread-out intervals to a neutral segment (mean near zero)
        for (int i = 0; i < 10; i++) {
            final SimpleInterval interval = new SimpleInterval("chr", 401 + i, 401 + i);
            intervals.add(interval);
            testData.add(ParamUtils.log2(0.1 * (i - 5) + 1));
        }

        final RealMatrix denoisedCopyRatioValues = new Array2DRowRealMatrix(1, intervals.size());
        denoisedCopyRatioValues.setRow(0, testData.stream().mapToDouble(x -> x).toArray());
        final CopyRatioCollection denoisedCopyRatios = new CopyRatioCollection(
                sampleMetadata,
                IntStream.range(0, intervals.size())
                        .mapToObj(i -> new CopyRatio(intervals.get(i), denoisedCopyRatioValues.getEntry(0, i)))
                        .collect(Collectors.toList()));

        final CopyRatioSegmentCollection copyRatioSegments = new CopyRatioSegmentCollection(sampleMetadata,
                Arrays.asList(
                        new CopyRatioSegment(new SimpleInterval("chr", 101, 110), 10, ParamUtils.log2(2.0)),   //amplification
                        new CopyRatioSegment(new SimpleInterval("chr", 201, 210), 10, ParamUtils.log2(0.5)),   //deletion
                        new CopyRatioSegment(new SimpleInterval("chr", 301, 310), 10, ParamUtils.log2(1)),     //neutral
                        new CopyRatioSegment(new SimpleInterval("chr", 401, 410), 10, ParamUtils.log2(1))));   //neutral

        final CalledCopyRatioSegmentCollection calledCopyRatioSegments = new ReCapSegCaller(denoisedCopyRatios, copyRatioSegments).makeCalls();

        Assert.assertEquals(copyRatioSegments.getSampleName(), calledCopyRatioSegments.getSampleName());
        Assert.assertEquals(
                copyRatioSegments.getIntervals(), calledCopyRatioSegments.getIntervals());
        Assert.assertEquals(
                copyRatioSegments.getRecords().stream().map(CopyRatioSegment::getMeanLog2CopyRatio).collect(Collectors.toList()),
                calledCopyRatioSegments.getRecords().stream().map(CopyRatioSegment::getMeanLog2CopyRatio).collect(Collectors.toList()));
        Assert.assertEquals(
                calledCopyRatioSegments.getRecords().stream().map(CalledCopyRatioSegment::getCall).collect(Collectors.toList()),
                Arrays.asList(AMPLIFICATION, DELETION, NEUTRAL, NEUTRAL));
    }
}