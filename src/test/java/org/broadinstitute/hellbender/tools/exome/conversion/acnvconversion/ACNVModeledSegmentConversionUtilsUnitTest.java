package org.broadinstitute.hellbender.tools.exome.conversion.acnvconversion;

import org.broadinstitute.hellbender.tools.exome.*;
import org.broadinstitute.hellbender.tools.pon.coverage.pca.PCATangentNormalizationUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.mcmc.PosteriorSummary;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class ACNVModeledSegmentConversionUtilsUnitTest extends BaseTest {

    @Test
    public void testSimpleConversionCannotYieldSegmentMeanOfZero() {
        final ACNVModeledSegment acnvModeledSegment = new ACNVModeledSegment(new SimpleInterval("1", 1000, 1500),
                new PosteriorSummary(-4000, -4001, -4002),
                new PosteriorSummary(-4000, -4001, -4002));
        final List<Target> targets = new ArrayList<>();
        targets.add(new Target("test", new SimpleInterval("1", 1300, 1302)));
        final double[] targetDummyValues = new double[targets.size()];
        final TargetCollection<ReadCountRecord.SingleSampleRecord> targetCollection = new HashedListTargetCollection<>(Collections.singletonList(new ReadCountRecord.SingleSampleRecord(targets.get(0), 0.0)));
        final ModeledSegment guess = ACNVModeledSegmentConversionUtils.convertACNVModeledSegmentToModeledSegment(acnvModeledSegment, targetCollection);
        Assert.assertTrue(guess.getSegmentMeanInCRSpace() > 0);
        Assert.assertEquals(guess.getSegmentMean(), ParamUtils.log2(PCATangentNormalizationUtils.EPSILON), 1e-10);
    }
}
