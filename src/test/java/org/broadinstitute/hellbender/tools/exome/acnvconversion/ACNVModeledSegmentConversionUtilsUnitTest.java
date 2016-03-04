package org.broadinstitute.hellbender.tools.exome.acnvconversion;

import org.broadinstitute.hellbender.tools.exome.*;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.mcmc.PosteriorSummary;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;

public class ACNVModeledSegmentConversionUtilsUnitTest extends BaseTest {

    @Test
    public void testSimpleConversionCannotYieldSegmentMeanOfZero() {
        final ACNVModeledSegment acnvModeledSegment = new ACNVModeledSegment(new SimpleInterval("1", 1000, 1500),
                new PosteriorSummary(-4000, -4001, -4002),
                new PosteriorSummary(-4000, -4001, -4002));
        final List<TargetCoverage> targets = new ArrayList<>();
        targets.add(new TargetCoverage("test", new SimpleInterval("1", 1300, 1302), 40));

        final TargetCollection<TargetCoverage> targetCollection = new HashedListTargetCollection<>(targets);
        final ModeledSegment guess = ACNVModeledSegmentConversionUtils.convertACNVModeledSegmentToModeledSegment(acnvModeledSegment, targetCollection);
        Assert.assertTrue(guess.getSegmentMeanInCRSpace() > 0);
    }
}
