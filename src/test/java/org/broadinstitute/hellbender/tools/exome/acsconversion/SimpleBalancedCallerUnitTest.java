package org.broadinstitute.hellbender.tools.exome.acsconversion;

import org.broadinstitute.hellbender.tools.exome.ACNVModeledSegment;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.mcmc.PosteriorSummary;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

public class SimpleBalancedCallerUnitTest extends BaseTest{

    @Test
    public void simpleCallerTest() {
        final SimpleBalancedSegmentCaller caller = new SimpleBalancedSegmentCaller();

        final PosteriorSummary segmentMeanPosteriorSummary = new PosteriorSummary(0.001, -0.001, 0.002);
        final PosteriorSummary minorAlleleFractionPosteriorSummary = new PosteriorSummary(.48, .47, .491);

        final ACNVModeledSegment testSeg = new ACNVModeledSegment(new SimpleInterval("1", 1000, 1010), segmentMeanPosteriorSummary,
                minorAlleleFractionPosteriorSummary);

        Assert.assertTrue(caller.isSegmentBalanced(testSeg));

        Assert.assertFalse(caller.isSegmentBalanced(new ACNVModeledSegment(new SimpleInterval("1", 1000, 1010),
                segmentMeanPosteriorSummary,
                new PosteriorSummary(.48, .47, .489))));
    }
}
