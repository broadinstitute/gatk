package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import org.apache.commons.math3.random.RandomDataGenerator;
import org.broadinstitute.hellbender.utils.Utils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;

/**
 * Created by tsato on 4/19/17.
 */
public class StrandArtifactFilterUnitTest {
    private static double DEFAULT_PRIOR = 0.001;

    // prior, forward count, reverse count, forward alt count, forward reverse count, is artifact
    @DataProvider(name = "ObviousCalls")
    public Object[][] makeObviousArtifactsData() {
        return new Object[][]{
                {DEFAULT_PRIOR, 100, 100, 20, 0, true},
                {DEFAULT_PRIOR, 100, 100, 50, 0, true},
                {DEFAULT_PRIOR, 100, 100, 100, 0, true},
                {DEFAULT_PRIOR, 100, 100, 100, 1, true},
                {DEFAULT_PRIOR, 100, 100, 100, 3, true},

                {DEFAULT_PRIOR, 100, 100, 50, 50, false},
                {DEFAULT_PRIOR, 100, 100, 70, 30, false},
                {DEFAULT_PRIOR, 100, 10, 50, 5, false},
                {DEFAULT_PRIOR, 100, 1, 50, 1, false},
        };
    }

    // obvious strand bias
    @Test(dataProvider = "ObviousCalls")
    public void testObviousCalls(final double prior, final int forwardCount, final int reverseCount, final int forwardAltCount, final int reverseAltCount, final boolean isArtifact) throws IOException {
        final double artifactProbability = new StrandArtifactFilter()
                .strandArtifactProbability(prior, forwardCount, reverseCount, forwardAltCount, reverseAltCount)
                .getArtifactProbability();
        Assert.assertEquals(artifactProbability, isArtifact ? 1.0 : 0.0, 1.0e-2);
    }

    @Test
    public void testSymmetry() {
        Utils.resetRandomGenerator();
        final RandomDataGenerator rdg = Utils.getRandomDataGenerator();

        for (int n = 0; n < 10; n++) {
            final double prior = rdg.nextUniform(0,1);
            final int forwardCount = rdg.nextInt(0, 100);
            final int reverseCount = rdg.nextInt(0, 100);
            final int forwardAltCount = rdg.nextInt(0, forwardCount);
            final int reverseAltCount = rdg.nextInt(0, reverseCount);

            final double prob = new StrandArtifactFilter()
                    .strandArtifactProbability(prior, forwardCount, reverseCount, forwardAltCount, reverseAltCount)
                    .getArtifactProbability();
            final double flipped = new StrandArtifactFilter()
                    .strandArtifactProbability(prior, reverseCount, forwardCount, reverseAltCount, forwardAltCount)
                    .getArtifactProbability();

            Assert.assertEquals(prob, flipped, 1.0e-8);
        }
    }
}