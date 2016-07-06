package org.broadinstitute.hellbender.tools.exome.alleliccount;

import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.Test;

/**
 * Tests for {@link AllelicCountWithPhasePosteriors}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public class AllelicCountWithPhasePosteriorsUnitTest {
    private static final double DELTA = 1E-10;

    @Test
    public void testNormalizationOfPosteriorProbabilities() {
        final AllelicCount count = new AllelicCount(new SimpleInterval("1", 1, 1), 10, 15);

        //create unnormalized log posterior probabilities
        final double normalizationFactor = 5;
        final double refMinorProb = 0.3;
        final double altMinorProb = 0.6;
        final double outlierProb = 0.1;
        final double refMinorUnnormalizedProb = refMinorProb * normalizationFactor;
        final double altMinorUnnormalizedProb = altMinorProb * normalizationFactor;
        final double outlierUnnormalizedProb = outlierProb * normalizationFactor;
        final double refMinorUnnormalizedLogProb = Math.log(refMinorUnnormalizedProb);
        final double altMinorUnnormalizedLogProb = Math.log(altMinorUnnormalizedProb);
        final double outlierUnnormalizedLogProb = Math.log(outlierUnnormalizedProb);

        //constructor should normalize posterior probabilities
        final AllelicCountWithPhasePosteriors countWithPhasePosteriors =
                new AllelicCountWithPhasePosteriors(count, refMinorUnnormalizedLogProb, altMinorUnnormalizedLogProb, outlierUnnormalizedLogProb);

        //getters should return un-logged posterior probabilities
        Assert.assertEquals(countWithPhasePosteriors.getRefMinorProb(), refMinorProb, DELTA);
        Assert.assertEquals(countWithPhasePosteriors.getAltMinorProb(), altMinorProb, DELTA);
        Assert.assertEquals(countWithPhasePosteriors.getOutlierProb(), outlierProb, DELTA);
    }
}