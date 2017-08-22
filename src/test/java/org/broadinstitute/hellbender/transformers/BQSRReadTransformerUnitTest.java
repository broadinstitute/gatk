package org.broadinstitute.hellbender.transformers;

import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.recalibration.RecalDatum;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;

public final class BQSRReadTransformerUnitTest extends GATKBaseTest {
    @Test
    public void basicHierarchicalBayesianQualityEstimateTest() {

        for( double epsilon = 15.0; epsilon <= 60.0; epsilon += 2.0 ) {
            double RG_Q = 45.0;
            RecalDatum RG = makeRecalDatum(100000000L, RG_Q);
            double Q = 30.0;
            RecalDatum QS = makeRecalDatum(100000000L, Q);
            RecalDatum COV = new RecalDatum(15L, 1.0, (byte)45.0); // no data here so Bayesian prior has a huge effect on the empirical quality

            // initial epsilon condition shouldn't matter when there are a lot of observations
            Assert.assertEquals(BQSRReadTransformer.hierarchicalBayesianQualityEstimate(epsilon, RG, QS, COV), Q, 1.0E-4);
        }

        for( double epsilon = 15.0; epsilon <= 60.0; epsilon += 2.0 ) {
            double RG_Q = 45.0;
            RecalDatum RG = makeRecalDatum(10L, RG_Q);
            double Q = 30.0;
            RecalDatum QS = makeRecalDatum(10L,Q);
            RecalDatum COV = new RecalDatum(15L, 1.0, (byte)45.0); // no data here so Bayesian prior has a huge effect on the empirical quality

            // initial epsilon condition dominates when there is no data
            Assert.assertEquals(BQSRReadTransformer.hierarchicalBayesianQualityEstimate(epsilon, RG, QS, COV), epsilon, 1.0E-4);
        }
    }

    private static  RecalDatum makeRecalDatum(final long count, final double qual){
        return new RecalDatum(count, count * 1.0 / (Math.pow(10.0, qual / 10.0)), (byte)qual);
    }

    @Test
    public void repeatedAndUnorderedFixedQualities() {
        // Test both repeated quals, and quals that aren't input in order
        List<Integer> quantizedQualsOrdered = Arrays.asList(11, 19);
        List<Integer> quantizedQualsUnordered = Arrays.asList(19, 11, 19, 19);

        // Unordered and Ordered qmapping should be identical
        byte[] qmappingUnordered = BQSRReadTransformer.constructStaticQuantizedMapping(quantizedQualsUnordered, true);
        byte[] qmappingOrdered = BQSRReadTransformer.constructStaticQuantizedMapping(quantizedQualsOrdered, true);
        Assert.assertEquals(qmappingOrdered.length, qmappingUnordered.length);
        for(int i = 0 ; i < qmappingUnordered.length ; i++) {
            Assert.assertEquals(qmappingOrdered[i], qmappingUnordered[i]);
        }
    }

    @Test
    public void nearestVsRoundDown() {
        List<Integer> fixedQuantizedQuals = Arrays.asList(10, 20, 30);

        byte[] qmappingRoundDown = BQSRReadTransformer.constructStaticQuantizedMapping(fixedQuantizedQuals, true);
        byte[] qmappingRoundNearest = BQSRReadTransformer.constructStaticQuantizedMapping(fixedQuantizedQuals, false);

        // Depending on rounding strategy, bin 19 should round to 10 or 20
        Assert.assertEquals(qmappingRoundDown[19], 10);
        Assert.assertEquals(qmappingRoundNearest[19], 20);

        // Regarless of rounding strategy, bin 21 should always round down to 20
        Assert.assertEquals(qmappingRoundDown[21], 20);
        Assert.assertEquals(qmappingRoundNearest[21], 20);
    }

    @Test
    public void onlyOneFixedQualUsed() {
        // Set all qualities to singleQual value (except for those below MIN_USABLE_Q_SCORE)
        int singleQual = 10;
        List<Integer> fixedQuantizedQuals = Arrays.asList(singleQual);

        byte[] qmapping = BQSRReadTransformer.constructStaticQuantizedMapping(fixedQuantizedQuals, true);

        for(int i = 0 ; i < qmapping.length ; i++) {
            if(i >= QualityUtils.MIN_USABLE_Q_SCORE) {
                Assert.assertEquals(qmapping[i], singleQual);
            }
            else {
                // Make sure that all values less than MIN_USABLE_Q_SCORE are preserved
                Assert.assertEquals(qmapping[i], i);
            }
        }
    }
}
