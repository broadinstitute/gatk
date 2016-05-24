package org.broadinstitute.hellbender.tools.exome.gcbias;

import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.fakedata.GCBiasSimulatedData;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;

/**
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public class GCCorrectorUnitTest {

    @Test
    public void testGCCorrection() {
        final int numSamples = 5;
        final int numTargets = 10000;

        final Pair<ReadCountCollection, double[]> data = GCBiasSimulatedData.simulatedData(numTargets, numSamples);
        final ReadCountCollection rcc1 = data.getLeft();
        final double[] gcContentByTarget = data.getRight();
        final ReadCountCollection correctedCounts1 = GCCorrector.correctCoverage(rcc1, gcContentByTarget);
        final double[] correctedNoiseBySample = GATKProtectedMathUtils.columnStdDevs(correctedCounts1.counts());
        Arrays.stream(correctedNoiseBySample).forEach(x -> Assert.assertTrue(x < 0.02));

        //check that GC correction is approximately idempotent -- if you correct again, very little should happen
        final ReadCountCollection correctedCounts2 = GCCorrector.correctCoverage(correctedCounts1, gcContentByTarget);
        final double change1 = correctedCounts1.counts().subtract(rcc1.counts()).getFrobeniusNorm();
        final double change2 = correctedCounts2.counts().subtract(correctedCounts1.counts()).getFrobeniusNorm();
        Assert.assertTrue(change2 < change1 / 10);
    }

}