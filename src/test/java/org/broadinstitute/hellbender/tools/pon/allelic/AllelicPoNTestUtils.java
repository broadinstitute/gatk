package org.broadinstitute.hellbender.tools.pon.allelic;

import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;

import java.util.List;
import java.util.Map;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class AllelicPoNTestUtils {
    private static final double DOUBLE_TOLERANCE = 1E-6;

    public static void assertAllelicPoNsEqual(final AllelicPanelOfNormals result, final AllelicPanelOfNormals expected) {
        Assert.assertEquals(result.getGlobalHyperparameterValues().getAlpha(), expected.getGlobalHyperparameterValues().getAlpha(), DOUBLE_TOLERANCE);
        Assert.assertEquals(result.getGlobalHyperparameterValues().getBeta(), expected.getGlobalHyperparameterValues().getBeta(), DOUBLE_TOLERANCE);
        final List<Map.Entry<SimpleInterval, AllelicPanelOfNormals.HyperparameterValues>> resultEntries = result.getSortedMapEntries();
        final List<Map.Entry<SimpleInterval, AllelicPanelOfNormals.HyperparameterValues>> expectedEntries = expected.getSortedMapEntries();
        Assert.assertEquals(resultEntries.size(), expectedEntries.size());
        for (int i = 0; i < resultEntries.size(); i++) {
            final Map.Entry<SimpleInterval, AllelicPanelOfNormals.HyperparameterValues> resultEntry = resultEntries.get(i);
            final Map.Entry<SimpleInterval, AllelicPanelOfNormals.HyperparameterValues> expectedEntry = expectedEntries.get(i);
            Assert.assertEquals(resultEntry.getKey(), expectedEntry.getKey());                                              //assert site SimpleIntervals are equal
            Assert.assertEquals(resultEntry.getValue().getAlpha(), expectedEntry.getValue().getAlpha(), DOUBLE_TOLERANCE);  //assert hyperparameter alphas are equal
            Assert.assertEquals(resultEntry.getValue().getBeta(), expectedEntry.getValue().getBeta(), DOUBLE_TOLERANCE);    //assert hyperparameter betas are equal
        }
    }
}
