package org.broadinstitute.hellbender.tools.exome.orientationbiasvariantfilter;

import htsjdk.samtools.metrics.MetricsFile;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.tools.picard.analysis.artifacts.SequencingArtifactMetrics;
import org.broadinstitute.hellbender.tools.picard.analysis.artifacts.Transition;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.FileReader;
import java.io.IOException;
import java.util.Map;

public class PreAdapterOrientationScorerUnitTest extends BaseTest {

    public static final String testPreAdapterDetailsMetrics = publicTestDir + "picard_metrics_test.pre_adapter_detail_metrics";

    /**
     * Note that (due to raw data), this test includes collapsing over libraries (not just the contexts).
     * @throws IOException
     */
    @Test
    public void testBasicScoring() throws IOException {
        final MetricsFile<SequencingArtifactMetrics.PreAdapterDetailMetrics, Comparable<?>> mf = new MetricsFile<>();
        mf.read(new FileReader(testPreAdapterDetailsMetrics));

        final Map<Transition, Double> scoreMap = PreAdapterOrientationScorer.scoreOrientationBiasMetricsOverContext(mf.getMetrics());

        Assert.assertNotNull(scoreMap);
        Assert.assertEquals(scoreMap.keySet().size(), 12);

        // Ground truth values painstakingly derived manually
        Assert.assertEquals(scoreMap.get(Transition.transitionOf('A', 'C')), 100.0, 1e-6);
        Assert.assertEquals(scoreMap.get(Transition.transitionOf('A', 'G')), 50.5788416297570, 1e-6);
        Assert.assertEquals(scoreMap.get(Transition.transitionOf('A', 'T')), 100.0, 1e-6);
        Assert.assertEquals(scoreMap.get(Transition.transitionOf('C', 'A')), 100.0, 1e-6);
        Assert.assertEquals(scoreMap.get(Transition.transitionOf('C', 'G')), 100.0, 1e-6);
        Assert.assertEquals(scoreMap.get(Transition.transitionOf('C', 'T')), 58.0641821538479, 1e-6);
    }

    /**
     * Note that (due to raw data), this test includes collapsing over libraries (not just the contexts).
     * @throws IOException
     */
    @Test
    public void testBasicCounting() throws IOException {
        final MetricsFile<SequencingArtifactMetrics.PreAdapterDetailMetrics, Comparable<?>> mf = new MetricsFile<>();
        mf.read(new FileReader(testPreAdapterDetailsMetrics));

        final Map<Transition, RealMatrix> countMap = PreAdapterOrientationScorer.countOrientationBiasMetricsOverContext(mf.getMetrics());

        Assert.assertNotNull(countMap);
        Assert.assertEquals(countMap.keySet().size(), 12);

        // Ground truth values painstakingly derived manually from both libraries (hence four values in the expression instead of two)
        // L1 norm can be used since all values are always positive.
        Assert.assertEquals(countMap.get(Transition.transitionOf('A', 'C')).getRowVector(0).getL1Norm(), 2836660 + 2631203 + 240 + 246, 1e-1);
        Assert.assertEquals(countMap.get(Transition.transitionOf('A', 'C')).getRowVector(1).getL1Norm(), 2852491 + 333 + 2646654 + 297, 1e-1);
        Assert.assertEquals(countMap.get(Transition.transitionOf('A', 'G')).getRowVector(0).getL1Norm(), 2631203 + 416 + 2836660 + 481, 1e-1);
        Assert.assertEquals(countMap.get(Transition.transitionOf('A', 'G')).getRowVector(1).getL1Norm(), 2852491 + 404 + 2646654 + 397, 1e-1);
    }
}
