package org.broadinstitute.hellbender.tools.exome.orientationbiasvariantfilter;

import com.google.cloud.dataflow.sdk.repackaged.com.google.common.annotations.VisibleForTesting;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.tools.picard.analysis.artifacts.SequencingArtifactMetrics;
import org.broadinstitute.hellbender.tools.picard.analysis.artifacts.Transition;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * This class handles the reading of a picard metrics file PreAdapterDetailMetrics.
 */
public class PreAdapterOrientationScorer {

    // rows {PRO, CON} x cols {ref, alt}
    final static int PRO = 0;
    final static int CON = 1;
    final static int REF = 0;
    final static int ALT = 1;

    /** Defined by convention that the Phred quality score is 100.0, so this value ensures that the Phred calculation will
     * be 100.0. */
    public static final double MAX_BASE_SCORE = Math.pow(10, -10);

    /**
     * Gets PreAdapterQ collapsed over contexts.
     *
     * rows {PRO, CON} x cols {ref, alt}
     * @param metrics metrics usually read from a picard preadapter detail file.
     * @return mapping to score all orientation bias artifact modes to a PreAdapterQ score.  This score can be used as a bam-file level score for level of artifacts.
     */
    @VisibleForTesting
    static Map<Transition, RealMatrix> countOrientationBiasMetricsOverContext(final List<SequencingArtifactMetrics.PreAdapterDetailMetrics> metrics) {
        Utils.nonNull(metrics, "Input metrics cannot be null");

        // Artifact mode to a matrix
        final Map<Transition, RealMatrix> result = new HashMap<>();

        // Collapse over context
        for (SequencingArtifactMetrics.PreAdapterDetailMetrics metric : metrics) {
            final Transition key = Transition.transitionOf(metric.REF_BASE, metric.ALT_BASE);
            result.putIfAbsent(key, new Array2DRowRealMatrix(2, 2));
            final RealMatrix preAdapterCountMatrix = result.get(key);
            preAdapterCountMatrix.addToEntry(PreAdapterOrientationScorer.PRO, PreAdapterOrientationScorer.ALT, metric.PRO_ALT_BASES);
            preAdapterCountMatrix.addToEntry(PreAdapterOrientationScorer.CON, PreAdapterOrientationScorer.ALT, metric.CON_ALT_BASES);
            preAdapterCountMatrix.addToEntry(PreAdapterOrientationScorer.PRO, PreAdapterOrientationScorer.REF, metric.PRO_REF_BASES);
            preAdapterCountMatrix.addToEntry(PreAdapterOrientationScorer.CON, PreAdapterOrientationScorer.REF, metric.CON_REF_BASES);
        }

        return result;
    }

    /** Scores all modes regardless of whether required in the input.
     *
     * @param metrics Can be read from a Picard PreAdapterDetail file.  Never {@code null}
     * @return mapping to score all orientation bias artifact modes to a PreAdapterQ score.  This score can be used as a bam-file level score for level of artifacts.
     */
    public static Map<Transition, Double> scoreOrientationBiasMetricsOverContext(final List<SequencingArtifactMetrics.PreAdapterDetailMetrics> metrics) {
        Utils.nonNull(metrics, "Input metrics cannot be null");

        // Artifact mode to a double
        final Map<Transition, Double> result = new HashMap<>();

        final Map<Transition, RealMatrix> counts = countOrientationBiasMetricsOverContext(metrics);

        for (final Transition transition : counts.keySet()) {
            final RealMatrix count = counts.get(transition);
            final double totalBases = count.getEntry(PRO, REF) + count.getEntry(PRO, ALT) +
                    count.getEntry(CON, REF) + count.getEntry(CON, ALT);
            final double rawFractionPro = count.getEntry(PRO, ALT)/totalBases;
            final double rawFractionCon = count.getEntry(CON, ALT)/totalBases;
            final double baseScore = rawFractionPro - rawFractionCon;
            final double score = QualityUtils.phredScaleErrorRate(Math.max(baseScore, MAX_BASE_SCORE));
            result.put(transition, score);
        }

        return result;
    }

    /** Do not allow instantiation of this class */
    private PreAdapterOrientationScorer() {};
}
