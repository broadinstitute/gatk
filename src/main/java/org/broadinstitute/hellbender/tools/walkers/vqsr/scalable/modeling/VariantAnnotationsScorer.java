package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling;

import java.io.File;
import java.util.function.Function;

public interface VariantAnnotationsScorer {
    void scoreSamples(final File inputAnnotationsFile,
                      final File outputScoresFile);

    // TODO
    static Function<Double, Double> createScoreToTruthSensitivityConverter(final double[] truthScores) {
        return x -> x;
    };
}
