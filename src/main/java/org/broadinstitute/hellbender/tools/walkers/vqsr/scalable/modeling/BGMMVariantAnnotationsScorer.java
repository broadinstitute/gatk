package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling;

import com.google.common.collect.ImmutableList;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.utils.clustering.BayesianGaussianMixtureModeller;

import java.io.File;
import java.io.Serializable;
import java.util.List;

public final class BGMMVariantAnnotationsScorer implements VariantAnnotationsScorer, Serializable {
    private static final long serialVersionUID = 1L;

    private final List<String> annotationNames;     // for validation
    private final BGMMVariantAnnotationsModel.Preprocesser preprocesser;
    private final BayesianGaussianMixtureModeller bgmm;

    public BGMMVariantAnnotationsScorer(final List<String> annotationNames,
                                        final BGMMVariantAnnotationsModel.Preprocesser preprocesser,
                                        final BayesianGaussianMixtureModeller bgmm) {
        this.annotationNames = ImmutableList.copyOf(annotationNames);
        this.preprocesser = preprocesser;
        this.bgmm = bgmm;
    }

    public static BGMMVariantAnnotationsScorer deserialize(final String pathPrefix) {
        return SerializationUtils.deserialize(
                new File(pathPrefix + BGMMVariantAnnotationsModel.BGMM_SCORER_SER_SUFFIX),
                BGMMVariantAnnotationsScorer.class);
    }

    private Pair<double[][], double[]> preprocessAndScoreSamples(final double[][] data) {
        final double[][] preprocessedData = preprocesser.transform(data);
        final double[] scores = bgmm.scoreSamples(preprocessedData);
        return Pair.of(preprocessedData, scores);
    }

    @Override
    public void scoreSamples(final File inputAnnotationsFile,
                             final File outputScoresFile) {

    }
}
