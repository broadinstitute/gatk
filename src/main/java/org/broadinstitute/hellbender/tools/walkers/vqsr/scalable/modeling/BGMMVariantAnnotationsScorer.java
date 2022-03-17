package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling;

import org.apache.commons.lang.NotImplementedException;
import org.broadinstitute.hellbender.utils.clustering.BayesianGaussianMixtureModeller;

import java.io.File;
import java.io.Serializable;
import java.util.List;

public final class BGMMVariantAnnotationsScorer implements VariantAnnotationsScorer, Serializable {
    private static final long serialVersionUID = 1L;

    public BGMMVariantAnnotationsScorer(final List<String> annotationNames,
                                        final BGMMVariantAnnotationsModel.Preprocesser preprocesser,
                                        final BayesianGaussianMixtureModeller bgmm) {
        throw new NotImplementedException("BGMM module implemented in separate PR.");
    }

    @Override
    public void score(final File inputAnnotationsFile,
                      final File outputScoresFile) {
        throw new NotImplementedException("BGMM module implemented in separate PR.");
    }

    public double[][] preprocess(final double[][] annotations) {
        throw new NotImplementedException("BGMM module implemented in separate PR.");
    }

    public void serialize(final File scorerFile) {
        throw new NotImplementedException("BGMM module implemented in separate PR.");
    }

    public static BGMMVariantAnnotationsScorer deserialize(final File scorerFile) {
        throw new NotImplementedException("BGMM module implemented in separate PR.");
    }
}
