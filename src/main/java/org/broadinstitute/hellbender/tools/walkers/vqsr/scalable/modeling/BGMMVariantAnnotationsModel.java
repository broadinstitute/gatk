package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling;

import org.apache.commons.lang.NotImplementedException;

import java.io.File;
import java.io.Serializable;

public final class BGMMVariantAnnotationsModel implements VariantAnnotationsModel {

    public static final String BGMM_SCORER_SER_SUFFIX = ".bgmmScorer.ser";

    public BGMMVariantAnnotationsModel(final File hyperparametersJSONFile) {
        throw new NotImplementedException("BGMM module implemented in separate PR.");
    }

    @Override
    public void trainAndSerialize(final File trainingAnnotationsFile,
                                  final String outputPrefix) {
        throw new NotImplementedException("BGMM module implemented in separate PR.");
    }

    static final class Preprocesser implements Serializable {

        Preprocesser() {
        }

        double[][] transform(final double[][] data) {
            throw new NotImplementedException("BGMM module implemented in separate PR.");
        }
    }
}