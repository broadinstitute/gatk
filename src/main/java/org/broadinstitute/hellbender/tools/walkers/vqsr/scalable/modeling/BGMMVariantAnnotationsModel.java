package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling;

import org.apache.commons.lang.NotImplementedException;

import java.io.File;
import java.io.Serializable;

public final class BGMMVariantAnnotationsModel implements VariantAnnotationsModel {

    public BGMMVariantAnnotationsModel(final File hyperparametersJSONFile) {
        throw new NotImplementedException("BGMM module implemented in separate PR.");
    }

    @Override
    public void trainAndSerialize(final File trainingAnnotationsFile,
                                  final String outputPrefix) {
        throw new NotImplementedException("BGMM module implemented in separate PR.");
    }

    static final class Preprocesser implements Serializable {
        private static final long serialVersionUID = 1L;

        Preprocesser() {
        }

        double[][] transform(final double[][] data) {
            throw new NotImplementedException("BGMM module implemented in separate PR.");
        }
    }
}