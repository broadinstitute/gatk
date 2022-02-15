package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling;

import java.io.File;

public interface VariantAnnotationsModel {

    /**
     * @param trainingAnnotationsFile Training data in HDF5 format. (TODO document paths)
     *                             A single model will be trained, so either
     *                             1) this data has already been subset to a single variant type (SNP or INDEL), or
     *                             2) we assume the model does not care about the variant type.
     *
     * @param outputPrefix
     */
    void trainAndSerialize(final File trainingAnnotationsFile,
                           final String outputPrefix);
}
