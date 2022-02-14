package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.model;

import java.io.File;
import java.util.List;

interface VariantAnnotationsModeller {
    /**
     *
     * @return
     */
    List<String> getAnnotationNames();

    /**
     * @param inputAnnotationsFile Training data in HDF5 format. (TODO document paths)
     *                             A single model will be trained, so either
     *                             1) this data has already been subset to a single variant type (SNP or INDEL), or
     *                             2) we assume the model does not care about the variant type.
     *
     * @param outputPrefix
     */
    void trainAndSerialize(final File inputAnnotationsFile,
                           final String outputPrefix);

    void scoreSamplesAndWriteScores(final String inputPrefix,
                                    final File inputAnnotationsFile,
                                    final File outputScoresFile);
}
