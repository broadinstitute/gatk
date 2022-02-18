package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling;

import java.io.File;

/**
 * File interface for passing annotations to a modeling backend and indicating a path prefix for resulting output.
 */
public interface VariantAnnotationsModel {

    /**
     * @param trainingAnnotationsFile   Training annotations in HDF5 format, containing at least the directory structure
     *
     *                                   ├── annotations
     *                                        ├── chunk_0
     *                                        ├── ...
     *                                        ├── chunk_{num_chunks - 1}
     *                                        ├── names
     *                                        ├── num_chunks
     *                                        ├── num_columns
     *                                        └── num_rows
     *
     *                                  Modeling backends are responsible for consuming annotations in this format.
     *                                  In current use, we assume that a single model will be trained, so either
     *                                    1) training annotations have already been subset to a single variant type (SNP or INDEL), or
     *                                    2) we assume the model does not care about the variant type.
     *                                  TODO we could also pass additional labels to be used in training,
     *                                       but all backends would have to likewise respect directory structure
     * @param outputPrefix              Path prefix for all output files
     */
    void trainAndSerialize(final File trainingAnnotationsFile,
                           final String outputPrefix);
}
