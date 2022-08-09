package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling;

import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.data.LabeledVariantAnnotationsData;

import java.io.File;

/**
 * File interface for passing annotations to a modeling backend and indicating a path prefix for resulting output.
 */
public interface VariantAnnotationsModel {

    /**
     * @param trainingAnnotationsFile   Training annotations in HDF5 format, containing at least the directory structure
     *
     *                                  <p>
     *                                    |--- annotations<br>
     *                                    |    |--- chunk_0<br>
     *                                    |    |--- ...<br>
     *                                    |    |--- chunk_{num_chunks - 1}<br>
     *                                    |    |--- names<br>
     *                                    |    |--- num_chunks<br>
     *                                    |    |--- num_columns<br>
     *                                    |    |--- num_rows<br>
     *                                  </p>
     *
     *                                  Here, each chunk is a double matrix, with dimensions given by
     *                                  (number of sites in the chunk) x (number of annotations).
     *                                  See {@link LabeledVariantAnnotationsData#writeHDF5}.
     *
     *                                  Modeling backends are responsible for consuming annotations in this format
     *                                  and outputting a {@link VariantAnnotationsScorer} for each variant type
     *                                  with the appropriate output names. This responsibility includes the
     *                                  implementation of functionality that allows validation of annotation names
     *                                  in downstream {@link VariantAnnotationsScorer} instances.
     *
     *                                  In current use, we assume that a single model will be trained, so either
     *                                    1) training annotations have already been subset to a single variant type (SNP or INDEL), or
     *                                    2) we assume the model does not care about the variant type.
     *                                  TODO we could also pass additional labels to be used in training,
     *                                       but all backends would have to likewise respect directory structure
     *
     * @param outputPrefix              Path prefix for all output files
     */
    void trainAndSerialize(final File trainingAnnotationsFile,
                           final String outputPrefix);
}
