package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling;

import com.google.common.collect.ImmutableList;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.data.LabeledVariantAnnotationsData;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.clustering.BayesianGaussianMixtureModeller;

import java.io.File;
import java.io.Serializable;
import java.util.List;

public final class PythonSklearnVariantAnnotationsScorer implements VariantAnnotationsScorer, Serializable {
    private static final long serialVersionUID = 1L;

    public PythonSklearnVariantAnnotationsScorer(final File pythonScriptFile,
                                                 final File hyperparametersJSONFile) {
    }

    @Override
    public void scoreSamples(final File inputAnnotationsFile,
                             final File outputScoresFile) {
        final List<String> inputAnnotationNames = LabeledVariantAnnotationsData.readAnnotationNames(inputAnnotationsFile);

    }
}
