package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling;

public enum VariantAnnotationsModelBackend {
    /**
     * Use {@link org.broadinstitute.hellbender.utils.clustering.BayesianGaussianMixtureModeller}.
     */
    JAVA_BGMM,

    /**
     * Use the script at org/broadinstitute/hellbender/tools/walkers/vqsr/scalable/isolation-forest.py
     */
    PYTHON_IFOREST,

    /**
     * Use a user-provided script.
     */
    PYTHON_SCRIPT
}
