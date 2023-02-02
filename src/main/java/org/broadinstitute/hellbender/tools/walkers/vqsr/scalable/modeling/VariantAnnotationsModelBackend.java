package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling;

public enum VariantAnnotationsModelBackend {
    // TODO will be added in a separate PR
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
