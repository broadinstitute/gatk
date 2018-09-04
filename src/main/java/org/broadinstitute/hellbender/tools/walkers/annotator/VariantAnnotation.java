package org.broadinstitute.hellbender.tools.walkers.annotator;

import java.util.List;

/**
 * Superclass of all variant annotations.
 */
public abstract class VariantAnnotation implements Annotation{

    /**
     * Return the keys
     */
    public abstract List<String> getKeyNames();

    @Override
    public String toString() {
        return getClass().getSimpleName();
    }
}
