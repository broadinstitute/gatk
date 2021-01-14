package org.broadinstitute.hellbender.tools.walkers.annotator;

import java.util.List;

/**
 * Interface of all variant annotations. See also InfoFieldAnnotation and GenotypeFieldAnnotation
 */
public interface VariantAnnotation extends Annotation {

    /**
     * Return the keys
     */
    List<String> getKeyNames();
}
