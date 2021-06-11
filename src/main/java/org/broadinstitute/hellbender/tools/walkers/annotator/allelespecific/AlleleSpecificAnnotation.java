package org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific;

import org.broadinstitute.hellbender.tools.walkers.annotator.Annotation;

/**
 * This is a marker interface used to indicate which annotations are allele-specific.
 */
public interface AlleleSpecificAnnotation extends Annotation {
    default String getEmptyRawValue() {
        return "0";
    }
}
