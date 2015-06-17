package org.broadinstitute.hellbender.utils.gene;

import org.broadinstitute.hellbender.exceptions.GATKException;

/**
 * Exception thrown when loading gene annotations. It is expected that there will be inconsistencies in annotation
 * files, such that these exceptions may be reported but not cause program termination.
 */
public final class GeneAnnotationException extends GATKException {
    private static final long serialVersionUID = 1L;
    public GeneAnnotationException(String message) {
        super(message);
    }

    public GeneAnnotationException(String message, Throwable throwable) {
        super(message, throwable);
    }
}
