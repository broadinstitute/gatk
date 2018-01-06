package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
import java.io.Serializable;
import java.nio.file.Path;

/**
 * An abstract ArgumentCollection for specifying a reference sequence file
 */
public abstract class ReferenceInputArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    /**
     * Get the name of the reference file specified at the command line.
     */
    public abstract String getReferenceFileName();

    /**
     * Get the Path to the reference, may be null
     */
    public Path getReferencePath() {
        return getReferenceFileName() != null ? IOUtils.getPath(getReferenceFileName()) : null;
    }
}
