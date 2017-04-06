package org.broadinstitute.hellbender.cmdline.argumentcollections;

import java.io.File;
import java.io.Serializable;

/**
 * An abstract ArgumentCollection for specifying a reference sequence file
 */
public abstract class ReferenceInputArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;
    /**
     * Get the reference file specified at the command line, creating the File object first if necessary.
     */
    public abstract File getReferenceFile();

    /**
     * Get the name of the reference file specified at the command line.
     */
    public abstract String getReferenceFileName();

    /**
     * Checks whether a reference was specified.
     * @return {@code true} iff a reference was specified.
     */
    public boolean hasReference() {
        return getReferenceFile() != null;
    }
}
