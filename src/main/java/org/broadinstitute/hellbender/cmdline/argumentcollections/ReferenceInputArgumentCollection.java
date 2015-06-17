package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.hellbender.cmdline.ArgumentCollectionDefinition;

import java.io.File;

/**
 * An abstract ArgumentCollection for specifying a reference sequence file
 */
public abstract class ReferenceInputArgumentCollection implements ArgumentCollectionDefinition {
    private static final long serialVersionUID = 1L;
    /**
     * Get the reference file specified at the command line, creating the File object first if necessary.
     */
    abstract public File getReferenceFile();

    /**
     * Get the name of the reference file specified at the command line.
     */
    abstract public String getReferenceFileName();
}
