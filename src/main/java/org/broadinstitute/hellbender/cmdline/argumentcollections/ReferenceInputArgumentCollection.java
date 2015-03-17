package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.hellbender.cmdline.ArgumentCollectionDefinition;

import java.io.File;

/**
 * An abstract ArgumentCollection for specifying a reference sequence file
 */
public abstract class ReferenceInputArgumentCollection implements ArgumentCollectionDefinition {
    /**
     * Get the reference file specified at the command line.
     */
    abstract public File getReferenceFile();
}
