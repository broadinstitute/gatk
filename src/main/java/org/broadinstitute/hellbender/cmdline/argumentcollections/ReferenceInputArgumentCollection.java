package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.hellbender.engine.GATKPathSpecifier;

import java.io.Serializable;
import java.nio.file.Path;

/**
 * An abstract ArgumentCollection for specifying a reference sequence file
 */
public abstract class ReferenceInputArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    /**
     * Get the name of the reference input specified at the command line.
     */
    public abstract GATKPathSpecifier getReferenceSpecifier();

    /**
     * Get the name of the reference file specified at the command line.
     */
    public String getReferenceFileName() {
        final GATKPathSpecifier inputPath = getReferenceSpecifier();
        return inputPath == null ? null : inputPath.getRawInputString();
    }

    /**
     * Get the Path to the reference, may be null
     */
    public Path getReferencePath() {
        return getReferenceSpecifier() != null ? getReferenceSpecifier().toPath() : null;
    }
}
