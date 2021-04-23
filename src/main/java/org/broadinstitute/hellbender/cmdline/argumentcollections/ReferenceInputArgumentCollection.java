package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.hellbender.engine.GATKPath;

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
    public abstract GATKPath getReferenceSpecifier();

    /**
     * Get the name of the reference file specified at the command line.
     */
    public String getReferenceFileName() {
        // We should remove this method completely once all the call sites that use it are updated
        // and the migration to GATKPath is complete, since its not not clear what format
        // the return value can have (i.e., can it have a URI scheme ? query params ?). Instead, all
        // consumers shoule either rely directly on GATKPath, or use a Path or URI, which
        // can both be obtained from a GATKPath.
        final GATKPath inputPathSpec = getReferenceSpecifier();
        return inputPathSpec == null ? null : inputPathSpec.getURI().getPath();
    }

    /**
     * Get the Path to the reference, may be null
     */
    public Path getReferencePath() {
        return getReferenceSpecifier() != null ? getReferenceSpecifier().toPath() : null;
    }
}
