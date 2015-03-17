package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;

import java.io.File;

/**
 * An argument collection for use with tools that optionally accept a reference file as input.
 */
public final class OptionalReferenceInputArgumentCollection extends ReferenceInputArgumentCollection {

    @Argument(fullName = StandardArgumentDefinitions.REFERENCE_LONG_NAME, shortName = StandardArgumentDefinitions.REFERENCE_SHORT_NAME, doc = "Reference sequence", optional = true)
    private File referenceFile;

    @Override
    public File getReferenceFile() {
        return referenceFile;
    }
}
