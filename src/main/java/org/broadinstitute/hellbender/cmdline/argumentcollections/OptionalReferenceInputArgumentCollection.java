package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollectionDefinition;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;

import java.io.File;

/**
 * An argument collection for use with tools that optionally accept a reference file as input.
 */
public class OptionalReferenceInputArgumentCollection implements ArgumentCollectionDefinition {

    @Argument(fullName = StandardArgumentDefinitions.REFERENCE_LONG_NAME, shortName = StandardArgumentDefinitions.REFERENCE_SHORT_NAME, doc = "Reference sequence", optional = true)
    public File referenceFile;

}
