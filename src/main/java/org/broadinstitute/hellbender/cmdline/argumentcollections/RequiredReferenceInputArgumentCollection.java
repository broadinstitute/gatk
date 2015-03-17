package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;

import java.io.File;

/**
 * An argument collection for use with tools that require a reference file as input.
 */
public class RequiredReferenceInputArgumentCollection extends ReferenceInputArgumentCollection {

    @Argument(fullName = StandardArgumentDefinitions.REFERENCE_LONG_NAME, shortName = StandardArgumentDefinitions.REFERENCE_SHORT_NAME, doc = "Reference sequence file", optional = false)
    private File referenceFile;

    public File getReferenceFile(){ return referenceFile;}
}
