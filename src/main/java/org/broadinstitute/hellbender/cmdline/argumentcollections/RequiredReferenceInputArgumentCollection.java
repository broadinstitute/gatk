package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;

import java.io.File;

/**
 * An argument collection for use with tools that require a reference file as input.
 */
public final class RequiredReferenceInputArgumentCollection extends ReferenceInputArgumentCollection {
    private static final long serialVersionUID = 1L;

    @Argument(fullName = StandardArgumentDefinitions.REFERENCE_LONG_NAME, shortName = StandardArgumentDefinitions.REFERENCE_SHORT_NAME, doc = "Reference sequence file", optional = false)
    private String referenceFileName;

    @Override
    public String getReferenceFileName() {
        return referenceFileName;
    }
}
