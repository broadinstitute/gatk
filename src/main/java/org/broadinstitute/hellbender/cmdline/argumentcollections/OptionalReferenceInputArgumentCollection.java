package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.WorkflowResource;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKPath;

/**
 * An argument collection for use with tools that optionally accept a reference file as input.
 */
public final class OptionalReferenceInputArgumentCollection extends ReferenceInputArgumentCollection {
    private static final long serialVersionUID = 1L;

    @WorkflowResource(input=true, output=false, companionResources = {"referenceDictionary", "referenceIndex"})
    @Argument(fullName = StandardArgumentDefinitions.REFERENCE_LONG_NAME, shortName = StandardArgumentDefinitions.REFERENCE_SHORT_NAME, doc = "Reference sequence", optional = true)
    private GATKPath referenceInputPathSpecifier;

    @Override
    public GATKPath getReferenceSpecifier() {
        return referenceInputPathSpecifier;
    }
}
