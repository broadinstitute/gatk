package org.broadinstitute.hellbender.tools.sv;

import htsjdk.tribble.Feature;

public abstract class SVEvidence implements Feature {

    // Returns record line for file output
    public abstract String encode();

}
