package org.broadinstitute.hellbender.engine.transformers;

import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.read.GATKRead;

public class IdentityTransformer implements ReadTransformer {
    public static final long serialVersionUID = 1L;

    @Override
    public GATKRead apply(GATKRead gatkRead) {
        return gatkRead;
    }
}
