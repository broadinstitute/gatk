package org.broadinstitute.hellbender.utils.codecs;

import htsjdk.tribble.AsciiFeatureCodec;
import org.broadinstitute.hellbender.tools.sv.SVEvidence;

public abstract class SVEvidenceCodec<T extends SVEvidence> extends AsciiFeatureCodec<T> {
    protected SVEvidenceCodec(final Class<T> myClass) {
        super(myClass);
    }

    public abstract String encode(T ev);
}
