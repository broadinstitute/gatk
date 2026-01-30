package org.broadinstitute.hellbender.tools.sv.concordance;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.cluster.CanonicalSVLinkage;

public class SVMatchingLinkage<T extends SVCallRecord> extends CanonicalSVLinkage<T> {

    public SVMatchingLinkage(final SAMSequenceDictionary dictionary) {
        super(dictionary, false);
    }



}
