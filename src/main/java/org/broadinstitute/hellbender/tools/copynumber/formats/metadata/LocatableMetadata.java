package org.broadinstitute.hellbender.tools.copynumber.formats.metadata;

import htsjdk.samtools.SAMSequenceDictionary;

/**
 * Interface for marking objects that contain metadata associated with a collection of locatables.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public interface LocatableMetadata extends Metadata {
    /**
     * @return  a non-null sequence dictionary
     */
    SAMSequenceDictionary getSequenceDictionary();
}
