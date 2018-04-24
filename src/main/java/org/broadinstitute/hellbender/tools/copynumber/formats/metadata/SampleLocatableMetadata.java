package org.broadinstitute.hellbender.tools.copynumber.formats.metadata;

import htsjdk.samtools.SAMSequenceDictionary;

/**
 * Interface for marking objects that contain metadata associated with a collection of locatables
 * associated with a single sample.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public interface SampleLocatableMetadata extends LocatableMetadata {
    /**
     * @return  a non-empty sample name
     */
    String getSampleName();

    /**
     * @return  a non-null sequence dictionary
     */
    SAMSequenceDictionary getSequenceDictionary();
}
