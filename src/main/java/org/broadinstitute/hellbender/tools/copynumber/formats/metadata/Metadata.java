package org.broadinstitute.hellbender.tools.copynumber.formats.metadata;

import htsjdk.samtools.SAMFileHeader;

/**
 * Interface for marking objects that contain metadata that can be represented as a {@link SAMFileHeader}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public interface Metadata {
    enum Type {
        SAMPLE,
        LOCATABLE,
        SAMPLE_LOCATABLE
    }

    /**
     * @return  the metadata represented as a {@link SAMFileHeader}
     */
    SAMFileHeader toHeader();
}
