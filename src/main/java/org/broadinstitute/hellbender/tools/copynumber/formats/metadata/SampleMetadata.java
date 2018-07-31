package org.broadinstitute.hellbender.tools.copynumber.formats.metadata;

/**
 * Interface for marking objects that contain metadata associated with a single sample.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public interface SampleMetadata extends Metadata {
    /**
     * @return  a non-empty sample name
     */
    String getSampleName();
}
