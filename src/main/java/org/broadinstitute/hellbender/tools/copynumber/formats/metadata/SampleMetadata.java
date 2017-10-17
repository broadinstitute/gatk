package org.broadinstitute.hellbender.tools.copynumber.formats.metadata;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public interface SampleMetadata {
    /**
     * @return  a non-empty sample name
     */
    String getSampleName();
}
