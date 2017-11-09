package org.broadinstitute.hellbender.tools.copynumber.formats.metadata;

import org.broadinstitute.hellbender.utils.Utils;

import java.io.Serializable;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public class SimpleSampleMetadata implements SampleMetadata, Serializable {

    private static final long serialVersionUID = 0L;

    private final String sampleName;

    public SimpleSampleMetadata(final String sampleName) {
        this.sampleName = Utils.nonEmpty(sampleName);
    }

    @Override
    public String getSampleName() {
        return sampleName;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }

        final SimpleSampleMetadata that = (SimpleSampleMetadata) o;
        return sampleName.equals(that.sampleName);
    }

    @Override
    public int hashCode() {
        return sampleName.hashCode();
    }
}
