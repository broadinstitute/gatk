package org.broadinstitute.hellbender.tools.exome.germlinehmm;

import org.broadinstitute.hellbender.utils.Utils;

/**
 * {@link CopyNumberTriStateSegment} file record.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class CopyNumberTriStateSegmentRecord {

    private final String sampleName;

    private final CopyNumberTriStateSegment segment;

    public CopyNumberTriStateSegmentRecord(final String sampleName, final CopyNumberTriStateSegment segment) {
        this.sampleName = Utils.nonNull(sampleName);
        this.segment = Utils.nonNull(segment);
    }

    public CopyNumberTriStateSegment getSegment() {
        return segment;
    }

    public String getSampleName() {
        return sampleName;
    }

    public String toString() {
        return String.format("{sample: %s, segment: %s}", sampleName, segment);
    }
}
