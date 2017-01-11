package org.broadinstitute.hellbender.utils.hmm.segmentation;

import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.utils.Utils;

/**
 * Represents a {@link HiddenStateSegment} file record.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class HiddenStateSegmentRecord<S, T extends Target> {

    private final String sampleName;

    private final HiddenStateSegment<S, T> segment;

    public HiddenStateSegmentRecord(final String sampleName, final HiddenStateSegment<S, T> segment) {
        this.sampleName = Utils.nonNull(sampleName);
        this.segment = Utils.nonNull(segment);
    }

    public HiddenStateSegment<S, T> getSegment() {
        return segment;
    }

    public String getSampleName() {
        return sampleName;
    }

    public String toString() {
        return String.format("{sample: %s, segment: %s}", sampleName, segment);
    }
}
