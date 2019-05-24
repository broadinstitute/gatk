package org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.segment;

/**
 * Class that represents the exon numbers overlapped by a genomic region.
 *
 * See {@link SegmentExonUtils}
 */
public class SegmentExonOverlaps {
    private final String segmentStartExonOverlap;
    private final String segmentEndExonOverlap;

    public SegmentExonOverlaps(final String segmentStartExonOverlap, final String segmentEndExonOverlap) {
        this.segmentStartExonOverlap = segmentStartExonOverlap;
        this.segmentEndExonOverlap = segmentEndExonOverlap;
    }

    public String getSegmentStartExonOverlap() {
        return segmentStartExonOverlap;
    }

    public String getSegmentEndExonOverlap() {
        return segmentEndExonOverlap;
    }
}