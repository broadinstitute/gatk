package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.positional;

import htsjdk.samtools.CigarElement;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReferenceConfidenceModel;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;

/**
 * A read prepared for positional assembly: its bases, qualities, the reference-window offset at which each base is
 * expected (from its alignment, with soft clips extrapolated), and its runs of bases usable in kmers.
 */
final class PositionalRead {
    final String fragment;
    final byte[] bases;
    final byte[] qualities;
    /** expected window offset of each base: aligned bases at their reference offset, inserted bases half-way after the
     * preceding aligned base, and clipped bases extrapolated from the nearest aligned base */
    final double[] projection;
    /** whether each base was soft-clipped in the original alignment */
    final boolean[] clipped;
    /** [start, end) runs of bases of quality at least the minimum and without N, at least a kmer long */
    final List<int[]> segments;
    /** for each segment, the graph node of each of its kmers */
    int[][] paths;

    private PositionalRead(final String fragment, final byte[] bases, final byte[] qualities, final double[] projection,
                           final boolean[] clipped, final List<int[]> segments) {
        this.fragment = fragment;
        this.bases = bases;
        this.qualities = qualities;
        this.projection = projection;
        this.clipped = clipped;
        this.segments = segments;
    }

    int length() {
        return bases.length;
    }

    String basesString(final int start, final int end) {
        return new String(bases, start, end - start, StandardCharsets.US_ASCII);
    }

    /**
     * @param read a read whose soft clips the assembly region has already reverted (recording the original alignment in
     *             the {@link ReferenceConfidenceModel#ORIGINAL_SOFTCLIP_START_TAG} and
     *             {@link ReferenceConfidenceModel#ORIGINAL_SOFTCLIP_END_TAG} tags) or removed
     * @param windowStart 1-based reference position of the window's first base
     * @return the prepared read, or {@code null} if it has no usable kmer
     */
    static PositionalRead create(final GATKRead read, final String sample, final int windowStart, final int kmerSize,
                                 final int minBaseQuality) {
        final byte[] bases = read.getBases();
        final byte[] qualities = read.getBaseQualities();
        final int length = bases.length;
        if (length < kmerSize || qualities.length != length) {
            return null;
        }
        final Integer originalStart = read.getAttributeAsInteger(ReferenceConfidenceModel.ORIGINAL_SOFTCLIP_START_TAG);
        final Integer originalEnd = read.getAttributeAsInteger(ReferenceConfidenceModel.ORIGINAL_SOFTCLIP_END_TAG);

        final double[] projection = new double[length];
        final boolean[] clipped = new boolean[length];
        int firstAligned = -1, lastAligned = -1;
        int offset = 0, refPosition = read.getStart();
        double lastAlignedProjection = refPosition - windowStart - 1;
        for (final CigarElement element : read.getCigarElements()) {
            final int n = element.getLength();
            switch (element.getOperator()) {
                case M, EQ, X -> {
                    for (int i = 0; i < n; i++, offset++, refPosition++) {
                        projection[offset] = refPosition - windowStart;
                        clipped[offset] = originalStart != null && originalEnd != null && (refPosition < originalStart || refPosition > originalEnd);
                        lastAlignedProjection = projection[offset];
                        if (firstAligned < 0) {
                            firstAligned = offset;
                        }
                        lastAligned = offset;
                    }
                }
                case I, S -> {
                    for (int i = 0; i < n; i++, offset++) {
                        projection[offset] = lastAlignedProjection + 0.5;
                    }
                }
                case D, N -> refPosition += n;
                default -> { }
            }
        }
        if (firstAligned < 0 || offset != length) {
            return null;
        }
        // bases before the first or after the last aligned base were clipped: extrapolate their positions
        for (int i = 0; i < firstAligned; i++) {
            projection[i] = projection[firstAligned] - (firstAligned - i);
            clipped[i] = true;
        }
        for (int i = lastAligned + 1; i < length; i++) {
            projection[i] = projection[lastAligned] + (i - lastAligned);
            clipped[i] = true;
        }

        final List<int[]> segments = new ArrayList<>();
        int segmentStart = -1;
        for (int i = 0; i <= length; i++) {
            final boolean usable = i < length && qualities[i] >= minBaseQuality && BaseUtils.simpleBaseToBaseIndex(bases[i]) >= 0;
            if (usable && segmentStart < 0) {
                segmentStart = i;
            } else if (!usable && segmentStart >= 0) {
                if (i - segmentStart >= kmerSize) {
                    segments.add(new int[] {segmentStart, i});
                }
                segmentStart = -1;
            }
        }
        return segments.isEmpty() ? null : new PositionalRead(sample + '\t' + read.getName(), bases, qualities, projection, clipped, segments);
    }
}
