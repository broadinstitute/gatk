package org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.segment;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.annotation.Strand;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.codecs.gencode.GencodeGtfExonFeature;
import org.broadinstitute.hellbender.utils.codecs.gencode.GencodeGtfTranscriptFeature;

import java.util.List;

public class SegmentExonUtils {

    private static final String AND_BELOW_STR = "-";
    private static final String AND_ABOVE_STR = "+";

    private SegmentExonUtils(){}

    // This should always be -1
    private final static int NO_EXON_OVERLAP = -1;

    /**
     * @param transcript Never {@code null}
     * @param segment Genomic region to determine which exons are covered in the transcript.  Never {@code null}
     * @return Instance of {@link SegmentExonOverlaps}.  A typical value will be a number appended to a "-" or a "+".
     *   The exon numbers are 0-based.  "-" if it covers upstream exons (considering coding direction -- i.e. previous exons
     *   are higher genomic coordinates on negative strand transcripts) and "+" for downstream exons.
     *   For example:
     *   - "6+" would be sixth exon and above.  This would indicate downstream exons (higher genomic
     *     coordinates on positive coding directions and lower genomic coordinates when negative coding) are covered by
     *     the segment.
     *   - "" indicates that the endpoints of the segment do not cover any transcript exons.  So either the entire
     *      transcript is covered or none of the transcript is covered.
     *
     *   Note that this will return the first exon covered.  So, this will yield a non-blank value when a segment
     *    breakpoint is in an intron.
     */
    public static SegmentExonOverlaps determineSegmentExonPosition(final GencodeGtfTranscriptFeature transcript, final Locatable segment) {

        Utils.nonNull(transcript);
        Utils.nonNull(segment);

        final String NOT_IN_TRANSCRIPT = "";

        // Internally, this method assumes that the exons are sorted in genomic order.  Which is NOT how these are
        //  stored in the GencodeGtfTranscriptFeature
        final List<GencodeGtfExonFeature> exons = transcript.getGenomicStrand() == Strand.NEGATIVE ?
                Lists.reverse(transcript.getExons()) : transcript.getExons();

        if (exons.size() == 0) {
            return new SegmentExonOverlaps(NOT_IN_TRANSCRIPT, NOT_IN_TRANSCRIPT);
        }

        final boolean[] isExonOverlappedMask = createExonOverlapMask(segment, exons);
        final SimpleInterval segmentStart = new SimpleInterval(segment.getContig(), segment.getStart(), segment.getStart());
        final SimpleInterval segmentEnd = new SimpleInterval(segment.getContig(), segment.getEnd(), segment.getEnd());

        // Find the proper index for the start.
        // If the start of the segment does not overlap the first exon, but the segment does, then the start index is -1 (no overlap)
        final int inclusiveIndexPositiveDirectionStart = findInclusiveIndexPositiveDirectionStart(transcript, exons, isExonOverlappedMask, segmentStart);

        // If the end of the segment does not overlap the last exon, but the segment does, then the end index is -1 (no overlap)
        final int inclusiveIndexPositiveDirectionEnd = findInclusiveIndexPositiveDirectionEnd(transcript, exons, isExonOverlappedMask, segmentEnd);

        // Construct the final strings
        final String startResult = inclusiveIndexPositiveDirectionStart != NO_EXON_OVERLAP ?
                inclusiveIndexPositiveDirectionStart + determineSegmentOverlapDirection(transcript.getGenomicStrand(), true)
                : NOT_IN_TRANSCRIPT;
        final String endResult = inclusiveIndexPositiveDirectionEnd != NO_EXON_OVERLAP ?
                inclusiveIndexPositiveDirectionEnd + determineSegmentOverlapDirection(transcript.getGenomicStrand(), false)
                : NOT_IN_TRANSCRIPT;

        return new SegmentExonOverlaps(startResult, endResult);
    }

    private static int findInclusiveIndexPositiveDirectionEnd(final GencodeGtfTranscriptFeature transcript, final List<GencodeGtfExonFeature> exons, final boolean[] isExonOverlappedMask, final SimpleInterval segmentEnd) {
        int inclusiveIndexPositiveDirectionEnd = NO_EXON_OVERLAP;
        final int lastExonIndex = exons.size() - 1;
        if (isExonOverlappedMask[lastExonIndex] && IntervalUtils.overlaps(segmentEnd, exons.get(lastExonIndex))) {
            inclusiveIndexPositiveDirectionEnd = convertInclusiveIndexForNegativeStrand(transcript, lastExonIndex);
        } else if (!isExonOverlappedMask[lastExonIndex]) {

            // Find the first exon that has overlap with the segment, going backwards in genomic space
            for (int i = lastExonIndex - 1; i >= 0; i --) {
                if (isExonOverlappedMask[i]) {
                    inclusiveIndexPositiveDirectionEnd =  convertInclusiveIndexForNegativeStrand(transcript, i);
                    break;
                }
            }
        } else {
            inclusiveIndexPositiveDirectionEnd = NO_EXON_OVERLAP;
        }
        return inclusiveIndexPositiveDirectionEnd;
    }

    private static int findInclusiveIndexPositiveDirectionStart(final GencodeGtfTranscriptFeature transcript, final List<GencodeGtfExonFeature> exons, final boolean[] isExonOverlappedMask, final SimpleInterval segmentStart) {
        int inclusiveIndexPositiveDirectionStart = NO_EXON_OVERLAP;
        if (isExonOverlappedMask[0] && IntervalUtils.overlaps(segmentStart, exons.get(0))) {
            inclusiveIndexPositiveDirectionStart = convertInclusiveIndexForNegativeStrand(transcript,0);
        } else if (!isExonOverlappedMask[0]) {
            // Find the first exon that has overlap with the segment
            for (int i = 1; i < exons.size(); i ++) {
                if (isExonOverlappedMask[i]) {
                    inclusiveIndexPositiveDirectionStart =  convertInclusiveIndexForNegativeStrand(transcript, i);
                    break;
                }
            }
        } else {
            inclusiveIndexPositiveDirectionStart = NO_EXON_OVERLAP;
        }
        return inclusiveIndexPositiveDirectionStart;
    }

    private static boolean[] createExonOverlapMask(final Locatable segment, final List<GencodeGtfExonFeature> exons) {
        final boolean[] isExonOverlappedMask = new boolean[exons.size()];

        for (int i = 0; i < exons.size(); i++) {
            final GencodeGtfExonFeature exon = exons.get(i);

            if (IntervalUtils.overlaps(segment, exon.getGenomicPosition())) {
                isExonOverlappedMask[i] = true;
            }
        }
        return isExonOverlappedMask;
    }

    @VisibleForTesting
    static String determineSegmentOverlapDirection(final Strand strand, final boolean isSegmentStart) {
        if (isSegmentStart ^ (strand == Strand.POSITIVE)) {
            return AND_BELOW_STR;
        } else {
            return AND_ABOVE_STR;
        }
    }

    @VisibleForTesting
    static int convertInclusiveIndexForNegativeStrand(final GencodeGtfTranscriptFeature transcript, final int initialInclusiveIndexPositiveDirection) {
        final List<GencodeGtfExonFeature> exons = transcript.getExons();

        if ((transcript.getGenomicStrand() == Strand.NEGATIVE) && (initialInclusiveIndexPositiveDirection != NO_EXON_OVERLAP)){
            return exons.size() - initialInclusiveIndexPositiveDirection - 1;
        }
        return initialInclusiveIndexPositiveDirection;
    }
}
