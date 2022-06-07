package org.broadinstitute.hellbender.tools.sv.cluster;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.utils.*;

import java.util.List;
import java.util.Map;
import java.util.TreeMap;

/**
 * CNV defragmenter for when the intervals used for coverage collection are available. The intervals are used to adjust
 * variant padding by rounding to the nearest bin boundary. This is particularly important for exomes, which typically
 * have a wide range of bin sizes.
 */
public class BinnedCNVLinkage<T extends SVCallRecord> extends CNVLinkage<T> {

    protected final TreeMap<GenomeLoc, Integer> genomicCoordinatesToBinIndexMap;
    protected final List<GenomeLoc> coverageIntervals;
    final GenomeLocParser parser;

    public BinnedCNVLinkage(final SAMSequenceDictionary dictionary, final double paddingFraction,
                            final double minSampleOverlap, final List<GenomeLoc> coverageIntervals) {
        super(dictionary, paddingFraction, minSampleOverlap);
        this.coverageIntervals = Utils.nonNull(coverageIntervals);
        genomicCoordinatesToBinIndexMap = new TreeMap<>();
        for (int i = 0; i < coverageIntervals.size(); i++) {
            genomicCoordinatesToBinIndexMap.put(coverageIntervals.get(i),i);
        }
        parser = new GenomeLocParser(this.dictionary);
    }

    // TODO: try to calculate this by sweeping through the bin intervals
    @Override
    public int getMaxClusterableStartingPosition(final T call) {
        Utils.nonNull(call);
        return dictionary.getSequence(call.getContigA()).getSequenceLength();
    }

    @Override
    protected SimpleInterval getPaddedRecordInterval(final String contig, final int start, final int end) {
        Utils.nonNull(contig);
        final GenomeLoc callStart = parser.createGenomeLoc(contig, start, start);
        final GenomeLoc callEnd = parser.createGenomeLoc(contig, end, end);

        //first interval that is equal to or "greater than" the call start, such that the start of the bin should match the call start, with a little wiggle room
        final Map.Entry<GenomeLoc, Integer> startBin = genomicCoordinatesToBinIndexMap.ceilingEntry(callStart);
        Utils.validate(startBin != null, "Call start " + callStart + " for  call at " + contig + ":" + start + "-" + end + " not found in model call intervals.");
        final int callStartIndex = startBin.getValue();

        //last interval that is equal to or "less than" the call end, such that the end of the bin should match the call end
        final Map.Entry<GenomeLoc, Integer> endBin = genomicCoordinatesToBinIndexMap.floorEntry(callEnd);
        Utils.validate(endBin != null, "Call end " + callEnd + " for call at " + contig + ":" + start + "-" + end + " not found in model call intervals.");
        final int callEndIndex = endBin.getValue();
        final int callLengthInBins = callEndIndex - callStartIndex + 1;
        Utils.validate(callLengthInBins > 0, "Copy number call at " + contig + ":" + start + "-"
                + end + " does not align with supplied model calling intervals. Use the filtered intervals input " +
                "from GermlineCNVCaller for this cohort/model.");

        final int paddedStartIndex = Math.max(callStartIndex - (int)Math.round(callLengthInBins * paddingFraction), 0);
        final int paddedCallStart;
        if (coverageIntervals.get(paddedStartIndex).getContig().equals(callStart.getContig())) {
            paddedCallStart = coverageIntervals.get(paddedStartIndex).getStart();
        } else {
            paddedCallStart = callStart.getStart();
        }

        final int paddedEndIndex = Math.min(callEndIndex + (int)Math.round(callLengthInBins * paddingFraction), genomicCoordinatesToBinIndexMap.size() - 1);
        final int paddedCallEnd;
        if (coverageIntervals.get(paddedEndIndex).getContig().equals(callEnd.getContig())) {
            paddedCallEnd = coverageIntervals.get(paddedEndIndex).getEnd();
        } else {
            paddedCallEnd = callEnd.getEnd();
        }

        final int contigLength = dictionary.getSequence(contig).getSequenceLength();
        return IntervalUtils.trimIntervalToContig(contig, paddedCallStart, paddedCallEnd, contigLength);
    }

}
