package org.broadinstitute.hellbender.tools.sv.aggregation;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.Strand;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SplitReadEvidence;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

public class SplitReadEvidenceAggregator extends SVEvidenceAggregator<SplitReadEvidence> {

    private final int window;
    private final boolean isStart; // Retrieve start position split reads, else end position

    /**
     * Constructor.
     * @param source        SR feature source
     * @param dictionary    sequence dictionary
     * @param window        evidence window around variant start/end (in bp)
     * @param isStart       true if aggregating start coordinate, else end coordinate
     */
    public SplitReadEvidenceAggregator(final FeatureDataSource<SplitReadEvidence> source,
                                       final SAMSequenceDictionary dictionary,
                                       final int window,
                                       final boolean isStart) {
        super(source, dictionary);
        Utils.validateArg(window >= 0, "Window cannot be negative");
        this.window = window;
        this.isStart = isStart;
    }

    public int getWindow() {
        return window;
    }

    @Override
    public SimpleInterval getEvidenceQueryInterval(final SVCallRecord call) {
        return isStart ? getStartEvidenceQueryInterval(call, window, dictionary) : getEndEvidenceQueryInterval(call, window, dictionary);
    }

    private static int getPositionMaybeAttribute(final int pos, final Object attribute) {
        if (attribute != null && !attribute.getClass().isAssignableFrom(String.class)) {
            throw new IllegalArgumentException("Expected integer type for variant attributes");
        }
        return attribute == null ? pos : Integer.parseInt((String) attribute);
    }

    public static SimpleInterval getStartEvidenceQueryInterval(final SVCallRecord call, final int window, final SAMSequenceDictionary dictionary) {
        final Object splitReadPos = call.getAttributes().get(GATKSVVCFConstants.FIRST_SPLIT_POSITION_ATTRIBUTE);
        final int pos = getPositionMaybeAttribute(call.getPositionA(), splitReadPos);
        final SimpleInterval result = shiftInterval(new SimpleInterval(call.getContigA(), pos, pos), 1).expandWithinContig(window, dictionary);
        Utils.nonNull(result, "Error generating padded interval for variant " + call.getId() + "; check that its coordinates are valid");
        return result;
    }

    public static SimpleInterval getEndEvidenceQueryInterval(final SVCallRecord call, final int window, final SAMSequenceDictionary dictionary) {
        final Object splitReadPos = call.getAttributes().get(GATKSVVCFConstants.SECOND_SPLIT_POSITION_ATTRIBUTE);
        final int pos = getPositionMaybeAttribute(call.getPositionB(), splitReadPos);
        final SimpleInterval result = shiftInterval(new SimpleInterval(call.getContigB(), pos, pos), 1).expandWithinContig(window, dictionary);
        Utils.nonNull(result, "Error generating padded interval for variant " + call.getId() + "; check that its coordinates are valid");
        return result;
    }

    private static SimpleInterval shiftInterval(final SimpleInterval interval, final int shiftBy) {
        return new SimpleInterval(interval.getContig(), interval.getStart() + shiftBy, interval.getEnd() + shiftBy);
    }

    @Override
    public boolean evidenceFilter(final SVCallRecord record, final SplitReadEvidence evidence) {
        Utils.validateArg(record.getStrandA() != null, "Attempted split read evidence filtering on " +
                "variant with null first strand");
        Utils.validateArg(record.getStrandB() != null, "Attempted split read evidence filtering on " +
                "variant with null second strand");
        if (isStart) {
            return evidence.getStrand() == record.getStrandA();
        } else {
            return evidence.getStrand() == record.getStrandB();
        }
    }
}
