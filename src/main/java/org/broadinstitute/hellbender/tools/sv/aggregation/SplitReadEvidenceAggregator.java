package org.broadinstitute.hellbender.tools.sv.aggregation;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SplitReadEvidence;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

public class SplitReadEvidenceAggregator extends SVEvidenceAggregator<SplitReadEvidence> {

    private final int window;
    private final boolean isStart; // Retrieve start position split reads, else end position

    public SplitReadEvidenceAggregator(final FeatureDataSource<SplitReadEvidence> source,
                                       final SAMSequenceDictionary dictionary,
                                       final int window,
                                       final boolean isStart) {
        super(source, window, dictionary);
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

    public static SimpleInterval getStartEvidenceQueryInterval(final SVCallRecord call, final int window, final SAMSequenceDictionary dictionary) {
        final SimpleInterval result = call.getPositionAInterval().expandWithinContig(window, dictionary);
        Utils.nonNull(result, "Error generating padded interval for variant " + call.getId() + "; check that its coordinates are valid");
        return result;
    }

    public static SimpleInterval getEndEvidenceQueryInterval(final SVCallRecord call, final int window, final SAMSequenceDictionary dictionary) {
        final SimpleInterval result = call.getPositionBInterval().expandWithinContig(window, dictionary);
        Utils.nonNull(result, "Error generating padded interval for variant " + call.getId() + "; check that its coordinates are valid");
        return result;
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
