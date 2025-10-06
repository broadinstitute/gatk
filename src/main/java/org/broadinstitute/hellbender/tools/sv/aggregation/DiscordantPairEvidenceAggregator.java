package org.broadinstitute.hellbender.tools.sv.aggregation;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.tools.sv.DiscordantPairEvidence;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;

public class DiscordantPairEvidenceAggregator extends SVEvidenceAggregator<DiscordantPairEvidence> {

    // PE evidence uses asymmetric windows for "outside" and "inside" the breakpoint
    private final int innerWindow;
    private final int outerWindow;

    /**
     * Constructor.
     * @param source        PE feature source
     * @param dictionary    sequence dictionary
     * @param innerWindow   window past start locus and before end locus
     * @param outerWindow   window before start locus and past end locus
     */
    public DiscordantPairEvidenceAggregator(final FeatureDataSource<DiscordantPairEvidence> source,
                                            final SAMSequenceDictionary dictionary,
                                            final int innerWindow,
                                            final int outerWindow) {
        super(source, dictionary);
        this.innerWindow = innerWindow;
        this.outerWindow = outerWindow;
    }

    @Override
    public SimpleInterval getEvidenceQueryInterval(final SVCallRecord record) {
        return getDiscordantPairStartInterval(record);
    }

    public int getInnerWindow() {
        return innerWindow;
    }

    public int getOuterWindow() {
        return outerWindow;
    }

    @Override
    public boolean evidenceFilter(final SVCallRecord call, final DiscordantPairEvidence evidence) {
        final SimpleInterval startInterval = getDiscordantPairStartInterval(call);
        final SimpleInterval endInterval = getDiscordantPairEndInterval(call);
        return discordantPairOverlapsInterval(evidence, startInterval, endInterval)
                && evidence.getStartStrand() == call.getStrandA()
                && evidence.getEndStrand() == call.getStrandB();
    }

    protected SimpleInterval getDiscordantPairStartInterval(final SVCallRecord call) {
        final String contig = call.getContigA();
        if (call.getStrandA() != null && call.getStrandA()) {
            return IntervalUtils.trimIntervalToContig(contig, call.getPositionA() - outerWindow, call.getPositionA() + innerWindow, dictionary.getSequence(contig).getSequenceLength());
        } else {
            return IntervalUtils.trimIntervalToContig(contig, call.getPositionA() - innerWindow, call.getPositionA() + outerWindow, dictionary.getSequence(contig).getSequenceLength());
        }
    }

    protected SimpleInterval getDiscordantPairEndInterval(final SVCallRecord call) {
        final String contig = call.getContigB();
        if (call.getStrandB() != null && call.getStrandB()) {
            return IntervalUtils.trimIntervalToContig(contig, call.getPositionB() - outerWindow, call.getPositionB() + innerWindow, dictionary.getSequence(contig).getSequenceLength());
        } else {
            return IntervalUtils.trimIntervalToContig(contig, call.getPositionB() - innerWindow, call.getPositionB() + outerWindow, dictionary.getSequence(contig).getSequenceLength());
        }
    }

    private boolean discordantPairOverlapsInterval(final DiscordantPairEvidence evidence,
                                                   final SimpleInterval startInterval,
                                                   final SimpleInterval endInterval) {
        return evidence.getContig().equals(startInterval.getContig())
                && evidence.getStart() >= startInterval.getStart()
                && evidence.getStart() <= startInterval.getEnd()
                && evidence.getEndContig().equals(endInterval.getContig())
                && evidence.getEndPosition() >= endInterval.getStart()
                && evidence.getEndPosition() <= endInterval.getEnd();
    }
}
