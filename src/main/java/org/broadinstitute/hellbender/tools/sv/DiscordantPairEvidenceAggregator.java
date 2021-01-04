package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Collections;
import java.util.List;

public class DiscordantPairEvidenceAggregator extends CachingSVEvidenceAggregator<DiscordantPairEvidence> {

    private final int innerWindow;
    private final int outerWindow;

    public DiscordantPairEvidenceAggregator(final FeatureDataSource<DiscordantPairEvidence> source,
                                            final SAMSequenceDictionary dictionary,
                                            final int innerWindow,
                                            final int outerWindow) {
        super(source, dictionary, "DiscordantPairSites");
        this.innerWindow = innerWindow;
        this.outerWindow = outerWindow;
    }

    @Override
    protected SimpleInterval getEvidenceQueryInterval(final SVCallRecordWithEvidence record) {
        return getDiscordantPairStartInterval(record);
    }

    @Override
    protected SVCallRecordWithEvidence assignEvidence(final SVCallRecordWithEvidence call, final List<DiscordantPairEvidence> evidence) {
        Utils.nonNull(call);
        final SVCallRecordWithEvidence callWithEvidence;
        if (call.isDepthOnly()) {
            callWithEvidence = new SVCallRecordWithEvidence(call, call.getStartSplitReadSites(), call.getEndSplitReadSites(), Collections.emptyList(), call.getCopyNumberDistribution());
        } else {
            callWithEvidence = new SVCallRecordWithEvidence(call, call.getStartSplitReadSites(), call.getEndSplitReadSites(), evidence, call.getCopyNumberDistribution());
        }
        return callWithEvidence;
    }

    public int getInnerWindow() {
        return innerWindow;
    }

    public int getOuterWindow() {
        return outerWindow;
    }

    @Override
    protected boolean evidenceFilter(final SVCallRecord record, final DiscordantPairEvidence evidence) {
        final SimpleInterval startInterval = getDiscordantPairStartInterval(record);
        final SimpleInterval endInterval = getDiscordantPairEndInterval(record);
        return discordantPairOverlapsInterval(evidence, startInterval, endInterval)
                && evidence.getStartStrand() == record.getStrandA()
                && evidence.getEndStrand() == record.getStrandB();
    }

    private SimpleInterval getDiscordantPairStartInterval(final SVCallRecord call) {
        final String contig = call.getContigA();
        if (call.getStrandA()) {
            return IntervalUtils.trimIntervalToContig(contig, call.getPositionA() - outerWindow, call.getPositionA() + innerWindow, dictionary.getSequence(contig).getSequenceLength());
        } else {
            return IntervalUtils.trimIntervalToContig(contig, call.getPositionA() - innerWindow, call.getPositionA() + outerWindow, dictionary.getSequence(contig).getSequenceLength());
        }
    }

    private SimpleInterval getDiscordantPairEndInterval(final SVCallRecord call) {
        final String contig = call.getContigB();
        if (call.getStrandB()) {
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
                && evidence.getStart() < startInterval.getEnd()
                && evidence.getEndContig().equals(endInterval.getContig())
                && evidence.getEnd() >= endInterval.getStart()
                && evidence.getEnd() < endInterval.getEnd();
    }
}
