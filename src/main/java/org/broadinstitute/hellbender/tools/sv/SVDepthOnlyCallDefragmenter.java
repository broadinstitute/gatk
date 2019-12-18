package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

public class SVDepthOnlyCallDefragmenter extends LocatableClusterEngine<SVCallRecordWithEvidence> {

    private final double MIN_SAMPLE_OVERLAP = 0.9;
    private final double PADDING_FRACTION = 0.2;

    public SVDepthOnlyCallDefragmenter(final SAMSequenceDictionary dictionary) {
        super(dictionary, CLUSTERING_TYPE.SINGLE_LINKAGE);
    }

    @Override
    protected SVCallRecordWithEvidence flattenCluster(final Collection<SVCallRecordWithEvidence> cluster) {
        final int newStart =  cluster.stream().mapToInt(SVCallRecordWithEvidence::getStart).min().getAsInt();
        final int newEnd = cluster.stream().mapToInt(SVCallRecordWithEvidence::getEnd).max().getAsInt();
        final SVCallRecordWithEvidence exampleCall = cluster.iterator().next();
        final int length = newEnd - newStart;
        final List<String> algorithms = cluster.stream().flatMap(v -> v.getAlgorithms().stream()).distinct().collect(Collectors.toList());
        final Set<String> clusterSamples = cluster.stream().flatMap(v -> v.getSamples().stream()).collect(Collectors.toSet());
        return new SVCallRecordWithEvidence(exampleCall.getContig(), newStart, exampleCall.getStartStrand(),
                exampleCall.getEndContig(), newEnd, exampleCall.getEndStrand(), exampleCall.getType(), length, algorithms, clusterSamples,
                exampleCall.getStartSplitReadSites(), exampleCall.getEndSplitReadSites(), exampleCall.getDiscordantPairs());
    }

    @Override
    protected boolean clusterTogether(final SVCallRecordWithEvidence a, final SVCallRecordWithEvidence b) {
        if (!isDepthOnlyCall(a) || !isDepthOnlyCall(b)) return false;
        Utils.validate(a.getContig().equals(a.getEndContig()), "Call A is depth-only but interchromosomal");
        Utils.validate(b.getContig().equals(b.getEndContig()), "Call B is depth-only but interchromosomal");
        if (!a.getType().equals(b.getType())) return false;
        final Set<String> sharedSamples = new HashSet<>(a.getSamples());
        sharedSamples.retainAll(b.getSamples());
        final double sampleOverlap = Math.min(sharedSamples.size() / (double) a.getSamples().size(), sharedSamples.size() / (double) b.getSamples().size());
        if (sampleOverlap < MIN_SAMPLE_OVERLAP) return false;
        return getClusteringInterval(a, null)
                .overlaps(getClusteringInterval(b, null));
    }

    @Override
    protected SimpleInterval getClusteringInterval(final SVCallRecordWithEvidence call, final SimpleInterval currentClusterInterval) {
        Utils.nonNull(call);
        final SimpleInterval callInterval = getCallInterval(call);
        final int paddedCallStart = (int) (callInterval.getStart() - PADDING_FRACTION * callInterval.getLengthOnReference());
        final int paddedCallEnd = (int) (callInterval.getEnd() + PADDING_FRACTION * callInterval.getLengthOnReference());
        final String currentContig = getCurrentContig();
        final int contigLength = dictionary.getSequence(currentContig).getSequenceLength();
        if (currentClusterInterval == null) {
            return IntervalUtils.trimIntervalToContig(currentContig, paddedCallStart, paddedCallEnd, contigLength);
        }
        final int newMinStart = Math.min(paddedCallStart, currentClusterInterval.getStart());
        final int newMaxEnd = Math.max(paddedCallEnd, currentClusterInterval.getEnd());
        return IntervalUtils.trimIntervalToContig(currentContig, newMinStart, newMaxEnd, contigLength);
    }

    // Not used for single-linkage clustering
    @Override
    protected boolean itemsAreIdentical(final SVCallRecordWithEvidence a, final SVCallRecordWithEvidence b) {
        return false;
    }

    // Not used for single-linkage clustering
    @Override
    protected SVCallRecordWithEvidence deduplicateIdenticalItems(final Collection<SVCallRecordWithEvidence> items) {
        return null;
    }

    private SimpleInterval getCallInterval(final SVCallRecordWithEvidence call) {
        return new SimpleInterval(call.getContig(), call.getStart(), call.getEnd());
    }

    public static boolean isDepthOnlyCall(final SVCallRecord call) {
        if (call.getAlgorithms().isEmpty()) return false;
        for (final String alg : call.getAlgorithms()) {
            if (!alg.equals(SVCluster.DEPTH_ALGORITHM)) return false;
        }
        return true;
    }
}
