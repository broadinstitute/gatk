package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.StructuralVariantType;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.Collection;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

public class SVClusterEngine extends LocatableClusterEngine<SVCallRecordWithEvidence> {

    private final double MIN_RECIPROCAL_OVERLAP_DEPTH = 0.8;
    private final double BREAKEND_CLUSTERING_WINDOW_FRACTION = 0.5;
    private final int MIN_BREAKEND_CLUSTERING_WINDOW = 50;
    private final int MAX_BREAKEND_CLUSTERING_WINDOW = 300;
    private final int MIXED_CLUSTERING_WINDOW = 2000;

    public SVClusterEngine(final SAMSequenceDictionary dictionary) {
        super(dictionary, CLUSTERING_TYPE.MAX_CLIQUE);
    }

    @Override
    protected SVCallRecordWithEvidence flattenCluster(final Collection<SVCallRecordWithEvidence> cluster) {
        final List<Integer> startPositions = cluster.stream().map(SVCallRecordWithEvidence::getStart).sorted().collect(Collectors.toList());
        final List<Integer> endPositions = cluster.stream().map(SVCallRecordWithEvidence::getEnd).sorted().collect(Collectors.toList());
        final int medianStart = startPositions.get(startPositions.size() / 2);
        final int medianEnd = endPositions.get(endPositions.size() / 2);
        final SVCallRecordWithEvidence exampleCall = cluster.iterator().next();
        final int length = exampleCall.getContig().equals(exampleCall.getEndContig()) && !exampleCall.getType().equals(StructuralVariantType.INS) ? medianEnd - medianStart : exampleCall.getLength();
        final List<String> algorithms = cluster.stream().flatMap(v -> v.getAlgorithms().stream()).distinct().collect(Collectors.toList());
        final Set<String> clusterSamples = cluster.stream().flatMap(v -> v.getSamples().stream()).collect(Collectors.toSet());

        final int newStart;
        final int newEnd;
        if (exampleCall.getType().equals(StructuralVariantType.INS)) {
            // Insertions should be a single locus; also fixes case where end-supporting split reads are to the
            // left of start-supporting split reads
            final int mean = (medianStart + medianEnd) / 2;
            newStart = mean;
            newEnd = mean + 1;
        } else {
            newStart = medianStart;
            newEnd = medianEnd;
        }
        return new SVCallRecordWithEvidence(exampleCall.getContig(), newStart, exampleCall.getStartStrand(),
                exampleCall.getEndContig(), newEnd, exampleCall.getEndStrand(), exampleCall.getType(), length, algorithms, clusterSamples,
                exampleCall.getStartSplitReadSites(), exampleCall.getEndSplitReadSites(), exampleCall.getDiscordantPairs());
    }

    @Override
    protected boolean clusterTogether(final SVCallRecordWithEvidence a, final SVCallRecordWithEvidence b) {
        if (!a.getType().equals(b.getType())) return false;
        final boolean depthOnlyA = isDepthOnlyCall(a);
        final boolean depthOnlyB = isDepthOnlyCall(b);
        if (depthOnlyA && depthOnlyB) {
            return clusterTogetherBothDepthOnly(a, b);
        } else if (depthOnlyA != depthOnlyB) {
            return clusterTogetherMixedEvidence(a, b);
        } else {
            return clusterTogetherBothWithEvidence(a, b);
        }
    }

    @Override
    protected SimpleInterval getClusteringInterval(final SVCallRecordWithEvidence call, final SimpleInterval clusterMinStartInterval) {
        final int minStart;
        final int maxStart;
        if (isDepthOnlyCall(call)) {
            minStart = (int) (call.getStart() - (1.0 - MIN_RECIPROCAL_OVERLAP_DEPTH) * call.getLength());
            maxStart = (int) (call.getStart() + (1.0 - MIN_RECIPROCAL_OVERLAP_DEPTH) * call.getLength());
        } else {
            minStart = call.getStart() - MAX_BREAKEND_CLUSTERING_WINDOW;
            maxStart = call.getStart() + MAX_BREAKEND_CLUSTERING_WINDOW;
        }
        final String currentContig = getCurrentContig();
        if (clusterMinStartInterval == null) {
            return IntervalUtils.trimIntervalToContig(currentContig, minStart, maxStart, dictionary.getSequence(currentContig).getSequenceLength());
        }
        final int newMinStart = Math.max(minStart, clusterMinStartInterval.getStart());
        final int newMaxStart = Math.min(maxStart, clusterMinStartInterval.getEnd());
        return IntervalUtils.trimIntervalToContig(currentContig, newMinStart, newMaxStart, dictionary.getSequence(currentContig).getSequenceLength());
    }

    @Override
    protected boolean itemsAreIdentical(final SVCallRecordWithEvidence a, final SVCallRecordWithEvidence b) {
        return a.getContig().equals(b.getContig())
                && a.getStart() == b.getStart()
                && a.getEndContig().equals(b.getEndContig())
                && a.getEnd() == b.getEnd()
                && a.getType().equals(b.getType())
                && a.getStartStrand() == b.getStartStrand()
                && a.getEndStrand() == b.getEndStrand();
    }

    @Override
    protected SVCallRecordWithEvidence deduplicateIdenticalItems(final Collection<SVCallRecordWithEvidence> items) {
        if (items.isEmpty()) return null;
        final Set<String> samples = items.stream()
                .map(SVCallRecordWithEvidence::getSamples)
                .flatMap(Collection::stream)
                .collect(Collectors.toSet());
        final List<String> algorithms = items.stream()
                .map(SVCallRecordWithEvidence::getAlgorithms)
                .flatMap(Collection::stream)
                .distinct()
                .collect(Collectors.toList());
        final SVCallRecordWithEvidence example = items.iterator().next();
        return new SVCallRecordWithEvidence(
                example.getContig(),
                example.getStart(),
                example.getStartStrand(),
                example.getEndContig(),
                example.getEnd(),
                example.getEndStrand(),
                example.getType(),
                example.getLength(),
                algorithms,
                samples,
                example.getStartSplitReadSites(),
                example.getEndSplitReadSites(),
                example.getDiscordantPairs());
    }

    private boolean clusterTogetherBothDepthOnly(final SVCallRecordWithEvidence a, final SVCallRecordWithEvidence b) {
        if (!a.getContig().equals(a.getEndContig()) || !b.getContig().equals(b.getEndContig())) {
            throw new IllegalArgumentException("Attempted to cluster depth-only calls with endpoints on different contigs");
        }
        final SimpleInterval intervalA = new SimpleInterval(a.getContig(), a.getStart(), a.getEnd());
        final SimpleInterval intervalB = new SimpleInterval(b.getContig(), b.getStart(), b.getEnd());
        return IntervalUtils.isReciprocalOverlap(intervalA, intervalB, MIN_RECIPROCAL_OVERLAP_DEPTH);
    }

    private boolean clusterTogetherBothWithEvidence(final SVCallRecordWithEvidence a, final SVCallRecordWithEvidence b) {
        // Reject if one is intrachromosomal and the other isn't
        final boolean intrachromosomalA = a.getContig().equals(a.getEndContig());
        final boolean intrachromosomalB = b.getContig().equals(b.getEndContig());
        if (intrachromosomalA != intrachromosomalB) return false;

        // Matching endpoints
        final SimpleInterval intervalAStart =  getStartClusteringInterval(a);
        final SimpleInterval intervalAEnd =  getEndClusteringInterval(a);
        final SimpleInterval intervalBStart =  getStartClusteringInterval(b);
        final SimpleInterval intervalBEnd =  getEndClusteringInterval(b);
        return intervalAStart.overlaps(intervalBStart) && intervalAEnd.overlaps(intervalBEnd);
    }

    private boolean clusterTogetherMixedEvidence(final SVCallRecordWithEvidence a, final SVCallRecordWithEvidence b) {
        final boolean intrachromosomalA = a.getContig().equals(a.getEndContig());
        final boolean intrachromosomalB = b.getContig().equals(b.getEndContig());
        if (!(intrachromosomalA && intrachromosomalB)) return false;
        if (!a.getContig().equals(b.getContig())) return false;
        return Math.abs(a.getStart() - b.getStart()) < MIXED_CLUSTERING_WINDOW
                && Math.abs(a.getEnd() - b.getEnd()) < MIXED_CLUSTERING_WINDOW;
    }

    private SimpleInterval getStartClusteringInterval(final SVCallRecordWithEvidence call) {
        final int padding = getEndpointClusteringPadding(call);
        return call.getStartAsInterval().expandWithinContig(padding, dictionary);
    }

    private SimpleInterval getEndClusteringInterval(final SVCallRecordWithEvidence call) {
        final int padding =  getEndpointClusteringPadding(call);
        return call.getEndAsInterval().expandWithinContig(padding, dictionary);
    }

    private int getEndpointClusteringPadding(final SVCallRecordWithEvidence call) {
        return (int) Math.min(MAX_BREAKEND_CLUSTERING_WINDOW, Math.max(MIN_BREAKEND_CLUSTERING_WINDOW, BREAKEND_CLUSTERING_WINDOW_FRACTION * call.getLength()));
    }

    public static boolean isDepthOnlyCall(final SVCallRecordWithEvidence call) {
        return call.getAlgorithms().size() == 1 && call.getAlgorithms().get(0).equals(SVCluster.DEPTH_ALGORITHM);
    }
}
