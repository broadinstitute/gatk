package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.StructuralVariantType;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
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

    /**
     * Find a single call representative of all the calls in the {@param cluster}
     * @param cluster   the events that are clustered together
     * @return  a call approximating the average event for the cluster and containing all the algorithms and genotypes
     */

    @Override
    protected SVCallRecordWithEvidence flattenCluster(final Collection<SVCallRecordWithEvidence> cluster) {
        final List<Integer> startPositions = cluster.stream().map(SVCallRecordWithEvidence::getStart).sorted().collect(Collectors.toList());
        final List<Integer> endPositions = cluster.stream().map(SVCallRecordWithEvidence::getEnd).sorted().collect(Collectors.toList());
        //use the mid value of the sorted list so the start and end represent real breakpoint observations
        final int medianStart = startPositions.get(startPositions.size() / 2);
        final int medianEnd = endPositions.get(endPositions.size() / 2);
        final SVCallRecordWithEvidence exampleCall = cluster.iterator().next();
        final int length = exampleCall.getContig().equals(exampleCall.getEndContig()) && !exampleCall.getType().equals(StructuralVariantType.INS) ? medianEnd - medianStart : exampleCall.getLength();
        final List<String> algorithms = cluster.stream().flatMap(v -> v.getAlgorithms().stream()).distinct().collect(Collectors.toList());
        final List<Genotype> clusterSamples = cluster.stream().flatMap(v -> v.getGenotypes().stream()).collect(Collectors.toList());

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
        //??? does evidence not need to be merged???
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

    /**
     * Determine an overlap interval for clustering using reciprocal overlap or breakend window, as applicable
     * Returned interval represents the interval in which the start position of a new event must fall in order to be added to the cluster (including {@param call})
     * @param call  new event to be clustered
     * @param clusterMinStartInterval    the cluster of interest, may be null
     * @return  an interval describing the cluster after {@param call} is added
     */
    @Override
    protected SimpleInterval getClusteringInterval(final SVCallRecordWithEvidence call, final SimpleInterval clusterMinStartInterval) {
        final int minStart;
        final int maxStart;
        if (isDepthOnlyCall(call)) {
            minStart = (int) (call.getEnd() - call.getLength() / MIN_RECIPROCAL_OVERLAP_DEPTH); //start of an overlapping event such that call represents (reciprocal overlap) of that event
            maxStart = (int) (call.getStart() + (1.0 - MIN_RECIPROCAL_OVERLAP_DEPTH) * call.getLength());
        } else {
            minStart = call.getStart() - MAX_BREAKEND_CLUSTERING_WINDOW;
            maxStart = call.getStart() + MAX_BREAKEND_CLUSTERING_WINDOW;
        }
        final String currentContig = getCurrentContig();
        if (clusterMinStartInterval == null) {
            return IntervalUtils.trimIntervalToContig(currentContig, minStart, maxStart, dictionary.getSequence(currentContig).getSequenceLength());
        }
        //NOTE: this is an approximation -- best method would back calculate cluster bounds, then rederive start and end based on call + cluster
        final int newMinStart = Math.min(minStart, clusterMinStartInterval.getStart());
        final int newMaxStart = Math.max(maxStart, clusterMinStartInterval.getEnd());
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

    /**
     *  Merge genotypes and algorithms for multiple calls describing the same event
     * @param items all entries are assumed to describe the same event, i.e. satisfy {@link #itemsAreIdentical}
     * @return  a single representative call
     */
    @Override
    protected SVCallRecordWithEvidence deduplicateIdenticalItems(final Collection<SVCallRecordWithEvidence> items) {
        if (items.isEmpty()) {
            return null;
        }
        final List<Genotype> genotypes = items.stream()
                .map(SVCallRecordWithEvidence::getGenotypes)
                .flatMap(Collection::stream)
                .collect(Collectors.toList());
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
                genotypes,
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
        return call.getAlgorithms().size() == 1 && call.getAlgorithms().get(0).equals(GATKSVVCFConstants.DEPTH_ALGORITHM);
    }

    public double getMinReciprocalOverlap() {
        return MIN_RECIPROCAL_OVERLAP_DEPTH;
    }
}
