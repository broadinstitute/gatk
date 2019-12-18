package org.broadinstitute.hellbender.tools.sv;

import com.google.common.primitives.Doubles;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.StructuralVariantType;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.Collection;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;

public class SVClusterEngine extends LocatableClusterEngine<SVCallRecord> {

    protected static final double MIN_RECIPROCAL_OVERLAP_DEPTH = 0.8;
    protected final double BREAKEND_CLUSTERING_WINDOW_FRACTION = 0.5;
    protected final int MIN_BREAKEND_CLUSTERING_WINDOW = 50;
    protected final int MAX_BREAKEND_CLUSTERING_WINDOW = 300;
    protected final int MIXED_CLUSTERING_WINDOW = 2000;
    protected BreakpointSummaryStrategy breakpointSummaryStrategy;

    public enum BreakpointSummaryStrategy {
        /**
         * Use the (first) middle value to summarize cluster starts and ends, such that the start and end were seen in the data
         */
        MEDIAN_START_MEDIAN_END,

        /**
         * A conservative strategy to summarize a cluster by its smallest extent
         */
        MIN_START_MAX_END,

        /**
         * A permissive strategy to summarize a cluster by it largest extent
         */
        MAX_START_MIN_END,

        /**
         * Summarize a cluster using the mean value for each end, even if that value was not represented in any sample
         */
        MEAN_START_MEAN_END

    }

    public SVClusterEngine(final SAMSequenceDictionary dictionary) {
        super(dictionary, CLUSTERING_TYPE.MAX_CLIQUE, null);
        this.breakpointSummaryStrategy = BreakpointSummaryStrategy.MEDIAN_START_MEDIAN_END;
    }

    public SVClusterEngine(final SAMSequenceDictionary dictionary, boolean depthOnly, BreakpointSummaryStrategy strategy) {
        super(dictionary, depthOnly ? CLUSTERING_TYPE.SINGLE_LINKAGE : CLUSTERING_TYPE.MAX_CLIQUE, null);
        this.breakpointSummaryStrategy = strategy;
    }

    /**
     * Find a single call representative of all the calls in the {@param cluster}
     * @param cluster   the events that are clustered together
     * @return  a call approximating the average event for the cluster and containing all the algorithms and genotypes
     */

    @Override
    protected SVCallRecord flattenCluster(final Collection<SVCallRecord> cluster) {
        final Collection<SVCallRecord> mostPreciseCalls;
        if (cluster.stream().allMatch(SVClusterEngine::isDepthOnlyCall)) {
            mostPreciseCalls = cluster;
        } else {
            mostPreciseCalls = cluster.stream().filter(call -> !SVDepthOnlyCallDefragmenter.isDepthOnlyCall(call)).collect(Collectors.toList());
        }
        final List<Integer> startPositions = mostPreciseCalls.stream().map(SVCallRecord::getPositionA).sorted().collect(Collectors.toList());
        final List<Integer> endPositions = mostPreciseCalls.stream().map(SVCallRecord::getPositionB).sorted().collect(Collectors.toList());
        //use the mid value of the sorted list so the start and end represent real breakpoint observations
        final int medianStart = startPositions.get(startPositions.size() / 2);
        final int medianEnd = endPositions.get(endPositions.size() / 2);
        final SVCallRecord exampleCall = mostPreciseCalls.iterator().next();
        final int length = exampleCall.getContigA().equals(exampleCall.getContigB()) && !exampleCall.getType().equals(StructuralVariantType.INS) ? medianEnd - medianStart + 1 : exampleCall.getLength();
        final List<String> algorithms = cluster.stream().flatMap(v -> v.getAlgorithms().stream()).distinct().collect(Collectors.toList());
        final List<Genotype> clusterSamples = cluster.stream().flatMap(v -> v.getGenotypes().stream()).collect(Collectors.toList());
        final List<StructuralVariantType> observedTypes = cluster.stream().map(SVCallRecord::getType).distinct().collect(Collectors.toList());
        final StructuralVariantType clusterType = observedTypes.size() == 1 ? observedTypes.get(0) : StructuralVariantType.CNV;

        final int newStart;
        final int newEnd;
        if (exampleCall.getType().equals(StructuralVariantType.INS)) {
            // Insertions should be a single locus; also fixes case where end-supporting split reads are to the
            // left of start-supporting split reads
            final int mean = (medianStart + medianEnd) / 2;
            newStart = mean;
            newEnd = mean;
        } else {
            switch (this.breakpointSummaryStrategy) {
                case MEDIAN_START_MEDIAN_END:
                    newStart = medianStart;
                    newEnd = medianEnd;
                    break;
                case MIN_START_MAX_END:
                    newStart = startPositions.stream().min(Integer::compareTo).orElse(startPositions.get(0));
                    newEnd = endPositions.stream().max(Integer::compareTo).orElse(endPositions.get(0));
                    break;
                case MAX_START_MIN_END:
                    newStart = startPositions.stream().max(Integer::compareTo).orElse(startPositions.get(0));
                    newEnd = endPositions.stream().min(Integer::compareTo).orElse(endPositions.get(0));
                    break;
                case MEAN_START_MEAN_END:
                    newStart = (int)Math.round(new Mean().evaluate(Doubles.toArray(startPositions)));
                    newEnd = (int)Math.round(new Mean().evaluate(Doubles.toArray(endPositions)));
                    break;
                default:
                    newStart = medianStart;
                    newEnd = medianEnd;
            }

        }

        //TODO: merge evidence for WGS data
        return new SVCallRecord(exampleCall.getId(), exampleCall.getContigA(), newStart, exampleCall.getStrandA(),
                exampleCall.getContigB(), newEnd, exampleCall.getStrandB(), clusterType, length, algorithms, clusterSamples);
    }

    @Override
    protected boolean clusterTogether(final SVCallRecord a, final SVCallRecord b) {
        //if (!a.getType().equals(b.getType())) return false;  //TODO: do we need to keep dels and dupes separate?
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
    protected SimpleInterval getClusteringInterval(final SVCallRecord call, final SimpleInterval clusterMinStartInterval) {
        final int minStart;
        final int maxStart;
        if (isDepthOnlyCall(call)) {
            minStart = (int) (call.getPositionB() - call.getLength() / MIN_RECIPROCAL_OVERLAP_DEPTH); //start of an overlapping event such that call represents (reciprocal overlap) of that event
            maxStart = (int) (call.getPositionA() + (1.0 - MIN_RECIPROCAL_OVERLAP_DEPTH) * call.getLength());
        } else {
            minStart = call.getPositionA() - MAX_BREAKEND_CLUSTERING_WINDOW;
            maxStart = call.getPositionA() + MAX_BREAKEND_CLUSTERING_WINDOW;
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
    protected SVDeduplicator<SVCallRecord> getDeduplicator() {
        final Function<Collection<SVCallRecord>,SVCallRecord> collapser = SVCallRecordUtils::deduplicateWithRawCallAttribute;
        return new SVCallRecordDeduplicator<>(collapser, dictionary);
    }

    protected boolean clusterTogetherBothDepthOnly(final SVCallRecord a, final SVCallRecord b) {
        if (!a.getContigA().equals(a.getContigB()) || !b.getContigA().equals(b.getContigB())) {
            throw new IllegalArgumentException("Attempted to cluster depth-only calls with endpoints on different contigs");
        }
        final SimpleInterval intervalA = new SimpleInterval(a.getContigA(), a.getPositionA(), a.getPositionB());
        final SimpleInterval intervalB = new SimpleInterval(b.getContigA(), b.getPositionA(), b.getPositionB());
        return IntervalUtils.isReciprocalOverlap(intervalA, intervalB, MIN_RECIPROCAL_OVERLAP_DEPTH);
    }

    protected boolean clusterTogetherBothWithEvidence(final SVCallRecord a, final SVCallRecord b) {
        // Reject if one is intrachromosomal and the other isn't
        final boolean intrachromosomalA = a.getContigA().equals(a.getContigB());
        final boolean intrachromosomalB = b.getContigB().equals(b.getContigB());
        if (intrachromosomalA != intrachromosomalB) return false;

        // Matching endpoints
        final SimpleInterval intervalAStart =  getStartClusteringInterval(a);
        final SimpleInterval intervalAEnd =  getEndClusteringInterval(a);
        final SimpleInterval intervalBStart =  getStartClusteringInterval(b);
        final SimpleInterval intervalBEnd =  getEndClusteringInterval(b);
        return intervalAStart.overlaps(intervalBStart) && intervalAEnd.overlaps(intervalBEnd);
    }

    protected boolean clusterTogetherMixedEvidence(final SVCallRecord a, final SVCallRecord b) {
        final boolean intrachromosomalA = a.getContigA().equals(a.getContigB());
        final boolean intrachromosomalB = b.getContigA().equals(b.getContigB());
        if (!(intrachromosomalA && intrachromosomalB)) return false;
        if (!a.getContigA().equals(b.getContigA())) return false;
        return Math.abs(a.getPositionA() - b.getPositionA()) < MIXED_CLUSTERING_WINDOW
                && Math.abs(a.getPositionB() - b.getPositionB()) < MIXED_CLUSTERING_WINDOW;
    }

    private SimpleInterval getStartClusteringInterval(final SVCallRecord call) {
        if (this.genomicToBinMap == null) {
            final int padding = getEndpointClusteringPadding(call);
            return call.getPositionAInterval().expandWithinContig(padding, dictionary);
        } else {
            return call.getPositionAInterval();
        }
    }

    protected SimpleInterval getEndClusteringInterval(final SVCallRecord call) {
        if (this.genomicToBinMap == null) {
            final int padding = getEndpointClusteringPadding(call);
            return call.getPositionBInterval().expandWithinContig(padding, dictionary);
        } else {
            return call.getPositionBInterval();
        }
    }

    protected int getEndpointClusteringPadding(final SVCallRecord call) {
        return (int) Math.min(MAX_BREAKEND_CLUSTERING_WINDOW, Math.max(MIN_BREAKEND_CLUSTERING_WINDOW, BREAKEND_CLUSTERING_WINDOW_FRACTION * call.getLength()));
    }

    public static boolean isDepthOnlyCall(final SVCallRecord call) {
        return call.getAlgorithms().size() == 1 && call.getAlgorithms().get(0).equals(GATKSVVCFConstants.DEPTH_ALGORITHM);
    }

    public static double getMinReciprocalOverlap() {
        return MIN_RECIPROCAL_OVERLAP_DEPTH;
    }
}
