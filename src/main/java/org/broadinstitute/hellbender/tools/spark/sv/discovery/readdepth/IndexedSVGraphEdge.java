package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalUtils;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.List;
import java.util.Objects;

/**
 * Edge with endpoint positions defined by graph node indices. This class is used for traversing SV graphs, as it
 * also defines whether it is inverted and which node is upstream. For strand-switch edges, the upstream node is
 * undefined by default but can be set during a traversal.
 */
public final class IndexedSVGraphEdge implements Comparable<IndexedSVGraphEdge> {

    final static int UNDFINED_SOURCE_NODE_INDEX = -1;
    private final int index;
    private final boolean strandA;
    private final boolean strandB;
    private final boolean isReference;
    private final boolean inverted;
    private final int nodeAIndex;
    private final int nodeBIndex;
    private final int upstreamNodeIndex;
    private final SVInterval interval;
    private final SVGraphEdgeEvidence evidence;
    private final SVGraphEdgePrior prior;

    private IndexedSVGraphEdge(final int index, final int nodeAIndex, final int nodeBIndex, final boolean strandA, final boolean strandB,
                               final boolean isReference, final boolean inverted, final int upstreamNodeIndex, final SVGraphEdgeEvidence evidence,
                               final SVGraphEdgePrior prior, final SVInterval interval, final SVGraph graph, final SAMSequenceDictionary dictionary,
                               final boolean isCopy) {
        Utils.nonNull(evidence, "Evidence cannot be null");
        if (!isCopy) {
            //Validate node indices
            Utils.nonNull(graph, "Graph cannot be null");
            Utils.nonNull(dictionary, "Dictionary cannot be null");
            final int nodesListSize = graph.getNodes().size();
            Utils.validateArg(nodeAIndex >= 0 && nodeAIndex < nodesListSize, "Invalid node A index: " + nodeAIndex);
            Utils.validateArg(nodeBIndex >= 0 && nodeBIndex < nodesListSize, "Invalid node B index: " + nodeBIndex);
            Utils.validateArg(upstreamNodeIndex == nodeAIndex || upstreamNodeIndex == nodeBIndex || upstreamNodeIndex == UNDFINED_SOURCE_NODE_INDEX, "Source node index must be undefined, node A, or node B");
            final SVGraphNode nodeA = graph.getNodes().get(nodeAIndex);
            final SVGraphNode nodeB = graph.getNodes().get(nodeBIndex);
            final SimpleInterval intervalA = SVIntervalUtils.createSimpleInterval(nodeA.getContig(), nodeA.getPosition(), nodeA.getPosition() + 1, dictionary);
            final SimpleInterval intervalB = SVIntervalUtils.createSimpleInterval(nodeB.getContig(), nodeB.getPosition(), nodeB.getPosition() + 1, dictionary);
            Utils.validateArg(IntervalUtils.compareLocatables(intervalA, intervalB, dictionary) <= 0, "Node A is not before node B");
            this.interval = generateInterval(nodeAIndex, nodeBIndex, graph.getNodes());
        } else {
            Utils.nonNull(interval, "Interval cannot be null");
            this.interval = interval;
        }
        this.index = index;
        this.nodeAIndex = nodeAIndex;
        this.nodeBIndex = nodeBIndex;
        this.strandA = strandA;
        this.strandB = strandB;
        this.isReference = isReference;
        this.inverted = inverted;
        this.upstreamNodeIndex = upstreamNodeIndex;
        this.evidence = evidence;
        this.prior = prior;
    }

    public IndexedSVGraphEdge(final int index, final int nodeAIndex, final int nodeBIndex, final boolean strandA, final boolean strandB,
                              final boolean isReference, final SVGraphEdgeEvidence evidence, final SVGraph graph, final SAMSequenceDictionary dictionary) {
        //TODO prior parameters
        this(index, nodeAIndex, nodeBIndex, strandA, strandB, isReference, false, defaultUpstreamNode(false, nodeAIndex, nodeBIndex, strandA, strandB), evidence, new SVGraphEdgePrior(evidence, 15, 6), null, graph, dictionary, false);
    }

    private static int defaultUpstreamNode(final boolean inverted, final int nodeAIndex, final int nodeBIndex, final boolean strandA, final boolean strandB) {
        //Undefined if strand-switch
        if (strandA == strandB) return UNDFINED_SOURCE_NODE_INDEX;
        //Negative strand node if inverted
        if (inverted) {
            return strandA ? nodeBIndex : nodeAIndex;
        }
        //Positive strand node by default
        return strandA ? nodeAIndex : nodeBIndex;
    }

    private static SVInterval generateInterval(final int nodeAIndex, final int nodeBIndex, final List<SVGraphNode> nodes) {
        final SVGraphNode nodeA = nodes.get(nodeAIndex);
        final SVGraphNode nodeB = nodes.get(nodeBIndex);
        if (nodeA.getContig() != nodeB.getContig()) return null;
        return new SVInterval(nodeA.getContig(), nodeA.getPosition(), nodeB.getPosition());
    }

    /**
     * Creates a copy in the non-inverted state
     */
    public IndexedSVGraphEdge nonInvertedCopy() {
        return new IndexedSVGraphEdge(index, nodeAIndex, nodeBIndex, strandA, strandB, isReference, false, getUpstreamNodeIndex(), evidence, prior, interval, null, null, true);
    }

    /**
     * Creates a copy in the inverted state
     */
    public IndexedSVGraphEdge invertedCopy() {
        return new IndexedSVGraphEdge(index, nodeAIndex, nodeBIndex, strandA, strandB, isReference, true, getDownstreamNodeIndex(), evidence, prior, interval, null, null, true);
    }

    /**
     * For strand switch edges, creates a copy in the current inversion state
     */
    public IndexedSVGraphEdge copyWithUpstreamNodeAs(final int upstreamNodeIndex) {
        Utils.validate(causesStrandSwitch(), "Cannot define source node for non-strand switch edge");
        return new IndexedSVGraphEdge(index, nodeAIndex, nodeBIndex, strandA, strandB, isReference, inverted, upstreamNodeIndex, evidence, prior, interval, null, null, true);
    }

    public SVInterval getInterval() {
        return interval;
    }

    public SVGraphEdgeEvidence getEvidence() {
        return evidence;
    }

    public SVGraphEdgePrior getPrior() { return prior; }

    public boolean causesStrandSwitch() {
        return strandA == strandB;
    }

    public boolean isInverted() {
        return inverted;
    }

    public int getIndex() {
        return index;
    }

    public int getNodeAIndex() {
        return nodeAIndex;
    }

    public int getUpstreamNodeIndex() {
        return upstreamNodeIndex;
    }

    public int getDownstreamNodeIndex() {
        if (upstreamNodeIndex == UNDFINED_SOURCE_NODE_INDEX) return UNDFINED_SOURCE_NODE_INDEX;
        return upstreamNodeIndex == nodeAIndex ? nodeBIndex : nodeAIndex;
    }

    public int getNodeBIndex() {
        return nodeBIndex;
    }

    public boolean isStrandA() {
        return strandA;
    }

    public boolean isStrandB() {
        return strandB;
    }

    public boolean isReference() {
        return isReference;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof IndexedSVGraphEdge)) return false;
        IndexedSVGraphEdge that = (IndexedSVGraphEdge) o;
        return index == that.index &&
                strandA == that.strandA &&
                strandB == that.strandB &&
                isReference == that.isReference &&
                inverted == that.inverted &&
                nodeAIndex == that.nodeAIndex &&
                nodeBIndex == that.nodeBIndex &&
                upstreamNodeIndex == that.upstreamNodeIndex &&
                Objects.equals(interval, that.interval) &&
                Objects.equals(evidence, that.evidence);
    }

    @Override
    public int hashCode() {
        return Objects.hash(index, strandA, strandB, isReference, inverted, nodeAIndex, nodeBIndex, upstreamNodeIndex, interval, evidence);
    }

    @Override
    public int compareTo(final IndexedSVGraphEdge other) {
        int result = Integer.compare(index, other.index);
        if (result == 0) {
            result = Boolean.compare(strandA, other.strandA);
            if (result == 0) {
                result = Boolean.compare(strandB, other.strandB);
                if (result == 0) {
                    result = Boolean.compare(isReference, other.isReference);
                    if (result == 0) {
                        result = Boolean.compare(inverted, other.inverted);
                        if (result == 0) {
                            result = Integer.compare(nodeAIndex, other.nodeAIndex);
                            if (result == 0) {
                                result = Integer.compare(nodeBIndex, other.nodeBIndex);
                                if (result == 0) {
                                    result = Integer.compare(upstreamNodeIndex, other.upstreamNodeIndex);
                                }
                            }
                        }
                    }
                }
            }
        }
        return result;
    }

}
