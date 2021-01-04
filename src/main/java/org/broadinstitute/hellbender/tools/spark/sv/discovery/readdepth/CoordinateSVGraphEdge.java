package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalUtils;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.List;

/**
 * Graph edge class with endpoint coordinates A and B, where A < B in dictionary sorted order
 */
public final class CoordinateSVGraphEdge {
    private final int contigA;
    private final int contigB;
    private final int nodeAPosition;
    private final int nodeBPosition;
    private final boolean strandA;
    private final boolean strandB;
    private final boolean inverted;
    private final boolean isReference;
    private final SVGraphEdgeEvidence evidence;
    //private final SVGraphEdgePrior prior;

    public CoordinateSVGraphEdge(final int contig1, final int position1, final boolean strand1,
                                 final int contig2, final int position2, final boolean strand2,
                                 final boolean isReference, final SVGraphEdgeEvidence evidence,
                                 final SAMSequenceDictionary dictionary) {
        Utils.nonNull(evidence, "Evidence cannot be null");
        Utils.nonNull(dictionary, "Sequence dictionary cannot be null");
        Utils.validateArg(contig1 == contig2, "Interchromosomal edges not currently supported");
        final SimpleInterval interval1 = SVIntervalUtils.createSimpleInterval(contig1, position1, position1 + 1, dictionary);
        final SimpleInterval interval2 = SVIntervalUtils.createSimpleInterval(contig2, position2, position2 + 1, dictionary);
        if (IntervalUtils.compareLocatables(interval1, interval2, dictionary) < 0) {
            contigA = contig1;
            contigB = contig2;
            nodeAPosition = position1;
            nodeBPosition = position2;
            strandA = strand1;
            strandB = strand2;
        } else {
            contigA = contig2;
            contigB = contig1;
            nodeAPosition = position2;
            nodeBPosition = position1;
            strandA = strand2;
            strandB = strand1;
        }
        this.isReference = isReference;
        this.inverted = false;
        this.evidence = evidence;
        //this.prior = new SVGraphEdgePrior(evidence, 15, 6); //TODO
    }

    public CoordinateSVGraphEdge(final IndexedSVGraphEdge indexedEdge, final List<SVGraphNode> nodes) {
        Utils.nonNull(indexedEdge, "Indexed edge cannot be null");
        Utils.nonNull(nodes, "Nodes list cannot be null");
        Utils.validateArg(indexedEdge.getNodeAIndex() < nodes.size() && indexedEdge.getNodeAIndex() >= 0, "Edge node A index must be non-negative and less than node list size");
        Utils.validateArg(indexedEdge.getNodeBIndex() < nodes.size() && indexedEdge.getNodeBIndex() >= 0, "Edge node B index must be non-negative and less than node list size");
        final SVGraphNode nodeA = nodes.get(indexedEdge.getNodeAIndex());
        final SVGraphNode nodeB = nodes.get(indexedEdge.getNodeBIndex());
        contigA = nodeA.getContig();
        contigB = nodeB.getContig();
        nodeAPosition = nodeA.getPosition();
        nodeBPosition = nodeB.getPosition();
        isReference = indexedEdge.isReference();
        strandA = indexedEdge.isStrandA();
        strandB = indexedEdge.isStrandB();
        inverted = indexedEdge.isInverted();
        evidence = indexedEdge.getEvidence();
        //prior = indexedEdge.getPrior();
    }

    public SVGraphEdgeEvidence getEvidence() {
        return evidence;
    }

    public SVInterval getInterval() {
        if (contigA != contigB) return null;
        if (nodeAPosition < nodeBPosition) return new SVInterval(contigA, nodeAPosition, nodeBPosition);
        return new SVInterval(contigA, nodeBPosition, nodeAPosition);
    }

    public SVInterval getIntervalA() {
        return new SVInterval(getContigA(), getNodeAPosition(), getNodeAPosition() + 1);
    }

    public SVInterval getIntervalB() {
        return new SVInterval(getContigB(), getNodeBPosition(), getNodeBPosition() + 1);
    }

    //public double getEdgeVisitPrior(int numVisits) {
    //    return prior.getLogPrior(numVisits);
    //}

    public boolean isIntrachromosomal() { return contigA == contigB; }

    public int getContigA() {
        return contigA;
    }

    public int getContigB() {
        return contigB;
    }

    public int getNodeAPosition() {
        return nodeAPosition;
    }

    public int getNodeBPosition() {
        return nodeBPosition;
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

    public boolean isInverted() {
        return inverted;
    }

    //public SVGraphEdgePrior getPrior() {
    //    return prior;
    //}
}
