package org.broadinstitute.hellbender.tools.longreads.graph;

import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.BaseEdge;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.SeqVertex;
import org.jgrapht.EdgeFactory;

public enum LabeledEdgeType {
    LABELED_EDGE(LabeledEdge.getFactory()),
    MULTI_LABELED_EDGE(MultiLabeledEdge.getFactory());

    private final EdgeFactory<SeqVertex, BaseEdge> edgeFactory;

    private LabeledEdgeType(final EdgeFactory<SeqVertex, BaseEdge> factory) {
        edgeFactory = factory;
    }

    public EdgeFactory<SeqVertex, BaseEdge> getEdgeFactory() {
        return edgeFactory;
    }
}
