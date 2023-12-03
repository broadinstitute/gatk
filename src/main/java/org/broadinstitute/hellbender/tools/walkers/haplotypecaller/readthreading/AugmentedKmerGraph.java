package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading;

import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.BaseGraph;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.MultiSampleEdge;
import org.jgrapht.EdgeFactory;

public class AugmentedKmerGraph extends BaseGraph<AugmentedVertex, MultiSampleEdge> {
    private static final long serialVersionUID = 1l;

    private static MyEdgeFactory EDGE_FACTORY = new MyEdgeFactory();

    protected AugmentedKmerGraph(int kmerSize) {
        super(kmerSize, EDGE_FACTORY);
    }

    private static final class MyEdgeFactory implements EdgeFactory<AugmentedVertex, MultiSampleEdge> {
        public MyEdgeFactory() { }

        @Override
        public MultiSampleEdge createEdge(final AugmentedVertex sourceVertex, final AugmentedVertex targetVertex) {
            return new MultiSampleEdge(false, 1, 1);
        }

        MultiSampleEdge createEdge(final boolean isRef, final int multiplicity) {
            return new MultiSampleEdge(isRef, multiplicity, 1);
        }
    }
}
