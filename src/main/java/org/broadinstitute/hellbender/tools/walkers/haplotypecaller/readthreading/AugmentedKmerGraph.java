package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading;

import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.BaseGraph;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.MultiSampleEdge;
import org.jgrapht.EdgeFactory;

public class AugmentedKmerGraph extends BaseGraph<MultiDeBruijnVertex, MultiSampleEdge> {

    private static MyEdgeFactory EDGE_FACTORY = new MyEdgeFactory();

    protected AugmentedKmerGraph(int kmerSize) {
        super(kmerSize, EDGE_FACTORY);
    }

    private static final class MyEdgeFactory implements EdgeFactory<MultiDeBruijnVertex, MultiSampleEdge> {
        public MyEdgeFactory() { }

        @Override
        public MultiSampleEdge createEdge(final MultiDeBruijnVertex sourceVertex, final MultiDeBruijnVertex targetVertex) {
            return new MultiSampleEdge(false, 1, 1);
        }

        MultiSampleEdge createEdge(final boolean isRef, final int multiplicity) {
            return new MultiSampleEdge(isRef, multiplicity, 1);
        }
    }
}
