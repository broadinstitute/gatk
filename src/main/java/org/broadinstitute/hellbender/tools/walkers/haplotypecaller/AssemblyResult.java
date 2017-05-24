package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.BaseEdge;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.BaseGraph;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.SeqGraph;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.SeqVertex;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ReadThreadingGraph;
import org.broadinstitute.hellbender.utils.Utils;

/**
 * Result of assembling, with the resulting graph and status
 */
public final class AssemblyResult {
    private final Status status;
    private final ReadThreadingGraph threadingGraph;
    private final SeqGraph graph;

    /**
     * Create a new assembly result
     * @param status the status, cannot be null
     * @param graph the resulting graph of the assembly, can only be null if result is failed
     */
    public AssemblyResult(final Status status, final SeqGraph graph, final ReadThreadingGraph threadingGraph) {
        Utils.nonNull(status, "status cannot be null");
        if ( status != Status.FAILED && graph == null ) {
            throw new IllegalArgumentException("graph is null but status is " + status);
        }

        this.status = status;
        this.graph = graph;
        this.threadingGraph = threadingGraph;
    }

    public ReadThreadingGraph getThreadingGraph() {
        return threadingGraph;
    }

    public Status getStatus() {
        return status;
    }

    public SeqGraph getGraph() {
        return graph;
    }

    public int getKmerSize() {
        return graph.getKmerSize();
    }

    /**
     * Status of the assembly result
     */
    public enum Status {
        /** Something went wrong, and we couldn't produce a meaningful graph */
        FAILED,
        /** Assembly succeeded, but graph degenerated into just the reference sequence */
        JUST_ASSEMBLED_REFERENCE,
        /** Assembly succeeded, and the graph has some meaningful structure */
        ASSEMBLED_SOME_VARIATION
    }
}
