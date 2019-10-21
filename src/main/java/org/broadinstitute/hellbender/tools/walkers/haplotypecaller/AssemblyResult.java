package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.SeqGraph;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.AbstractReadThreadingGraph;
import org.broadinstitute.hellbender.utils.Utils;

/**
 * Result of assembling, with the resulting graph and status
 */
public final class AssemblyResult {
    private final Status status;
    private final AbstractReadThreadingGraph threadingGraph;
    private final SeqGraph graph;

    /**
     * Create a new assembly result
     * @param status the status, cannot be null
     * @param graph the resulting graph of the assembly, can only be null if result is failed
     */
    public AssemblyResult(final Status status, final SeqGraph graph, final AbstractReadThreadingGraph threadingGraph) {
        Utils.nonNull(status, "status cannot be null");
        Utils.validateArg( status == Status.FAILED || (graph != null || threadingGraph != null) , "graph is null but status is " + status);

        this.status = status;
        this.graph = graph;
        this.threadingGraph = threadingGraph;
    }

    public AbstractReadThreadingGraph getThreadingGraph() {
        return threadingGraph;
    }

    public Status getStatus() {
        return status;
    }

    public SeqGraph getSeqGraph() {
        return graph;
    }

    public int getKmerSize() {
        return graph == null? threadingGraph.getKmerSize() : graph.getKmerSize();
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
