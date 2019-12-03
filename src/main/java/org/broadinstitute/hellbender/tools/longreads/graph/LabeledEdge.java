package org.broadinstitute.hellbender.tools.longreads.graph;

import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.BaseEdge;

public class LabeledEdge extends BaseEdge {

    private String label;

    public LabeledEdge(final String label) {
        super(false, 1);
        this.label = label;
    }

    public String getLabel() {
        return label;
    }

}
