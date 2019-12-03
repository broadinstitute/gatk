package org.broadinstitute.hellbender.tools.longreads.graph;

import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.BaseEdge;

public class ReadTypedEdge extends BaseEdge {

    private String readType;

    public ReadTypedEdge(final String readType) {
        super(false, 1);
        this.readType = readType;
    }

    public String getReadType() {
        return readType;
    }

    public String getDotExtraInfo() {
        return "readType=\"" + readType + "\"";
    }

}
