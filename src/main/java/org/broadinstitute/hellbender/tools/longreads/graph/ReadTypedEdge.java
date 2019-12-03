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

    @Override
    public String getDotExtraInfo() {
        return "readType=\"" + readType + "\"";
    }

    @Override
    public BaseEdge copy() {
        return new ReadTypedEdge( readType );
    }

    @Override
    public String getGexfAttributesString() {
        return "<attvalue for=\"0\" value=\"" + readType + "\" />";
    }

}
