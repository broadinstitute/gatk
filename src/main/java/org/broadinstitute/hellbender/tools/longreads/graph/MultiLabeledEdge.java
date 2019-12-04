package org.broadinstitute.hellbender.tools.longreads.graph;

import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.BaseEdge;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.SeqVertex;
import org.jgrapht.EdgeFactory;

import java.io.Serializable;
import java.util.LinkedHashSet;
import java.util.List;

public class MultiLabeledEdge extends LabeledEdge {

    private static final long     serialVersionUID = 0x1337L;
    private LinkedHashSet<String> labelSet         = new LinkedHashSet<>();

    public MultiLabeledEdge(final String readType) {
        super(readType);
        labelSet.add(readType);
    }

    private MultiLabeledEdge(final List<String> labelList) {
        super(String.join(",", labelList));
        this.labelSet = new LinkedHashSet<>(labelList);
    }

    private MultiLabeledEdge(final LinkedHashSet<String> labelList) {
        super(String.join(",", labelList));
        this.labelSet = new LinkedHashSet<>(labelList);
    }

    @Override
    public String getLabel() {
        return String.join(",", labelSet);
    }

    @Override
    public void updateLabel(final String label) {
        labelSet.add(label);
    }

    public void setLabel(final String label) {
        labelSet = new LinkedHashSet<>();
        labelSet.add(label);
    }

    @Override
    public String getDotExtraInfo() {
        return "readType=\"" + String.join(",", labelSet) + "\"";
    }

    @Override
    public BaseEdge copy() {
        return new MultiLabeledEdge(labelSet);
    }

    @Override
    public String getGexfAttributesString() {
        return "<attvalue for=\"0\" value=\"" + String.join(",", labelSet) + "\" />";
    }

    public static EdgeFactory<SeqVertex, BaseEdge> getFactory() {
        return new MultiLabeledEdgeFactory();
    }

    /**
     * Edge factory that creates labeled edges for this graph.
     */
    public static class MultiLabeledEdgeFactory implements EdgeFactory<SeqVertex, BaseEdge>, Serializable {

        private static final long serialVersionUID = 0x1337L;

        @Override
        /**
         * {@inheritDoc}
         *
         * Because we add edges to new nodes from existing nodes, we only need to check the target to know what kind of
         * read the edge came from.
         * See {@link AlignedBaseGraphCollection#initializeGraphWithNodes} and {@link AlignedBaseGraphCollection#mergeNodesIntoGraph}.
         */
        public BaseEdge createEdge(final SeqVertex sourceVertex, final SeqVertex targetVertex) {
            final String readName = ((AlignedBaseVertex)targetVertex).getReadName();
            final String label = createLabel(readName);
            return new MultiLabeledEdge(label);
        }
    }
}
