package org.broadinstitute.hellbender.tools.longreads.graph;

import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.BaseEdge;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.SeqVertex;
import org.jgrapht.EdgeFactory;

import java.io.Serializable;
import java.util.regex.Pattern;

public class LabeledEdge extends BaseEdge {

    private static final long serialVersionUID = 0x1337L;

    protected static final Pattern CCS_READ_PATTERN       = Pattern.compile("/ccs[ \t]*$");
    protected static final Pattern RECLAIMED_READ_PATTERN = Pattern.compile("^m[0-9]*.*/[0-9]*/[0-9]*_[0-9]*_[0-9]*_[0-9]*_[0-9]*[ \t]*$");

    private static final String CCS_LABEL = "CCS";
    private static final String RECLAIMED_LABEL = "RECLAIMED";
    private static final String RAW_READ_LABEL = "RAW";

    private String label;

    public LabeledEdge(final String label) {
        super(false, 1);
        this.label = label;
    }

    public String getLabel() {
        return label;
    }

    public void updateLabel(final String label) {
        this.label = label;
    }

    /**
     * Creates a read label string based on the given read name.
     * @param vertex {@link AlignedBaseVertex} containing the name of the read to label.
     * @return A {@link String} containing a read label.
     */
    public static String createLabel(final AlignedBaseVertex vertex) {
        return createLabel(vertex.getReadName());
    }

    /**
     * Creates a read label string based on the given read name.
     * @param readName {@link String} containing the name of the read to label.
     * @return A {@link String} containing a read label.
     */
    public static String createLabel(final String readName) {
//            // TODO: DEBUGGING TO BE REMOVED AT CODE REVIEW TIME
//            System.out.println("ReadName = " + readName);
//            System.out.println("CCS_PATTERN_FIND = " + CCS_READ_PATTERN.matcher(readName).find());
//            System.out.println("RECLAIMED_READ_PATTERN_FIND = " + RECLAIMED_READ_PATTERN.matcher(readName).find());

        if ( CCS_READ_PATTERN.matcher(readName).find() ) {
            // This is a CCS read!
            return CCS_LABEL;
        }
        else if ( RECLAIMED_READ_PATTERN.matcher(readName).find() ) {
            // This is a RECLAIMED read!
            return RECLAIMED_LABEL;
        }
        else {
            // This is a NORMAL read!
            return RAW_READ_LABEL;
        }
    }

    @Override
    public String getDotExtraInfo() {
        return "readType=\"" + label + "\"";
    }

    @Override
    public BaseEdge copy() {
        return new LabeledEdge(label);
    }

    @Override
    public String getGexfAttributesString() {
        return "<attvalue for=\"0\" value=\"" + label + "\" />";
    }

    public static EdgeFactory<SeqVertex, BaseEdge> getFactory() {
        return new LabeledEdgeFactory();
    }

    /**
     * Edge factory that creates labeled edges for this graph.
     */
    public static class LabeledEdgeFactory implements EdgeFactory<SeqVertex, BaseEdge>, Serializable {

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
            final String label = LabeledEdge.createLabel(readName);
            return new LabeledEdge(label);
        }
    }
}
