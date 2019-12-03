package org.broadinstitute.hellbender.tools.longreads.graph;

import com.google.common.primitives.Bytes;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.BaseEdge;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.SeqGraph;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.SeqVertex;
import org.jgrapht.EdgeFactory;

import java.util.*;
import java.util.regex.Pattern;

public class AlignedBaseGraph extends SeqGraph {

    private static final long serialVersionUID = 0x1337;

    private final Map<String, AlignedBaseVertex> contigToFirstNodeMap = new HashMap<>();
    private final Map<String, AlignedBaseVertex> contigToLastNodeMap = new HashMap<>();

    private final Set<String> contigSet = new LinkedHashSet<>();

    public AlignedBaseGraph() {
        super(1, new LabeledEdgeFactory());
    }

    /**
     * Edge factory that creates labeled edges for this graph.
     */
    private static class LabeledEdgeFactory implements EdgeFactory<SeqVertex, BaseEdge> {

        private static final Pattern CCS_READ_PATTERN = Pattern.compile("/ccs[ \t]*$");
        private static final Pattern RECLAIMED_READ_PATTERN = Pattern.compile("^m[0-9]*.*/[0-9]*/[0-9]*_[0-9]*_[0-9]*_[0-9]*_[0-9]*[ \t]*$");

        private static final String CCS_LABEL = "CCS";
        private static final String RECLAIMED_LABEL = "RECLAIMED";
        private static final String RAW_READ_LABEL = "RAW";

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

            if ( CCS_READ_PATTERN.matcher(readName).find() ) {
                // This is a CCS read!
                return new ReadTypedEdge(CCS_LABEL);
            }
            else if ( RECLAIMED_READ_PATTERN.matcher(readName).find() ) {
                // This is a RECLAIMED read!
                return new ReadTypedEdge(RECLAIMED_LABEL);
            }
            else {
                // This is a NORMAL read!
                return new ReadTypedEdge(RAW_READ_LABEL);
            }
        }
    }

    @Override
    public String getGexfEdgeAttributesDefinition() {
        return "  <attributes class=\"edge\"><attribute id=\"0\" title=\"readType\" type=\"string\"/></attributes>";
    }

    /**
     * Gets the first node in this graph for the first contig seen in this graph.
     * @return The first {@link AlignedBaseVertex} in the graph if the graph is not empty.  {@code null} otherwise.
     */
    public AlignedBaseVertex getFirstNode() {
        return getFirstNode(contigSet.iterator().next());
    }

    /**
     * Gets the last node in this graph for the last contig seen in this graph.
     * @return The last {@link AlignedBaseVertex} in the graph if the graph is not empty.  {@code null} otherwise.
     */
    public AlignedBaseVertex getLastNode() {
        return getLastNode(contigSet.iterator().next());
    }

    /**
     * Gets the first node in this graph for the given contig.
     * @param contig Contig in which to find the first node.
     * @return If the given contig exists in the graph, the first {@link AlignedBaseVertex} in the graph for the given contig.  {@code null} otherwise.
     */
    public AlignedBaseVertex getFirstNode(final String contig) {
        return contigToFirstNodeMap.get(contig);
    }

    /**
     * Gets the last node in this graph for the given contig.
     * @param contig Contig in which to find the last node.
     * @return If the given contig exists in the graph, the last {@link AlignedBaseVertex} in the graph for the given contig.  {@code null} otherwise.
     */
    public AlignedBaseVertex getLastNode(final String contig) {
        return contigToFirstNodeMap.get(contig);
    }

    @Override
    public boolean addVertex(final SeqVertex v)
    {
        final AlignedBaseVertex abv = (AlignedBaseVertex) v;

        final boolean newlyAdded = super.addVertex(abv);
        if (newlyAdded) {

            // Track this contig:
            if ( !contigSet.contains(abv.getPos().getContig()) ) {
                contigSet.add(abv.getPos().getContig());
            }

            // Check if we must track this node for our start:
            if ( contigToFirstNodeMap.containsKey(abv.getPos().getContig()) ) {
                contigToFirstNodeMap.replace(abv.getPos().getContig(), abv);
            }
            else {
                contigToFirstNodeMap.put(abv.getPos().getContig(), abv);
            }

            // Check if we must track this node for our end:
            if ( contigToLastNodeMap.containsKey(abv.getPos().getContig()) ) {
                if (abv.getPos().compareTo(contigToLastNodeMap.get(abv.getPos().getContig()).getPos()) == -1) {
                    contigToLastNodeMap.replace(abv.getPos().getContig(), abv);
                }
            }
            else {
                contigToLastNodeMap.put(abv.getPos().getContig(), abv);
            }
        }
        return newlyAdded;
    }

    @Override
    public void simplifyGraph() {
        zipLinearChains();
    }

    /**
     * Zip up all of the simple linear chains present in this graph which occur before the given position.
     *
     * Functions exactly as {@link #zipLinearChains()} except will only zip chains for nodes
     * that occur before the given position.
     *
     * @return true if any such pair of vertices could be found, false otherwise
     */
    public boolean zipLinearChainsBefore(final GenomicAndInsertionPosition position) {

        // create the list of start sites [doesn't modify graph yet]
        final Collection<SeqVertex> zipStarts = new LinkedList<>();
        for ( final SeqVertex source : vertexSet() ) {

            final AlignedBaseVertex abv = (AlignedBaseVertex) source;

            if ( (abv.getPos().getContig().equals(position.getContig())) &&
                    (abv.getPos().compareTo(position) < 0) &&
                    isLinearChainStart(source) ) {
                zipStarts.add(source);
            }
        }

        if ( zipStarts.isEmpty() ) // nothing to do, as nothing could start a chain
        {
            return false;
        }

        // At this point, zipStarts contains all of the vertices in this graph that might start some linear
        // chain of vertices.  We walk through each start, building up the linear chain of vertices and then
        //        // zipping them up with mergeLinearChain, if possible
        boolean mergedOne = false;
        for ( final SeqVertex zipStart : zipStarts ) {
            final LinkedList<SeqVertex> linearChain = traceLinearChainBefore((AlignedBaseVertex)zipStart, position);

            // merge the linearized chain, recording if we actually did some useful work
            mergedOne |= mergeLinearChain(linearChain);
        }

        return mergedOne;
    }

    /**
     * {@inheritDoc}
     *
     * Requires that nodes are at most distance 1 away from each other for zipping.
     */
    protected LinkedList<SeqVertex> traceLinearChain(final SeqVertex zipStart) {
        final LinkedList<SeqVertex> linearChain = new LinkedList<>();
        linearChain.add(zipStart);

        boolean lastIsRef = isReferenceNode(zipStart); // remember because this calculation is expensive
        SeqVertex last = zipStart;
        while (true) {
            if ( outDegreeOf(last) != 1 )
            // cannot extend a chain from last if last has multiple outgoing branches
            {
                break;
            }

            // there can only be one (outgoing edge of last) by contract
            final SeqVertex target = getEdgeTarget(outgoingEdgeOf(last));

            if ( inDegreeOf(target) != 1 || last.equals(target) )
            // cannot zip up a target that has multiple incoming nodes or that's a cycle to the last node
            {
                break;
            }

            final boolean targetIsRef = isReferenceNode(target);
            if ( lastIsRef != targetIsRef ) // both our isRef states must be equal
            {
                break;
            }

            // Make sure we only zip adjacent nodes.
            if ( !((AlignedBaseVertex)last).isAdjacentTo((AlignedBaseVertex)target) ) {
                break;
            }

            linearChain.add(target); // extend our chain by one

            // update our last state to be the current state, and continue
            last = target;
            lastIsRef = targetIsRef;
        }

        return linearChain;
    }

    /**
     * Traces a linear chain starting at the given vertex through all nodes with position less than the given pos.
     * 
     * This method is almost identical to {@link #traceLinearChain(SeqVertex)}, with the addition of another break
     * condition.
     * 
     * @param vertex {@link AlignedBaseVertex} at which to start the trace.  This node is included in the resulting chain.
     * @param pos {@link GenomicAndInsertionPosition} before which to include nodes in the chain.
     * @return A {@link LinkedList<SeqVertex>} containing a linear chain of nodes starting at the given vertex, each of which has a position less than the given pos.
     */
    private LinkedList<SeqVertex> traceLinearChainBefore(final AlignedBaseVertex vertex, final GenomicAndInsertionPosition pos) {
        final LinkedList<SeqVertex> linearChain = new LinkedList<>();
        linearChain.add(vertex);

        boolean lastIsRef = isReferenceNode(vertex); // remember because this calculation is expensive
        AlignedBaseVertex last = vertex;
        while (true) {
            if ( outDegreeOf(last) != 1 )
            // cannot extend a chain from last if last has multiple outgoing branches
            {
                break;
            }

            // there can only be one (outgoing edge of last) by contract
            final AlignedBaseVertex target = (AlignedBaseVertex)getEdgeTarget(outgoingEdgeOf(last));

            if ( inDegreeOf(target) != 1 || last.equals(target) )
            // cannot zip up a target that has multiple incoming nodes or that's a cycle to the last node
            {
                break;
            }

            final boolean targetIsRef = isReferenceNode(target);
            if ( lastIsRef != targetIsRef ) // both our isRef states must be equal
            {
                break;
            }

            if ( target.getPos().compareTo(pos) >= 0 ) {
                // Cannot add a target that occurs at or after the given position!
                break;
            }

            linearChain.add(target); // extend our chain by one

            // update our last state to be the current state, and continue
            last = target;
            lastIsRef = targetIsRef;
        }

        return linearChain;
    }

    @Override
    protected SeqVertex mergeLinearChainVertices(final Iterable<SeqVertex> vertices) {
        final List<byte[]> seqs = new LinkedList<>();

        boolean isFirst = true;
        GenomicAndInsertionPosition pos = null;

        int numMerged = 0;
        boolean mustReplaceFirstNode = false;
        boolean mustReplaceLastNode = false;
        for ( final SeqVertex v : vertices ) {
            seqs.add(v.getSequence());

            // Make sure we keep our first and last nodes correct:
            mustReplaceFirstNode |= contigToFirstNodeMap.get(((AlignedBaseVertex)v).getPos().getContig()).equals(v);
            mustReplaceLastNode |= contigToLastNodeMap.get(((AlignedBaseVertex)v).getPos().getContig()).equals(v);

            ++numMerged;
            if ( isFirst ) {
                pos = ((AlignedBaseVertex)v).getPos();
                isFirst = false;
            }
        }
        final byte[] seqsCat = Bytes.concat(seqs.toArray(new byte[][]{}));

        final AlignedBaseVertex abv = new AlignedBaseVertex(seqsCat, pos, "Merged_" + numMerged + "_Seqs");

        // Replace our nodes if we must:
        if ( mustReplaceFirstNode ) {
            contigToFirstNodeMap.replace( abv.getPos().getContig(), abv );
        }
        if ( mustReplaceLastNode ) {
            contigToLastNodeMap.replace( abv.getPos().getContig(), abv );
        }

        return abv;
    }
}
