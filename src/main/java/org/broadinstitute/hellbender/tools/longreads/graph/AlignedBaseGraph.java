package org.broadinstitute.hellbender.tools.longreads.graph;

import com.google.common.primitives.Bytes;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.SeqGraph;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.SeqVertex;

import java.util.*;

public class AlignedBaseGraph extends SeqGraph {

    private static final long serialVersionUID = 0x1337;

    private final Map<String, AlignedBaseVertex> contigToFirstNodeMap = new HashMap<>();
    private final Map<String, AlignedBaseVertex> contigToLastNodeMap = new HashMap<>();

    private final Set<String> contigSet = new LinkedHashSet<>();

    public AlignedBaseGraph() {
        super(1);
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
