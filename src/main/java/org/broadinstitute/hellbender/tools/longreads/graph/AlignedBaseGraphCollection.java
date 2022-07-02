package org.broadinstitute.hellbender.tools.longreads.graph;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.BaseEdge;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.SeqVertex;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.File;
import java.io.Serializable;
import java.util.*;

public class AlignedBaseGraphCollection implements Serializable {
    transient private static final Logger logger = LogManager.getLogger(AlignedBaseGraphCollection.class);

    //==================================================================================================================
    // Private Members:

    private static final long serialVersionUID = 0x1337L;

    private final HashMap<String, AlignedBaseGraph> contigSubGraphMap = new HashMap<>();
    private final HashMap<String, TreeMap<GenomicAndInsertionPosition, Set<AlignedBaseVertex>>> contigPositionVertexMap
            = new HashMap<>();

    private final HashMap<String, GenomicAndInsertionPosition> contigUncollapsedPositionMap = new HashMap<>();

    private boolean isGraphCollapsed = false;

    private boolean relabelEdgeTypes = false;

    private int numSequencesAdded = 0;

    private int periodicMergeDistance = 100;

    private LabeledEdgeType edgeType = LabeledEdgeType.LABELED_EDGE;

    //==================================================================================================================
    // Constructors:

    public AlignedBaseGraphCollection() {}
    public AlignedBaseGraphCollection(final LabeledEdgeType edgeType) { this.edgeType = edgeType; }

    //==================================================================================================================
    // Override Methods:

    //==================================================================================================================
    // Static Methods:

    //==================================================================================================================
    // Instance Methods:

    /**
     * Dumps all graphs in this aligned base graph to disk as DOT files.
     * The files will be named according to their contigs and the given {@code baseName}.
     *
     * @param baseName The common base for the names of all output files.
     */
    public void serializeToDotFiles(final String baseName) {
        for ( final Map.Entry<String, AlignedBaseGraph> entry : contigSubGraphMap.entrySet() ) {
            // Print the graph to a file.  We don't want to prune anything, so give the pruner a factor of zero.
            entry.getValue().printGraph(new File(baseName + '.' + entry.getKey() + ".dot"), 0 );
        }
    }

    /**
     * Dumps all graphs in this aligned base graph to disk as GEXF files.
     * The files will be named according to their contigs and the given {@code baseName}.
     *
     * @param baseName The common base for the names of all output files.
     */
    public void serializeToGexfFiles(final String baseName) {
        for ( final Map.Entry<String, AlignedBaseGraph> entry : contigSubGraphMap.entrySet() ) {
            // Print the graph to a file.  We don't want to prune anything, so give the pruner a factor of zero.
            entry.getValue().serializeToGexf(
                    new File(baseName + '.' + entry.getKey() + ".gexf"),
                    0,
                    "GATK - " + this.getClass().getCanonicalName(),
                    "Graph representing genomic data on " + entry.getKey());
        }
    }

    /**
     * Dumps all graphs in this aligned base graph to disk as GFA 1.0 files.
     * The files will be named according to their contigs and the given {@code baseName}.
     *
     * GFA 1 spec is located here: http://gfa-spec.github.io/GFA-spec/GFA1.html
     *
     * @param baseName The common base for the names of all output files.
     */
    public void serializeToGfa1Files(final String baseName) {
        helpSerializeToGfaFile(baseName, false);
    }

    public boolean isGraphCollapsed() {
        return isGraphCollapsed;
    }

    public int getNumSequencesAdded() {
        return numSequencesAdded;
    }

    /**
     * Dumps all graphs in this aligned base graph to disk as GFA 2.0 files.
     * The files will be named according to their contigs and the given {@code baseName}.
     *
     * GFA 2 spec is located here: https://github.com/GFA-spec/GFA-spec/blob/master/GFA2.md
     *
     * @param baseName The common base for the names of all output files.
     */
    public void serializeToGfa2Files(final String baseName) {
        helpSerializeToGfaFile(baseName, true);
    }

    private void helpSerializeToGfaFile(final String baseName, final boolean isGfa2) {

        for ( final Map.Entry<String, AlignedBaseGraph> entry : contigSubGraphMap.entrySet() ) {

            final String contig = entry.getKey();
            final AlignedBaseGraph graph = entry.getValue();

            final String fileBaseName = baseName + '.' + contig;

            if ( isGfa2 ) {
                graph.serializeToGfa2Files(fileBaseName);
            }
            else {
                graph.serializeToGfa1Files(fileBaseName);
            }
        }
    }

    /**
     * Collapse all linear chains of nodes into one large node containing all bases in the collapsed nodes.
     * This should be done before writing the graph to GFA format so you can stay sane.
     */
    public void collapseAdjacentNodes() {

        logger.info("Starting to collapse all graphs.  Total sequences added: " + numSequencesAdded);

        // Zip all linear chains:
        for ( final Map.Entry<String, AlignedBaseGraph> entry : contigSubGraphMap.entrySet() ) {

            final String contig = entry.getKey();
            final AlignedBaseGraph graph = entry.getValue();

            logger.info("  Collapsing graph (" + graph.vertexSet().size() + " nodes) for contig: " + contig);

            // Before we get started, clear the existing position map to make some room in memory.
            // When we finish, we'll rebuild it and there will be fewer nodes in there so it should take less space.
            contigPositionVertexMap.get( contig ).clear();

            // First zip chains to make the next step faster.
            logger.info("    Zipping linear chains (1 / 2)");
            graph.zipLinearChains();

            // the new ZipLinearChains does these steps as well:

//            // Here we must check for adjacent positions in the graph which have no links between them and add such links.
//            // We should do this because we know they must be adjacent based on the alignment information we were given
//            // in the supporting reads.
//            logger.info("    Linking adjacent nodes (2 / 4)");
//            linkAdjacentNodes(graph);
//
//            // Because we've modified our graphs, we must Zip all linear chains again.
//            // It will be faster this time around because there are fewer nodes in the graph.
//            logger.info("    Zipping linear chains again (3 / 4)");
//            graph.zipLinearChains();

            //NOTE: If we want to add more nodes after this it will be hard because
            //      we don't have single base nodes anymore, so the pileups will be hard.
            //      For now, we do not allow it.
            // Now we rebuild the contigPositionVertexMap for this contig so we can write out the data later:
            logger.info("    Rebuilding position node map (2 / 2)");
            rebuildContigPositionVertexMap(graph);

            logger.info("  Collapsed graph (" + contig + ") size: " + graph.vertexSet().size());
        }

        isGraphCollapsed = true;
    }

    private void rebuildContigPositionVertexMap(final AlignedBaseGraph graph) {
        final AlignedBaseVertex firstNode = graph.getFirstNode();
        final String contig = firstNode.getPos().getContig();

        // Clear the existing map:
        contigPositionVertexMap.get( contig ).clear();

        // Iterate through all nodes in the graph and add them in:
        for ( final SeqVertex v : graph.vertexSet() ) {
            final AlignedBaseVertex vertex = (AlignedBaseVertex) v;

            // Add in our holding set (if needed) and any nodes at that position:
            if ( !contigPositionVertexMap.get( contig ).containsKey(vertex.getPos()) ) {
                contigPositionVertexMap.get( contig ).put( vertex.getPos(), new HashSet<>() );
            }
            contigPositionVertexMap.get( contig ).get( vertex.getPos() ).add( vertex );
        }
    }

    /**
     * Link nodes that are at adjacent positions in the graph according to their position information.
     * We have already resolved the positions from the input alignment information, so if two nodes are exactly
     * adjacent, they can be chained together.  This can occur in any data set.
     * @param graph An {@link AlignedBaseGraph} contaning nodes to link together.
     */
    private static void linkAdjacentNodes(final AlignedBaseGraph graph) {

        final HashSet<AlignedBaseVertex> startVertices = new HashSet<>();
        final HashSet<AlignedBaseVertex> endVertices = new HashSet<>();

        // Get all vertices in the graph that have in- or out-degree of 0:
        for ( final SeqVertex vertex : graph.vertexSet() ) {
            if ( graph.inDegreeOf(vertex) == 0 ) {
                startVertices.add((AlignedBaseVertex) vertex);
            }
            else if ( graph.outDegreeOf(vertex) == 0 ) {
                endVertices.add((AlignedBaseVertex) vertex);
            }
        }

        // Now we can check the starts / ends to see if they're next to eachother.
        // If so, we can add a link between them:
        for ( final AlignedBaseVertex endVertex : endVertices ) {
            for ( final AlignedBaseVertex startVertex : startVertices ) {
                if ( endVertex.isAdjacentTo(startVertex) ) {
                    graph.addEdge(endVertex, startVertex);
                }
            }
        }
    }

    /**
     * Link nodes that are at adjacent positions in the graph according to their position information if and only if
     * they have genomic positions that occur before the given {@code pos}.
     * We have already resolved the positions from the input alignment information, so if two nodes are exactly
     * adjacent, they can be chained together.  This can occur in any data set.
     * @param graph An {@link AlignedBaseGraph} contaning nodes to link together.
     * @param pos A {@link GenomicAndInsertionPosition} to use as a boundary for linking nodes in the given {@code graph}.
     */
    private static void linkAdjacentNodesBefore(final AlignedBaseGraph graph, final GenomicAndInsertionPosition pos) {
        final HashSet<AlignedBaseVertex> startVertices = new HashSet<>();
        final HashSet<AlignedBaseVertex> endVertices = new HashSet<>();

        // Get all vertices in the graph that have in- or out-degree of 0:
        for ( final SeqVertex vertex : graph.vertexSet() ) {

            // Only look at nodes before the given pos:
            if ( (((AlignedBaseVertex)vertex).getPos().compareTo(pos) < 0) ) {
                if ( graph.inDegreeOf(vertex) == 0 ) {
                    startVertices.add((AlignedBaseVertex) vertex);
                }
                else if ( graph.outDegreeOf(vertex) == 0 ) {
                    endVertices.add((AlignedBaseVertex) vertex);
                }
            }
        }

        // Now we can check the starts / ends to see if they're next to eachother.
        // If so, we can add a link between them:
        for ( final AlignedBaseVertex endVertex : endVertices ) {
            for ( final AlignedBaseVertex startVertex : startVertices ) {
                if ( endVertex.isAdjacentTo(startVertex) ) {
                    graph.addEdge(endVertex, startVertex);
                }
            }
        }
    }

    /**
     * Add in the given aligned sequence information to this {@link AlignedBaseGraphCollection}.
     * @param read An aligned {@link GATKRead} object to add to this {@link AlignedBaseGraphCollection}.
     */
    public void addSequence(final GATKRead read) {

        final byte[]             bases    = read.getBasesNoCopy();
        final List<CigarElement> cigar    = read.getCigarElements();
        final String             contig   = read.getContig();
        final int                startPos = read.getStart();
        final String             readName = read.getName();

        addSequence(bases, cigar, contig, startPos, readName);
    }

    /**
     * Add in the given aligned sequence information to this {@link AlignedBaseGraphCollection}.
     * This information is assumed to be pulled directly from a {@link GATKRead} object.  This has implications for
     * certain parameters (the {@code startPos} in particular).
     *
     * Data can only be added to this graph before it is collapsed.  After that, it is an error to add anything else.
     *
     * @param bases Bases to be added in the order in which they appear in the read.
     * @param cigar {@link List<CigarElement>} containing alignment information for each base in {@code bases}.
     * @param contig The contig to which the given base sequence is aligned.
     * @param startPos The start position of the first aligned base in the given base sequence.
     * @param readName The name of the parent read from which this aligned sequence information is derived.
     */
    public void addSequence(final byte[] bases,
                            final List<CigarElement> cigar,
                            final String contig,
                            final int startPos,
                            final String readName) {

        // Make sure we're still OK to add data:
        if ( isGraphCollapsed ) {
            throw new GATKException("Attempted to addSequence to " + this.getClass().getName() + " after collapsing it.  This is a logic error.");
        }

        // Create a queue of nodes that represent the linear graph of all nodes in the given data:
        final Queue<AlignedBaseVertex> nodeQueue = createAlignedNodes(bases, cigar, contig, startPos, readName);

        if ( nodeQueue.peek() != null ) {
            final GenomicAndInsertionPosition firstNewNodePos = nodeQueue.peek().getPos();

            // Special handling for empty starting graphs:
            if ( !contigSubGraphMap.containsKey(contig) ) {
                // Initialize the graph:
                initializeGraphWithNodes(contig, nodeQueue);
            }
            else {
                // Merge the given node queue into the existing graph:
                mergeNodesIntoGraph(contig, nodeQueue);
            }

            // Collapse nodes in the graph that occur before this sequence if it's time:
//            if ( contigUncollapsedPositionMap.get(contig).getDistanceTo(firstNewNodePos) > periodicMergeDistance ) {
//                collapseEarlierNodes(contig, firstNewNodePos);
//            }

            ++numSequencesAdded;
        }
        else {
            logger.info( "No nodes to add in given read: " + readName + " (" + contig + ":" + startPos + ")" );
        }
    }

    /**
     * Collapses the graph for all nodes that occur earlier than the given position.
     * This prevents memory usage from going too high by consolidating nodes and removing them from the positional map.
     * @param contig Contig on which to collapse nodes (must not be {@code null}).
     * @param pos The start position of the read before which all nodes should be collapsed (can be {@code null}).
     */
    private void collapseEarlierNodes(final String contig, final GenomicAndInsertionPosition pos) {

        // If we have a null position, do nothing:
        if ( pos == null ) {
            return;
        }

        // This operation is the same as collapsing the whole graph, but we must do it for a range instead.
        // We do this to save memory, because this thing eats more memory than a humpback whale in a krill factory.
        final Map.Entry<GenomicAndInsertionPosition, Set<AlignedBaseVertex>> entry =
                contigPositionVertexMap.get(contig).floorEntry(pos);

        // Check to see if we have data to collapse:
        if ( entry.getKey().compareTo(pos) < 0 ) {
            // OK, we have to collapse some nodes...
            logger.info("Collapsing nodes before position: " + pos.getContig() + ":" + pos.getStart() + "_" + pos.getInsertionOffset());

            // 1 - Zip linear chains before the given pos:
            contigSubGraphMap.get(contig).zipLinearChainsBefore(pos);

            // 2 - Link adjacent nodes before the position:
            linkAdjacentNodesBefore(contigSubGraphMap.get(contig), pos);

            // 3 - Zip linear chains before the given pos (again):
            contigSubGraphMap.get(contig).zipLinearChainsBefore(pos);

            // 4 - Adjust the contig / position vertex map with data before the given position:
            updateContigPositionVertexMap(contig, pos);

            // 5 - Update the position of our first uncollapsed node:
            contigUncollapsedPositionMap.put(contig, pos);
        }
    }

    private void updateContigPositionVertexMap(final String contig, final GenomicAndInsertionPosition pos) {
        // 1 - remove all nodes from the tree that occur before pos:
        for ( final Map.Entry<GenomicAndInsertionPosition, Set<AlignedBaseVertex>> entry :
                contigPositionVertexMap.get(contig).entrySet()) {
            if ( entry.getKey().compareTo(pos) < 0 ) {
                contigPositionVertexMap.get(contig).remove(entry.getKey(), entry.getValue());
            }
        }

        // 2 - rebuild entries in contigPositionVertexMap for the newly created combined nodes:
        // TODO: FINISHME
        throw new GATKException("Not implemented yet.");
    }

    /**
     * Initializes the {@link #contigSubGraphMap} and {@link #contigPositionVertexMap} with the given node information.
     * @param contig Contig from which the given nodes are derived.
     * @param nodes {@link Queue<AlignedBaseVertex>} of nodes representing the aligned linear graph of nodes from a given read that should be added to our graph.
     */
    private void initializeGraphWithNodes(final String contig, final Queue<AlignedBaseVertex> nodes) {
        final AlignedBaseGraph graph = new AlignedBaseGraph(edgeType);

        // Set up our node position map:
        final TreeMap<GenomicAndInsertionPosition, Set<AlignedBaseVertex>> nodePositionMap = new TreeMap<>();
        contigPositionVertexMap.put(contig, nodePositionMap);

        AlignedBaseVertex lastVertex = null;

        if ( nodes.peek() != null ) {
            // Add the first position to our uncollapsed position map:
            contigUncollapsedPositionMap.put(contig, nodes.peek().getPos());
        }

        for ( final AlignedBaseVertex vertex : nodes ) {

            // Add the vertex to our graph
            graph.addVertex(vertex);

            // Add the link in the graph from the last vertex to this one:
            if ( lastVertex != null ) {
                graph.addEdge(lastVertex, vertex);
            }

            // Add the nodes to our node store:
            final Map<GenomicAndInsertionPosition, Set<AlignedBaseVertex>> posVtxMap =
                    contigPositionVertexMap.get(contig);
            posVtxMap.computeIfAbsent(vertex.getPos(), v -> new HashSet<>());
            posVtxMap.get(vertex.getPos()).add(vertex);

            // Store the last vertex:
            lastVertex = vertex;
        }

        // Add the graph to our contig graph map:
        contigSubGraphMap.put(contig, graph);
    }

    /**
     * Adds all given {@code nodes} into {@link #contigSubGraphMap}.
     * Assumes that {@code contig} exists in {@link #contigSubGraphMap} and {@link #contigPositionVertexMap}.
     * @param contig Contig from which the given nodes originate.
     * @param nodes {@link Queue<AlignedBaseVertex>} representing a linear graph of bases in a given read.
     */
    private void mergeNodesIntoGraph(final String contig, final Queue<AlignedBaseVertex> nodes) {
        // Get the graph we'll be working with:
        final AlignedBaseGraph graph = contigSubGraphMap.get(contig);

        // Get the position map we'll be dealing with:
        final TreeMap<GenomicAndInsertionPosition, Set<AlignedBaseVertex>> positionVertexMap =
                contigPositionVertexMap.get(contig);

        AlignedBaseVertex lastVertex = null;
        AlignedBaseVertex vertex;

        boolean lastNodeWasAddedToGraph = false;

        // Iterate through our nodes and find the best place to insert them:
        for ( final AlignedBaseVertex rawVertex : nodes ) {
            vertex = rawVertex;

            // Get the nodes in the graph that aligns to this one:
            final Set<AlignedBaseVertex> graphVertices = positionVertexMap.get(vertex.getPos());

            // Does any graph node have the same allele info as the current node?
            // If so, we don't add a new vertex to the graph and we set the vertex to the graph vertex to
            // preserve link continuity later:
            boolean vertexIsInGraph = false;
            if ( graphVertices != null ) {
                for ( final AlignedBaseVertex graphVertex : graphVertices ) {
                    if ( vertex.isSequenceEqual(graphVertex.getSequence()) ) {
                        // Setting vertex here will cause lastVertex to be set later.
                        vertex = graphVertex;
                        vertexIsInGraph = true;

                        // If we added our last node, we need to make a link from the newly-added last node to
                        // the current graph node.
                        if ( lastNodeWasAddedToGraph ) {
                            // This if statement can probably be removed:
                            if (lastVertex.equals(vertex)) {
                                final String message = "Equal Vertices detected: " + vertex.toString() + " @ " + vertex.getPos().toString();
                                logger.error(message);
                                throw new GATKException(message);
                            }
                            graph.addEdge(lastVertex, vertex);
                        }
                        lastNodeWasAddedToGraph = false;
                        break;
                    }
                }
            }

            // We must add in the new vertex and its edges:
            if ( !vertexIsInGraph ) {
                addNewNodeToGraph(graph, positionVertexMap, lastVertex, vertex);
                lastNodeWasAddedToGraph = true;
            }
            else if ( relabelEdgeTypes ) {
                relabelEdges(graph, lastVertex, vertex, rawVertex);
            }

            lastVertex = vertex;
        }
    }

    private void relabelEdges(final AlignedBaseGraph graph, final AlignedBaseVertex lastVertex, final AlignedBaseVertex vertex, final AlignedBaseVertex rawVertex) {
        // If we need to relable edges, we should do so:
        final BaseEdge edge = graph.getEdge(lastVertex, vertex);
        if ( edge != null ) {
            final LabeledEdge labeledEdge = (LabeledEdge)edge;
            final String      newLabel    = LabeledEdge.createLabel(rawVertex);
            if ( !labeledEdge.getLabel().equals(newLabel) ) {
                logger.debug("Relabeling edge: " + labeledEdge.getLabel() + " -> " + newLabel);
                labeledEdge.updateLabel(newLabel);
            }
        }
    }

    private void addNewNodeToGraph(final AlignedBaseGraph graph, final TreeMap<GenomicAndInsertionPosition, Set<AlignedBaseVertex>> positionVertexMap, final AlignedBaseVertex lastVertex, final AlignedBaseVertex vertex) {
        // Add our vertex to the graph:
        graph.addVertex(vertex);

        // Add the vertex to our positional data store:
        positionVertexMap.computeIfAbsent(vertex.getPos(), v -> new HashSet<>());
        positionVertexMap.get(vertex.getPos()).add(vertex);

        // Add an edge from the previous position if we must:
        if (lastVertex != null) {
            graph.addEdge(lastVertex, vertex);
        }
    }

    /**
     * Create a linear sequence of node objects based on aligned sequence information.
     * @param bases Bases to be added in the order in which they appear in the read.
     * @param cigar {@link List<CigarElement>} containing alignment information for each base in {@code bases}.
     * @param contig The contig to which the given base sequence is aligned.
     * @param startPos The start position of the first aligned base in the given base sequence.
     * @param readName The name of the parent read from which this aligned sequence information is derived.
     * @return A {@link Queue<AlignedBaseVertex>} representing a linear graph of bases in a given read.
     */
    private Queue<AlignedBaseVertex> createAlignedNodes(final byte[] bases, final List<CigarElement> cigar, final String contig, final int startPos, final String readName) {

        final Queue<AlignedBaseVertex> nodeQueue = new ArrayDeque<>();

        int readBasePos = 0;
        // Start at start-1 just in case we begin with Soft Clips or Insertions:
        int genomePos = startPos - 1;

        for ( final CigarElement cigarElement : cigar ) {
            final CigarOperator cigarOperator = cigarElement.getOperator();

            // Ignore hard clips and padding:
            if ( cigarOperator.equals(CigarOperator.HARD_CLIP) || cigarOperator.equals(CigarOperator.PADDING) ) {
                continue;
            }

            // Set up insertion position if we need to track inserted bases:
            int insertionOffset = 0;
            insertionOffset = updateInsertionOffset(insertionOffset, cigarOperator);

            for ( int i = 0; i < cigarElement.getLength(); ++i ) {
                if ( cigarOperator.consumesReferenceBases() ) {
                    ++genomePos;
                }

                // Add our vertex to the node queue:
                nodeQueue.add(
                        new AlignedBaseVertex(
                                new byte[] {bases[readBasePos]},
                                contig,
                                genomePos,
                                insertionOffset,
                                readName)
                );

                if ( cigarOperator.consumesReadBases() ) {
                    ++readBasePos;
                    insertionOffset = updateInsertionOffset(insertionOffset, cigarOperator);
                }
            }
        }

        return nodeQueue;
    }

    /**
     * Updates the given insertion offset based on the cigar operator.
     *
     * {@link CigarOperator#INSERTION} and {@link CigarOperator#SOFT_CLIP} should result in the returned value being
     * incremented by one.
     *
     * @param insertionOffset Initial insertion offset to update.
     * @param cigarOperator {@link CigarOperator} to check for updating the insertion offset.
     * @return An updated insertion offset based on the given {@link CigarOperator}.
     */
    private int updateInsertionOffset(final int insertionOffset, final CigarOperator cigarOperator) {
        if ( cigarOperator.equals(CigarOperator.INSERTION) ||
                cigarOperator.equals(CigarOperator.SOFT_CLIP) ) {
            return insertionOffset + 1;
        }
        return insertionOffset;
    }

    /**
     * Get the {@link CigarOperator} at a given position in a read based on a cigar represented by a {@link List<CigarElement>}.
     * @param cigarElementList {@link List<CigarElement>} representing a cigar alignment.
     * @param readBasePos Position in a read to query for a {@link CigarOperator}.
     * @return The {@link CigarOperator} at the given read position.
     */
    private CigarOperator getCigarOperatorAtReadBasePosition(final List<CigarElement> cigarElementList, final int readBasePos) {
        int cumulativePos = 0;
        for (final CigarElement element: cigarElementList) {
            if ( element.getOperator().consumesReadBases() ) {
                final int cigarElementReadEnd = element.getLength() + cumulativePos;

                if ( readBasePos < cigarElementReadEnd ) {
                    return element.getOperator();
                }
                cumulativePos += element.getLength();
            }
        }
        return cigarElementList.get(cigarElementList.size()-1).getOperator();
    }

    public boolean isRelabelEdgeTypes() {
        return relabelEdgeTypes;
    }

    public void setRelabelEdgeTypes(final boolean relabelEdgeTypes) {
        this.relabelEdgeTypes = relabelEdgeTypes;
    }

    //==================================================================================================================
    // Helper Data Types:

}
