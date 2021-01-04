package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import scala.Tuple2;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Graph representation for resolving structural variants.
 *
 * In this graph, nodes represent breakpoint positions in the reference. There are two classes of edges: breakpoint
 * edges and reference edges. A breakpoint edge is added for each piece of evidence. Reference edges, representing
 * reference sequence, are added between adjacent nodes.
 *
 * For example, consider the reference sequence
 *
 * ATTCGAGACGAAGG
 *
 * and with a deletion breakpoint at positions 3 and 10:
 *
 * ATT^CGAGACG^AAGG
 *
 * where the ^ indicate the breakpoint positions. Then there would be nodes at these positions, as well as the beginning and
 * end of the sequence (positions 0, 3, 10, and 14). A breakpoint edge b connects the nodes at 3 and 10:
 *
 *            ___b____
 *            |      |
 *            |      V
 *    (0)    (3)    (10)    (14)
 *
 * Note the deletion breakpoint edge is directed, starting at (3) and ending at (10), since the deletion carries the
 * sequence in that direction (left + strand, right - strand). In general, the direction of the edge is determined by
 * the strand of each end, flowing from + to -. For inversions (+/+ and -/-), either direction is possible. In that case,
 * the direction is undefined (and later determined during path enumeration).
 *
 * Next, reference edges are added between nodes 0 and 3 (representing r1 = ATT), 3 and 10 (r2 = CGAGACG), and 10 and 14 (r3 = AAGG):
 *
 *            ___b____
 *            |      |
 *            |      V
 *    (0)--->(3)--->(10)--->(14)
 *        r1     r2      r3
 *
 * By traversing this graph, one can generate sequences of r1, r2, and r3 corresponding to possible haplotypes:
 *
 *  r1,r2,r3 (ref)
 *  r1,r3    (del of r2)
 *
 * which can later by filtered using read depth information.
 */
public final class SVGraph {

    final List<IndexedSVGraphEdge> edges;
    final List<SVGraphNode> nodes;
    final Set<Integer> startingNodes;
    final Set<Integer> endingNodes;
    final SAMSequenceDictionary dictionary;
    final List<SVInterval> intervals;

    public SVGraph(final Collection<CoordinateSVGraphEdge> breakpointEdges, final SAMSequenceDictionary dictionary) {
        this.dictionary = dictionary;

        //Compute breakpoint positions
        final Map<Integer, List<Integer>> contigToBreakpointPositionsMap = getContigToBreakpointPositionsMap(breakpointEdges);

        //Compute reference edges between each pair of positions
        final Collection<CoordinateSVGraphEdge> referenceEdges = getReferenceEdges(contigToBreakpointPositionsMap);

        //Combine breakpoint and reference edges that form the entire graph
        final Collection<CoordinateSVGraphEdge> inputEdges = new ArrayList<>(breakpointEdges.size() + referenceEdges.size());
        inputEdges.addAll(breakpointEdges);
        inputEdges.addAll(referenceEdges);

        //Compute nodes and indexed edges
        nodes = generateNodes(inputEdges, dictionary);
        edges = generateIndexedEdgesAndUpdateNodes(inputEdges);

        if (hasInterchromosomalEdge()) {
            throw new IllegalArgumentException("Interchromosomal edges not currently supported");
        }

        //Find start and end nodes
        startingNodes = generateStartNodes();
        endingNodes = generateEndNodes();
        intervals = generateIntervals();
    }

    private boolean hasInterchromosomalEdge() {
        return edges.stream().anyMatch(edge -> nodes.get(edge.getNodeAIndex()).getContig() != nodes.get(edge.getNodeBIndex()).getContig());
    }

    public List<SVGraphNode> getNodes() {
        return nodes;
    }

    private static Stream<Tuple2<Integer, Integer>> getEdgeContigsAndPositionsStream(final CoordinateSVGraphEdge edge) {
        return Stream.of(new Tuple2<>(edge.getContigA(), edge.getNodeAPosition()), new Tuple2<>(edge.getContigB(), edge.getNodeBPosition()));
    }

    private static List<Integer> flattenAndSortUniqueValues(List<Tuple2<Integer, Integer>> tuplesList) {
        return tuplesList.stream().map(pair -> pair._2).distinct().sorted().collect(Collectors.toList());
    }

    private static Tuple2<Integer, List<Integer>> flattenAndSortContigPositionPairs(final Integer contig, final List<Tuple2<Integer, Integer>> contigPositionPairsList) {
        return new Tuple2<>(contig, flattenAndSortUniqueValues(contigPositionPairsList));
    }

    private static Map<Integer, List<Integer>> getContigToBreakpointPositionsMap(final Collection<CoordinateSVGraphEdge> edges) {
        return edges.stream()
                .flatMap(edge -> getEdgeContigsAndPositionsStream(edge)) //Get the two positions as contig-position pairs
                .collect(Collectors.groupingBy(Tuple2::_1))     //Group pairs by contig
                .entrySet().stream()
                .map(entry -> flattenAndSortContigPositionPairs(entry.getKey(), entry.getValue())) //Convert contig-position pairs to sorted lists of positions keyed by contig
                .collect(Collectors.toMap(Tuple2::_1, Tuple2::_2));
    }

    private Collection<CoordinateSVGraphEdge> getReferenceEdges(final Map<Integer, List<Integer>> contigToBreakpointPositionsMap) {
        final Collection<CoordinateSVGraphEdge> referenceEdges = new ArrayList<>();
        for (final Integer contig : contigToBreakpointPositionsMap.keySet()) {
            final List<Integer> positions = contigToBreakpointPositionsMap.get(contig);
            for (int i = 0; i < positions.size() - 1; i++) {
                final int position1 = contigToBreakpointPositionsMap.get(contig).get(i);
                final int position2 = contigToBreakpointPositionsMap.get(contig).get(i + 1);
                final CoordinateSVGraphEdge edge = new CoordinateSVGraphEdge(contig, position1, true, contig, position2, false, true, new SVGraphEdgeEvidence(), dictionary);
                referenceEdges.add(edge);
            }
        }
        return referenceEdges;
    }

    private static List<SVGraphNode> generateNodes(final Collection<CoordinateSVGraphEdge> edges, final SAMSequenceDictionary dictionary) {
        return edges.stream()
                .flatMap(edge -> Stream.of(new SVGraphNode(edge.getContigA(), edge.getNodeAPosition()), new SVGraphNode(edge.getContigB(), edge.getNodeBPosition())))
                .distinct()
                .sorted((a, b) -> SVIntervalUtils.compareIntervals(new SVInterval(a.getContig(), a.getPosition(), a.getPosition()), new SVInterval(b.getContig(), b.getPosition(), b.getPosition())))
                .collect(Collectors.toList());
    }

    private List<IndexedSVGraphEdge> generateIndexedEdgesAndUpdateNodes(final Collection<CoordinateSVGraphEdge> edges) {

        //Needed for converting edge types
        final Map<Tuple2<Integer,Integer>, Integer> contigAndPositionToNodeIndexMap = new HashMap<>(SVUtils.hashMapCapacity(nodes.size()));
        int nodeIndex = 0;
        for (final SVGraphNode node : nodes) {
            contigAndPositionToNodeIndexMap.put(new Tuple2<>(node.getContig(), node.getPosition()), nodeIndex);
            nodeIndex++;
        }

        //Generate indexed edges list from coordinate edges
        final List<IndexedSVGraphEdge> indexedEdges = new ArrayList<>(edges.size());
        for (final CoordinateSVGraphEdge edge : edges) {
            final int nodeAIndex = contigAndPositionToNodeIndexMap.get(new Tuple2<>(edge.getContigA(), edge.getNodeAPosition()));
            final int nodeBIndex = contigAndPositionToNodeIndexMap.get(new Tuple2<>(edge.getContigB(), edge.getNodeBPosition()));
            final int edgeIndex = indexedEdges.size();
            final IndexedSVGraphEdge newEdge = new IndexedSVGraphEdge(edgeIndex, nodeAIndex, nodeBIndex, edge.isStrandA(), edge.isStrandB(), edge.isReference(), edge.getEvidence(), this, dictionary);
            indexedEdges.add(newEdge);
            addIndexedEdgeToNodeEdgesList(newEdge, nodes.get(nodeAIndex), true);
            addIndexedEdgeToNodeEdgesList(newEdge, nodes.get(nodeBIndex), false);
        }
        return indexedEdges;
    }

    private static void addIndexedEdgeToNodeEdgesList(final IndexedSVGraphEdge edge, final SVGraphNode node, final boolean isNodeA) {
        if ((isNodeA && edge.isStrandA()) || (!isNodeA && edge.isStrandB())) {
            node.getOutEdges().add(edge);
        } else {
            node.getInEdges().add(edge);
        }
    }

    private Set<Integer> generateStartNodes() {
        return IntStream.range(0, nodes.size()).filter(i -> nodes.get(i).countReferenceInEdges() == 0).boxed().collect(Collectors.toSet());
    }

    private Set<Integer> generateEndNodes() {
        return IntStream.range(0, nodes.size()).filter(i -> nodes.get(i).countReferenceOutEdges() == 0).boxed().collect(Collectors.toSet());
    }

    public SAMSequenceDictionary getDictionary() {
        return dictionary;
    }

    public List<IndexedSVGraphEdge> getReferenceEdges() {
        return edges.stream().filter(edge -> edge.isReference()).distinct().collect(Collectors.toList());
    }

    private List<SVInterval> generateIntervals() {
        return getReferenceEdges().stream().map(edge -> edge.getInterval()).collect(Collectors.groupingBy(interval -> interval.getContig())).entrySet().stream()
                .map(entry -> new SVInterval(entry.getKey(), entry.getValue().stream().mapToInt(SVInterval::getStart).min().getAsInt(), entry.getValue().stream().mapToInt(SVInterval::getEnd).max().getAsInt()))
                .collect(Collectors.toList());
    }

    public List<SVInterval> getContigIntervals() {
        return intervals;
    }

    public List<IndexedSVGraphEdge> getEdges() {
        return edges;
    }

    public List<SVGraphNode> generateNodes() {
        return nodes;
    }

    public Set<Integer> getStartingNodes() {
        return startingNodes;
    }

    public Set<Integer> getEndingNodes() {
        return endingNodes;
    }


}