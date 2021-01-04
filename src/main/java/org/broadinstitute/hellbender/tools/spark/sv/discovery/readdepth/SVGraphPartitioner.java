package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.function.BiFunction;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Partitions SV graph into subgraphs using a provided edge dependence function
 */
public final class SVGraphPartitioner {

    private final List<IndexedSVGraphEdge> edges;
    private final List<SVGraphNode> nodes;
    private final SAMSequenceDictionary dictionary;

    public SVGraphPartitioner(final SVGraph graph) {
        edges = graph.getEdges();
        nodes = graph.generateNodes();
        dictionary = graph.getDictionary();
    }

    public List<SVGraph> getIndependentSubgraphs(final BiFunction<SVInterval, SVInterval, Boolean> edgeDependenceFunction) {
        //Tree must store collection of edges in case multiple edges have identical intervals
        final SVIntervalTree<Collection<IndexedSVGraphEdge>> breakpointEdgeTree = new SVIntervalTree<>();
        for (final IndexedSVGraphEdge edge : edges) {
            addEdgeToCollectionTree(edge, breakpointEdgeTree);
        }

        //Visit each edge and cluster with dependent edges
        final Collection<Set<IndexedSVGraphEdge>> subgraphEdgesCollection = new ArrayList<>();
        final Set<Integer> visitedEdges = new HashSet<>(SVUtils.hashMapCapacity(breakpointEdgeTree.size()));
        for (int i = 0; i < edges.size(); i++) {
            //Seed a new cluster at each unvisited breakpoint edge
            if (!visitedEdges.contains(i)) {
                final IndexedSVGraphEdge seedEdge = edges.get(i);
                if (!seedEdge.isReference()) {
                    subgraphEdgesCollection.add(generateSubgraphEdges(seedEdge, visitedEdges, breakpointEdgeTree, edgeDependenceFunction));
                }
            }
        }
        return convertEdgeCollectionsToGraphs(subgraphEdgesCollection);
    }

    private Stream<IndexedSVGraphEdge> getDependentEdgesStream(final SVInterval interval,
                                                               final SVIntervalTree<Collection<IndexedSVGraphEdge>> tree,
                                                               final BiFunction<SVInterval, SVInterval, Boolean> edgeDependenceFunction) {
        return Utils.stream(tree.overlappers(interval))
                .flatMap(entry -> entry.getValue().stream())
                .filter(overlappingEdge -> edgeDependenceFunction.apply(overlappingEdge.getInterval(), interval));
    }

    private Collection<IndexedSVGraphEdge> getNewDependentEdges(final Collection<IndexedSVGraphEdge> edges,
                                                                final SVIntervalTree<Collection<IndexedSVGraphEdge>> tree,
                                                                final BiFunction<SVInterval, SVInterval, Boolean> edgeDependenceFunction,
                                                                final Set<IndexedSVGraphEdge> currentEdges) {
        return edges.stream().map(edge -> edge.getInterval())
                .flatMap(interval -> getDependentEdgesStream(interval, tree, edgeDependenceFunction))
                .distinct()
                .filter(edge -> !currentEdges.contains(edge))
                .collect(Collectors.toList());
    }

    private Collection<CoordinateSVGraphEdge> convertIndexedEdges(final Collection<IndexedSVGraphEdge> edges) {
        return edges.stream().map(edge -> new CoordinateSVGraphEdge(edge, nodes)).collect(Collectors.toList());
    }

    private List<SVGraph> convertEdgeCollectionsToGraphs(final Collection<Set<IndexedSVGraphEdge>> edgesCollection) {
        return edgesCollection.stream().map(edgeSet -> convertIndexedEdges(edgeSet))
                .map(edgeSet -> new SVGraph(edgeSet, dictionary))
                .collect(Collectors.toList());
    }

    private final Set<IndexedSVGraphEdge> generateSubgraphEdges(final IndexedSVGraphEdge seedEdge,
                                                                final Set<Integer> visitedEdges,
                                                                final SVIntervalTree<Collection<IndexedSVGraphEdge>> tree,
                                                                final BiFunction<SVInterval, SVInterval, Boolean> edgeDependenceFunction) {
        visitedEdges.add(seedEdge.getIndex());
        final Set<IndexedSVGraphEdge> subgraphEdges = new HashSet<>();
        final Set<IndexedSVGraphEdge> newEdges = new HashSet<>();
        subgraphEdges.add(seedEdge);
        newEdges.add(seedEdge);

        //Agglomerate intersecting edges
        while (!newEdges.isEmpty()) {
            final Collection<IndexedSVGraphEdge> newDependentEdges = getNewDependentEdges(newEdges, tree, edgeDependenceFunction, subgraphEdges);
            newEdges.clear();
            newEdges.addAll(newDependentEdges);
            subgraphEdges.addAll(newEdges);
            final Collection<Integer> newEdgeIndices = newEdges.stream().map(IndexedSVGraphEdge::getIndex).collect(Collectors.toList());
            visitedEdges.addAll(newEdgeIndices);
        }
        return subgraphEdges;
    }

    private void addEdgeToCollectionTree(final IndexedSVGraphEdge edge, final SVIntervalTree<Collection<IndexedSVGraphEdge>> tree) {
        if (!edge.isReference()) {
            final SVInterval interval = edge.getInterval();
            final SVIntervalTree.Entry<Collection<IndexedSVGraphEdge>> entry = tree.find(interval);
            if (entry == null) {
                final Collection<IndexedSVGraphEdge> newEdgeCollection = new LinkedList<>();
                newEdgeCollection.add(edge);
                tree.put(interval, newEdgeCollection);
            } else {
                entry.getValue().add(edge);
            }
        }
    }
}
