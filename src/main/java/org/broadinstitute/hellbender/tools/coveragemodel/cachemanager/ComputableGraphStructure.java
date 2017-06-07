package org.broadinstitute.hellbender.tools.coveragemodel.cachemanager;

import avro.shaded.com.google.common.collect.Sets;
import org.broadinstitute.hellbender.utils.Utils;

import javax.annotation.Nonnull;
import java.io.Serializable;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * This class performs consistency checks and computes several structural properties for the DAG specified by a
 * set of {@link CacheNode}s. These include:
 *
 * - assertion for existence of no cycles
 * - construction of maps of nodes to their descendants and ancestors
 * - propagation of tags from descendants to ancestors (tags are instances of {@link CacheNode.NodeTag} that are used
 *   for marking one or more cache nodes that are required for a specific computation; for example, see the javadoc
 *   of {@link ImmutableComputableGraph} for a concrete use case)
 * - topological order for evaluating a computable node
 * - topological order for mutating a primitive/externally-computed node and updating the caches of all involved nodes
 * - topological order for evaluating all nodes associated to a tag (see {@link ImmutableComputableGraph})
 * - topological order for evaluating all nodes in the graph
 *
 * @implNote the implementation of some of the algorithms in this class, while being fast for small graphs, are not
 * optimized for large graphs and can run into {@link StackOverflowError}.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
final class ComputableGraphStructure implements Serializable {

    private static final long serialVersionUID = -3124293279477371159L;

    private final Set<CacheNode.NodeKey> nodeKeysSet;
    private final Set<CacheNode.NodeTag> nodeTagsSet;
    private final Map<CacheNode.NodeKey, Set<CacheNode.NodeTag>> inducedTagsMap;
    private final Map<CacheNode.NodeKey, Set<CacheNode.NodeKey>> descendantsMap;
    private final Map<CacheNode.NodeKey, Set<CacheNode.NodeKey>> childrenMap;
    private final Map<CacheNode.NodeKey, Set<CacheNode.NodeKey>> parentsMap;
    private final Map<CacheNode.NodeKey, Set<CacheNode.NodeKey>> ancestorsMap;
    private final Map<CacheNode.NodeKey, Integer> topologicalOrderMap;
    private final Map<CacheNode.NodeKey, List<CacheNode.NodeKey>> topologicalOrderForNodeEvaluation;
    private final Map<CacheNode.NodeKey, List<CacheNode.NodeKey>> topologicalOrderForNodeMutation;
    private final Map<CacheNode.NodeTag, List<CacheNode.NodeKey>> topologicalOrderForTagEvaluation;
    private final List<CacheNode.NodeKey> topologicalOrderForCompleteEvaluation;

    /**
     * An arbitrary negative number to denote the to-be-determined topological order a node
     */
    private static final int UNDEFINED_TOPOLOGICAL_ORDER = -1;

    /**
     * Package-private constructor from a set of {@link CacheNode}s. The graph is specified by the immediate
     * parents and descendants of each node. A {@link CyclicGraphException} is thrown of the graph has a cycle.
     *
     * @param nodeSet a set of {@link CacheNode}s
     */
    ComputableGraphStructure(@Nonnull final Set<CacheNode> nodeSet) {
        Utils.nonNull(nodeSet, "The given set of nodes must be non-null");
        nodeKeysSet = extractKeys(nodeSet);
        nodeTagsSet = extractTags(nodeSet);
        assertParentKeysExist(nodeSet, nodeKeysSet);

        parentsMap = getParentsMap(nodeSet);
        childrenMap = getChildrenMap(nodeSet);
        topologicalOrderMap = getTopologicalOrderMap(parentsMap, nodeKeysSet);
        final Map<Integer, Set<CacheNode.NodeKey>> nodesByTopologicalOrderMap =
                getNodesByTopologicalOrderMap(topologicalOrderMap, nodeKeysSet);
        descendantsMap = getDescendantsMap(childrenMap, nodesByTopologicalOrderMap);
        ancestorsMap = getAncestorsMap(parentsMap, nodesByTopologicalOrderMap);
        inducedTagsMap = getInducedTagsMap(nodeSet, nodesByTopologicalOrderMap, ancestorsMap);
        final Map<CacheNode.NodeTag, Set<CacheNode.NodeKey>> nodesByInducedTagMap =
                getNodesByInducedTagMap(inducedTagsMap, nodeKeysSet, nodeTagsSet);
        topologicalOrderForNodeEvaluation = getTopologicalOrderForNodeEvaluation(nodeKeysSet, topologicalOrderMap,
                ancestorsMap);
        topologicalOrderForTagEvaluation = getTopologicalOrderForTagEvaluation(nodeTagsSet, topologicalOrderMap,
                nodesByInducedTagMap, ancestorsMap);
        topologicalOrderForCompleteEvaluation = getTopologicalOrderForCompleteEvaluation(nodeKeysSet,
                topologicalOrderMap);
        topologicalOrderForNodeMutation = getTopologicalOrderForNodeMutation(nodeKeysSet, ancestorsMap,
                descendantsMap, topologicalOrderMap);
    }

    private static void assertParentKeysExist(@Nonnull final Set<CacheNode> nodeSet,
                                              @Nonnull final Set<CacheNode.NodeKey> nodeKeysSet) {
        for (final CacheNode node : nodeSet) {
            if (!nodeKeysSet.containsAll(node.getParents())) {
                final Set<CacheNode.NodeKey> undefinedParents = Sets.difference(new HashSet<>(node.getParents()), nodeKeysSet);
                throw new NonexistentParentNodeKey("Node " + ImmutableComputableGraphUtils.quote(node.getKey().toString()) +
                        " depends on undefined parent(s): " + undefinedParents.stream()
                        .map(CacheNode.NodeKey::toString)
                        .map(ImmutableComputableGraphUtils::quote).collect(Collectors.joining(", ")));
            }
        }
    }

    private static Set<CacheNode.NodeTag> extractTags(@Nonnull Set<CacheNode> nodeSet) {
        return nodeSet.stream().map(CacheNode::getTags).flatMap(Collection::stream).collect(Collectors.toSet());
    }

    private static Set<CacheNode.NodeKey> extractKeys(@Nonnull Set<CacheNode> nodeSet) {
        return nodeSet.stream().map(CacheNode::getKey).collect(Collectors.toSet());
    }

    /**
     * Sorts nodes by topological order using depth-first search. The output is a map from node keys to
     * their depth (root nodes have 0 depth).
     */
    private static Map<CacheNode.NodeKey, Integer> getTopologicalOrderMap(
            @Nonnull final Map<CacheNode.NodeKey, Set<CacheNode.NodeKey>> parentsMap,
            @Nonnull final Set<CacheNode.NodeKey> nodeKeysSet) {
        final Map<CacheNode.NodeKey, Integer> topologicalOrderMap = new HashMap<>();
        nodeKeysSet.forEach(key -> topologicalOrderMap.put(key, UNDEFINED_TOPOLOGICAL_ORDER));
        nodeKeysSet.forEach(nodeKey -> updateDepth(nodeKey, 0, nodeKeysSet, parentsMap, topologicalOrderMap));
        return topologicalOrderMap;
    }

    /**
     * Generates a (topological order -> set of node keys) map. The map clusters the nodes at the
     * same depth (root nodes have 0 depth).
     */
    private static Map<Integer, Set<CacheNode.NodeKey>> getNodesByTopologicalOrderMap(
            @Nonnull final Map<CacheNode.NodeKey, Integer> topologicalOrderMap,
            @Nonnull final Set<CacheNode.NodeKey> nodeKeysSet) {
        final int maxDepth = Collections.max(topologicalOrderMap.values());
        final Map<Integer, Set<CacheNode.NodeKey>> nodesByTopologicalOrderMap = new HashMap<>();
        IntStream.range(0, maxDepth + 1)
                .forEach(depth -> nodesByTopologicalOrderMap.put(depth,
                        nodeKeysSet.stream()
                                .filter(node -> topologicalOrderMap.get(node) == depth)
                                .collect(Collectors.toSet())));
        return nodesByTopologicalOrderMap;
    }

    /**
     * Generates a (nodeKey -> set of parents' keys) map.
     */
    private static Map<CacheNode.NodeKey, Set<CacheNode.NodeKey>> getParentsMap(@Nonnull final Set<CacheNode> nodeSet) {
        return nodeSet.stream()
                .collect(Collectors.toMap(CacheNode::getKey, node -> new HashSet<>(node.getParents())));
    }

    /**
     * Generates a (nodeKey -> set of children's keys) map.
     */
    private static Map<CacheNode.NodeKey, Set<CacheNode.NodeKey>> getChildrenMap(@Nonnull final Set<CacheNode> nodeSet) {
        final Map<CacheNode.NodeKey, Set<CacheNode.NodeKey>> childrenMap = nodeSet.stream()
                .collect(Collectors.toMap(CacheNode::getKey, node -> new HashSet<CacheNode.NodeKey>()));
        nodeSet.forEach(node -> node.getParents().forEach(parentKey -> childrenMap.get(parentKey)
                .add(node.getKey())));
        return childrenMap;
    }

    /**
     * Generates a (nodeKey -> set of upward-propagated tags) map.
     */
    private static Map<CacheNode.NodeKey, Set<CacheNode.NodeTag>> getInducedTagsMap(
            @Nonnull final Set<CacheNode> nodeSet,
            @Nonnull final Map<Integer, Set<CacheNode.NodeKey>> nodesByTopologicalOrderMap,
            @Nonnull final Map<CacheNode.NodeKey, Set<CacheNode.NodeKey>> allParentsMap) {
        final int maxDepth = Collections.max(nodesByTopologicalOrderMap.keySet());
        /* initialize with given tags */
        final Map<CacheNode.NodeKey, Set<CacheNode.NodeTag>> allTagsMap = nodeSet.stream()
                .collect(Collectors.toMap(CacheNode::getKey, node -> new HashSet<>(node.getTags())));
        /* propagate tags to all parents */
        for (int depth = maxDepth; depth >= 0; depth--) {
            nodesByTopologicalOrderMap.get(depth)
                    .forEach(nodeKey -> allParentsMap.get(nodeKey)
                            .forEach(parentKey -> allTagsMap.get(parentKey).addAll(allTagsMap.get(nodeKey))));
        }
        return allTagsMap;
    }

    /**
     * Generates a (tag -> set of tagged nodes, including upward-propagation) map.
     */
    private static Map<CacheNode.NodeTag, Set<CacheNode.NodeKey>> getNodesByInducedTagMap(
            @Nonnull final Map<CacheNode.NodeKey, Set<CacheNode.NodeTag>> allTagsMap,
            @Nonnull final Set<CacheNode.NodeKey> nodeKeysSet,
            @Nonnull final Set<CacheNode.NodeTag> nodeTagsSet) {
        final Map<CacheNode.NodeTag, Set<CacheNode.NodeKey>> nodesByTagMap = nodeTagsSet.stream()
                .collect(Collectors.toMap(Function.identity(), tag -> new HashSet<CacheNode.NodeKey>()));
        nodeKeysSet.forEach(nodeKey ->
                allTagsMap.get(nodeKey).forEach(tag ->
                        nodesByTagMap.get(tag).add(nodeKey)));
        return nodesByTagMap;
    }

    /**
     * Generates a (nodeKey -> set of descendants' keys) map.
     */
    private static Map<CacheNode.NodeKey, Set<CacheNode.NodeKey>> getDescendantsMap(
            @Nonnull final Map<CacheNode.NodeKey, Set<CacheNode.NodeKey>> childrenMap,
            @Nonnull final Map<Integer, Set<CacheNode.NodeKey>> nodesByTopologicalOrderMap) {
        final Map<CacheNode.NodeKey, Set<CacheNode.NodeKey>> descendantsMap = new HashMap<>();
        final int maxDepth = Collections.max(nodesByTopologicalOrderMap.keySet());
        /* deepest nodes have no descendants */
        nodesByTopologicalOrderMap.get(maxDepth).forEach(node -> descendantsMap.put(node, new HashSet<>()));
        /* get all descendants by ascending the tree */
        for (int depth = maxDepth - 1; depth >= 0; depth -= 1) {
            for (final CacheNode.NodeKey node : nodesByTopologicalOrderMap.get(depth)) {
                final Set<CacheNode.NodeKey> nodeDescendants = new HashSet<>();
                nodeDescendants.addAll(childrenMap.get(node));
                for (final CacheNode.NodeKey child : childrenMap.get(node)) {
                    nodeDescendants.addAll(descendantsMap.get(child));
                }
                descendantsMap.put(node, nodeDescendants);
            }
        }
        return descendantsMap;
    }

    /**
     * Generates a (nodeKey -> set of ancestors' keys) map.
     */
    private static Map<CacheNode.NodeKey, Set<CacheNode.NodeKey>> getAncestorsMap(
            @Nonnull final Map<CacheNode.NodeKey, Set<CacheNode.NodeKey>> parentsMap,
            @Nonnull final Map<Integer, Set<CacheNode.NodeKey>> nodesByTopologicalOrderMap) {
        final Map<CacheNode.NodeKey, Set<CacheNode.NodeKey>> ancestorsMap = new HashMap<>();
        final int maxDepth = Collections.max(nodesByTopologicalOrderMap.keySet());
        nodesByTopologicalOrderMap.get(0).forEach(node -> ancestorsMap.put(node, new HashSet<>()));
        for (int depth = 1; depth <= maxDepth; depth += 1) {
            for (final CacheNode.NodeKey node : nodesByTopologicalOrderMap.get(depth)) {
                final Set<CacheNode.NodeKey> nodeAncestors = new HashSet<>();
                nodeAncestors.addAll(parentsMap.get(node));
                for (final CacheNode.NodeKey parent : parentsMap.get(node)) {
                    nodeAncestors.addAll(ancestorsMap.get(parent));
                }
                ancestorsMap.put(node, nodeAncestors);
            }
        }
        return ancestorsMap;
    }

    /**
     * Generates the topological order for evaluating a single node as a (nodeKey -> list of topologically-ordered
     * node keys) map.
     */
    private static Map<CacheNode.NodeKey, List<CacheNode.NodeKey>> getTopologicalOrderForNodeEvaluation(
            @Nonnull final Set<CacheNode.NodeKey> nodeKeysSet,
            @Nonnull final Map<CacheNode.NodeKey, Integer> topologicalOrderMap,
            @Nonnull final Map<CacheNode.NodeKey, Set<CacheNode.NodeKey>> ancestorsMap) {
        final Map<CacheNode.NodeKey, List<CacheNode.NodeKey>> topologicalOrderForNodeEvaluation = new HashMap<>();
        for (final CacheNode.NodeKey nodeKey : nodeKeysSet) {
            final List<CacheNode.NodeKey> allParentsIncludingTheNode = new ArrayList<>();
            allParentsIncludingTheNode.addAll(ancestorsMap.get(nodeKey));
            allParentsIncludingTheNode.add(nodeKey);
            /* sort by depth */
            allParentsIncludingTheNode.sort(Comparator.comparingInt(topologicalOrderMap::get));
            topologicalOrderForNodeEvaluation.put(nodeKey, allParentsIncludingTheNode);
        }
        return topologicalOrderForNodeEvaluation;
    }

    /**
     * Generates the topological order for evaluating all nodes associated to a tag as a (nodeKey -> list of
     * topologically-ordered node keys) map.
     */
    private static Map<CacheNode.NodeTag, List<CacheNode.NodeKey>> getTopologicalOrderForTagEvaluation(
            @Nonnull final Set<CacheNode.NodeTag> nodeTagsSet,
            @Nonnull final Map<CacheNode.NodeKey, Integer> topologicalOrderMap,
            @Nonnull final Map<CacheNode.NodeTag, Set<CacheNode.NodeKey>> nodesByTagMap,
            @Nonnull final Map<CacheNode.NodeKey, Set<CacheNode.NodeKey>> ancestorsMap) {
        final Map<CacheNode.NodeTag, List<CacheNode.NodeKey>> topologicalOrderForTagEvaluation = new HashMap<>();
        for (final CacheNode.NodeTag tag : nodeTagsSet) {
            final Set<CacheNode.NodeKey> ancestorsIncludingTheTaggedNodesSet = new HashSet<>();
            for (final CacheNode.NodeKey node : nodesByTagMap.get(tag)) {
                ancestorsIncludingTheTaggedNodesSet.addAll(ancestorsMap.get(node));
                ancestorsIncludingTheTaggedNodesSet.add(node);
            }
            final List<CacheNode.NodeKey> ancestorsIncludingTheTaggedNodesList = new ArrayList<>();
            ancestorsIncludingTheTaggedNodesList.addAll(ancestorsIncludingTheTaggedNodesSet);
            /* sort by depth */
            ancestorsIncludingTheTaggedNodesList.sort(Comparator.comparingInt(topologicalOrderMap::get));
            topologicalOrderForTagEvaluation.put(tag, ancestorsIncludingTheTaggedNodesList);
        }
        return topologicalOrderForTagEvaluation;
    }

    /**
     * Generates topological order for evaluating all nodes in the graph as a (nodeKey -> list of topologically-ordered
     * node keys) map.
     */
    private static List<CacheNode.NodeKey> getTopologicalOrderForCompleteEvaluation(
            @Nonnull final Set<CacheNode.NodeKey> nodeKeysSet,
            @Nonnull final Map<CacheNode.NodeKey, Integer> topologicalOrderMap) {
        return new ArrayList<>(nodeKeysSet).stream()
                .sorted(Comparator.comparingInt(topologicalOrderMap::get))
                .collect(Collectors.toList());
    }

    /**
     * Generate topological order for mutating a primitive/externally-computed node and updating the caches of the
     * involved nodes. These include mutated node, its descendants, and the ancestors of all of its descendants.
     * The latter is required for updating the caches of the descendants.
     */
    private static Map<CacheNode.NodeKey, List<CacheNode.NodeKey>> getTopologicalOrderForNodeMutation(
            @Nonnull final Set<CacheNode.NodeKey> nodeKeysSet,
            @Nonnull final Map<CacheNode.NodeKey, Set<CacheNode.NodeKey>> ancestorsMap,
            @Nonnull final Map<CacheNode.NodeKey, Set<CacheNode.NodeKey>> descendantsMap,
            @Nonnull final Map<CacheNode.NodeKey, Integer> topologicalOrderMap) {
        final Map<CacheNode.NodeKey, List<CacheNode.NodeKey>> topologicalOrderForNodeMutation = new HashMap<>();
        for (final CacheNode.NodeKey nodeKey : nodeKeysSet) {
            final Set<CacheNode.NodeKey> allInvolvedNodesSet = new HashSet<>();
            allInvolvedNodesSet.add(nodeKey);
            allInvolvedNodesSet.addAll(descendantsMap.get(nodeKey));
            for (final CacheNode.NodeKey descendantNodeKey : descendantsMap.get(nodeKey)) {
                allInvolvedNodesSet.add(descendantNodeKey);
                allInvolvedNodesSet.addAll(ancestorsMap.get(descendantNodeKey));
            }
            final List<CacheNode.NodeKey> allInvolvedNodesList = new ArrayList<>();
            allInvolvedNodesList.addAll(allInvolvedNodesSet);
            allInvolvedNodesList.sort(Comparator.comparingInt(topologicalOrderMap::get));
            topologicalOrderForNodeMutation.put(nodeKey, allInvolvedNodesList);
        }
        return topologicalOrderForNodeMutation;
    }

    /**
     * Updates the depth of a node recursively
     */
    private static void updateDepth(@Nonnull final CacheNode.NodeKey nodeKey,
                                    final int recursion,
                                    @Nonnull Set<CacheNode.NodeKey> nodeKeysSet,
                                    @Nonnull Map<CacheNode.NodeKey, Set<CacheNode.NodeKey>> parentsMap,
                                    @Nonnull Map<CacheNode.NodeKey, Integer> topologicalOrderMap) {
        if (recursion > nodeKeysSet.size()) {
            throw new CyclicGraphException("The graph is not acyclic");
        }
        if (parentsMap.get(nodeKey).isEmpty()) {
            topologicalOrderMap.put(nodeKey, 0);
        } else if (topologicalOrderMap.get(nodeKey) == UNDEFINED_TOPOLOGICAL_ORDER) {
            parentsMap.get(nodeKey).forEach(parentNodeKey -> updateDepth(parentNodeKey,
                    recursion + 1, nodeKeysSet, parentsMap, topologicalOrderMap));
            final int maxParentDepth = parentsMap.get(nodeKey).stream()
                    .map(topologicalOrderMap::get)
                    .max(Integer::compareTo)
                    .get(); /* guaranteed to have a value */
            topologicalOrderMap.put(nodeKey, maxParentDepth + 1);
        }
        /* do nothing otherwise -- we already have the order for this node */
    }

    Set<CacheNode.NodeKey> getNodeKeysSet() { return nodeKeysSet; }

    Set<CacheNode.NodeTag> getNodeTagsSet() { return nodeTagsSet; }

    Set<CacheNode.NodeTag> getInducedTagsForNode(final CacheNode.NodeKey nodeKey) {
        return inducedTagsMap.get(nodeKey);
    }

    /**
     * Return the topological order for a node. Note: nodes at the same depth have the same topological order.
     * @param nodeKey node key
     * @return (integer) topological order
     */
    int getTopologicalOrder(@Nonnull final CacheNode.NodeKey nodeKey) {
        return topologicalOrderMap.get(nodeKey);
    }

    Set<CacheNode.NodeKey> getChildren(@Nonnull final CacheNode.NodeKey nodeKey) {
        return childrenMap.get(nodeKey);
    }

    Set<CacheNode.NodeKey> getParents(@Nonnull final CacheNode.NodeKey nodeKey) {
        return parentsMap.get(nodeKey);
    }

    Set<CacheNode.NodeKey> getAncestors(@Nonnull final CacheNode.NodeKey nodeKey) {
        return ancestorsMap.get(nodeKey);
    }

    Set<CacheNode.NodeKey> getDescendants(@Nonnull final CacheNode.NodeKey nodeKey) {
        return descendantsMap.get(nodeKey);
    }

    List<CacheNode.NodeKey> getTopologicalOrderForNodeEvaluation(final CacheNode.NodeKey nodeKey) {
        return topologicalOrderForNodeEvaluation.get(nodeKey);
    }

    List<CacheNode.NodeKey> getTopologicalOrderForNodeMutation(final CacheNode.NodeKey nodeKey) {
        return topologicalOrderForNodeMutation.get(nodeKey);
    }

    List<CacheNode.NodeKey> getTopologicalOrderForTagEvaluation(final CacheNode.NodeTag tagKey) {
        return topologicalOrderForTagEvaluation.get(tagKey);
    }

    List<CacheNode.NodeKey> getTopologicalOrderForCompleteEvaluation() {
        return topologicalOrderForCompleteEvaluation;
    }

    /**
     * This exception will be thrown if the graph has loops
     */
    static final class CyclicGraphException extends RuntimeException {
        private static final long serialVersionUID = 5887360871425098163L;

        CyclicGraphException(String s) {
            super(s);
        }
    }

    /**
     * This exception will be thrown if an alleged parent node key is missing
     */
    static final class NonexistentParentNodeKey extends RuntimeException {
        private static final long serialVersionUID = 586676245552229897L;

        NonexistentParentNodeKey(String s) {
            super(s);
        }
    }
}