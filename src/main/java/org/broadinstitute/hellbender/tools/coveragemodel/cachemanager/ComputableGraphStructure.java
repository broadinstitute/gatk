package org.broadinstitute.hellbender.tools.coveragemodel.cachemanager;

import org.broadinstitute.hellbender.utils.Utils;

import javax.annotation.Nonnull;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * This class computes a number of useful structural properties for DAG specified by the
 * immediate parents and immediate children of each node.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class ComputableGraphStructure {

    private final Set<String> nodeKeysSet;
    private final Set<String> nodeTagsSet;
    private final Map<String, Set<String>> allTagsMap;
    private final Map<String, Set<String>> immediateDescendentsMap;
    private final Map<String, Set<String>> immediateParentsMap;
    private final Map<String, Set<String>> allDescendentsMap;
    private final Map<String, Set<String>> allParentsMap;
    private final Map<String, Integer> depthsMap;
    private final Map<Integer, Set<String>> nodesByDepthMap;
    private final Map<String, Set<String>> nodesByTagMap;
    private final Map<String, List<String>> evalStrategyForNodeEvaluation;
    private final Map<String, List<String>> evalStrategyForNodeMutation;
    private final Map<String, List<String>> evalStrategyForTagEvaluation;
    private final List<String> evalStrategyForAllNodes;

    public ComputableGraphStructure(@Nonnull final Set<String> nodeKeysSet,
                                    @Nonnull final Map<String, Set<String>> initialTagsMap,
                                    @Nonnull final Map<String, Set<String>> immediateParentsMap,
                                    @Nonnull final Map<String, Set<String>> immediateDescendentsMap) {
        Utils.nonNull(nodeKeysSet, "The node key collection can not be null");
        Utils.nonNull(immediateParentsMap, "The immediate parents map can not be null");
        Utils.nonNull(initialTagsMap, "The initial tags map can not be null");
        Utils.nonNull(immediateDescendentsMap, "The immediate descendents map can not be null");
        /* check the descendents and parents */
        for (final String node : nodeKeysSet) {
            Utils.validateArg(immediateDescendentsMap.get(node) != null, "The immediate descendents node map" +
                    " has null values");
            Utils.validateArg(nodeKeysSet.containsAll(immediateDescendentsMap.get(node)), "The immediate descendents" +
                    " node map refers to unknown nodes.");
            Utils.validateArg(immediateParentsMap.get(node) != null, "The immediate parents node map has null values");
            Utils.validateArg(nodeKeysSet.containsAll(immediateParentsMap.get(node)), "The immediate parents node map" +
                    " refers to unknown nodes.");
        }

        this.nodeKeysSet = nodeKeysSet;
        this.immediateParentsMap = immediateParentsMap;
        this.immediateDescendentsMap = immediateDescendentsMap;

        /* create maps for descendents, parents, and depthsMap */
        allDescendentsMap = new HashMap<>();
        allParentsMap = new HashMap<>();
        depthsMap = new HashMap<>();
        allTagsMap = new HashMap<>();
        nodeKeysSet.forEach(key -> {
            allDescendentsMap.put(key, new HashSet<>());
            allParentsMap.put(key, new HashSet<>());
            allTagsMap.put(key, new HashSet<>());
            depthsMap.put(key, null);
        });

        /* calculate the depth of each node; primitive nodes or nodes with no immediateParentsMap have depth 0 */
        immediateParentsMap.keySet().stream()
                .filter(node -> immediateParentsMap.get(node).size() == 0)
                .forEach(rootNodeKey -> depthsMap.replace(rootNodeKey, 0));
        nodeKeysSet.forEach(this::updateDepth);

        final int maxDepth = Collections.max(depthsMap.values());

        /* list of nodes by their depthsMap */
        nodesByDepthMap = new HashMap<>();
        IntStream.range(0, maxDepth + 1).forEach(depth -> nodesByDepthMap.put(depth,
                nodeKeysSet.stream().filter(node -> depthsMap.get(node) == depth).collect(Collectors.toSet())));

        /* all descendents of the deepest nodes (empty set) */
        nodesByDepthMap.get(maxDepth).forEach(node -> allDescendentsMap.put(node, new HashSet<>()));
        /* get all descendents by descending the tree */
        for (int depth = maxDepth - 1; depth >= 0; depth -= 1) {
            for (final String node : nodesByDepthMap.get(depth)) {
                final Set<String> nodeAllDescendents = new HashSet<>();
                nodeAllDescendents.addAll(immediateDescendentsMap.get(node));
                for (final String child : immediateDescendentsMap.get(node)) {
                    nodeAllDescendents.addAll(allDescendentsMap.get(child));
                }
                allDescendentsMap.put(node, nodeAllDescendents);
            }
        }

        /* all parents of the primitive nodes (empty set) */
        nodesByDepthMap.get(0).forEach(node -> allParentsMap.put(node, new HashSet<>()));
        for (int depth = 1; depth <= maxDepth; depth += 1) {
            for (final String node : nodesByDepthMap.get(depth)) {
                final Set<String> nodeAllParents = new HashSet<>();
                nodeAllParents.addAll(immediateParentsMap.get(node));
                for (final String parent : immediateParentsMap.get(node)) {
                    nodeAllParents.addAll(allParentsMap.get(parent));
                }
                allParentsMap.put(node, nodeAllParents);
            }
        }

        /* build the full tags map; the parents inherit the tags of the descendents */
        nodesByDepthMap.get(maxDepth).forEach(node -> allTagsMap.get(node).addAll(initialTagsMap.get(node)));
        for (int depth = maxDepth - 1; depth >= 0; depth -= 1) {
            nodesByDepthMap.get(depth).forEach(node -> {
                allTagsMap.get(node).addAll(initialTagsMap.get(node));
                immediateDescendentsMap.get(node).forEach(desc -> allTagsMap.get(node).addAll(initialTagsMap.get(desc)));
            });
        }

        /* build a nodes-by-tag map and nodeTagsSet */
        nodeTagsSet = initialTagsMap.values().stream().flatMap(Set::stream).collect(Collectors.toSet());
        nodesByTagMap = new HashMap<>();
        nodeTagsSet.forEach(tag ->
                nodesByTagMap.put(tag, new HashSet<>()));
        nodeKeysSet.forEach(node ->
                allTagsMap.get(node).forEach(tag ->
                        nodesByTagMap.get(tag).add(node)));

        /* evaluation strategy for a single node */
        evalStrategyForNodeEvaluation = new HashMap<>();
        for (final String node : nodeKeysSet) {
            final List<String> allParentsIncludingTheNode = new ArrayList<>();
            allParentsIncludingTheNode.addAll(allParentsMap.get(node));
            allParentsIncludingTheNode.add(node);
            /* sort by depth */
            allParentsIncludingTheNode.sort(Comparator.comparingInt(depthsMap::get));
            evalStrategyForNodeEvaluation.put(node, allParentsIncludingTheNode);
        }

        /* evaluation strategy for all nodes associated to a tag */
        evalStrategyForTagEvaluation = new HashMap<>();
        for (final String tag : nodeTagsSet) {
            final Set<String> allParentsIncludingTheNodesSet = new HashSet<>();
            for (final String node : nodesByTagMap.get(tag)) {
                    allParentsIncludingTheNodesSet.addAll(allParentsMap.get(node));
                    allParentsIncludingTheNodesSet.add(node);
            }
            final List<String> allParentsIncludingTheNodesList = new ArrayList<>();
            allParentsIncludingTheNodesList.addAll(allParentsIncludingTheNodesSet);
            /* sort by depth */
            allParentsIncludingTheNodesList.sort(Comparator.comparingInt(depthsMap::get));
            evalStrategyForTagEvaluation.put(tag, allParentsIncludingTheNodesList);
        }

        /* evaluation strategy for all nodes */
        evalStrategyForAllNodes = new ArrayList<>();
        evalStrategyForAllNodes.addAll(nodeKeysSet);
        evalStrategyForAllNodes.sort(Comparator.comparingInt(depthsMap::get));

        /* evaluation strategy for updating the descendents of a mutated node */
        evalStrategyForNodeMutation = new HashMap<>();
        for (final String mutNode : nodeKeysSet) {
            final Set<String> allInvolvedSet = new HashSet<>();
            allInvolvedSet.add(mutNode);
            for (final String desc : allDescendentsMap.get(mutNode)) {
                allInvolvedSet.add(desc);
                allInvolvedSet.addAll(allParentsMap.get(desc));
            }
            final List<String> allInvolvedList = new ArrayList<>();
            allInvolvedList.addAll(allInvolvedSet);
            /* sort by depth */
            allInvolvedList.sort(Comparator.comparingInt(depthsMap::get));
            evalStrategyForNodeMutation.put(mutNode, allInvolvedList);
        }
    }

    /**
     * Updates the depth of a node recursively
     *
     * @param nodeKey the key of the node to update
     */
    private void updateDepth(final String nodeKey) {
        if (depthsMap.get(nodeKey) == null) {
            immediateParentsMap.get(nodeKey).forEach(this::updateDepth);
            final int depth = Collections.max(immediateParentsMap.get(nodeKey).stream().map(depthsMap::get)
                    .collect(Collectors.toList())) + 1;
            if (depth > nodeKeysSet.size() - 1) {
                throw new CyclicGraphException("The graph has cycles");
            }
            depthsMap.replace(nodeKey, depth);
        }
    }

    public Set<String> getNodeKeysSet() { return nodeKeysSet; }

    public Set<String> getNodeTagsSet() { return nodeTagsSet; }

    public Map<String, Set<String>> getAllTagsMap() { return allTagsMap; }

    public Map<String, Set<String>> getImmediateDescendentsMap() { return immediateDescendentsMap; }

    public Map<String, Set<String>> getImmediateParentsMap() { return immediateParentsMap; }

    public Map<String, Set<String>> getAllDescendentsMap() { return allDescendentsMap; }

    public Map<String, Set<String>> getAllParentsMap() { return allParentsMap; }

    public Map<String, Integer> getDepthsMap() { return depthsMap; }

    public Map<Integer, Set<String>> getNodesByDepthMap() { return nodesByDepthMap; }

    public Map<String, Set<String>> getNodesByTagMap() { return  nodesByTagMap; }

    public List<String> getEvalStrategyByEvaluatedNode(final String nodeKey) {
        return evalStrategyForNodeEvaluation.get(nodeKey);
    }

    public List<String> getEvalStrategyByMutatedNode(final String nodeKey) {
        return evalStrategyForNodeMutation.get(nodeKey);
    }

    public List<String> getEvalStrategyByTag(final String tagKey) {
        return evalStrategyForTagEvaluation.get(tagKey);
    }

    public List<String> getEvalStrategyForAllNodes() {
        return evalStrategyForAllNodes;
    }

    public String statusToString() {
        String status = "";
        for (final String nodeKey : nodeKeysSet) {
            status += "node: " + nodeKey + "\n" +
                    "\tdepth: " + depthsMap.get(nodeKey) + "\n" +
                    "\timmediate parents: " +
                    immediateParentsMap.get(nodeKey).stream().collect(Collectors.joining(", ", "[", "]\n")) +
                    "\timmediate descendents: " +
                    immediateDescendentsMap.get(nodeKey).stream().collect(Collectors.joining(", ", "[", "]\n")) +
                    "\tall parents: " +
                    allParentsMap.get(nodeKey).stream().collect(Collectors.joining(", ", "[", "]\n")) +
                    "\tall descendents: " +
                    allDescendentsMap.get(nodeKey).stream().collect(Collectors.joining(", ", "[", "]\n")) +
                    "\tall tags: " +
                    allTagsMap.get(nodeKey).stream().collect(Collectors.joining(", ", "[", "]\n"));
        }

        status += "\n";
        for (final String tag : nodeTagsSet) {
            status += "tag: " + tag + ", nodes: " +
                    nodesByTagMap.get(tag).stream().collect(Collectors.joining(", ", "[", "]\n"));
        }

        status += "\n";
        for (final String tag : nodeTagsSet) {
            status += "eval strategy per tag: " + tag + ", nodes in depth order: " +
                    evalStrategyForTagEvaluation.get(tag).stream().
                            map(nodeKey -> nodeKey + "(" + depthsMap.get(nodeKey) + ")").
                            collect(Collectors.joining(", ", "[", "]\n"));
        }

        status += "\n";
        for (final String node : nodeKeysSet) {
            status += "eval strategy for node evaluation: " + node + ", nodes in depth order: " +
                    evalStrategyForNodeEvaluation.get(node).stream().
                            map(nodeKey -> nodeKey + "(" + depthsMap.get(nodeKey) + ")").
                            collect(Collectors.joining(", ", "[", "]\n"));
        }

        status += "\n";
        for (final String node : nodesByDepthMap.get(0)) {
            status += "eval strategy for node mutation: " + node + ", nodes in depth order: " +
                    evalStrategyForNodeMutation.get(node).stream().
                            map(nodeKey -> nodeKey + "(" + depthsMap.get(nodeKey) + ")").
                            collect(Collectors.joining(", ", "[", "]\n"));
        }
        return status;
    }

    /**
     * This exception will be thrown if the graph has loops
     */
    public final class CyclicGraphException extends RuntimeException {
        private static final long serialVersionUID = 5887360871425098163L;

        public CyclicGraphException(String s) {
            super(s);
        }
    }
}