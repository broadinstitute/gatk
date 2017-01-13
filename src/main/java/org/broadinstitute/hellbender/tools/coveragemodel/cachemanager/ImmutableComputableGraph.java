package org.broadinstitute.hellbender.tools.coveragemodel.cachemanager;

import org.broadinstitute.hellbender.utils.Utils;

import javax.annotation.Nonnull;
import javax.annotation.Nullable;
import java.io.Serializable;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * This class provides a general purpose framework for evaluating functions on directed acyclic graphs consisting of
 * "primitive" and "computable" nodes. Primitive nodes {@link PrimitiveCacheNode} are value placeholders while computable
 * nodes {@link ComputableCacheNode} evaluate a function on the graph. Importantly, computable nodes may cache their
 * evaluations to avoid redundant expensive computations.
 *
 * This class implements a number of automatic bookkeeping strategies for caching the results of computable nodes
 * and updating the status of the cached values after mutating the primitive nodes.
 *
 * The nodes may store any object that implements {@link Duplicable}. These objects provide a recipe for
 * making a deep copy of the value(s) that they hold by implementing {@link Duplicable#deepCopy()}, and
 * {@link Duplicable#isNull()} to indicate whether the object holds any null pointers.
 *
 * Typical use case: evaluating computationally expensive expressions with common subexpressions.
 *
 * Example: let X, Y and Z be three mutable values (primitive nodes) and we want to calculate f(X,Y) and g(f(X,Y),Z).
 * To this end, we may first calculate Q_1 = f(X,Y) and proceed to calculate Q_2 = g(f(X,Y),Z). By caching the value of
 * Q_1, we save ourselves recomputing the subexpression f(X,Y) everytime we mutate the primitive value Z. Graphically,
 * the computation can be represented as a level-ordered directed acylic graph (the edges are assumed to have downward
 * arrows):
 *
 *                X      Y   Z    (depth 0)
 *                 \    /   /
 *                  \  /   /
 *                   Q_1  /       (depth 1)
 *                    \  /
 *                     \/
 *                    Q_2         (depth 2)
 *
 * One must manually identify the common subexpressions and assign a label to each. In the future, this can be
 * streamlined using a CAS library. The evaluation scheme can be coveninently set up using
 * {@link ImmutableComputableGraphBuilder} and it's two main methods:
 * {@link ImmutableComputableGraphBuilder#addPrimitiveNode(String, String[], Duplicable)}, and
 * {@link ImmutableComputableGraphBuilder#addComputableNode(String, String[], String[], Function, boolean)}.
 *
 * For computable nodes, one must provide the function, a list of immediate parent nodes, a list of tags,
 * and whether or not the value is to be cached. Computable nodes come in 3 species depending on the way
 * they are constructed:
 *
 *      (1) Caching: these nodes may store the values they evaluate for future lookup
 *      (2) Non-caching: these nodes are compute-on-demand
 *      (3) Externally set: these nodes are constructed using a {@code null} evaluation function. These user is
 *          responsible for evaluating these nodes when required and using {@link #setValue(String, Duplicable)}
 *          to update them. This class performs the bookkeeping of the status of their values (up to date or
 *          out of date) based on the provided parents list.
 *
 * One require the caching computable nodes to be updated and cached automatically after each mutation of a primitive
 * or externally set computable node by invoking {@link ImmutableComputableGraphBuilder#enableCacheAutoUpdate()}.
 * This is feature, however, is NOT recommended as one may not need all the caches to be up-to-date at all times
 * (see below for updating caches selectively). This class with throw an {@link IllegalStateException} if an old
 * cache is invoked in order to notify the user to update the cache manually.
 *
 * Updating caches:
 * ================
 *
 * If cache auto update is not enabled, the user is responsible for updating the caches by calling
 * either {@link #updateAllCaches()}, {@link #updateCachesForNode(String)}, or {@link #updateCachesForTag(String)}.
 *
 * Querying:
 * =========
 *
 * The graph is queried either by calling {@link #getValueDirect(String)} or
 * {@link #getValueWithRequiredEvaluations(String)}. The former only fetches the values and throws an exception
 * if the some of the required caches are out of date, or a non-caching computable node is encountered.
 * The latter performs the required evaluations along the way.
 *
 * IMPORTANT NOTE: the queried values are return by reference. It is the user's responsibility not to mutate
 * them. Otherwise, the functional structure will be broken.
 *
 * Mutation:
 * =========
 *
 * The primitive values can be mutated by calling {@link ImmutableComputableGraph#setValue(String, Duplicable)}.
 * The mutations do not occur in place; rather, a new instance of {@link ImmutableComputableGraph} is created along
 * with new instances for the updated nodes. If cache auto update is enabled, the affected nodes will be evaluated.
 * Otherwise, only the cache status will go out of date. Unchanged nodes are passed as reference to the new instance.
 * JVM's garbage collector will free up the memory for old cached nodes in sequential computations.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class ImmutableComputableGraph implements Serializable {

    private static final long serialVersionUID = -1162776031416105027L;

    private static Map<String, Duplicable> EMPTY_MAP = new HashMap<>();

    /**
     * A simple builder class for {@link ImmutableComputableGraph}
     */
    public static class ImmutableComputableGraphBuilder {
        private final Collection<CacheNode> nodes;
        private boolean cacheAutoUpdate;

        ImmutableComputableGraphBuilder() {
            nodes = new ArrayList<>();
            cacheAutoUpdate = false;
        }

        public ImmutableComputableGraphBuilder addPrimitiveNode(@Nonnull final String key,
                                                                @Nonnull final String[] tags,
                                                                @Nonnull Duplicable value) {
            nodes.add(new PrimitiveCacheNode(key, Arrays.stream(tags).collect(Collectors.toList()), value));
            return this;
        }

        public ImmutableComputableGraphBuilder addComputableNode(@Nonnull final String key,
                                                                 @Nonnull final String[] tags,
                                                                 @Nonnull final String[] parents,
                                                                 @Nullable final Function<Map<String, ? extends Duplicable>, ? extends Duplicable> func,
                                                                 final boolean cacheEvals) {
            nodes.add(new ComputableCacheNode(key,
                    Arrays.stream(tags).collect(Collectors.toList()),
                    Arrays.stream(parents).collect(Collectors.toList()),
                    func, cacheEvals));
            return this;
        }

        public ImmutableComputableGraphBuilder enableCacheAutoUpdate() {
            cacheAutoUpdate = true;
            return this;
        }

        public ImmutableComputableGraphBuilder disableCacheAutoUpdate() {
            cacheAutoUpdate = false;
            return this;
        }

        public ImmutableComputableGraph build() {
            if (nodes.size() == 0) {
                throw new IllegalStateException("Can not make an empty cache node collection");
            } else {
                return new ImmutableComputableGraph(nodes, cacheAutoUpdate);
            }
        }
    }

    public static ImmutableComputableGraphBuilder builder() {
        return new ImmutableComputableGraphBuilder();
    }

    /******************************************************************************************************************/

    private final Map<String, CacheNode> nodesMap;
    private final boolean cacheAutoUpdate;
    private final ComputableGraphStructure cgs;

    /**
     * Constructor
     * @param nodesCollection a collection of {@link CacheNode} subclasses
     */
    public ImmutableComputableGraph(@Nonnull final Collection<? extends CacheNode> nodesCollection,
                                    final boolean cacheAutoUpdate) {
        Utils.nonNull(nodesCollection, "The nodes collection must be non null.");
        this.cacheAutoUpdate = cacheAutoUpdate;
        final Set<String> nodeKeys = nodesCollection.stream()
                .map(CacheNode::getKey)
                .collect(Collectors.toSet());
        if (nodeKeys.size() != nodesCollection.size()) {
            final Set<String> badNodes = nodesCollection.stream()
                    .map(CacheNode::getKey)
                    .collect(Collectors.groupingBy(Function.identity(), Collectors.counting()))
                    .entrySet().stream()
                    .filter(entry -> entry.getValue() > 1)
                    .map(Map.Entry::getKey)
                    .collect(Collectors.toSet());

            throw new IllegalArgumentException("The node list contains nodes with repeated keys: " +
                badNodes.stream().collect(Collectors.joining(", ")));
        }
        nodesMap = nodesCollection.stream().collect(Collectors.toMap(CacheNode::getKey, Function.identity()));

        /* create maps for descendents, parents, and tags */
        final Map<String, Set<String>> immediateDescendents = new HashMap<>();
        final Map<String, Set<String>> immediateParents = new HashMap<>();
        final Map<String, Set<String>> initialTags = new HashMap<>();
        nodeKeys.forEach(key -> {
            immediateDescendents.put(key, new HashSet<>());
            immediateParents.put(key, new HashSet<>());
            initialTags.put(key, new HashSet<>());
        });

        /* immediate parents, descendents and tags */
        nodeKeys.stream().forEach(node -> {
            nodesMap.get(node).getParents().stream().forEach(parent -> immediateDescendents.get(parent).add(node));
            immediateParents.get(node).addAll(nodesMap.get(node).getParents());
            initialTags.get(node).addAll(nodesMap.get(node).getTags());
        });

        /* build the graph */
        cgs = new ComputableGraphStructure(nodeKeys, initialTags, immediateParents, immediateDescendents);
    }

    /**
     * A private constructor (used by duplicators)
     * @param nodesMap a previously constructed key -> {@link CacheNode} map
     * @param cgs a previously constructed {@link ComputableGraphStructure}
     */
    private ImmutableComputableGraph(@Nonnull final Map<String, CacheNode> nodesMap,
                                     @Nonnull final ComputableGraphStructure cgs,
                                     final boolean cacheAutoUpdate) {
        this.nodesMap = nodesMap;
        this.cgs = cgs;
        this.cacheAutoUpdate = cacheAutoUpdate;
    }

    /**
     * Sets the value of a primitive cache node
     * @param nodeKey the key of the node
     * @param newValue the new value of the node (note: it will not be duplicated; the user is responsible for duplicating)
     * @return a new instace of {@link ImmutableComputableGraph} with refernce to unchanged nodes and duplicated
     *         changed nodes
     * @throws IllegalArgumentException if the node does not exist
     * @throws UnsupportedOperationException if the node is non-primitive
     */
    public ImmutableComputableGraph setValue(@Nonnull final String nodeKey,
                                             @Nonnull final Duplicable newValue)
            throws IllegalArgumentException, UnsupportedOperationException {
        assertNodeExists(Utils.nonNull(nodeKey, "The key of the node can not be null."));
        CacheNode node = nodesMap.get(nodeKey);
        if (!node.isExternallyMutable()) {
            throw new UnsupportedOperationException("Can not explicitly set the value of a non-primitive cache node.");
        }
        final Map<String, CacheNode> updatedNodesMap = new HashMap<>();
        updatedNodesMap.put(nodeKey, node.duplicateWithUpdatedValue(newValue));
        final ImmutableComputableGraph out = duplicateWithUpdatedNodes(
                addDuplicateOfOutdatedDescendents(nodeKey, updatedNodesMap));
        if (cacheAutoUpdate) {
            Map<String, Duplicable> accumulatedValues = out.accumulateEvaluationsForDepthSortedNodeList(
                    cgs.getEvalStrategyByMutatedNode(nodeKey), true);
            return out.updateCachesFromAccumulatedValues(accumulatedValues);
        } else {
            return out;
        }
    }

    /**
     * Make a key -> node map containing new instances of the nodes which become out of date as a matter of
     * updating node {@code key}
     * @param key key of the updated note
     * @param updatedNodesMap a key -> node map
     */
    private Map<String, CacheNode> addDuplicateOfOutdatedDescendents(@Nonnull final String key,
                                                                     @Nonnull final Map<String, CacheNode> updatedNodesMap) {
        for (final String child : cgs.getAllDescendentsMap().get(key)) {
            CacheNode oldChild = nodesMap.get(child);
            /* all of the desendents are computable nodes and can be safely upcasted */
            updatedNodesMap.put(child, ((ComputableCacheNode)oldChild).duplicateWithUpdatedCacheStatus(false));
        }
        return updatedNodesMap;
    }

    /**
     * Returns a reference to the value of a given node
     *
     * Note: this function is purposefully meant to be light:
     *
     * (1) it does not update out-of-date cacheable computable nodes, and
     * (2) it does not evaluate non-cacheable computable nodes.
     *
     * @param nodeKey the key of the node
     * @return value of the node
     * @throws IllegalStateException if a cached node is out of date, or if the value of a primitive node is
     *                               not initialized
     * @throws IllegalArgumentException if the node does not exist
     */
    public Duplicable getValueDirect(@Nonnull final String nodeKey)
            throws IllegalStateException, IllegalArgumentException {
        assertNodeExists(Utils.nonNull(nodeKey, "The key of the node can not be null."));
        return nodesMap.get(nodeKey).getValue(EMPTY_MAP);
    }

    /**
     *
     * @param nodeKey key of the node
     * @return value of the node
     * @throws IllegalArgumentException if the node does not exist
     */
    public Duplicable getValueWithRequiredEvaluations(@Nonnull final String nodeKey)
            throws IllegalStateException, IllegalArgumentException {
        assertNodeExists(Utils.nonNull(nodeKey, "The key of the node can not be null."));
        return accumulateEvaluationsForDepthSortedNodeList(cgs.getEvalStrategyByEvaluatedNode(nodeKey), true).get(nodeKey);
    }

    /**
     * This method calculates and accumulates the values of a list of depth-sorted nodes (from low depth to high depth).
     * Compute-on-demand nodes are always evaluated. Caching nodes are only evaluated if {@code evalOnOutOfDateCachingNodes}
     * is true.
     *
     * @param depthSortedNodes depth-sorted list of nodes
     * @param evalOutdatedCachingNodes evaluate
     * @return
     */
    private Map<String, Duplicable> accumulateEvaluationsForDepthSortedNodeList(@Nonnull final List<String> depthSortedNodes,
                                                                                final boolean evalOutdatedCachingNodes) {
        final Map<String, Duplicable> accumulatedValues = new HashMap<>();
        for (final String nodeKey : depthSortedNodes) {
            if (nodesMap.get(nodeKey).isPrimitive()) {
                Duplicable value;
                try {
                    value = nodesMap.get(nodeKey).getValue(accumulatedValues);
                } catch (IllegalStateException e) {
                    value = null;
                }
                accumulatedValues.put(nodeKey, value);
            } else {
                final ComputableCacheNode computableNode = ((ComputableCacheNode)nodesMap.get(nodeKey));
                Duplicable value;
                if (computableNode.cacheEvals()) {
                    if (!evalOutdatedCachingNodes) { /* directly fetch the cache; may throw an exception */
                        value = computableNode.getValue(accumulatedValues);
                    } else if (computableNode.isStoredValueAvailableAndCurrent()) { /* cache is current; fetch */
                        value = computableNode.getValue(EMPTY_MAP);
                    } else { /* caches are old; try to evaluate */
                        try {
                            value = computableNode.getFunction().apply(accumulatedValues);
                        } catch (NullPointerException e) { /* may occur of some of parents are null */
                            value = null;
                        }
                    }
                } else { /* compute-on-demand node */
                    try {
                        value = computableNode.getFunction().apply(accumulatedValues);
                    } catch (RuntimeException e) { /* may occur of some of the parents are null */
                        value = null;
                    }
                }
                accumulatedValues.put(nodeKey, value);
            }
        }
        return accumulatedValues;
    }

    /**
     * Check {@code accumulatedValues} for truly affected nodes, create new node instances for affected nodes,
     * and make a new instance of {@link ImmutableComputableGraph} with reference to unaffected nodes, and new nodes
     *
     * @param accumulatedValues a nodekey -> duplicable map for possibly affected nodes
     * @return a new instance of {@link ImmutableComputableGraph} with new instances of updated nodes
     */
    private ImmutableComputableGraph updateCachesFromAccumulatedValues(final Map<String, Duplicable> accumulatedValues) {
        /* since accumulatedValues may contain unchanged values (by reference), we filter and only update
         * the affected nodes */
        return duplicateWithUpdatedNodes(
                accumulatedValues.keySet().stream()
                        /* filter out primitives and caching nodes that are current */
                        .filter(node -> !nodesMap.get(node).isPrimitive() &&
                                ((ComputableCacheNode)nodesMap.get(node)).cacheEvals() &&
                                !((ComputableCacheNode)nodesMap.get(node)).isCacheCurrent())
                        /* collect to a map: key -> duplicated node with updated value */
                        .collect(Collectors.toMap(Function.identity(), node ->
                                ((ComputableCacheNode)nodesMap.get(node))
                                        .duplicateWithUpdatedValue(accumulatedValues.get(node)))));
    }

    /**
     * Update the cached values by node key
     * @return a new instance of {@link ImmutableComputableGraph} with new instances of updated nodes
     */
    public ImmutableComputableGraph updateCachesForNode(@Nonnull final String nodeKey) {
        assertNodeExists(Utils.nonNull(nodeKey, "The key of the node can not be null."));
        final Map<String, Duplicable> accumulatedValues = accumulateEvaluationsForDepthSortedNodeList(
                cgs.getEvalStrategyByEvaluatedNode(nodeKey), true);
        return updateCachesFromAccumulatedValues(accumulatedValues);
    }

    /**
     * Update the cached values by tag
     * @return a new instance of {@link ImmutableComputableGraph} with new instances of updated nodes
     */
    public ImmutableComputableGraph updateCachesForTag(final String tagKey) {
        assertTagExists(Utils.nonNull(tagKey, "The key of the tag can not be null."));
        final Map<String, Duplicable> accumulatedValues = accumulateEvaluationsForDepthSortedNodeList(
                cgs.getEvalStrategyByTag(tagKey), true);
        return updateCachesFromAccumulatedValues(accumulatedValues);
    }

    /**
     * Update all caches values
     * @return a new instance of {@link ImmutableComputableGraph} with new instances of updated nodes
     */
    public ImmutableComputableGraph updateAllCaches() {
        final Map<String, Duplicable> accumulatedValues = accumulateEvaluationsForDepthSortedNodeList(
                cgs.getEvalStrategyForAllNodes(), true);
        return updateCachesFromAccumulatedValues(accumulatedValues);
    }

    /**
     * Make a new instance of {@link ImmutableComputableGraph} by replacing the nodes from {@code updatedNodesMap}
     *
     * @param updatedNodesMap nodes to be replaced and their new values
     * @return a new instance of {@link ImmutableComputableGraph} with new instances of updated nodes
     */
    private ImmutableComputableGraph duplicateWithUpdatedNodes(final Map<String, CacheNode> updatedNodesMap) {
        final Map<String, CacheNode> newNodesMap = new HashMap<>();
        /* nodes that are present in the new map */
        cgs.getNodeKeysSet().stream()
                .filter(node -> !updatedNodesMap.keySet().contains(node))
                .forEach(node -> newNodesMap.put(node, nodesMap.get(node)));
        newNodesMap.putAll(updatedNodesMap);
        return new ImmutableComputableGraph(newNodesMap, cgs, cacheAutoUpdate);
    }

    private void assertNodeExists(final String nodeKey) {
        Utils.validateArg(cgs.getNodeKeysSet().contains(nodeKey), "The node \"" + nodeKey + "\" does not exist.");
    }

    private void assertTagExists(final String tagKey) {
        Utils.validateArg(cgs.getNodeTagsSet().contains(tagKey), "The tag \"" + tagKey + "\" does not exist.");
    }

    /**
     * A text output informing about the status of the nodes
     * @return a string
     */
    public String statusToString(final boolean showValues) {
        String out = "";
        for (final String key : cgs.getNodeKeysSet()) {
            out += "node: " + key + "\n";
            out += "\ttype: " + (nodesMap.get(key).isPrimitive() ? "primitive" : "computable") + "\n";
            String value;
            try {
                Duplicable dup = getValueDirect(key);
                if (dup instanceof DuplicableNDArray) {
                    value = ((DuplicableNDArray)dup).value().toString();
                } else if (dup instanceof DuplicableNumber) {
                    value = ((DuplicableNumber)dup).value().toString();
                } else {
                    value = dup.toString();
                }
                if (!showValues) {
                    value = "is available";
                }
            } catch (IllegalStateException | NullPointerException e) {
                value = "not available";
            }
            out += "\tvalue: " + value + "\n";
            if (!nodesMap.get(key).isPrimitive()) {
                out += "\tup to date: " + ((ComputableCacheNode)nodesMap.get(key)).isStoredValueAvailableAndCurrent() + "\n";
                out += "\tcaches values: " + ((ComputableCacheNode)nodesMap.get(key)).cacheEvals() + "\n";
                out += "\thas a value: " + nodesMap.get(key).isStoredValueAvailable() + "\n";
            }
        }
        return out;
    }

    /**
     * A text output informing about the graph statucture of the cache nodes
     * @return a string
     */
    public String graphToString() {
        return cgs.statusToString();
    }
}
