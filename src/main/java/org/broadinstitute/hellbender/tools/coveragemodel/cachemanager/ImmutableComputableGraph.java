package org.broadinstitute.hellbender.tools.coveragemodel.cachemanager;

import com.google.common.annotations.VisibleForTesting;
import org.broadinstitute.hellbender.tools.coveragemodel.cachemanager.ImmutableComputableGraphUtils.ImmutableComputableGraphBuilder;
import org.broadinstitute.hellbender.utils.Utils;

import javax.annotation.Nonnull;
import java.io.Serializable;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * This class provides a general purpose framework for evaluating functions on directed acyclic graphs (DAG) consisting of
 * "primitive" and "computable" nodes. Primitive nodes {@link PrimitiveCacheNode} are value placeholders while computable
 * nodes {@link ComputableCacheNode} evaluate a function on the graph. Importantly, computable nodes may cache their
 * values in order to avoid redundant/expensive computations. This class implements automatic bookkeeping strategies
 * for caching the values of computable nodes and updating the status of the cached values after mutating the primitive nodes.
 *
 * The {@link CacheNode}s can store any object that implements {@link Duplicable}. These objects provide a recipe for
 * making a deep copy of the stored value(s) via {@link Duplicable#duplicate()}. If {@link Duplicable#hasValue()}
 * is true, it indicates that the {@link Duplicable} currently stores a non-null Object.
 *
 * Typical use case: evaluating computationally expensive expressions with common subexpressions.
 *
 * Example: let X, Y and Z be three immutable values (primitive nodes) and we want to calculate f(X,Y) and g(f(X,Y),Z).
 * To this end, we may first calculate Q_1 = f(X,Y) and proceed to calculate Q_2 = g(f(X,Y),Z). By caching the value of
 * Q_1, we save on recomputing the subexpression f(X,Y) every time we mutate the primitive value Z and we need g.
 * Graphically, the computation can be represented as a topologically-ordered DAG (the edges are assumed to have downward
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
 * streamlined using a CAS library. The {@link ImmutableComputableGraph} can be conveniently constructed using
 * the builder pattern provided by {@link ImmutableComputableGraphBuilder}.
 *
 * Primitive nodes live at the top of the DAG and are specified by their key, initial value, and a set of tags.
 *
 * For computable nodes, one must provide a {@link ComputableNodeFunction}, a list of parent nodes, a list of
 * tags, and whether or not the value is to be cached. Computable nodes come in three different species depending on
 * the way they are constructed:
 *
 *      (1) Caching: these nodes may store the values they evaluate for future lookup
 *      (2) Non-caching: these nodes are compute-on-demand
 *      (3) Externally-computed: these nodes are constructed using a {@code null} evaluation function. The user is
 *          responsible for evaluating these nodes when required and using {@link #setValue(CacheNode.NodeKey, Duplicable)}
 *          to update them. This class performs the bookkeeping of the status of their values (up to date or
 *          out of date) based on the provided parents list.
 *
 * One can require the caching computable nodes to be updated and cached automatically after each mutation of a
 * primitive or externally-computed node by invoking {@link ImmutableComputableGraphBuilder#withCacheAutoUpdate()}.
 * This feature, however, is <b>not</b> recommended as one may not need all the caches to be up-to-date at all times
 * (see below for updating caches selectively). This class throws exceptions if an out-of-date cached value is queried
 * in order to notify the user to update the cache manually.
 *
 * Updating caches:
 * ================
 *
 * If cache auto update is not enabled, the user is responsible for updating the caches by calling
 * either {@link #updateAllCaches()}, {@link #updateCachesForNode(CacheNode.NodeKey)}, or {@link #updateCachesForTag(CacheNode.NodeTag)}.
 * These methods also come with counterparts {@link #updateAllCachesIfPossible()},
 * {@link #updateCachesForNodeIfPossible(CacheNode.NodeKey)}, and {@link #updateCachesForTagIfPossible(CacheNode.NodeTag)}. These methods
 * are safeguarded against throwing exceptions (e.g. if a required primitive or externally-computed value is not
 * available) and try to update as many cache nodes as possible.
 *
 * Tags:
 * =====
 *
 * Tags are arbitrary string identifiers used for grouping nodes that are semantically related. The user may want to update the
 * nodes associated to the same tag simultaneously, for example, when performing an operation that requires updated
 * values for all nodes that denote subexpressions of a larger expression. In this case, the user will tag all nodes
 * that appear in the larger expression with a common name and calls {@link #updateCachesForTag(CacheNode.NodeTag)}.
 *
 * Tags are inherited from descendents to ancestors. In the above example, if Q_2 is tagged with {"FOO"} and Q_1 is tagged
 * with {"BAR"}, then X and Y both inherit {"FOO", "BAR"} tags whereas Z only inherits {"FOO"}.
 *
 * Querying:
 * =========
 *
 * The graph is queried either by calling {@link #fetchDirectly(CacheNode.NodeKey)} or
 * {@link #fetchWithRequiredEvaluations(CacheNode.NodeKey)}. The former only fetches the values and throws an exception
 * if the some of the required nodes are out of date, or a non-caching computable node is queried. The latter performs
 * the required <b>necessary</b> evaluations along the way.
 *
 * IMPORTANT NOTE: the queried values are returned by reference. It is the user's responsibility not to mutate
 * them. Otherwise, immutability will be broken.
 *
 * Mutation:
 * =========
 *
 * The primitive values can be mutated by calling {@link #setValue(CacheNode.NodeKey, Duplicable)}. Mutations do not occur in place;
 * rather, a new instance of {@link ImmutableComputableGraph} is created along with new instances for the updated nodes.
 * Immutability is desired, for instance, if this class is used as elements of a {@link org.apache.spark.api.java.JavaRDD}.
 * If cache auto update is enabled, the affected nodes will be automatically computed and cached. Otherwise, only the
 * cache status go out of date and the old stored values are {@code null}ed.
 *
 * Note: the new {@link ImmutableComputableGraph} instance returned by {@link #setValue(CacheNode.NodeKey, Duplicable)} is <b>not</b> a
 * deep copy and may hold references to {@link CacheNode}s contained the previous instance(s).
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class ImmutableComputableGraph implements Serializable {

    private static final long serialVersionUID = -1162776031416105027L;

    private static Map<CacheNode.NodeKey, Duplicable> EMPTY_NODE_KEY_VALUE_MAP = new HashMap<>();

    private final Map<CacheNode.NodeKey, CacheNode> nodesMap;
    private final boolean cacheAutoUpdate;
    private final ComputableGraphStructure cgs;

    public static ImmutableComputableGraphBuilder builder() {
        return new ImmutableComputableGraphBuilder();
    }

    /**
     * Package-private constructor from a node collection (used by the builder).
     *
     * @param nodeSet a collection of {@link CacheNode}s
     */
    ImmutableComputableGraph(@Nonnull final Set<CacheNode> nodeSet,
                             final boolean cacheAutoUpdate) {
        Utils.nonNull(nodeSet, "The nodes collection must be non-null.");
        this.cacheAutoUpdate = cacheAutoUpdate;
        nodesMap = nodeSet.stream().collect(Collectors.toMap(CacheNode::getKey, Function.identity()));
        cgs = new ComputableGraphStructure(nodeSet);
    }

    /**
     * A private constructor (used by duplicators).
     *
     * @param nodesMap a previously constructed key -> {@link CacheNode} map
     * @param cgs a previously constructed {@link ComputableGraphStructure}
     */
    private ImmutableComputableGraph(@Nonnull final Map<CacheNode.NodeKey, CacheNode> nodesMap,
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
     * @return a new instance of {@link ImmutableComputableGraph} with reference to unchanged nodes and duplicated
     *         changed nodes
     * @throws IllegalArgumentException if the node does not exist
     * @throws UnsupportedOperationException if the node is non-primitive
     */
    public ImmutableComputableGraph setValue(@Nonnull final CacheNode.NodeKey nodeKey,
                                             @Nonnull final Duplicable newValue)
            throws IllegalArgumentException, UnsupportedOperationException {
        CacheNode node = nodesMap.get(assertNodeExists(nodeKey));
        if (!node.isExternallyComputed()) {
            throw new UnsupportedOperationException("Can not explicitly set the value of a non-primitive cache node.");
        }
        final Map<CacheNode.NodeKey, CacheNode> updatedNodesMap = new HashMap<>();
        updatedNodesMap.put(nodeKey, node.duplicateWithUpdatedValue(newValue));
        final ImmutableComputableGraph out = duplicateWithUpdatedNodes(
                addDuplicateOfOutdatedDescendants(nodeKey, updatedNodesMap));
        if (cacheAutoUpdate) {
            try { /* try to update caches; it is not guaranteed if some of the nodes are not initialized */
                final Map<CacheNode.NodeKey, Duplicable> accumulatedValues = out.evaluateInTopologicalOrder(
                        cgs.getTopologicalOrderForNodeMutation(nodeKey));
                return out.updateCachesFromAccumulatedValues(accumulatedValues);
            } catch (final PrimitiveCacheNode.PrimitiveValueNotInitializedException |
                    ComputableCacheNode.ExternallyComputableNodeValueUnavailableException ex) {
                /* cache auto-update failed; will return "out" = ICG with updated node and outdated descendents */
            }
        }
        return out;
    }

    /**
     * Make a key -> node map containing new instances of the nodes that go out of date as a result of
     * updating node {@code key}
     *
     * @param key key of the updated node
     * @param updatedNodesMap a key -> node map
     */
    private Map<CacheNode.NodeKey, CacheNode> addDuplicateOfOutdatedDescendants(
            @Nonnull final CacheNode.NodeKey key,
            @Nonnull final Map<CacheNode.NodeKey, CacheNode> updatedNodesMap) {
        for (final CacheNode.NodeKey descendant : cgs.getDescendants(key)) {
            CacheNode oldDescendant = nodesMap.get(descendant);
            /* all of the descendants are computable nodes and can be safely up-casted */
            updatedNodesMap.put(descendant, ((ComputableCacheNode)oldDescendant).duplicateWithOutdatedCacheStatus());
        }
        return updatedNodesMap;
    }

    /**
     * Returns a reference to the value of a given node
     *
     * Note: this function is purposefully meant to be <i>light</i> in the following sense:
     *
     * (1) it does not update out-of-date caching computable nodes, and
     * (2) it does not evaluate non-caching computable nodes.
     *
     * @param nodeKey the key of the node
     * @return value of the node
     * @throws IllegalStateException if a cached node is out of date, or if the value of a primitive node is
     *                               not initialized
     * @throws IllegalArgumentException if the node does not exist
     */
    public Duplicable fetchDirectly(@Nonnull final CacheNode.NodeKey nodeKey) throws IllegalArgumentException {
        return nodesMap.get(assertNodeExists(nodeKey)).get(EMPTY_NODE_KEY_VALUE_MAP);
    }

    /**
     * Returns the value of a node. If the node is computable and its value is not available, all of the required
     * intermediate calculations will be done.
     *
     * Note: the result of intermediate calculations are <b>not</b> stored back into the (possibly out-of-date) ancestor
     * nodes. This may result in redundant calculations. In generic situations, the most efficient approach is to
     * update the required cache nodes first, and then fetch the up-to-date cached values using
     * {@link #fetchDirectly(CacheNode.NodeKey)} instead.
     *
     * @param nodeKey key of the node
     * @return value of the node
     * @throws IllegalArgumentException if the node does not exist
     */
    public Duplicable fetchWithRequiredEvaluations(@Nonnull final CacheNode.NodeKey nodeKey)
            throws IllegalStateException, IllegalArgumentException {
        final CacheNode node = nodesMap.get(assertNodeExists(nodeKey));
        if (node.hasValue()) {
            return node.get(EMPTY_NODE_KEY_VALUE_MAP);
        } else {
            return evaluateInTopologicalOrder(cgs.getTopologicalOrderForNodeEvaluation(nodeKey)).get(nodeKey);
        }
    }

    /**
     * This method computes and accumulates the values on a list of topologically ordered set of nodes (here meaning
     * from low depth to high depth).
     *
     * Along the way and en route to the output map,
     *
     *   - The value of a primitive node is passed by reference
     *   - A non-caching computable node is always computed
     *   - A caching computable node is only computed if it has no cached value or its stored valued is not current
     *
     * Note: this method does not check whether {@code topologicallyOrderedNodeKeys} is actually topologically ordered.
     *
     * @param topologicallyOrderedNodeKeys topologically sorted list of nodes
     * @throws ComputableNodeFunction.ParentValueNotFoundException if a parent value required for a computation function
     *         is not found; it can be thrown if {@code topologicallyOrderedNodeKeys} is not truly topologically ordered
     * @throws ComputableCacheNode.ExternallyComputableNodeValueUnavailableException if the value of an externally computable
     *         node is out of date or not initialized
     * @throws PrimitiveCacheNode.PrimitiveValueNotInitializedException if the value of a required primitive node is not
     *         initialized
     * @return a map from node keys to their values accumulated during computation
     */
    private Map<CacheNode.NodeKey, Duplicable> evaluateInTopologicalOrder(
            @Nonnull final List<CacheNode.NodeKey> topologicallyOrderedNodeKeys) {
        final Map<CacheNode.NodeKey, Duplicable> accumulatedValues = new HashMap<>();
        for (final CacheNode.NodeKey nodeKey : topologicallyOrderedNodeKeys) {
            accumulatedValues.put(nodeKey, nodesMap.get(nodeKey).get(accumulatedValues));
        }
        return accumulatedValues;
    }

    /**
     * This method is the same as {@link #evaluateInTopologicalOrder(List)} except for it catches all exceptions
     * that would otherwise be thrown by {@link #evaluateInTopologicalOrder(List)}, and performs a partial evaluation
     * to the possible extent
     *
     * @param topologicallyOrderedNodeKeys topologically sorted list of nodes
     * @return a map from node keys to their values accumulated during computation
     */
    private Map<CacheNode.NodeKey, Duplicable> evaluateInTopologicalOrderIfPossible(
            @Nonnull final List<CacheNode.NodeKey> topologicallyOrderedNodeKeys) {
        final Map<CacheNode.NodeKey, Duplicable> accumulatedValues = new HashMap<>();
        for (final CacheNode.NodeKey nodeKey : topologicallyOrderedNodeKeys) {
            Duplicable value = null;
            try {
                value = nodesMap.get(nodeKey).get(accumulatedValues);
            } catch (final Exception ex) {
                /* do nothing */
            }
            if (value != null) {
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
    private ImmutableComputableGraph updateCachesFromAccumulatedValues(
            @Nonnull final Map<CacheNode.NodeKey, Duplicable> accumulatedValues) {
        /* since accumulatedValues may contain unchanged values (by reference), we filter and only update
         * the affected nodes */
        return duplicateWithUpdatedNodes(
                accumulatedValues.keySet().stream()
                        /* filter out primitives and caching nodes that are current */
                        .filter(node -> !(nodesMap.get(node).isPrimitive() ||
                                nodesMap.get(node).hasValue() ||
                                !((ComputableCacheNode)nodesMap.get(node)).isCaching()))
                        /* collect to a map: key -> duplicated node with updated value */
                        .collect(Collectors.toMap(Function.identity(), node ->
                                ((ComputableCacheNode)nodesMap.get(node))
                                        .duplicateWithUpdatedValue(accumulatedValues.get(node)))));
    }

    /**
     * Update the cached values by node key
     *
     * @return a new instance of {@link ImmutableComputableGraph} with new instances of updated nodes
     */
    public ImmutableComputableGraph updateCachesForNode(@Nonnull final CacheNode.NodeKey nodeKey) {
        return evaluateAndUpdateCaches(this::evaluateInTopologicalOrder,
                cgs.getTopologicalOrderForNodeEvaluation(assertNodeExists(nodeKey)));
    }

    /**
     * Update the cached values by node key
     *
     * @return a new instance of {@link ImmutableComputableGraph} with new instances of updated nodes
     */
    public ImmutableComputableGraph updateCachesForNodeIfPossible(@Nonnull final CacheNode.NodeKey nodeKey) {
        return evaluateAndUpdateCaches(this::evaluateInTopologicalOrderIfPossible,
                cgs.getTopologicalOrderForNodeEvaluation(assertNodeExists(nodeKey)));
    }

    /**
     * Update the cached values by tag
     *
     * @return a new instance of {@link ImmutableComputableGraph} with new instances of updated nodes
     */
    public ImmutableComputableGraph updateCachesForTag(final CacheNode.NodeTag tagKey) {
        return evaluateAndUpdateCaches(this::evaluateInTopologicalOrder,
                cgs.getTopologicalOrderForTagEvaluation(assertTagExists(tagKey)));
    }

    /**
     * Update all possibly updatable tagged caching nodes
     *
     * @return a new instance of {@link ImmutableComputableGraph} with new instances of updated nodes
     */
    public ImmutableComputableGraph updateCachesForTagIfPossible(final CacheNode.NodeTag tagKey) {
        return evaluateAndUpdateCaches(this::evaluateInTopologicalOrderIfPossible,
                cgs.getTopologicalOrderForTagEvaluation(assertTagExists(tagKey)));
    }

    /**
     * Update all caching nodes
     *
     * @return a new instance of {@link ImmutableComputableGraph} with new instances of updated nodes
     */
    public ImmutableComputableGraph updateAllCaches() {
        return evaluateAndUpdateCaches(this::evaluateInTopologicalOrder, cgs.getTopologicalOrderForCompleteEvaluation());
    }

    /**
     * Update all possibly updatable caches
     *
     * @return a new instance of {@link ImmutableComputableGraph} with new instances of updated nodes
     */
    public ImmutableComputableGraph updateAllCachesIfPossible() {
        return evaluateAndUpdateCaches(this::evaluateInTopologicalOrderIfPossible, cgs.getTopologicalOrderForCompleteEvaluation());
    }

    private ImmutableComputableGraph evaluateAndUpdateCaches(
            @Nonnull final Function<List<CacheNode.NodeKey>, Map<CacheNode.NodeKey, Duplicable>> topologicalEvaluator,
            @Nonnull final List<CacheNode.NodeKey> topologicallyOrderedNodeKeys) {
        return updateCachesFromAccumulatedValues(topologicalEvaluator.apply(topologicallyOrderedNodeKeys));
    }

    /**
     * Make a new instance of {@link ImmutableComputableGraph} by replacing the nodes from {@code updatedNodesMap}
     *
     * @param updatedNodesMap nodes to be replaced and their new values
     * @return a new instance of {@link ImmutableComputableGraph} with new instances of updated nodes
     */
    private ImmutableComputableGraph duplicateWithUpdatedNodes(final Map<CacheNode.NodeKey, CacheNode> updatedNodesMap) {
        final Map<CacheNode.NodeKey, CacheNode> newNodesMap = new HashMap<>();
        final Set<CacheNode.NodeKey> updatedNodeKeys = updatedNodesMap.keySet();
        /* intact nodes */
        cgs.getNodeKeysSet().stream()
                .filter(node -> !updatedNodeKeys.contains(node))
                .forEach(node -> newNodesMap.put(node, nodesMap.get(node)));
        /* updated nodes */
        newNodesMap.putAll(updatedNodesMap);
        return new ImmutableComputableGraph(newNodesMap, cgs, cacheAutoUpdate);
    }

    public boolean isValueDirectlyAvailable(final CacheNode.NodeKey nodeKey) {
        return nodesMap.get(assertNodeExists(nodeKey)).hasValue();
    }

    private CacheNode.NodeKey assertNodeExists(final CacheNode.NodeKey nodeKey) {
        Utils.nonNull(nodeKey, "The node key must be non-null");
        Utils.validateArg(cgs.getNodeKeysSet().contains(nodeKey), "The node" +
                ImmutableComputableGraphUtils.quote(nodeKey.toString()) + " does not exist.");
        return nodeKey;
    }

    private CacheNode.NodeTag assertTagExists(final CacheNode.NodeTag tagKey) {
        Utils.nonNull(tagKey, "The tag key must be non-null");
        Utils.validateArg(cgs.getNodeTagsSet().contains(tagKey), "The tag " +
                ImmutableComputableGraphUtils.quote(tagKey.toString()) + " does not exist.");
        return tagKey;
    }

    @VisibleForTesting
    ComputableGraphStructure getComputableGraphStructure() {
        return cgs;
    }

    @VisibleForTesting
    CacheNode getCacheNode(@Nonnull final CacheNode.NodeKey nodeKey) {
        return nodesMap.get(assertNodeExists(nodeKey));
    }

    @Override
    public String toString() {
        String out = "";
        for (final CacheNode.NodeKey key : cgs.getNodeKeysSet()) {
            out += "node: " + key + "\n";
            out += "\ttype: " + (nodesMap.get(key).isPrimitive() ? "primitive" : "computable") + "\n";
            String value;
            try {
                value = fetchDirectly(key).toString();
            } catch (final Exception ex) {
                value = "not available";
            }
            out += "\tvalue: " + value + "\n";
            if (!nodesMap.get(key).isPrimitive()) {
                out += "\tcaches values: " + ((ComputableCacheNode)nodesMap.get(key)).isCaching() + "\n";
                out += "\thas a value: " + nodesMap.get(key).hasValue() + "\n";
            }
        }
        return out;
    }
}
