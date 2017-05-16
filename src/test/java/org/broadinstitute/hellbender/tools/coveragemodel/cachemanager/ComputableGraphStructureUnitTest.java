package org.broadinstitute.hellbender.tools.coveragemodel.cachemanager;

import avro.shaded.com.google.common.collect.ImmutableMap;
import avro.shaded.com.google.common.collect.Sets;
import org.apache.commons.lang.RandomStringUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Unit tests for {@link ComputableGraphStructure}.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class ComputableGraphStructureUnitTest extends BaseTest {

    private static final Random rng = new Random(1984);
    private static final int MAX_DAG_DEPTH = 10;
    private static final int MAX_NODES_PER_LAYER = 10;
    private static final int MAX_TAGS_PER_NODE = 10;
    private static final int MAX_PARENTS_PER_NODE = 10;
    private static final int NUM_TRIALS = 5;

    private static final CacheNode.NodeKey X_KEY = new CacheNode.NodeKey("x");
    private static final CacheNode.NodeKey Y_KEY = new CacheNode.NodeKey("y");
    private static final CacheNode.NodeKey Z_KEY = new CacheNode.NodeKey("z");
    private static final CacheNode.NodeKey F_KEY = new CacheNode.NodeKey("f");
    private static final CacheNode.NodeKey G_KEY = new CacheNode.NodeKey("g");
    private static final CacheNode.NodeKey H_KEY = new CacheNode.NodeKey("h");

    @Test(expectedExceptions = ComputableGraphStructure.NonexistentParentNodeKey.class)
    public void testMissingParents() {
        final CacheNode.NodeKey Q_KEY = new CacheNode.NodeKey("q");
        ImmutableComputableGraph.builder()
                .primitiveNode(X_KEY, new CacheNode.NodeTag[] {}, new DuplicableNDArray())
                .primitiveNode(Y_KEY, new CacheNode.NodeTag[] {}, new DuplicableNumber<Double>())
                .primitiveNode(Z_KEY, new CacheNode.NodeTag[] {}, new DuplicableNDArray())
                .computableNode(F_KEY, new CacheNode.NodeTag[] {}, new CacheNode.NodeKey[] {X_KEY, Y_KEY, Z_KEY, Q_KEY}, null, true) /* q is undefined */
                .computableNode(G_KEY, new CacheNode.NodeTag[] {}, new CacheNode.NodeKey[] {Y_KEY, Z_KEY}, null, true)
                .computableNode(H_KEY, new CacheNode.NodeTag[] {}, new CacheNode.NodeKey[] {F_KEY, G_KEY}, null, true)
                .build();
    }

    @Test(expectedExceptions = ComputableGraphStructure.CyclicGraphException.class)
    public void testCyclicGraphException_1() {
        final CacheNode.NodeKey W_KEY = new CacheNode.NodeKey("w");
        ImmutableComputableGraph.builder()
                .primitiveNode(X_KEY, new CacheNode.NodeTag[] {}, new DuplicableNDArray())
                .computableNode(Y_KEY, new CacheNode.NodeTag[] {}, new CacheNode.NodeKey[] {X_KEY, W_KEY}, null, true) /* cycle */
                .computableNode(Z_KEY, new CacheNode.NodeTag[] {}, new CacheNode.NodeKey[] {Y_KEY}, null, true)
                .computableNode(W_KEY, new CacheNode.NodeTag[] {}, new CacheNode.NodeKey[] {Z_KEY}, null, true)
                .build();
    }

    @Test(expectedExceptions = ComputableGraphStructure.CyclicGraphException.class)
    public void testCyclicGraphException_2() {
        ImmutableComputableGraph.builder()
                .primitiveNode(X_KEY, new CacheNode.NodeTag[] {}, new DuplicableNDArray())
                .primitiveNode(Y_KEY, new CacheNode.NodeTag[] {}, new DuplicableNDArray())
                .primitiveNode(Z_KEY, new CacheNode.NodeTag[] {}, new DuplicableNDArray())
                .computableNode(F_KEY, new CacheNode.NodeTag[] {}, new CacheNode.NodeKey[] {X_KEY, Y_KEY, H_KEY}, null, true) /* cycle */
                .computableNode(G_KEY, new CacheNode.NodeTag[] {}, new CacheNode.NodeKey[] {Y_KEY, Z_KEY}, null, true)
                .computableNode(H_KEY, new CacheNode.NodeTag[] {}, new CacheNode.NodeKey[] {F_KEY, G_KEY}, null, true)
                .build();
    }

    @Test(invocationCount = NUM_TRIALS)
    public void testNodeTagsAndKeysInitializationAndAccessors() {
        final RandomDAG dag = RandomDAG.getRandomDAG();
        final ComputableGraphStructure cgs = new ComputableGraphStructure(dag.getEquivalentCacheNodeSet());
        final Set<CacheNode.NodeTag> cgsNodeTagsSet = cgs.getNodeTagsSet();
        final Set<CacheNode.NodeTag> dagNodeTagsSet = dag.tagsSet;
        final Set<CacheNode.NodeKey> cgsNodeKeysSet = cgs.getNodeKeysSet();
        final Set<CacheNode.NodeKey> dagNodeKeysSet = dag.nodeKeysSet;
        Assert.assertTrue(cgsNodeKeysSet.equals(dagNodeKeysSet));
        Assert.assertTrue(cgsNodeTagsSet.equals(dagNodeTagsSet));
    }

    @Test(invocationCount = NUM_TRIALS)
    public void testTopologicalOrder() {
        final RandomDAG dag = RandomDAG.getRandomDAG();
        final ComputableGraphStructure cgs = new ComputableGraphStructure(dag.getEquivalentCacheNodeSet());
        dag.nodeKeysSet.forEach(nodeKey -> Assert.assertTrue(dag.topologicalOrderMap.get(nodeKey) ==
                cgs.getTopologicalOrder(nodeKey)));
    }

    @Test(invocationCount = NUM_TRIALS)
    public void testAncestors() {
        final RandomDAG dag = RandomDAG.getRandomDAG();
        final ComputableGraphStructure cgs = new ComputableGraphStructure(dag.getEquivalentCacheNodeSet());
        dag.nodeKeysSet.forEach(nodeKey -> Assert.assertTrue(dag.getAncestors(nodeKey).equals(cgs.getAncestors(nodeKey))));
    }

    @Test(invocationCount = NUM_TRIALS)
    public void testDescendants() {
        final RandomDAG dag = RandomDAG.getRandomDAG();
        final ComputableGraphStructure cgs = new ComputableGraphStructure(dag.getEquivalentCacheNodeSet());
        dag.nodeKeysSet.forEach(nodeKey -> Assert.assertTrue(dag.getDescendents(nodeKey).equals(cgs.getDescendants(nodeKey))));
    }

    @Test(invocationCount = NUM_TRIALS)
    public void testParents() {
        final RandomDAG dag = RandomDAG.getRandomDAG();
        final ComputableGraphStructure cgs = new ComputableGraphStructure(dag.getEquivalentCacheNodeSet());
        dag.nodeKeysSet.forEach(nodeKey -> Assert.assertTrue(dag.getParents(nodeKey).equals(cgs.getParents(nodeKey))));
    }

    @Test(invocationCount = NUM_TRIALS)
    public void testChildren() {
        final RandomDAG dag = RandomDAG.getRandomDAG();
        final ComputableGraphStructure cgs = new ComputableGraphStructure(dag.getEquivalentCacheNodeSet());
        dag.nodeKeysSet.forEach(nodeKey -> Assert.assertTrue(dag.getChildren(nodeKey).equals(cgs.getChildren(nodeKey))));
    }

    @Test(invocationCount = NUM_TRIALS)
    public void testInducedTags() {
        final RandomDAG dag = RandomDAG.getRandomDAG();
        final ComputableGraphStructure cgs = new ComputableGraphStructure(dag.getEquivalentCacheNodeSet());
        dag.nodeKeysSet.forEach(nodeKey -> Assert.assertTrue(dag.getInducedTags(nodeKey).equals(cgs.getInducedTagsForNode(nodeKey))));
    }

    @Test(invocationCount = NUM_TRIALS)
    public void testTopologicalOrderForNodeEvaluation() {
        final RandomDAG dag = RandomDAG.getRandomDAG();
        final ComputableGraphStructure cgs = new ComputableGraphStructure(dag.getEquivalentCacheNodeSet());
        dag.nodeKeysSet.forEach(nodeKey -> Assert.assertTrue(isTopologicallyEquivalent(dag.getTopologicalOrderForNodeEvaluation(nodeKey),
                cgs.getTopologicalOrderForNodeEvaluation(nodeKey), dag.topologicalOrderMap)));
    }

    @Test(invocationCount = NUM_TRIALS)
    public void testTopologicalOrderForNodeMutation() {
        final RandomDAG dag = RandomDAG.getRandomDAG();
        final ComputableGraphStructure cgs = new ComputableGraphStructure(dag.getEquivalentCacheNodeSet());
        dag.nodeKeysSet.forEach(nodeKey -> Assert.assertTrue(isTopologicallyEquivalent(dag.getTopologicalOrderForNodeMutation(nodeKey),
                cgs.getTopologicalOrderForNodeMutation(nodeKey), dag.topologicalOrderMap)));
    }

    @Test(invocationCount = NUM_TRIALS)
    public void testTopologicalOrderForTagEvaluation() {
        final RandomDAG dag = RandomDAG.getRandomDAG();
        final ComputableGraphStructure cgs = new ComputableGraphStructure(dag.getEquivalentCacheNodeSet());
        dag.tagsSet.forEach(tag -> Assert.assertTrue(isTopologicallyEquivalent(dag.getTopologicalOrderForTagEvaluation(tag),
                cgs.getTopologicalOrderForTagEvaluation(tag), dag.topologicalOrderMap)));
    }

    @Test(invocationCount = NUM_TRIALS)
    public void testTopologicalOrderForCompleteEvaluation() {
        final RandomDAG dag = RandomDAG.getRandomDAG();
        final ComputableGraphStructure cgs = new ComputableGraphStructure(dag.getEquivalentCacheNodeSet());
        final List<CacheNode.NodeKey> orderedNodes = new ArrayList<>(dag.nodeKeysSet);
        orderedNodes.sort(Comparator.comparingInt(dag.topologicalOrderMap::get));
        Assert.assertTrue(isTopologicallyEquivalent(cgs.getTopologicalOrderForCompleteEvaluation(), orderedNodes,
                dag.topologicalOrderMap));
    }

    @Test
    public void testAssertTopologicallyEquivalentLists() {
        final Map<String, Integer> topologicalOrderMap = ImmutableMap.<String, Integer>builder()
                .put("a0", 0).put("b0", 0).put("c0", 0)
                .put("a1", 1).put("b1", 1).put("c1", 1)
                .put("a2", 2).put("b2", 2).put("c2", 2).put("d2", 2).build();
        Assert.assertTrue(isTopologicallyEquivalent(
                Arrays.asList("a0", "c0", "a2", "d2"),
                Arrays.asList("c0", "a0", "a2", "d2"),
                topologicalOrderMap));
        Assert.assertTrue(isTopologicallyEquivalent(
                Arrays.asList("a0", "c0", "b1", "a1", "c1", "a2", "d2"),
                Arrays.asList("c0", "a0", "b1", "c1", "a1", "a2", "d2"),
                topologicalOrderMap));
        Assert.assertTrue(!isTopologicallyEquivalent(
                Arrays.asList("c0", "a2", "d2"),
                Arrays.asList("c0", "a0", "a2", "d2"),
                topologicalOrderMap));
        Assert.assertTrue(!isTopologicallyEquivalent(
                Arrays.asList("a0", "c0", "b1", "a1", "c1", "a2", "c2"),
                Arrays.asList("c0", "a0", "b1", "c1", "a1", "a2", "d2"),
                topologicalOrderMap));
    }

    /**
     * This test helper class creates a random DAG starting from a topological order. All helper methods are
     * implemented in a brute-force manner.
     */
    private static final class RandomDAG {
        private static final int TAG_LENGTH = 32;
        private static final int NODE_KEY_LENGTH = 32;

        final Map<Integer, Set<CacheNode.NodeKey>> nodesByTopologicalOrder;
        final Map<CacheNode.NodeKey, Integer> topologicalOrderMap;
        final Map<CacheNode.NodeKey, Set<CacheNode.NodeKey>> parentsMap;
        final Map<CacheNode.NodeKey, Set<CacheNode.NodeTag>> tagsMap;
        final Set<CacheNode.NodeKey> nodeKeysSet;
        final Set<CacheNode.NodeTag> tagsSet;

        private RandomDAG(final int depth, final int maxNodesPerLayer, final int maxParentsPerNode, final int maxTagsPerNode) {
            Utils.validateArg(depth >= 0, "DAG depth must be  >= 0");
            Utils.validateArg(maxNodesPerLayer > 0, "Max nodes per layer must be positive");
            Utils.validateArg(maxParentsPerNode > 0, "Max parents per node must be positive");
            Utils.validateArg(maxTagsPerNode > 0, "Max tags per node must be positive");

            nodesByTopologicalOrder = new HashMap<>();
            parentsMap = new HashMap<>();
            tagsMap = new HashMap<>();
            nodeKeysSet = new HashSet<>();

            for (int d = 0; d <= depth; d++) {
                final Set<CacheNode.NodeKey> randomNodes = getUniqueRandomNodeKeys(maxNodesPerLayer, d);
                nodeKeysSet.addAll(randomNodes);
                nodesByTopologicalOrder.put(d, randomNodes);
                randomNodes.forEach(nodeKey -> tagsMap.put(nodeKey, getRandomTags(maxTagsPerNode)));
                randomNodes.forEach(nodeKey -> parentsMap.put(nodeKey, new HashSet<>()));
                if (d > 0) {
                    final Set<CacheNode.NodeKey> possibleAncestors = IntStream.range(0, d)
                            .mapToObj(nodesByTopologicalOrder::get)
                            .flatMap(Set::stream)
                            .collect(Collectors.toSet());
                    for (final CacheNode.NodeKey nodeKey : randomNodes) {
                        final CacheNode.NodeKey randomParent = getRandomElement(nodesByTopologicalOrder.get(d - 1));
                        final int numParents = rng.nextInt(maxParentsPerNode);
                        final Set<CacheNode.NodeKey> randomAncestors = IntStream.range(0, numParents)
                                .mapToObj(i -> getRandomElement(possibleAncestors))
                                .collect(Collectors.toSet());
                        parentsMap.put(nodeKey, new HashSet<>());
                        parentsMap.get(nodeKey).add(randomParent);
                        parentsMap.get(nodeKey).addAll(randomAncestors);
                    }
                }
            }
            topologicalOrderMap = new HashMap<>();
            nodesByTopologicalOrder.entrySet().forEach(entry -> entry.getValue()
                    .forEach(nodeKey -> topologicalOrderMap.put(nodeKey, entry.getKey())));
            tagsSet = tagsMap.values().stream().flatMap(Set::stream).collect(Collectors.toSet());
        }

        static RandomDAG getRandomDAG() {
            return new RandomDAG(rng.nextInt(MAX_DAG_DEPTH),
                    1 + rng.nextInt(MAX_NODES_PER_LAYER),
                    1 + rng.nextInt(MAX_PARENTS_PER_NODE),
                    1 + rng.nextInt(MAX_TAGS_PER_NODE));
        }

        Set<CacheNode.NodeKey> getParents(final CacheNode.NodeKey nodeKey) {
            return parentsMap.get(nodeKey);
        }

        /**
         * Brute-force method
         */
        Set<CacheNode.NodeKey> getAncestors(final CacheNode.NodeKey nodeKey) {
            final Set<CacheNode.NodeKey> parents = getParents(nodeKey);
            return Sets.union(parents, parents.stream()
                    .map(this::getAncestors)
                    .flatMap(Set::stream)
                    .collect(Collectors.toSet()));
        }

        /**
         * Brute-force method
         */
        Set<CacheNode.NodeKey> getChildren(final CacheNode.NodeKey nodeKey) {
            return nodeKeysSet.stream()
                    .filter(key -> getParents(key).contains(nodeKey))
                    .collect(Collectors.toSet());
        }

        /**
         * Brute-force method
         */
        Set<CacheNode.NodeKey> getDescendents(final CacheNode.NodeKey nodeKey) {
            final Set<CacheNode.NodeKey> children = getChildren(nodeKey);
            return Sets.union(children, children.stream()
                    .map(this::getDescendents)
                    .flatMap(Set::stream)
                    .collect(Collectors.toSet()));
        }

        Set<CacheNode.NodeTag> getInducedTags(final CacheNode.NodeKey nodeKey) {
            final Set<CacheNode.NodeKey> allNodes = Sets.union(getDescendents(nodeKey), Collections.singleton(nodeKey));
            return allNodes.stream()
                    .map(tagsMap::get)
                    .flatMap(Set::stream)
                    .collect(Collectors.toSet());
        }

        List<CacheNode.NodeKey> getTopologicalOrderForNodeEvaluation(final CacheNode.NodeKey nodeKey) {
            final List<CacheNode.NodeKey> sortedNodes = new ArrayList<>(Sets.union(getAncestors(nodeKey),
                    Collections.singleton(nodeKey)));
            sortedNodes.sort(Comparator.comparingInt(topologicalOrderMap::get));
            return sortedNodes;
        }

        /**
         * Brute-force method
         */
        List<CacheNode.NodeKey> getTopologicalOrderForNodeMutation(final CacheNode.NodeKey nodeKey) {
            final Set<CacheNode.NodeKey> mutatedNodeAndDescendants = Sets.union(getDescendents(nodeKey),
                    Collections.singleton(nodeKey));
            final Set<CacheNode.NodeKey> ancestorsOfDescendants = getDescendents(nodeKey)
                    .stream()
                    .map(this::getAncestors)
                    .flatMap(Set::stream)
                    .collect(Collectors.toSet());
            final Set<CacheNode.NodeKey> involvedNodes = Sets.union(mutatedNodeAndDescendants, ancestorsOfDescendants);
            final List<CacheNode.NodeKey> topologicallySortedInvolvedNodes = new ArrayList<>(involvedNodes);
            topologicallySortedInvolvedNodes.sort(Comparator.comparingInt(topologicalOrderMap::get));
            return topologicallySortedInvolvedNodes;
        }

        /**
         * Brute-force method
         */
        List<CacheNode.NodeKey> getTopologicalOrderForTagEvaluation(final CacheNode.NodeTag tag) {
            final Set<CacheNode.NodeKey> taggedNodes = nodeKeysSet.stream()
                    .filter(nodeKey -> getInducedTags(nodeKey).contains(tag))
                    .collect(Collectors.toSet());
            final Set<CacheNode.NodeKey> taggedNodesAndTheirAncestors = Sets.union(taggedNodes,
                    taggedNodes.stream()
                            .map(this::getAncestors)
                            .flatMap(Set::stream)
                            .collect(Collectors.toSet()));
            final List<CacheNode.NodeKey> topologicallySortedtaggedNodesAndTheirAncestors =
                    new ArrayList<>(taggedNodesAndTheirAncestors);
            topologicallySortedtaggedNodesAndTheirAncestors.sort(Comparator.comparingInt(topologicalOrderMap::get));
            return topologicallySortedtaggedNodesAndTheirAncestors;
        }

        Set<CacheNode> getEquivalentCacheNodeSet() {
            final List<CacheNode.NodeKey> shuffledNodeList = new ArrayList<>(nodeKeysSet);
            Collections.shuffle(shuffledNodeList, rng);
            return shuffledNodeList.stream()
                    .map(nodeKey -> topologicalOrderMap.get(nodeKey) == 0
                            ? new PrimitiveCacheNode(nodeKey, tagsMap.get(nodeKey), null)
                            : new ComputableCacheNode(nodeKey, tagsMap.get(nodeKey), getParents(nodeKey), null, true))
                    .collect(Collectors.toSet());
        }

        private static Set<CacheNode.NodeKey> getUniqueRandomNodeKeys(final int maxNodes, final int depth) {
            Set<CacheNode.NodeKey> randomNodeKeySet = new HashSet<>();
            final int count = rng.nextInt(maxNodes) + 1;
            while (randomNodeKeySet.stream().distinct().count() != count) {
             randomNodeKeySet = IntStream.range(0, count)
                     .mapToObj(i -> new CacheNode.NodeKey("NODE_KEY_" + depth + "_" + RandomStringUtils.randomAlphanumeric(NODE_KEY_LENGTH)))
                     .collect(Collectors.toSet());
            }
            return randomNodeKeySet;
        }

        private static Set<CacheNode.NodeTag> getRandomTags(final int maxTags) {
            return IntStream.range(0, rng.nextInt(maxTags))
                    .mapToObj(i -> new CacheNode.NodeTag("TAG_" + RandomStringUtils.randomAlphanumeric(TAG_LENGTH)))
                    .collect(Collectors.toSet());
        }

        private static <T> T getRandomElement(final Set<T> set) {
            final List<T> list = new ArrayList<>(set);
            return list.get(rng.nextInt(list.size()));
        }
    }

    private static <T> boolean isTopologicallyEquivalent(final List<T> actual, final List<T> expected,
                                              final Map<T, Integer> topologicalOrderMap) {
        if (actual == null || expected == null || !(new HashSet<>(actual).equals(new HashSet<>(expected)))) {
            return false;
        }
        Utils.validateArg(topologicalOrderMap.keySet().containsAll(actual), "Some strings have unknown topological order");
        return IntStream.range(0, expected.size()).allMatch(i ->
                topologicalOrderMap.get(actual.get(i)).equals(topologicalOrderMap.get(expected.get(i))));
    }
}
