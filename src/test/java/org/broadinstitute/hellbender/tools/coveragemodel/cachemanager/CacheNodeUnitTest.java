package org.broadinstitute.hellbender.tools.coveragemodel.cachemanager;

import org.broadinstitute.hellbender.utils.MathObjectAsserts;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.*;

/**
 * Unit tests for {@link CacheNode}.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class CacheNodeUnitTest extends BaseTest {
    private static final List<CacheNode.NodeTag> EMPTY_NODE_TAG_LIST = Collections.emptyList();
    private static final List<CacheNode.NodeKey> EMPTY_NODE_KEY_LIST = Collections.emptyList();
    private static final Map<CacheNode.NodeKey, Duplicable> EMPTY_PARENTS = Collections.emptyMap();

    /**
     * Asserts that the equality comparison of two {@link CacheNode}s is done just based on their key
     */
    @Test
    public void testEquality() {
        final List<CacheNode> nodesWithOneKey = getRandomCollectionOfNodesWithTheSameKey(new CacheNode.NodeKey("ONE_KEY"));
        final List<CacheNode> nodesWithAnotherKey = getRandomCollectionOfNodesWithTheSameKey(new CacheNode.NodeKey("ANOTHER_KEY"));

        for (final CacheNode node_0 : nodesWithOneKey) {
            for (final CacheNode node_1 : nodesWithOneKey) {
                Assert.assertTrue(node_0.equals(node_1) || (node_0.getClass() != node_1.getClass()));
            }
        }

        for (final CacheNode node_0 : nodesWithAnotherKey) {
            for (final CacheNode node_1 : nodesWithAnotherKey) {
                Assert.assertTrue(node_0.equals(node_1) || (node_0.getClass() != node_1.getClass()));
            }
        }

        for (final CacheNode node_0 : nodesWithOneKey) {
            for (final CacheNode node_1 : nodesWithAnotherKey) {
                Assert.assertTrue(!node_0.equals(node_1));
            }
        }
    }

    @Test
    public void testToString() {
        final List<CacheNode> nodesWithOneKey = getRandomCollectionOfNodesWithTheSameKey(new CacheNode.NodeKey("ONE_KEY"));
        for (final CacheNode node : nodesWithOneKey) {
            Assert.assertTrue(node.toString().equals("ONE_KEY"));
        }
    }

    private List<CacheNode> getRandomCollectionOfNodesWithTheSameKey(final CacheNode.NodeKey key) {
        final List<CacheNode> collection = new ArrayList<>();
        collection.add(new PrimitiveCacheNode(key, EMPTY_NODE_TAG_LIST, null));
        collection.add(new PrimitiveCacheNode(key,
                Arrays.asList(new CacheNode.NodeTag("a"), new CacheNode.NodeTag("b"), new CacheNode.NodeTag("c")),
                new DuplicableNumber<>(1.0)));
        collection.add(new ComputableCacheNode(key, EMPTY_NODE_TAG_LIST,
                Arrays.asList(new CacheNode.NodeKey("d"), new CacheNode.NodeKey("e")),
                ImmutableComputableGraphUnitTest.f_computation_function, false));
        collection.add(new ComputableCacheNode(key,
                EMPTY_NODE_TAG_LIST,
                EMPTY_NODE_KEY_LIST,
                ImmutableComputableGraphUnitTest.f_computation_function, false));
        collection.add(new ComputableCacheNode(key,
                Arrays.asList(new CacheNode.NodeTag("f")),
                Arrays.asList(new CacheNode.NodeKey("g")), null, true));
        return collection;
    }

    @Test(expectedExceptions = UnsupportedOperationException.class)
    public void testSetValueOfAutomaticallyComputableNode() {
        new ComputableCacheNode(new CacheNode.NodeKey("TEST"),
                EMPTY_NODE_TAG_LIST, EMPTY_NODE_KEY_LIST, ImmutableComputableGraphUnitTest.h_computation_function, false)
                .set(new DuplicableNDArray(ImmutableComputableGraphUnitTest.getRandomINDArray()));
    }

    @Test
    public void testSetValueOfExternallyComputableNode() {
        final ComputableCacheNode node = new ComputableCacheNode(new CacheNode.NodeKey("TEST"),
                EMPTY_NODE_TAG_LIST, EMPTY_NODE_KEY_LIST, null, true);
        final INDArray arr = ImmutableComputableGraphUnitTest.getRandomINDArray();
        node.set(new DuplicableNDArray(arr));
        MathObjectAsserts.assertNDArrayEquals((INDArray)node.get(EMPTY_PARENTS).value(), arr);
    }

    @Test
    public void testSetValueOfPrimitiveNode() {
        final PrimitiveCacheNode node = new PrimitiveCacheNode(new CacheNode.NodeKey("TEST"), EMPTY_NODE_TAG_LIST, null);
        final INDArray arr = ImmutableComputableGraphUnitTest.getRandomINDArray();
        node.set(new DuplicableNDArray(arr));
        MathObjectAsserts.assertNDArrayEquals((INDArray)node.get(EMPTY_PARENTS).value(), arr);
    }

    @Test
    public void testPrimitiveNodeDuplication() {
        final PrimitiveCacheNode node = new PrimitiveCacheNode(new CacheNode.NodeKey("TEST"), EMPTY_NODE_TAG_LIST,
                new DuplicableNDArray(ImmutableComputableGraphUnitTest.getRandomINDArray()));
        final PrimitiveCacheNode dupNode = node.duplicate();
        MathObjectAsserts.assertNDArrayEquals((INDArray)node.get(EMPTY_PARENTS).value(),
                (INDArray)dupNode.get(EMPTY_PARENTS).value());
        Assert.assertTrue(dupNode.hasValue());
        Assert.assertTrue(dupNode.getKey().equals(new CacheNode.NodeKey("TEST")));
    }

    @Test
    public void testCachingComputableNodeDuplication() {
        final INDArray testArray = ImmutableComputableGraphUnitTest.getRandomINDArray();
        final Duplicable testDuplicable = new DuplicableNDArray(testArray);
        final ComputableNodeFunction trivialFunction = parents -> testDuplicable;

        final ComputableCacheNode cachingAutoNodeUncached = new ComputableCacheNode(new CacheNode.NodeKey("TEST"),
                EMPTY_NODE_TAG_LIST, EMPTY_NODE_KEY_LIST, trivialFunction, true);
        final ComputableCacheNode cachingAutoNodeUncachedDup = cachingAutoNodeUncached.duplicate();

        final ComputableCacheNode cachingAutoNodeCached = cachingAutoNodeUncached.duplicateWithUpdatedValue(testDuplicable);
        final ComputableCacheNode cachingAutoNodeCachedDup = cachingAutoNodeCached.duplicate();

        final ComputableCacheNode cachingAutoNodeCachedOutdated = cachingAutoNodeCached.duplicateWithOutdatedCacheStatus();
        final ComputableCacheNode cachingAutoNodeCachedOutdatedDup = cachingAutoNodeCachedOutdated.duplicate();

        Assert.assertTrue(cachingAutoNodeUncached.isCaching());
        Assert.assertTrue(cachingAutoNodeUncachedDup.isCaching());
        Assert.assertTrue(!cachingAutoNodeUncached.isExternallyComputed());
        Assert.assertTrue(!cachingAutoNodeUncachedDup.isExternallyComputed());
        Assert.assertTrue(!cachingAutoNodeUncached.hasValue());
        Assert.assertTrue(!cachingAutoNodeUncachedDup.hasValue());

        Assert.assertTrue(cachingAutoNodeCached.isCaching());
        Assert.assertTrue(cachingAutoNodeCachedDup.isCaching());
        Assert.assertTrue(!cachingAutoNodeCached.isExternallyComputed());
        Assert.assertTrue(!cachingAutoNodeCachedDup.isExternallyComputed());
        Assert.assertTrue(cachingAutoNodeCached.hasValue());
        Assert.assertTrue(cachingAutoNodeCachedDup.hasValue());
        MathObjectAsserts.assertNDArrayEquals((INDArray)cachingAutoNodeCached.get(EMPTY_PARENTS).value(),
                (INDArray)cachingAutoNodeCachedDup.get(EMPTY_PARENTS).value());

        Assert.assertTrue(cachingAutoNodeCachedOutdated.isCaching());
        Assert.assertTrue(cachingAutoNodeCachedOutdatedDup.isCaching());
        Assert.assertTrue(!cachingAutoNodeCachedOutdated.isExternallyComputed());
        Assert.assertTrue(!cachingAutoNodeCachedOutdatedDup.isExternallyComputed());
        Assert.assertTrue(!cachingAutoNodeCachedOutdated.hasValue()); /* outdated caches must drop out */
        Assert.assertTrue(!cachingAutoNodeCachedOutdatedDup.hasValue()); /* outdated caches must drop out */
    }
}
