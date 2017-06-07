package org.broadinstitute.hellbender.tools.coveragemodel.cachemanager;

import avro.shaded.com.google.common.collect.ImmutableMap;
import avro.shaded.com.google.common.collect.Sets;
import org.apache.commons.lang.RandomStringUtils;
import org.broadinstitute.hellbender.utils.MathObjectAsserts;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;
import org.nd4j.linalg.ops.transforms.Transforms;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import javax.annotation.Nonnull;
import java.util.*;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Unit tests for {@link ImmutableComputableGraph}
 *
 * Most of the tests are done on the following, fairly generic, graph:
 *
 *     x    y    z
 *     /\  / \  /
 *     \ \/   \/
 *      \ f   g
 *       \ \  /
 *        \ \/
 *          h
 *
 *  x stores {@code DuplicableNDArray}
 *  y stores {@code DuplicableNumber<Double>}
 *  z stores {@code DuplicableNDArray}
 *
 *  f = f(x, y)
 *  g = g(y, z)
 *  h = h(f, g, x)
 *
 *  f, g, and h will be randomly composed from sine, cosine, identity, addition, subtraction, and
 *  multiplication every time by calling {@link #generateNewRandomFunctionalComposition()}.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class ImmutableComputableGraphUnitTest extends BaseTest {

    private static final Random rng = new Random(1984);

    /**
     * Number of trials for tests that involve random numbers
     */
    private static final int NUM_TRIALS = 5;

    private static final CacheNode.NodeKey X_KEY = new CacheNode.NodeKey("x");
    private static final CacheNode.NodeKey Y_KEY = new CacheNode.NodeKey("y");
    private static final CacheNode.NodeKey Z_KEY = new CacheNode.NodeKey("z");
    private static final CacheNode.NodeKey F_KEY = new CacheNode.NodeKey("f");
    private static final CacheNode.NodeKey G_KEY = new CacheNode.NodeKey("g");
    private static final CacheNode.NodeKey H_KEY = new CacheNode.NodeKey("h");

    private static final Set<CacheNode.NodeKey> ALL_NODES = new HashSet<>(Arrays.asList(X_KEY, Y_KEY, Z_KEY, F_KEY, G_KEY, H_KEY));
    private static final Set<CacheNode.NodeKey> ALL_PRIMITIVE_NODES = new HashSet<>(Arrays.asList(X_KEY, Y_KEY, Z_KEY));
    private static final Set<CacheNode.NodeKey> ALL_COMPUTABLE_NODES = new HashSet<>(Arrays.asList(F_KEY, G_KEY, H_KEY));

    /**
     * A static counter to keep track of the number of function evaluations
     */
    private static final Counter counter = new Counter(F_KEY, G_KEY, H_KEY);

    /**
     * Shape of NDArrays in the ICG
     */
    private static final int[] TEST_NDARRAY_SHAPE = new int[] {2, 3};

    private static final Function<INDArray, INDArray> UNARY_FUNCTION_COSINE = x -> Transforms.cos(x, true);
    private static final Function<INDArray, INDArray> UNARY_FUNCTION_SINE = x -> Transforms.sin(x, true);
    private static final Function<INDArray, INDArray> UNARY_FUNCTION_IDENTITY = x -> x;
    private static final BiFunction<INDArray, INDArray, INDArray> BINARY_FUNCTION_ADD = INDArray::add;
    private static final BiFunction<INDArray, INDArray, INDArray> BINARY_FUNCTION_SUB = INDArray::sub;
    private static final BiFunction<INDArray, INDArray, INDArray> BINARY_FUNCTION_MUL = INDArray::mul;

    private static final Map<String, BiFunction<INDArray, INDArray, INDArray>> TEST_BINARY_FUNCTIONS = ImmutableMap.of(
            "add", BINARY_FUNCTION_ADD,
            "sub", BINARY_FUNCTION_SUB,
            "mul", BINARY_FUNCTION_MUL);

    private static final Map<String, Function<INDArray, INDArray>> TEST_UNARY_FUNCTIONS = ImmutableMap.of(
            "sin", UNARY_FUNCTION_SINE,
            "cos", UNARY_FUNCTION_COSINE,
            "id", UNARY_FUNCTION_IDENTITY);

    /**
     * Instructions for computing f(x, y) = F[2] ( F[0](x), F[1](y) )
     *
     * The first two strings describe one of the unary functions in {@link #TEST_UNARY_FUNCTIONS}
     * The last string describes a binary function in {@link #TEST_BINARY_FUNCTIONS}
     */
    private static final List<String> F_COMPUTATION_INSTRUCTIONS = new ArrayList<>(Arrays.asList("id", "id", "add"));

    /**
     * Instructions for computing g(y, z) = F[2] ( F[0](y), F[1](z) )
     *
     * The first two strings describe one of the unary functions in {@link #TEST_UNARY_FUNCTIONS}
     * The last string describes a binary function in {@link #TEST_BINARY_FUNCTIONS}
     */
    private static final List<String> G_COMPUTATION_INSTRUCTIONS = new ArrayList<>(Arrays.asList("id", "id", "mul"));

    /**
     * Instructions for computing h(f, g, x) = F[4]( F[3]( F[0](f), F[1](g) ), F[2](x) )
     *
     * The first three strings describe one of the unary functions in {@link #TEST_UNARY_FUNCTIONS}
     * The last two strings describe one of the binary functions in {@link #TEST_BINARY_FUNCTIONS}
     */
    private static final List<String> H_COMPUTATION_INSTRUCTIONS = new ArrayList<>(Arrays.asList("id", "id", "id", "mul", "sub"));

    /**
     * Generates new functions for f, g, and h by updating {@link #F_COMPUTATION_INSTRUCTIONS},
     * {@link #G_COMPUTATION_INSTRUCTIONS}, and {@link #H_COMPUTATION_INSTRUCTIONS}
     */
    private static void generateNewRandomFunctionalComposition() {
        F_COMPUTATION_INSTRUCTIONS.clear();
        F_COMPUTATION_INSTRUCTIONS.add(getRandomChoice(TEST_UNARY_FUNCTIONS.keySet()));
        F_COMPUTATION_INSTRUCTIONS.add(getRandomChoice(TEST_UNARY_FUNCTIONS.keySet()));
        F_COMPUTATION_INSTRUCTIONS.add(getRandomChoice(TEST_BINARY_FUNCTIONS.keySet()));

        G_COMPUTATION_INSTRUCTIONS.clear();
        G_COMPUTATION_INSTRUCTIONS.add(getRandomChoice(TEST_UNARY_FUNCTIONS.keySet()));
        G_COMPUTATION_INSTRUCTIONS.add(getRandomChoice(TEST_UNARY_FUNCTIONS.keySet()));
        G_COMPUTATION_INSTRUCTIONS.add(getRandomChoice(TEST_BINARY_FUNCTIONS.keySet()));

        H_COMPUTATION_INSTRUCTIONS.clear();
        H_COMPUTATION_INSTRUCTIONS.add(getRandomChoice(TEST_UNARY_FUNCTIONS.keySet()));
        H_COMPUTATION_INSTRUCTIONS.add(getRandomChoice(TEST_UNARY_FUNCTIONS.keySet()));
        H_COMPUTATION_INSTRUCTIONS.add(getRandomChoice(TEST_UNARY_FUNCTIONS.keySet()));
        H_COMPUTATION_INSTRUCTIONS.add(getRandomChoice(TEST_BINARY_FUNCTIONS.keySet()));
        H_COMPUTATION_INSTRUCTIONS.add(getRandomChoice(TEST_BINARY_FUNCTIONS.keySet()));
    }

    static INDArray getRandomINDArray() {
        return Nd4j.rand(TEST_NDARRAY_SHAPE);
    }

    private static double getRandomDouble() {
        return rng.nextDouble();
    }

    private static <T> T getRandomChoice(final Set<T> collection) {
        return getRandomChoice(new ArrayList<>(collection));
    }

    private static <T> T getRandomChoice(final List<T> collection) {
        return collection.get(rng.nextInt(collection.size()));
    }

    private Set<CacheNode.NodeTag> getRandomSetOfTags() {
        final int MAX_NUM_TAGS = 5;
        final int TAG_LENGTH = 12;
        return IntStream.range(0, rng.nextInt(MAX_NUM_TAGS))
                .mapToObj(n -> new CacheNode.NodeTag(RandomStringUtils.randomAlphanumeric(TAG_LENGTH)))
                .collect(Collectors.toSet());
    }

    private static Counter getCounterInstance() {
        return counter.copy();
    }

    /**
     * Computes F_KEY from X_KEY and Y_KEY according to the instructions in {@link #F_COMPUTATION_INSTRUCTIONS}
     */
    private static INDArray f_computer(final INDArray x, final INDArray y) {
        final INDArray xTrans = TEST_UNARY_FUNCTIONS.get(F_COMPUTATION_INSTRUCTIONS.get(0)).apply(x);
        final INDArray yTrans = TEST_UNARY_FUNCTIONS.get(F_COMPUTATION_INSTRUCTIONS.get(1)).apply(y);
        return TEST_BINARY_FUNCTIONS.get(F_COMPUTATION_INSTRUCTIONS.get(2)).apply(xTrans, yTrans);
    }

    /**
     * Computes G_KEY from Y_KEY and Z_KEY according to the instructions in {@link #G_COMPUTATION_INSTRUCTIONS}
     */
    private static INDArray g_computer(final INDArray y, final INDArray z) {
        final INDArray yTrans = TEST_UNARY_FUNCTIONS.get(G_COMPUTATION_INSTRUCTIONS.get(0)).apply(y);
        final INDArray zTrans = TEST_UNARY_FUNCTIONS.get(G_COMPUTATION_INSTRUCTIONS.get(1)).apply(z);
        return TEST_BINARY_FUNCTIONS.get(G_COMPUTATION_INSTRUCTIONS.get(2)).apply(yTrans, zTrans);
    }

    /**
     * Computes H_KEY from F_KEY, G_KEY and X_KEY according to the instructions in {@link #H_COMPUTATION_INSTRUCTIONS}
     */
    private static INDArray h_computer(final INDArray f, final INDArray g, final INDArray x) {
        final INDArray fTrans = TEST_UNARY_FUNCTIONS.get(H_COMPUTATION_INSTRUCTIONS.get(0)).apply(f);
        final INDArray gTrans = TEST_UNARY_FUNCTIONS.get(H_COMPUTATION_INSTRUCTIONS.get(1)).apply(g);
        final INDArray xTrans = TEST_UNARY_FUNCTIONS.get(H_COMPUTATION_INSTRUCTIONS.get(2)).apply(x);
        final INDArray fgResult = TEST_BINARY_FUNCTIONS.get(H_COMPUTATION_INSTRUCTIONS.get(3)).apply(fTrans, gTrans);
        return TEST_BINARY_FUNCTIONS.get(H_COMPUTATION_INSTRUCTIONS.get(4)).apply(fgResult, xTrans);
    }

    /**
     * An instance of {@link ComputableNodeFunction} for calculating f(x, y) automatically in the
     * {@link ImmutableComputableGraph} representation of the problem. It computes F_KEY by calling
     * {@link #f_computer(INDArray, INDArray)} and increments the static F_KEY-function evaluation counter.
     */
    static ComputableNodeFunction f_computation_function = new ComputableNodeFunction() {
        @Override
        public Duplicable apply(Map<CacheNode.NodeKey, Duplicable> parents) throws ParentValueNotFoundException {
            final INDArray x = fetchINDArray(X_KEY, parents);
            final INDArray y = Nd4j.zeros(x.shape()).add(fetchDouble(Y_KEY, parents));
            final INDArray result = f_computer(x, y);
            counter.increment(F_KEY);
            return new DuplicableNDArray(result);
        }
    };

    /**
     * An instance of {@link ComputableNodeFunction} for calculating g(y, z) automatically in the
     * {@link ImmutableComputableGraph} representation of the problem. It computes G_KEY by calling
     * {@link #g_computer(INDArray, INDArray)} and increments the static G_KEY-function evaluation counter.
     *
     * Note: Y_KEY will be casted into an {@link INDArray}
     */
    static ComputableNodeFunction g_computation_function = new ComputableNodeFunction() {
        @Override
        public Duplicable apply(Map<CacheNode.NodeKey, Duplicable> parents) throws ParentValueNotFoundException {
            final INDArray z = fetchINDArray(Z_KEY, parents);
            final INDArray y = Nd4j.zeros(z.shape()).add(fetchDouble(Y_KEY, parents));
            final INDArray result = g_computer(y, z);
            counter.increment(G_KEY);
            return new DuplicableNDArray(result);
        }
    };

    /**
     * An instance of {@link ComputableNodeFunction} for calculating h(f, g, x) automatically in the
     * {@link ImmutableComputableGraph} representation of the problem. It computes H_KEY by calling
     * {@link #h_computer(INDArray, INDArray, INDArray)} and increments the static H_KEY-function evaluation
     * counter.
     */
    static ComputableNodeFunction h_computation_function = new ComputableNodeFunction() {
        @Override
        public Duplicable apply(Map<CacheNode.NodeKey, Duplicable> parents) throws ParentValueNotFoundException {
            final INDArray f = fetchINDArray(F_KEY, parents);
            final INDArray g = fetchINDArray(G_KEY, parents);
            final INDArray x = fetchINDArray(X_KEY, parents);
            final INDArray result = h_computer(f, g, x);
            counter.increment(H_KEY);
            return new DuplicableNDArray(result);
        }
    };

    private static ImmutableComputableGraphUtils.ImmutableComputableGraphBuilder getTestICGBuilder(
            final boolean f_caching, final boolean f_external,
            final boolean g_caching, final boolean g_external,
            final boolean h_caching, final boolean h_external,
            final CacheNode.NodeTag[] x_tags, final CacheNode.NodeTag[] y_tags, final CacheNode.NodeTag[] z_tags,
            final CacheNode.NodeTag[] f_tags, final CacheNode.NodeTag[] g_tags, final CacheNode.NodeTag[] h_tags) {
        return ImmutableComputableGraph.builder()
                .primitiveNode(X_KEY, x_tags, new DuplicableNDArray())
                .primitiveNode(Y_KEY, y_tags, new DuplicableNumber<Double>())
                .primitiveNode(Z_KEY, z_tags, new DuplicableNDArray())
                .computableNode(F_KEY, f_tags, new CacheNode.NodeKey[] {X_KEY, Y_KEY},
                        f_external ? null : f_computation_function, f_caching)
                .computableNode(G_KEY, g_tags, new CacheNode.NodeKey[] {Y_KEY, Z_KEY},
                        g_external ? null : g_computation_function, g_caching)
                .computableNode(H_KEY, h_tags, new CacheNode.NodeKey[] {F_KEY, G_KEY, X_KEY},
                        h_external ? null : h_computation_function, h_caching);
    }

    private static ImmutableComputableGraphUtils.ImmutableComputableGraphBuilder getTestICGBuilder(
            final boolean f_caching, final boolean f_external,
            final boolean g_caching, final boolean g_external,
            final boolean h_caching, final boolean h_external) {
        return getTestICGBuilder(f_caching, f_external, g_caching, g_external, h_caching, h_external,
                new CacheNode.NodeTag[] {}, new CacheNode.NodeTag[] {}, new CacheNode.NodeTag[] {},
                new CacheNode.NodeTag[] {}, new CacheNode.NodeTag[] {}, new CacheNode.NodeTag[] {});
    }

    /**
     * Calculates f, g, and h directly and asserts correctness
     */
    private static void assertCorrectness(final INDArray xExpected, final INDArray yExpected, final INDArray zExpected,
                                          final INDArray xActual, final INDArray yActual, final INDArray zActual,
                                          final INDArray fActual, final INDArray gActual, final INDArray hActual) {
        final INDArray fExpected = f_computer(xExpected, yExpected);
        final INDArray gExpected = g_computer(yExpected, zExpected);
        final INDArray hExpected = h_computer(fExpected, gExpected, xExpected);
        MathObjectAsserts.assertNDArrayEquals(xActual, xExpected);
        MathObjectAsserts.assertNDArrayEquals(yActual, yExpected);
        MathObjectAsserts.assertNDArrayEquals(zActual, zExpected);
        MathObjectAsserts.assertNDArrayEquals(fActual, fExpected);
        MathObjectAsserts.assertNDArrayEquals(gActual, gExpected);
        MathObjectAsserts.assertNDArrayEquals(hActual, hExpected);
    }

    private boolean assertIntactReferences(@Nonnull final ImmutableComputableGraph original,
                                           @Nonnull final ImmutableComputableGraph other,
                                           @Nonnull final Set<CacheNode.NodeKey> unaffectedNodeKeys) {
        final Set<CacheNode.NodeKey> affectedNodeKeys = unaffectedNodeKeys.stream()
                .filter(nodeKey -> original.getCacheNode(nodeKey) != other.getCacheNode(nodeKey))
                .collect(Collectors.toSet());
        if (!affectedNodeKeys.isEmpty()) {
            throw new AssertionError("Some of the node references have changed but they were supposed to remain intact: "
                    + affectedNodeKeys.stream().map(CacheNode.NodeKey::toString).collect(Collectors.joining(", ")));
        }
        return true;
    }

    private boolean assertChangedReferences(@Nonnull final ImmutableComputableGraph original,
                                            @Nonnull final ImmutableComputableGraph other,
                                            @Nonnull final Set<CacheNode.NodeKey> affectedNodeKeys) {
        final Set<CacheNode.NodeKey> unaffectedNodeKeys = affectedNodeKeys.stream()
                .filter(nodeKey -> original.getCacheNode(nodeKey) == other.getCacheNode(nodeKey))
                .collect(Collectors.toSet());
        if (!unaffectedNodeKeys.isEmpty()) {
            throw new AssertionError("Some of the node references have not changed but they were supposed to change: " +
                    unaffectedNodeKeys.stream().map(CacheNode.NodeKey::toString).collect(Collectors.joining(", ")));
        }
        return true;
    }

    private boolean assertIntactReferences(@Nonnull final ImmutableComputableGraph original,
                                           @Nonnull final ImmutableComputableGraph other,
                                           @Nonnull final CacheNode.NodeKey... unaffectedNodeKeys) {
        return assertIntactReferences(original, other, Arrays.stream(unaffectedNodeKeys).collect(Collectors.toSet()));
    }

    private boolean assertChangedReferences(@Nonnull final ImmutableComputableGraph original,
                                            @Nonnull final ImmutableComputableGraph other,
                                            @Nonnull final CacheNode.NodeKey... affectedNodeKeys) {
        return assertChangedReferences(original, other, Arrays.stream(affectedNodeKeys).collect(Collectors.toSet()));
    }

    /**
     * Tests a fully automated auto-updating {@link ImmutableComputableGraph}
     */
    @Test(invocationCount = NUM_TRIALS)
    public void testAutoUpdateCache() {
        final ImmutableComputableGraph icg_0 = getTestICGBuilder(true, false, true, false, true, false)
                .withCacheAutoUpdate().build();
        generateNewRandomFunctionalComposition();
        final INDArray x = getRandomINDArray();
        final double y = getRandomDouble();
        final INDArray z = getRandomINDArray();

        Counter startCounts = getCounterInstance();
        ImmutableComputableGraph icg_1 = icg_0
                .setValue(X_KEY, new DuplicableNDArray(x))
                .setValue(Y_KEY, new DuplicableNumber<>(y))
                .setValue(Z_KEY, new DuplicableNDArray(z));
        final INDArray xICG = (INDArray)icg_1.fetchDirectly(X_KEY).value();
        final double yICG = (Double)icg_1.fetchDirectly(Y_KEY).value();
        final INDArray zICG = (INDArray)icg_1.fetchDirectly(Z_KEY).value();
        final INDArray fICG = (INDArray)icg_1.fetchDirectly(F_KEY).value();
        final INDArray gICG = (INDArray)icg_1.fetchDirectly(G_KEY).value();
        final INDArray hICG = (INDArray)icg_1.fetchDirectly(H_KEY).value();
        Counter diffCounts = getCounterInstance().diff(startCounts);

        assertCorrectness(x, Nd4j.zeros(TEST_NDARRAY_SHAPE).add(y), z,
                xICG, Nd4j.zeros(TEST_NDARRAY_SHAPE).add(yICG), zICG,
                fICG, gICG, hICG);

        /* each function must be calculated only once; otherwise, ICG is doing redundant computations */
        Assert.assertEquals(diffCounts.getCount(F_KEY), 1);
        Assert.assertEquals(diffCounts.getCount(G_KEY), 1);
        Assert.assertEquals(diffCounts.getCount(H_KEY), 1);

        /* if we update all caches again, nothing should change */
        startCounts = getCounterInstance();
        ImmutableComputableGraph icg_2 = icg_1.updateAllCaches();
        diffCounts = getCounterInstance().diff(startCounts);
        assertIntactReferences(icg_1, icg_2, ALL_NODES);
        diffCounts.assertZero();
    }

    @DataProvider(name = "allPossibleNodeFlags")
    public Object[][] getAllPossibleNodeFlags() {
        final List<Object[]> data = new ArrayList<>();
        for (final boolean f_caching : new boolean[] {true, false})
            for (final boolean f_external : f_caching ? new boolean[] {true, false} : new boolean[] {false})
                for (final boolean g_caching : new boolean[] {true, false})
                    for (final boolean g_external : g_caching ? new boolean[] {true, false} : new boolean[] {false})
                        for (final boolean h_caching : new boolean[] {true, false})
                            for (final boolean h_external : h_caching ? new boolean[] {true, false} : new boolean[] {false})
                                data.add(new Object[] {f_caching, f_external, g_caching, g_external, h_caching, h_external});
        return data.toArray(new Object[data.size()][6]);
    }

    /**
     * Tests bookkeeping of outdated nodes
     */
    @Test(dataProvider = "allPossibleNodeFlags", invocationCount = NUM_TRIALS)
    public void testBookkeeping(final boolean f_caching, final boolean f_external,
                                final boolean g_caching, final boolean g_external,
                                final boolean h_caching, final boolean h_external) {
        generateNewRandomFunctionalComposition();
        final ImmutableComputableGraph icg_0 = getTestICGBuilder(f_caching, f_external, g_caching, g_external,
                h_caching, h_external).build();

        Assert.assertTrue(!icg_0.isValueDirectlyAvailable(X_KEY));
        Assert.assertTrue(!icg_0.isValueDirectlyAvailable(Y_KEY));
        Assert.assertTrue(!icg_0.isValueDirectlyAvailable(Z_KEY));
        Assert.assertTrue(!icg_0.isValueDirectlyAvailable(F_KEY));
        Assert.assertTrue(!icg_0.isValueDirectlyAvailable(G_KEY));
        Assert.assertTrue(!icg_0.isValueDirectlyAvailable(H_KEY));

        ImmutableComputableGraph icg_tmp = icg_0;
        icg_tmp = icg_tmp.setValue(X_KEY, new DuplicableNDArray(getRandomINDArray()));
        Assert.assertTrue(icg_tmp.isValueDirectlyAvailable(X_KEY));
        Assert.assertTrue(!icg_tmp.isValueDirectlyAvailable(Y_KEY));
        Assert.assertTrue(!icg_tmp.isValueDirectlyAvailable(Z_KEY));
        Assert.assertTrue(!icg_tmp.isValueDirectlyAvailable(F_KEY));
        Assert.assertTrue(!icg_tmp.isValueDirectlyAvailable(G_KEY));
        Assert.assertTrue(!icg_tmp.isValueDirectlyAvailable(H_KEY));
        assertIntactReferences(icg_0, icg_tmp, Y_KEY, Z_KEY, G_KEY);

        ImmutableComputableGraph icg_tmp_old = icg_tmp;
        icg_tmp = icg_tmp.setValue(Y_KEY, new DuplicableNumber<>(getRandomDouble()));
        Assert.assertTrue(icg_tmp.isValueDirectlyAvailable(X_KEY));
        Assert.assertTrue(icg_tmp.isValueDirectlyAvailable(Y_KEY));
        Assert.assertTrue(!icg_tmp.isValueDirectlyAvailable(Z_KEY));
        Assert.assertTrue(!icg_tmp.isValueDirectlyAvailable(F_KEY));
        Assert.assertTrue(!icg_tmp.isValueDirectlyAvailable(G_KEY));
        Assert.assertTrue(!icg_tmp.isValueDirectlyAvailable(H_KEY));
        assertIntactReferences(icg_tmp_old, icg_tmp, X_KEY, Z_KEY);

        icg_tmp_old = icg_tmp;
        icg_tmp = icg_tmp.setValue(Z_KEY, new DuplicableNDArray(getRandomINDArray()));
        Assert.assertTrue(icg_tmp.isValueDirectlyAvailable(X_KEY));
        Assert.assertTrue(icg_tmp.isValueDirectlyAvailable(Y_KEY));
        Assert.assertTrue(icg_tmp.isValueDirectlyAvailable(Z_KEY));
        Assert.assertTrue(!icg_tmp.isValueDirectlyAvailable(F_KEY));
        Assert.assertTrue(!icg_tmp.isValueDirectlyAvailable(G_KEY));
        Assert.assertTrue(!icg_tmp.isValueDirectlyAvailable(H_KEY));
        assertIntactReferences(icg_tmp_old, icg_tmp, X_KEY, Y_KEY, F_KEY);

        icg_tmp_old = icg_tmp;
        try {
            icg_tmp = icg_tmp.updateAllCaches();
        } catch (final Exception ex) {
            if (!f_external && !g_external && !h_external) {
                throw new AssertionError("Could not update all caches but it should have been possible");
            } else {
                icg_tmp = icg_tmp.updateAllCachesIfPossible(); /* this will not throw exception by design */
            }
        }
        assertIntactReferences(icg_tmp_old, icg_tmp, X_KEY, Y_KEY, Z_KEY);

        Assert.assertTrue((!f_caching && assertIntactReferences(icg_tmp_old, icg_tmp, F_KEY)) ||
                (f_external && !icg_tmp.isValueDirectlyAvailable(F_KEY) && assertIntactReferences(icg_tmp_old, icg_tmp, F_KEY)) ||
                (!f_external && icg_tmp.isValueDirectlyAvailable(F_KEY) && assertChangedReferences(icg_tmp_old, icg_tmp, F_KEY)));

        Assert.assertTrue((!g_caching && assertIntactReferences(icg_tmp_old, icg_tmp, G_KEY)) ||
                (g_external && !icg_tmp.isValueDirectlyAvailable(G_KEY) && assertIntactReferences(icg_tmp_old, icg_tmp, G_KEY)) ||
                (!g_external && icg_tmp.isValueDirectlyAvailable(G_KEY) && assertChangedReferences(icg_tmp_old, icg_tmp, G_KEY)));

        if (!f_external && !g_external) {
            Assert.assertTrue((!h_caching && assertIntactReferences(icg_tmp_old, icg_tmp, H_KEY)) ||
                    (h_external && !icg_tmp.isValueDirectlyAvailable(H_KEY) && assertIntactReferences(icg_tmp_old, icg_tmp, H_KEY)) ||
                    (!h_external && icg_tmp.isValueDirectlyAvailable(H_KEY) && assertChangedReferences(icg_tmp_old, icg_tmp, H_KEY)));
        } else {
            Assert.assertTrue(!icg_tmp.isValueDirectlyAvailable(H_KEY) && assertIntactReferences(icg_tmp_old, icg_tmp, H_KEY));
        }

        /* fill in the external values */
        if (f_external) {
            icg_tmp_old = icg_tmp;
            icg_tmp = icg_tmp.setValue(F_KEY, f_computation_function.apply(
                    ImmutableMap.of(X_KEY, icg_tmp.fetchDirectly(X_KEY), Y_KEY, icg_tmp.fetchDirectly(Y_KEY))));
            Assert.assertTrue(icg_tmp.isValueDirectlyAvailable(F_KEY));
            assertIntactReferences(icg_tmp_old, icg_tmp, X_KEY, Y_KEY, Z_KEY, G_KEY);
            assertChangedReferences(icg_tmp_old, icg_tmp, F_KEY, H_KEY);
        }

        if (g_external) {
            icg_tmp_old = icg_tmp;
            icg_tmp = icg_tmp.setValue(G_KEY, g_computation_function.apply(
                    ImmutableMap.of(Y_KEY, icg_tmp.fetchDirectly(Y_KEY), Z_KEY, icg_tmp.fetchDirectly(Z_KEY))));
            Assert.assertTrue(icg_tmp.isValueDirectlyAvailable(G_KEY));
            assertIntactReferences(icg_tmp_old, icg_tmp, X_KEY, Y_KEY, Z_KEY, F_KEY);
            assertChangedReferences(icg_tmp_old, icg_tmp, G_KEY, H_KEY);
        }

        if (h_external) {
            icg_tmp_old = icg_tmp;
            icg_tmp = icg_tmp.setValue(H_KEY, h_computation_function.apply(ImmutableMap.of(
                    F_KEY, icg_tmp.fetchWithRequiredEvaluations(F_KEY),
                    G_KEY, icg_tmp.fetchWithRequiredEvaluations(G_KEY),
                    X_KEY, icg_tmp.fetchDirectly(X_KEY))));
            Assert.assertTrue(icg_tmp.isValueDirectlyAvailable(H_KEY));
            assertIntactReferences(icg_tmp_old, icg_tmp, X_KEY, Y_KEY, Z_KEY, F_KEY, G_KEY);
            assertChangedReferences(icg_tmp_old, icg_tmp, H_KEY);
        }

        /* since all externally computed nodes are initialized, a call to updateAllCaches() must succeed */
        icg_tmp = icg_tmp.updateAllCaches();

        /* at this point, every caching node must be up-to-date */
        Assert.assertTrue(icg_tmp.isValueDirectlyAvailable(X_KEY));
        Assert.assertTrue(icg_tmp.isValueDirectlyAvailable(Y_KEY));
        Assert.assertTrue(icg_tmp.isValueDirectlyAvailable(Z_KEY));
        Assert.assertTrue(!f_caching || icg_tmp.isValueDirectlyAvailable(F_KEY));
        Assert.assertTrue(!g_caching || icg_tmp.isValueDirectlyAvailable(G_KEY));
        Assert.assertTrue(!h_caching || icg_tmp.isValueDirectlyAvailable(H_KEY));

        /* update x -- f and h must go out of date */
        ImmutableComputableGraph icg_tmp_x = icg_tmp.setValue(X_KEY, new DuplicableNDArray(getRandomINDArray()));
        Assert.assertTrue(!icg_tmp_x.isValueDirectlyAvailable(F_KEY));
        Assert.assertTrue(!g_caching || icg_tmp_x.isValueDirectlyAvailable(G_KEY));
        Assert.assertTrue(!icg_tmp_x.isValueDirectlyAvailable(H_KEY));

        /* update y -- f, g and h must go out of date */
        ImmutableComputableGraph icg_tmp_y = icg_tmp.setValue(Y_KEY, new DuplicableNumber<>(getRandomDouble()));
        Assert.assertTrue(!icg_tmp_y.isValueDirectlyAvailable(F_KEY));
        Assert.assertTrue(!icg_tmp_y.isValueDirectlyAvailable(G_KEY));
        Assert.assertTrue(!icg_tmp_y.isValueDirectlyAvailable(H_KEY));

        /* update z -- g and h must go out of date */
        ImmutableComputableGraph icg_tmp_z = icg_tmp.setValue(Z_KEY, new DuplicableNDArray(getRandomINDArray()));
        Assert.assertTrue(!f_caching || icg_tmp_z.isValueDirectlyAvailable(F_KEY));
        Assert.assertTrue(!icg_tmp_z.isValueDirectlyAvailable(G_KEY));
        Assert.assertTrue(!icg_tmp_z.isValueDirectlyAvailable(H_KEY));

        /* update x and y -- f, g and h must go out of date */
        ImmutableComputableGraph icg_tmp_xy = icg_tmp
                .setValue(X_KEY, new DuplicableNDArray(getRandomINDArray()))
                .setValue(Y_KEY, new DuplicableNumber<>(getRandomDouble()));
        Assert.assertTrue(!icg_tmp_xy.isValueDirectlyAvailable(F_KEY));
        Assert.assertTrue(!icg_tmp_xy.isValueDirectlyAvailable(G_KEY));
        Assert.assertTrue(!icg_tmp_xy.isValueDirectlyAvailable(H_KEY));

        /* update x and z -- f, g and h must go out of date */
        ImmutableComputableGraph icg_tmp_xz = icg_tmp
                .setValue(X_KEY, new DuplicableNDArray(getRandomINDArray()))
                .setValue(Z_KEY, new DuplicableNDArray(getRandomINDArray()));
        Assert.assertTrue(!icg_tmp_xz.isValueDirectlyAvailable(F_KEY));
        Assert.assertTrue(!icg_tmp_xz.isValueDirectlyAvailable(G_KEY));
        Assert.assertTrue(!icg_tmp_xz.isValueDirectlyAvailable(H_KEY));

        /* update x and z -- f, g and h must go out of date */
        ImmutableComputableGraph icg_tmp_xyz = icg_tmp
                .setValue(X_KEY, new DuplicableNDArray(getRandomINDArray()))
                .setValue(Y_KEY, new DuplicableNumber<>(getRandomDouble()))
                .setValue(Z_KEY, new DuplicableNDArray(getRandomINDArray()));
        Assert.assertTrue(!icg_tmp_xyz.isValueDirectlyAvailable(F_KEY));
        Assert.assertTrue(!icg_tmp_xyz.isValueDirectlyAvailable(G_KEY));
        Assert.assertTrue(!icg_tmp_xyz.isValueDirectlyAvailable(H_KEY));

        if (f_external) {
            /* update f -- h must go out of date */
            ImmutableComputableGraph icg_tmp_f = icg_tmp
                    .setValue(F_KEY, new DuplicableNDArray(getRandomINDArray()));
            Assert.assertTrue(!g_caching || icg_tmp_f.isValueDirectlyAvailable(G_KEY));
            Assert.assertTrue(!icg_tmp_f.isValueDirectlyAvailable(H_KEY));
        }

        if (g_external) {
            /* update g -- h must go out of date */
            ImmutableComputableGraph icg_tmp_g = icg_tmp
                    .setValue(G_KEY, new DuplicableNDArray(getRandomINDArray()));
            Assert.assertTrue(!f_caching || icg_tmp_g.isValueDirectlyAvailable(F_KEY));
            Assert.assertTrue(!icg_tmp_g.isValueDirectlyAvailable(H_KEY));
        }

        if (f_external && g_external) {
            /* update f and g -- h must go out of date */
            ImmutableComputableGraph icg_tmp_fg = icg_tmp
                    .setValue(F_KEY, new DuplicableNDArray(getRandomINDArray()))
                    .setValue(G_KEY, new DuplicableNDArray(getRandomINDArray()));
            Assert.assertTrue(!icg_tmp_fg.isValueDirectlyAvailable(H_KEY));
        }
    }

    /**
     * Tests propagation of tags from descendents to parents
     */
    @Test(dataProvider = "allPossibleNodeFlags", invocationCount = NUM_TRIALS)
    public void testTagPropagation(final boolean f_caching, final boolean f_external,
                                   final boolean g_caching, final boolean g_external,
                                   final boolean h_caching, final boolean h_external) {
        final Set<CacheNode.NodeTag> x_tags = getRandomSetOfTags();
        final Set<CacheNode.NodeTag> y_tags = getRandomSetOfTags();
        final Set<CacheNode.NodeTag> z_tags = getRandomSetOfTags();
        final Set<CacheNode.NodeTag> f_tags = getRandomSetOfTags();
        final Set<CacheNode.NodeTag> g_tags = getRandomSetOfTags();
        final Set<CacheNode.NodeTag> h_tags = getRandomSetOfTags();
        final ImmutableComputableGraph icg = getTestICGBuilder(
                f_caching, f_external, g_caching, g_external, h_caching, h_external,
                x_tags.toArray(new CacheNode.NodeTag[0]), y_tags.toArray(new CacheNode.NodeTag[0]),
                z_tags.toArray(new CacheNode.NodeTag[0]), f_tags.toArray(new CacheNode.NodeTag[0]),
                g_tags.toArray(new CacheNode.NodeTag[0]), h_tags.toArray(new CacheNode.NodeTag[0])).build();

        final Set<CacheNode.NodeTag> all_x_tags = Sets.union(Sets.union(x_tags, f_tags), h_tags);
        final Set<CacheNode.NodeTag> all_y_tags = Sets.union(Sets.union(Sets.union(y_tags, f_tags), g_tags), h_tags);
        final Set<CacheNode.NodeTag> all_z_tags = Sets.union(Sets.union(z_tags, g_tags), h_tags);
        final Set<CacheNode.NodeTag> all_f_tags = Sets.union(f_tags, h_tags);
        final Set<CacheNode.NodeTag> all_g_tags = Sets.union(g_tags, h_tags);
        final Set<CacheNode.NodeTag> all_h_tags = h_tags;

        final Set<CacheNode.NodeTag> all_x_tags_actual = icg.getComputableGraphStructure().getInducedTagsForNode(X_KEY);
        final Set<CacheNode.NodeTag> all_y_tags_actual = icg.getComputableGraphStructure().getInducedTagsForNode(Y_KEY);
        final Set<CacheNode.NodeTag> all_z_tags_actual = icg.getComputableGraphStructure().getInducedTagsForNode(Z_KEY);
        final Set<CacheNode.NodeTag> all_f_tags_actual = icg.getComputableGraphStructure().getInducedTagsForNode(F_KEY);
        final Set<CacheNode.NodeTag> all_g_tags_actual = icg.getComputableGraphStructure().getInducedTagsForNode(G_KEY);
        final Set<CacheNode.NodeTag> all_h_tags_actual = icg.getComputableGraphStructure().getInducedTagsForNode(H_KEY);

        Assert.assertTrue(all_x_tags.equals(all_x_tags_actual));
        Assert.assertTrue(all_y_tags.equals(all_y_tags_actual));
        Assert.assertTrue(all_z_tags.equals(all_z_tags_actual));
        Assert.assertTrue(all_f_tags.equals(all_f_tags_actual));
        Assert.assertTrue(all_g_tags.equals(all_g_tags_actual));
        Assert.assertTrue(all_h_tags.equals(all_h_tags_actual));
    }

    private Map<CacheNode.NodeKey, INDArray> getExpectedComputableNodeValues(final Duplicable x, final Duplicable y, final Duplicable z) {
        final INDArray xVal = (INDArray)x.value();
        final INDArray yVal = Nd4j.zeros(TEST_NDARRAY_SHAPE).add((Double)y.value());
        final INDArray zVal = (INDArray)z.value();
        final INDArray fExpected = f_computer(xVal, yVal);
        final INDArray gExpected = g_computer(yVal, zVal);
        final INDArray hExpected = h_computer(fExpected, gExpected, xVal);
        return ImmutableMap.of(F_KEY, fExpected, G_KEY, gExpected, H_KEY, hExpected);
    }

    /**
     * Tests {@link ImmutableComputableGraph#updateCachesForTag(CacheNode.NodeTag)}}
     */
    @Test(dataProvider = "allPossibleNodeFlags", invocationCount = NUM_TRIALS)
    public void testUpdateCachesByTag(final boolean f_caching, final boolean f_external,
                                      final boolean g_caching, final boolean g_external,
                                      final boolean h_caching, final boolean h_external) {
        generateNewRandomFunctionalComposition();
        final ImmutableComputableGraph icg_empty = getTestICGBuilder(
                f_caching, f_external, g_caching, g_external, h_caching, h_external,
                getRandomSetOfTags().toArray(new CacheNode.NodeTag[0]), getRandomSetOfTags().toArray(new CacheNode.NodeTag[0]),
                getRandomSetOfTags().toArray(new CacheNode.NodeTag[0]), getRandomSetOfTags().toArray(new CacheNode.NodeTag[0]),
                getRandomSetOfTags().toArray(new CacheNode.NodeTag[0]), getRandomSetOfTags().toArray(new CacheNode.NodeTag[0])).build();

        final Set<CacheNode.NodeTag> all_x_tags = icg_empty.getComputableGraphStructure().getInducedTagsForNode(X_KEY);
        final Set<CacheNode.NodeTag> all_y_tags = icg_empty.getComputableGraphStructure().getInducedTagsForNode(Y_KEY);
        final Set<CacheNode.NodeTag> all_z_tags = icg_empty.getComputableGraphStructure().getInducedTagsForNode(Z_KEY);
        final Set<CacheNode.NodeTag> all_f_tags = icg_empty.getComputableGraphStructure().getInducedTagsForNode(F_KEY);
        final Set<CacheNode.NodeTag> all_g_tags = icg_empty.getComputableGraphStructure().getInducedTagsForNode(G_KEY);
        final Set<CacheNode.NodeTag> all_h_tags = icg_empty.getComputableGraphStructure().getInducedTagsForNode(H_KEY);
        final Set<CacheNode.NodeTag> all_tags = new HashSet<>();
        all_tags.addAll(all_x_tags); all_tags.addAll(all_y_tags); all_tags.addAll(all_z_tags);
        all_tags.addAll(all_f_tags); all_tags.addAll(all_g_tags); all_tags.addAll(all_h_tags);

        final INDArray x = getRandomINDArray();
        final double y = getRandomDouble();
        final INDArray z = getRandomINDArray();
        final ImmutableComputableGraph icg_0 = icg_empty
                .setValue(X_KEY, new DuplicableNDArray(x))
                .setValue(Y_KEY, new DuplicableNumber<>(y))
                .setValue(Z_KEY, new DuplicableNDArray(z));
        final Map<CacheNode.NodeKey, INDArray> expectedComputableNodeValues = getExpectedComputableNodeValues(
                icg_0.fetchDirectly(X_KEY), icg_0.fetchDirectly(Y_KEY), icg_0.fetchDirectly(Z_KEY));

        for (final CacheNode.NodeTag tag : all_tags) {
            ImmutableComputableGraph icg_1;
            Counter startCounter;
            try {
                startCounter = getCounterInstance();
                icg_1 = icg_0.updateCachesForTag(tag);
            } catch (final Exception ex) { /* should fail only if some of the tagged nodes are external */
                if (!f_external && !g_external && !h_external) {
                    throw new AssertionError("Could not update tagged nodes but it should have been possible");
                }
                /* perform a partial update and continue */
                startCounter = getCounterInstance();
                icg_1 = icg_0.updateCachesForTagIfPossible(tag);
            }
            final Counter evalCounts = getCounterInstance().diff(startCounter);

            /* check updated caches */
            final Set<CacheNode.NodeKey> updatedNodesExpected = new HashSet<>();
            if (!f_external && f_caching && all_f_tags.contains(tag)) {
                updatedNodesExpected.add(F_KEY);
            }
            if (!g_external && g_caching && all_g_tags.contains(tag)) {
                updatedNodesExpected.add(G_KEY);
            }
            if (!h_external && !f_external && !g_external && h_caching && all_h_tags.contains(tag)) {
                updatedNodesExpected.add(H_KEY);
            }
            assertChangedReferences(icg_0, icg_1, updatedNodesExpected);
            assertIntactReferences(icg_0, icg_1, Sets.difference(ALL_NODES, updatedNodesExpected));

            for (final CacheNode.NodeKey nodeKey : updatedNodesExpected) {
                Assert.assertTrue(icg_1.isValueDirectlyAvailable(nodeKey));
                MathObjectAsserts.assertNDArrayEquals((INDArray)icg_1.fetchDirectly(nodeKey).value(),
                        expectedComputableNodeValues.get(nodeKey));
            }
            for (final CacheNode.NodeKey nodeKey : Sets.difference(ALL_COMPUTABLE_NODES, updatedNodesExpected)) {
                Assert.assertTrue(!icg_1.isValueDirectlyAvailable(nodeKey));
            }

            /* check function evaluation counts */
            if ((!f_external && all_f_tags.contains(tag)) /* f is computable and caching */ ||
                    (all_h_tags.contains(tag) && !f_external && !g_external && !h_external) /* h, as a descendant, is computable */) {
                Assert.assertEquals(evalCounts.getCount(F_KEY), 1);
            } else {
                Assert.assertEquals(evalCounts.getCount(F_KEY), 0);
            }
            if ((!g_external && all_g_tags.contains(tag)) /* g is computable and caching */ ||
                    (all_h_tags.contains(tag) && !g_external && !f_external && !h_external) /* h, as a descendant, is computable */) {
                Assert.assertEquals(evalCounts.getCount(G_KEY), 1);
            } else {
                Assert.assertEquals(evalCounts.getCount(G_KEY), 0);
            }
            if (all_h_tags.contains(tag) && !f_external && !g_external && !h_external) {
                Assert.assertEquals(evalCounts.getCount(H_KEY), 1);
            } else {
                Assert.assertEquals(evalCounts.getCount(H_KEY), 0);
            }
        }
    }

    /**
     * Tests {@link ImmutableComputableGraph#updateCachesForNode(CacheNode.NodeKey)}}
     */
    @Test(dataProvider = "allPossibleNodeFlags", invocationCount = NUM_TRIALS)
    public void testUpdateCacheByNode(final boolean f_caching, final boolean f_external,
                                      final boolean g_caching, final boolean g_external,
                                      final boolean h_caching, final boolean h_external) {
        generateNewRandomFunctionalComposition();
        final ImmutableComputableGraph icg_empty = getTestICGBuilder(
                f_caching, f_external, g_caching, g_external, h_caching, h_external).build();

        final INDArray x = getRandomINDArray();
        final double y = getRandomDouble();
        final INDArray z = getRandomINDArray();
        final ImmutableComputableGraph icg_0 = icg_empty
                .setValue(X_KEY, new DuplicableNDArray(x))
                .setValue(Y_KEY, new DuplicableNumber<>(y))
                .setValue(Z_KEY, new DuplicableNDArray(z));
        final Map<CacheNode.NodeKey, INDArray> expectedComputableNodeValues = getExpectedComputableNodeValues(
                icg_0.fetchDirectly(X_KEY), icg_0.fetchDirectly(Y_KEY), icg_0.fetchDirectly(Z_KEY));

        for (final CacheNode.NodeKey nodeKey : ALL_PRIMITIVE_NODES) {
            Counter startCounter = getCounterInstance();
            ImmutableComputableGraph icg_1 = icg_0.updateCachesForNode(nodeKey);
            final Counter evalCounts = getCounterInstance().diff(startCounter);
            assertIntactReferences(icg_0, icg_1, ALL_NODES);
            evalCounts.assertZero();
        }

        ImmutableComputableGraph icg_1;
        Counter startCounter;
        Counter diff;

        /* tests for F_KEY */
        try {
            startCounter = getCounterInstance();
            icg_1 = icg_0.updateCachesForNode(F_KEY);
        } catch (final Exception ex) { /* should fail only if some of the tagged nodes are external */
            if (!f_external && !g_external && !h_external) {
                throw new AssertionError("Could not update tagged nodes but it should have been possible");
            }
            startCounter = getCounterInstance();
            icg_1 = icg_0.updateCachesForNodeIfPossible(F_KEY);
        }
        diff = getCounterInstance().diff(startCounter);

        assertIntactReferences(icg_0, icg_1, X_KEY, Y_KEY, Z_KEY, G_KEY);
        if (f_external) {
            assertIntactReferences(icg_0, icg_1, ALL_NODES);
            diff.assertZero();
        } else {
            Assert.assertEquals(diff.getCount(F_KEY), 1);
            Assert.assertEquals(diff.getCount(G_KEY), 0);
            Assert.assertEquals(diff.getCount(H_KEY), 0);
            if (f_caching) {
                Assert.assertTrue(icg_1.isValueDirectlyAvailable(F_KEY));
                MathObjectAsserts.assertNDArrayEquals((INDArray)icg_1.fetchDirectly(F_KEY).value(),
                        expectedComputableNodeValues.get(F_KEY));
            } else {
                final Counter before = getCounterInstance();
                Assert.assertTrue(!icg_1.isValueDirectlyAvailable(F_KEY));
                MathObjectAsserts.assertNDArrayEquals((INDArray)icg_1.fetchWithRequiredEvaluations(F_KEY).value(),
                        expectedComputableNodeValues.get(F_KEY));
                final Counter diff2 = getCounterInstance().diff(before);
                Assert.assertEquals(diff2.getCount(F_KEY), 1);
                Assert.assertEquals(diff2.getCount(G_KEY), 0);
                Assert.assertEquals(diff2.getCount(H_KEY), 0);
            }
        }

        /* tests for G_KEY */
        try {
            startCounter = getCounterInstance();
            icg_1 = icg_0.updateCachesForNode(G_KEY);
        } catch (final Exception ex) { /* should fail only if some of the tagged nodes are external */
            if (!f_external && !g_external && !h_external) {
                throw new AssertionError("Could not update tagged nodes but it should have been possible");
            }
            startCounter = getCounterInstance();
            icg_1 = icg_0.updateCachesForNodeIfPossible(G_KEY);
        }
        diff = getCounterInstance().diff(startCounter);

        assertIntactReferences(icg_0, icg_1, X_KEY, Y_KEY, Z_KEY, F_KEY);
        if (g_external) {
            assertIntactReferences(icg_0, icg_1, ALL_NODES);
            diff.assertZero();
        } else {
            Assert.assertEquals(diff.getCount(F_KEY), 0);
            Assert.assertEquals(diff.getCount(G_KEY), 1);
            Assert.assertEquals(diff.getCount(H_KEY), 0);
            if (g_caching) {
                Assert.assertTrue(icg_1.isValueDirectlyAvailable(G_KEY));
                MathObjectAsserts.assertNDArrayEquals((INDArray)icg_1.fetchDirectly(G_KEY).value(),
                        expectedComputableNodeValues.get(G_KEY));
            } else {
                Assert.assertTrue(!icg_1.isValueDirectlyAvailable(G_KEY));
                final Counter before = getCounterInstance();
                MathObjectAsserts.assertNDArrayEquals((INDArray)icg_1.fetchWithRequiredEvaluations(G_KEY).value(),
                        expectedComputableNodeValues.get(G_KEY));
                final Counter diff2 = getCounterInstance().diff(before);
                Assert.assertEquals(diff2.getCount(F_KEY), 0);
                Assert.assertEquals(diff2.getCount(G_KEY), 1);
                Assert.assertEquals(diff2.getCount(H_KEY), 0);
            }
        }

        /* tests for H_KEY */
        try {
            startCounter = getCounterInstance();
            icg_1 = icg_0.updateCachesForNode(H_KEY);
        } catch (final Exception ex) { /* should fail only if some of the tagged nodes are external */
            if (!f_external && !g_external && !h_external) {
                throw new AssertionError("Could not update tagged nodes but it should have been possible");
            }
            startCounter = getCounterInstance();
            icg_1 = icg_0.updateCachesForNodeIfPossible(H_KEY);
        }
        diff = getCounterInstance().diff(startCounter);
        assertIntactReferences(icg_0, icg_1, X_KEY, Y_KEY, Z_KEY);
        if (h_external && f_external && g_external) {
            assertIntactReferences(icg_0, icg_1, ALL_NODES);
            diff.assertZero();
        } else if (!h_external && !f_external && !g_external) {
            Assert.assertEquals(diff.getCount(F_KEY), 1);
            Assert.assertEquals(diff.getCount(G_KEY), 1);
            Assert.assertEquals(diff.getCount(H_KEY), 1);
            if (h_caching) {
                Assert.assertTrue(icg_1.isValueDirectlyAvailable(H_KEY));
                MathObjectAsserts.assertNDArrayEquals(
                        (INDArray)icg_1.fetchDirectly(H_KEY).value(),
                        expectedComputableNodeValues.get(H_KEY));
            } else {
                Assert.assertTrue(!icg_1.isValueDirectlyAvailable(H_KEY));
                final Counter before = getCounterInstance();
                MathObjectAsserts.assertNDArrayEquals(
                        (INDArray)icg_1.fetchWithRequiredEvaluations(H_KEY).value(),
                        expectedComputableNodeValues.get(H_KEY));
                final Counter diff2 = getCounterInstance().diff(before);
                Assert.assertEquals(diff2.getCount(F_KEY), f_caching ? 0 : 1);
                Assert.assertEquals(diff2.getCount(G_KEY), g_caching ? 0 : 1);
                Assert.assertEquals(diff2.getCount(H_KEY), 1);
            }
        }
    }

    @Test
    public void testUninitializedPrimitiveNode() {
        final ImmutableComputableGraph icg = getTestICGBuilder(true, false, true, false, true, false).build()
                .setValue(X_KEY, new DuplicableNDArray(getRandomINDArray()))
                .setValue(Y_KEY, new DuplicableNumber<>(getRandomDouble()));
        boolean failed = false;
        try {
            icg.updateAllCaches();
        } catch (final PrimitiveCacheNode.PrimitiveValueNotInitializedException ex) {
            failed = true;
        }
        if (!failed) {
            throw new AssertionError("Expected PrimitiveValueNotInitializedException but it was not thrown");
        }

        icg.updateCachesForNode(F_KEY); /* should not fail */

        failed = false;
        try {
            icg.updateCachesForNode(G_KEY);
        } catch (final PrimitiveCacheNode.PrimitiveValueNotInitializedException ex) {
            failed = true;
        }
        if (!failed) {
            throw new AssertionError("Expected PrimitiveValueNotInitializedException but it was not thrown");
        }

        failed = false;
        try {
            icg.updateCachesForNode(H_KEY);
        } catch (final PrimitiveCacheNode.PrimitiveValueNotInitializedException ex) {
            failed = true;
        }
        if (!failed) {
            throw new AssertionError("Expected PrimitiveValueNotInitializedException but it was not thrown");
        }
    }

    @Test
    public void testExternallyComputedNode() {
        final ImmutableComputableGraph icg = getTestICGBuilder(true, true, true, false, true, false).build()
                .setValue(X_KEY, new DuplicableNDArray(getRandomINDArray()))
                .setValue(Y_KEY, new DuplicableNumber<>(getRandomDouble()))
                .setValue(Z_KEY, new DuplicableNDArray(getRandomINDArray()));
        boolean failed = false;
        try {
            icg.updateAllCaches();
        } catch (final ComputableCacheNode.ExternallyComputableNodeValueUnavailableException ex) {
            failed = true;
        }
        if (!failed) {
            throw new AssertionError("Expected ExternallyComputableNodeValueUnavailableException but it was not thrown");
        }

        icg.updateCachesForNode(G_KEY); /* should not fail */

        failed = false;
        try {
            icg.updateCachesForNode(F_KEY);
        } catch (final ComputableCacheNode.ExternallyComputableNodeValueUnavailableException ex) {
            failed = true;
        }
        if (!failed) {
            throw new AssertionError("Expected ExternallyComputableNodeValueUnavailableException but it was not thrown");
        }

        failed = false;
        try {
            icg.updateCachesForNode(H_KEY);
        } catch (final ComputableCacheNode.ExternallyComputableNodeValueUnavailableException ex) {
            failed = true;
        }
        if (!failed) {
            throw new AssertionError("Expected ExternallyComputableNodeValueUnavailableException but it was not thrown");
        }

        /* supply f */
        ImmutableComputableGraph icg_1 = icg.setValue(F_KEY, f_computation_function.apply(
                ImmutableMap.of(X_KEY, icg.fetchDirectly(X_KEY), Y_KEY, icg.fetchDirectly(Y_KEY))));
        Assert.assertTrue(icg_1.isValueDirectlyAvailable(F_KEY));

        /* cache g */
        Assert.assertTrue(!icg_1.isValueDirectlyAvailable(G_KEY));
        Counter before = getCounterInstance();
        icg_1 = icg_1.updateCachesForNode(G_KEY);
        Assert.assertTrue(icg_1.isValueDirectlyAvailable(G_KEY));
        Counter diff = getCounterInstance().diff(before);
        Assert.assertEquals(diff.getCount(F_KEY), 0);
        Assert.assertEquals(diff.getCount(G_KEY), 1);
        Assert.assertEquals(diff.getCount(H_KEY), 0);

        /* cache h -- now, it is computable */
        Assert.assertTrue(!icg_1.isValueDirectlyAvailable(H_KEY));
        before = getCounterInstance();
        icg_1 = icg_1.updateCachesForNode(H_KEY);
        Assert.assertTrue(icg_1.isValueDirectlyAvailable(H_KEY));
        diff = getCounterInstance().diff(before);
        Assert.assertEquals(diff.getCount(F_KEY), 0);
        Assert.assertEquals(diff.getCount(G_KEY), 0);
        Assert.assertEquals(diff.getCount(H_KEY), 1);

        /* updating all caches must have no effect */
        before = getCounterInstance();
        ImmutableComputableGraph icg_2 = icg_1.updateAllCaches();
        getCounterInstance().diff(before).assertZero();
        Assert.assertTrue(icg_2.isValueDirectlyAvailable(F_KEY));
        Assert.assertTrue(icg_2.isValueDirectlyAvailable(G_KEY));
        Assert.assertTrue(icg_2.isValueDirectlyAvailable(H_KEY));
        assertIntactReferences(icg_1, icg_2, ALL_NODES);
    }

    /**
     * A simple helper class for keeping track of ICG function evaluations
     */
    private static class Counter {
        final Map<CacheNode.NodeKey, Integer> counts;

        Counter(CacheNode.NodeKey ... keys) {
            counts = new HashMap<>();
            for (final CacheNode.NodeKey key : keys) {
                counts.put(key, 0);
            }
        }

        private Counter(final Map<CacheNode.NodeKey, Integer> otherCounts) {
            counts = new HashMap<>(otherCounts.size());
            counts.putAll(otherCounts);
        }

        void increment(final CacheNode.NodeKey key) {
            counts.put(key, getCount(key) + 1);
        }

        int getCount(final CacheNode.NodeKey key) {
            return counts.get(key);
        }

        Set<CacheNode.NodeKey> getKeys() {
            return counts.keySet();
        }

        Counter copy() {
            return new Counter(counts);
        }

        Counter diff(final Counter oldCounter) {
            Utils.validateArg(Sets.symmetricDifference(oldCounter.getKeys(), getKeys()).isEmpty(),
                    "the counters must have the same keys");
            final Map<CacheNode.NodeKey, Integer> diffMap = new HashMap<>(getKeys().size());
            getKeys().forEach(key -> diffMap.put(key, getCount(key) - oldCounter.getCount(key)));
            return new Counter(diffMap);
        }

        void assertZero() {
            Assert.assertTrue(counts.values().stream().allMatch(val -> val == 0));
        }
    }
}
