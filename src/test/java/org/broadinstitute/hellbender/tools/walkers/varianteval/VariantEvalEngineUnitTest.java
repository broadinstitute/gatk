package org.broadinstitute.hellbender.tools.walkers.varianteval;


// the imports for unit testing.

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.testutils.BaseTest;
import org.broadinstitute.hellbender.tools.walkers.varianteval.evaluators.VariantEvaluator;
import org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications.VariantStratifier;
import org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications.manager.StratificationManager;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.EvaluationContext;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.VariantEvalContext;
import org.broadinstitute.hellbender.utils.Utils;
import org.testng.Assert;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;


public class VariantEvalEngineUnitTest extends BaseTest {
    VariantEvalEngine engine;
    VariantContext eval;


    @BeforeMethod
    public void init() {
        engine = new VariantEvalEngine();
        eval = new VariantContextBuilder("x", "chr1", 1, 1, Collections.singleton(Allele.create("A", true))).make();
    }

    // --------------------------------------------------------------------------------
    //
    // Test stratifications / evaluations
    //
    // --------------------------------------------------------------------------------

    private class StratifiedEvalTestProvider extends BaseTest.TestDataProvider {
        final List<VariantStratifier> stratificationObjects = new ArrayList<VariantStratifier>();
        final Set<Class<? extends VariantEvaluator>> evaluationObjects = new HashSet<Class<? extends VariantEvaluator>>();
        final List<Integer> expectedCounts;
        final int maxI;

        /**
         *
         * @param maxI test integers from 1 ... maxI
         * @param expectedCounts the expected number of integers from 1 ... maxI divisible by each combination, in order, of allStates
         * @param allStates all stratification tests, in order
         */
        @SuppressWarnings("unchecked")
        public StratifiedEvalTestProvider(int maxI,
                                          final List<Integer> expectedCounts,
                                          final List<Integer> ... allStates) {
            super(StratifiedEvalTestProvider.class);

            this.maxI = maxI;
            this.expectedCounts = expectedCounts;
            this.evaluationObjects.add(CounterEval.class);

            String stateName = "";
            for ( List<Integer> states : allStates ) {
                stratificationObjects.add(new IntegerStratifier(states));
                stateName = stateName + Utils.join(",", states) + " ";
            }

            setName(String.format("maxI=%d expectedCounts=%s states=%s", maxI, Utils.join(",", expectedCounts), stateName));
        }
    }

    /**
     * Test stratifier -> holds a list of integers, and the states are if the integer value of evalName is divisable
     * by that number
     */
    public static class IntegerStratifier extends VariantStratifier {
        List<Integer> integers;

        public IntegerStratifier(VariantEvalEngine engine) {
            super(null);
        }

        private IntegerStratifier(final List<Integer> integers) {
            super(null);
            this.integers = integers;
            states.addAll(integers);
        }

        @Override
        public List<Object> getRelevantStates(final VariantEvalContext context, final VariantContext comp, final String compName, final VariantContext eval, final String evalName, final String sampleName, final String familyName) {
            int i = Integer.valueOf(evalName); // a terrible hack, but we can now provide accessible states
            List<Object> states = new ArrayList<>();
            for ( int state : integers )
                if ( i % state == 0 )
                    states.add(state);
            return states;
        }
    }

    /**
     * Test evaluator -> just counts the number of calls to update1
     */
    public static class CounterEval extends VariantEvaluator {
        public CounterEval(VariantEvalEngine engine) {
            super(engine);
        }

        public int count = 0;

        public CounterEval() {
            super(null);
        }

        @Override public int getComparisonOrder() { return 1; }

        @Override
        public void update1(final VariantContext eval, final VariantEvalContext context) {
            count++;
        }

        @Override
        public boolean supportsCombine() {
            return true;
        }

        @Override
        public void combine(final VariantEvaluator other) {
            this.count += ((CounterEval)other).count;
        }
    }

    private void initialize(StratifiedEvalTestProvider cfg) {
        engine.createStratificationStates(cfg.stratificationObjects, cfg.evaluationObjects);

        final VariantEvalContext vec = null;
        final VariantContext comp = null;
        final String compName = null, sampleName = null, familyName = null;

        // increment eval counts for each stratification of divisors of i from from 1...maxI
        for ( int i = 1; i <= cfg.maxI; i++ ) {
            final String evalName = String.valueOf(i); // terrible hack to stratify by divisor
            for ( EvaluationContext nec : engine.getEvaluationContexts(vec, eval, evalName, comp, compName, sampleName, familyName) ) {
                synchronized (nec) {
                    nec.apply(vec, comp, eval);
                }
            }
        }
    }

    @SuppressWarnings("unchecked")
    @DataProvider(name = "StratifiedEvalTestProvider")
    public Object[][] makeStratifiedEvalTestProvider() {

        new StratifiedEvalTestProvider(4, // test 1, 2, 3, 4
                Arrays.asList(4, 2), //  4 divisible by 1, 2 by 2
                Arrays.asList(1, 2));

        new StratifiedEvalTestProvider(6, // test 1, 2, 3, 4, 5, 6
                Arrays.asList(6, 3, 2), //  6 divisible by 1, 3 by 2, 2 by 3
                Arrays.asList(1, 2, 3));

        // test that some states can be empty -- does this work in VE?
        new StratifiedEvalTestProvider(6,
                Arrays.asList(3, 2),
                Arrays.asList(2, 3));

        // test a single stratification
        new StratifiedEvalTestProvider(6,
                Arrays.asList(3),
                Arrays.asList(2));

        // test a meaningless state
        new StratifiedEvalTestProvider(4, // test 1, 2, 3, 4
                Arrays.asList(4, 2), //  4 divisible by 1, 2 by 2
                Arrays.asList(1, 2), Arrays.asList(1));

        // test a adding a state that divides space in half
        new StratifiedEvalTestProvider(4,
                Arrays.asList(2, 2),
                Arrays.asList(1, 2), Arrays.asList(2));

        // test pairs of strats
        new StratifiedEvalTestProvider(12,
                Arrays.asList(3, 4, 3, 2),
                Arrays.asList(1, 2), Arrays.asList(3, 4));

        return StratifiedEvalTestProvider.getTests(StratifiedEvalTestProvider.class);
    }

    /**
     * Ensures that counting and stratifications all are working properly by iterating
     * over integers 1...cfg.N and stratify according to cfg, and that the counts in
     * each bin are as expected.
     *
     * @param cfg
     */
    @Test(dataProvider = "StratifiedEvalTestProvider")
    public void testBasicOperation(StratifiedEvalTestProvider cfg) {
        initialize(cfg);
        checkStratificationCountsAreExpected(engine.getStratManager(), cfg.expectedCounts);
    }

    private final void checkStratificationCountsAreExpected(final StratificationManager<VariantStratifier, EvaluationContext> manager,
                                                            final List<Integer> expectedCounts) {
        for ( int key = 0; key < manager.size(); key++ ) {
            final String stratStateString = manager.getStratsAndStatesStringForKey(key);
            final EvaluationContext nec = manager.get(key);

            for ( final VariantEvaluator ve : nec.getVariantEvaluators() ) {
                // test for count here
                final CounterEval counterEval = (CounterEval)ve;
                final int expected = expectedCounts.get(key);
                Assert.assertEquals(counterEval.count, expected, "Count seen of " + counterEval.count + " not expected " + expected + " at " + stratStateString);
            }
        }
    }

    /**
     * A derived test on testBasicOperation that checks that combining stratifications
     * works as expected by ensuring the results are the same when the remapped
     * strats are the identity map (A -> A, B -> B, etc)
     */
    @Test(dataProvider = "StratifiedEvalTestProvider", dependsOnMethods = {"testBasicOperation"})
    public void testIdentityCombine(StratifiedEvalTestProvider cfg) {
        for ( int i = 0; i < cfg.stratificationObjects.size(); i++ ) {
            initialize(cfg);
            final VariantStratifier toReplace = cfg.stratificationObjects.get(i);
            final VariantStratifier newStrat = cfg.stratificationObjects.get(i);
            final Map<Object, Object> remappedStates = makeIdentityFunctionMap(newStrat.getAllStates());
            StratificationManager<VariantStratifier, EvaluationContext> combined =
                    engine.getStratManager().combineStrats(toReplace, newStrat, EvaluationContext.COMBINER, remappedStates);
            checkStratificationCountsAreExpected(combined, cfg.expectedCounts);
        }
    }

    private <T> Map<T, T> makeIdentityFunctionMap(Collection<T> values) {
        Map<T,T> map = new HashMap<T, T>(values.size());
        for ( final T value : values )
            map.put(value, value);
        return Collections.unmodifiableMap(map);
    }

    /**
     * A derived test on testBasicOperation that checks that combining stratifications
     * works as expected. We look into cfg, and if there are multiple states we create
     * dynamically create a combinations of the stratifications, and ensure that the
     * combined results are as we expected.
     */
    @Test(dataProvider = "StratifiedEvalTestProvider", dependsOnMethods = {"testBasicOperation"})
    public void testCombinedEachStrat(StratifiedEvalTestProvider cfg) {
        for ( int i = 0; i < cfg.stratificationObjects.size(); i++ ) {
            initialize(cfg);
            final VariantStratifier toReplace = cfg.stratificationObjects.get(i);

            // TODO -- replace this code with something that combines values in strat
            final VariantStratifier newStrat = cfg.stratificationObjects.get(i);
            final Map<Object, Object> remappedStates = makeIdentityFunctionMap(newStrat.getAllStates());
            final List<Integer> expected = cfg.expectedCounts;

            StratificationManager<VariantStratifier, EvaluationContext> combined =
                    engine.getStratManager().combineStrats(toReplace, newStrat, EvaluationContext.COMBINER, remappedStates);
            checkStratificationCountsAreExpected(combined, expected);
        }
    }
}