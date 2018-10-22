package org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications.manager;


// the imports for unit testing.


import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.testutils.BaseTest;
import org.broadinstitute.hellbender.utils.Utils;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.FileNotFoundException;
import java.util.*;


public class StratificationManagerUnitTest extends BaseTest {
    @BeforeClass
    public void init() throws FileNotFoundException {
    }

    // --------------------------------------------------------------------------------
    //
    // Basic tests Provider
    //
    // --------------------------------------------------------------------------------

    private class StratificationStatesTestProvider extends BaseTest.TestDataProvider {
        final List<List<Object>> allStates = new ArrayList<>();
        final List<IntegerStratifier> asSetOfStates = new ArrayList<>();
        final int nStates;

        @SuppressWarnings("unchecked")
        public StratificationStatesTestProvider(final List<Integer> ... allStates) {
            super(StratificationStatesTestProvider.class);
            
            for ( List<Integer> states : allStates ) {
                this.allStates.add(new ArrayList<>(states));
            }

            for ( List<Object> states : this.allStates ) { 
                asSetOfStates.add(new IntegerStratifier(states));
            }
            this.nStates = nCombinations(allStates);

            setName(getName());
        }

        private String getName() {
            StringBuilder b = new StringBuilder();
            int c = 1;
            for ( List<Object> state : allStates )
                b.append(String.format("%d = [%s] ", c++, Utils.join(",", state)));
            return b.toString();
        }
        
        public List<IntegerStratifier> getStateSpaceList() {
            return asSetOfStates;
        }
        
        public ArrayList<Integer> values() {
            final ArrayList<Integer> l = new ArrayList<Integer>();
            for ( int i = 0; i < nStates; i++ )
                l.add(i);
            return l;
        }
            
        public Queue<List<Object>> getAllCombinations() {
            return getAllCombinations(new LinkedList<List<Object>>(allStates));
        }

        private Queue<List<Object>> getAllCombinations(Queue<List<Object>> states) {
            if ( states.isEmpty() ) 
                return new LinkedList<List<Object>>();
            else {
                List<Object> head = states.poll();
                Queue<List<Object>> substates = getAllCombinations(states);
                Queue<List<Object>> newStates = new LinkedList<List<Object>>();
                for ( final Object e : head) {
                    if ( substates.isEmpty() ) {
                        newStates.add(new LinkedList<Object>(Collections.singleton(e)));
                    } else {
                        for ( final List<Object> state : substates ) {
                            List<Object> newState = new LinkedList<Object>();
                            newState.add(e);
                            newState.addAll(state);
                            newStates.add(newState);
                        }
                    }
                }
                return newStates;
            }
        }
    }

    private class IntegerStratifier implements Stratifier<Object> {
        final List<Object> integers;

        private IntegerStratifier(final List<Object> integers) {
            this.integers = integers;
        }
        
        @Override
        public List<Object> getAllStates() {
            return integers;
        }
    }

    @DataProvider(name = "StratificationStatesTestProvider")
    @SuppressWarnings("unchecked")
    public Object[][] makeStratificationStatesTestProvider() {
        new StratificationStatesTestProvider(Arrays.asList(0));
        new StratificationStatesTestProvider(Arrays.asList(0, 1));
        new StratificationStatesTestProvider(Arrays.asList(0, 1), Arrays.asList(2, 3));
        new StratificationStatesTestProvider(Arrays.asList(0, 1), Arrays.asList(2, 3), Arrays.asList(4, 5));
        new StratificationStatesTestProvider(Arrays.asList(0, 1), Arrays.asList(2, 3, 4), Arrays.asList(5, 6));
        new StratificationStatesTestProvider(Arrays.asList(0, 1), Arrays.asList(2, 3, 4, 5), Arrays.asList(6));
        new StratificationStatesTestProvider(Arrays.asList(0, 1), Arrays.asList(2, 3, 4, 5), Arrays.asList(6, 7));
        new StratificationStatesTestProvider(Arrays.asList(0, 1), Arrays.asList(2, 3), Arrays.asList(4, 5), Arrays.asList(6, 7));
        return StratificationStatesTestProvider.getTests(StratificationStatesTestProvider.class);
    }
    
    private final StratificationManager<IntegerStratifier, Integer> createManager(StratificationStatesTestProvider cfg) {
        final StratificationManager<IntegerStratifier, Integer> manager = new StratificationManager<>(cfg.getStateSpaceList());
        List<Integer> values = cfg.values();
        for ( int i = 0; i < cfg.nStates; i++ )
            manager.set(i, values.get(i));
        
        Assert.assertEquals(manager.values(), values, "Values not equal");
        
        return manager;
    }

    @Test(dataProvider = "StratificationStatesTestProvider")
    public void testLeafCount(StratificationStatesTestProvider cfg) {
        final StratificationManager<IntegerStratifier, Integer> stratificationManager = createManager(cfg);
        
        Assert.assertEquals(stratificationManager.size(), cfg.nStates);
        
        int nLeafs = 0;
        for ( final StratNode<StratificationManagerUnitTest.IntegerStratifier> node : stratificationManager.getRoot() ) {
            if ( node.isLeaf() )
                nLeafs++;
        }
        Assert.assertEquals(nLeafs, cfg.nStates, "Unexpected number of leaves");
    }

    @Test(dataProvider = "StratificationStatesTestProvider")
    public void testKeys(StratificationStatesTestProvider cfg) {
        final StratificationManager<IntegerStratifier, Integer> stratificationManager = createManager(cfg);
        final Set<Integer> seenKeys = new HashSet<Integer>(cfg.nStates);
        for ( final StratNode<StratificationManagerUnitTest.IntegerStratifier> node : stratificationManager.getRoot() ) {
            if ( node.isLeaf() ) {
                Assert.assertFalse(seenKeys.contains(node.getKey()), "Already seen the key");
                seenKeys.add(node.getKey());
            }
        }
    }

    @Test(dataProvider = "StratificationStatesTestProvider")
    public void testFindSingleKeys(StratificationStatesTestProvider cfg) {
        final StratificationManager<IntegerStratifier, Integer> stratificationManager = createManager(cfg);
        final Set<Integer> seenKeys = new HashSet<Integer>(cfg.nStates);
        for ( List<Object> state : cfg.getAllCombinations() ) {
            final int key = stratificationManager.getKey(state);
            Assert.assertFalse(seenKeys.contains(key), "Already saw state mapping to this key");
            Assert.assertTrue(stratificationManager.containsKey(state));
            seenKeys.add(key);

            // test value
            Assert.assertEquals(stratificationManager.get(key), cfg.values().get(key));
            Assert.assertEquals(stratificationManager.get(state), cfg.values().get(key));

            state.set(0, 12345); // not present
            Assert.assertEquals(stratificationManager.getKey(state), -1);
            Assert.assertFalse(stratificationManager.containsKey(state));
        }
    }

    @Test(dataProvider = "StratificationStatesTestProvider")
    public void testFindMultipleKeys(StratificationStatesTestProvider cfg) {
        final StratificationManager<IntegerStratifier, Integer> stratificationManager = createManager(cfg);
        final List<List<Object>> states = new ArrayList<List<Object>>(cfg.allStates);
        final Set<Integer> keys = stratificationManager.getKeys(states);
        Assert.assertEquals(keys.size(), cfg.nStates, "Find all states didn't find all of the expected unique keys");

        final Queue<List<Object>> combinations = cfg.getAllCombinations();
        while ( ! combinations.isEmpty() ) {
            List<Object> first = combinations.poll();
            List<Object> second = combinations.peek();
            if ( second != null ) {
                List<List<Object>> combined = StratificationManager.combineStates(first, second);
                int nExpectedKeys = nCombinations(combined);

                final int key1 = stratificationManager.getKey(first);
                final int key2 = stratificationManager.getKey(second);
                final Set<Integer> keysCombined = stratificationManager.getKeys(combined);
            
                Assert.assertTrue(keysCombined.contains(key1), "couldn't find key in data set");
                Assert.assertTrue(keysCombined.contains(key2), "couldn't find key in data set");
                
                Assert.assertEquals(keysCombined.size(), nExpectedKeys);
            }
        }
    }

    @Test(dataProvider = "StratificationStatesTestProvider")
    public void testMapSet(StratificationStatesTestProvider cfg) {
        final StratificationManager<IntegerStratifier, Integer> stratificationManager = createManager(cfg);
        stratificationManager.set(0, -1);
        Assert.assertEquals((int)stratificationManager.get(0), -1);
    }

    @Test(dataProvider = "StratificationStatesTestProvider")
    public void testStratifierByKey(StratificationStatesTestProvider cfg) {
        final StratificationManager<IntegerStratifier, Integer> manager = createManager(cfg);
        for ( int key = 0; key < cfg.nStates; key++ ) {
            List<Pair<IntegerStratifier, Object>> stratsAndStates = manager.getStratsAndStatesForKey(key);
            final List<Object> strats = manager.getStatesForKey(key);
            Assert.assertEquals((int)manager.get(strats), key, "Key -> strats -> key failed to return same key");

            for ( int i = 0; i < strats.size(); i++ ) {
                Assert.assertEquals(stratsAndStates.get(i).getRight(), strats.get(i), "Strats and StratsAndStates differ");
            }
        }
    }

    private <T> int nCombinations(final Collection<T>[] options) {
        int nStates = 1;
        for ( Collection<T> states : options ) {
            nStates *= states.size();
        }
        return nStates;
    }

    private <T> int nCombinations(final List<List<T>> options) {
        int nStates = 1;
        for ( Collection<T> states : options ) {
            nStates *= states.size();
        }
        return nStates;
    }
}