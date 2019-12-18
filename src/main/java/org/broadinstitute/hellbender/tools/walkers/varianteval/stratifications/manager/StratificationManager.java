package org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications.manager;

import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * Represents the full state space of all stratification combinations
 *
 * @author Mark DePristo
 * @since 3/27/12
 */
public class StratificationManager<K extends Stratifier<Object>, V> implements Map<List<Object>, V> {
    private final StratNode<K> root;
    private final int size;

    private final ArrayList<K> stratifiers;

    // values associated with each key
    private final ArrayList<V> valuesByKey;
    private final ArrayList<List<Object>> stratifierValuesByKey;
    private final ArrayList<String> keyStrings;

    // -------------------------------------------------------------------------------------
    //
    // creating the manager
    //
    // -------------------------------------------------------------------------------------

    /**
     * Create a new StratificationManager with nodes to store data for all combinations
     * of the ordered list of strats
     *
     * @param strats ordered list of stratifications to representation
     */
    public StratificationManager(final List<K> strats) {
        this.stratifiers = new ArrayList<K>(strats);

        // construct and store the full tree of strats
        this.root = buildStratificationTree(new LinkedList<K>(strats));
        // assign the linear key ordering to the leafs
        assignKeys(root);

        // cache the size, and check for a bad state
        this.size = root.size();
        if ( this.size == 0 )
            throw new GATKException("Size == 0 in StratificationManager");

        // prepare the assocated data vectors mapping from key -> data
        this.valuesByKey = new ArrayList<V>(size());
        this.stratifierValuesByKey = new ArrayList<List<Object>>(size());
        this.keyStrings = new ArrayList<String>(size());
        for ( int i = 0; i < size(); i++ ) {
            this.valuesByKey.add(null);
            this.stratifierValuesByKey.add(null);
            this.keyStrings.add(null);
        }

        assignStratifierValuesByKey(root);
    }

    /**
     * Recursive construction helper for main constructor that fills into the
     * complete tree of StratNodes.  This function returns the complete tree
     * suitable for associating data with each combinatino of keys.  Note
     * that the tree is not fully complete as the keys are not yet set for
     * each note (see assignStratifierValuesByKey)
     *
     * @param strats
     * @return
     */
    private StratNode<K> buildStratificationTree(final Queue<K> strats) {
        final K first = strats.poll();
        if ( first == null ) {
            // we are at a leaf
            return new StratNode<K>();
        } else {
            // we are in the middle of the tree
            final Collection<Object> states = first.getAllStates();
            
            if ( states.isEmpty() )
                throw new GATKException("State " + first + " is empty!");
            
            final LinkedHashMap<Object, StratNode<K>> subNodes = new LinkedHashMap<Object, StratNode<K>>(states.size());
            for ( final Object state : states ) {
                // have to copy because poll modifies the queue
                final Queue<K> copy = new LinkedList<K>(strats);
                subNodes.put(state, buildStratificationTree(copy));
            }
            return new StratNode<K>(first, subNodes);
        }
    }

    /**
     * Set the key for each leaf from root, in order from 0 to N - 1 for N leaves in the tree
     * @param root
     */
    private void assignKeys(final StratNode<K> root) {
        int key = 0;
        for ( final StratNode<K> node : root ) {
            if ( node.isLeaf() )
                node.setKey(key++);
        }
    }

    /**
     * Entry point to recursive tool that fills in the list of state values corresponding
     * to each key.  After this function is called you can map from key -> List of StateValues
     * instead of walking the tree to find the key and reading the list of state values
     *
     * @param root
     */
    private void assignStratifierValuesByKey(final StratNode<K> root) {
        assignStratifierValuesByKey(root, new LinkedList<>());

        // do a last sanity check that no key has null value after assigning
        for ( List<Object> stateValues : stratifierValuesByKey )
            if ( stateValues == null )
                throw new GATKException("Found a null state value set that's null");
    }

    private void assignStratifierValuesByKey(final StratNode<K> node, final LinkedList<Object> states) {
        if ( node.isLeaf() ) { // we're here!
            if ( states.isEmpty() )
                throw new GATKException("Found a leaf node with an empty state values vector");
            stratifierValuesByKey.set(node.getKey(), Collections.unmodifiableList(new ArrayList<Object>(states)));
        } else {
            for ( Map.Entry<Object, StratNode<K>> entry : node.getSubnodes().entrySet() ) {
                final LinkedList<Object> newStates = new LinkedList<Object>(states);
                newStates.addLast(entry.getKey());
                assignStratifierValuesByKey(entry.getValue(), newStates);
            }
        }
    }
    
    // -------------------------------------------------------------------------------------
    //
    // simple accessors
    //
    // -------------------------------------------------------------------------------------

    /**
     * How many states are held in this stratification manager?
     * @return
     */
    public int size() {
        return size;
    }

    protected StratNode<K> getRoot() {
        return root;
    }

    public List<K> getStratifiers() {
        return stratifiers;
    }

    // -------------------------------------------------------------------------------------
    //
    // mapping from states -> keys
    //
    // -------------------------------------------------------------------------------------

    public int getKey(final List<Object> states) {
        return root.find(states, 0);
    }

    public Set<Integer> getKeys(final List<List<Object>> allStates) {
        final HashSet<Integer> keys = new HashSet<Integer>();
        root.find(allStates, 0, keys);
        return keys;
    }

    public List<Object> getStatesForKey(final int key) {
        final List<Object> states = new ArrayList<Object>(stratifiers.size());
        for ( int i = 0; i < stratifiers.size(); i++ ) {
            final Object stratValue = stratifierValuesByKey.get(key).get(i);
            states.add(stratValue);
        }
        return states;
    }

    public List<Pair<K, Object>> getStratsAndStatesForKey(final int key) {
        final List<Pair<K, Object>> states = new ArrayList<Pair<K, Object>>(stratifiers.size());
        for ( int i = 0; i < stratifiers.size(); i++ ) {
            final K strat = stratifiers.get(i);
            final Object stratValue = stratifierValuesByKey.get(key).get(i);
            states.add(Pair.of(strat, stratValue));
        }
        return states;
    }

    public String getStratsAndStatesStringForKey(final int key) {
        if ( keyStrings.get(key) == null ) {
            StringBuilder b = new StringBuilder();
            for ( int i = 0; i < stratifiers.size(); i++ ) {
                final K strat = stratifiers.get(i);
                final Object stratValue = stratifierValuesByKey.get(key).get(i);
                b.append(strat.toString()).append(":").append(stratValue.toString());
            }
            keyStrings.set(key, b.toString());
        }
        
        return keyStrings.get(key);
    }

    // -------------------------------------------------------------------------------------
    //
    // valuesByKey
    //
    // -------------------------------------------------------------------------------------

    @Override
    public ArrayList<V> values() {
        return valuesByKey;
    }
    
    public Collection<V> values(List<List<Object>> states) {
        // TODO -- SHOULD BE INLINE TO AVOID CREATING LIST OF KEYS JUST TO ITERATE OVER IT
        Collection<V> vals = new LinkedList<V>();
        for ( int key : getKeys(states) ) 
            vals.add(get(key));
        return vals;
    }

    public void set(final int key, final V value) {
        valuesByKey.set(key, value);
    }

    public V get(final int key) {
        return valuesByKey.get(key);
    }

    public V get(final List<Object> states) {
        return get(getKey(states));
    }

    @Override
    public V get(final Object o) {
        return get(cast((List)o));
    }

    @Override
    public boolean isEmpty() {
        return false;
    }

    public boolean containsKey(final List<Object> o) {
        return getKey(o) != -1;
    }

    @Override
    public boolean containsKey(final Object o) {
        return containsKey(cast((List)o));
    }

    //TODO: improve this
    private List<Object> cast(List<?> list){
        return list.stream()
                .map(Object.class::cast)
                .collect(Collectors.toList());
    }

    @Override
    public boolean containsValue(final Object o) {
        throw new GATKException("containsValue() not implemented for StratificationManager");
    }

    @Override
    public V put(final List<Object> objects, final V v) {
        throw new GATKException("put() not implemented for StratificationManager");
    }

    @Override
    public V remove(final Object o) {
        throw new GATKException("remove() not implemented for StratificationManager");
    }

    @Override
    public void putAll(final Map<? extends List<Object>, ? extends V> map) {
        throw new GATKException("clear() not implemented for StratificationManager");
    }

    @Override
    public void clear() {
        throw new GATKException("clear() not implemented for StratificationManager");
    }

    @Override
    public Set<List<Object>> keySet() {
        throw new GATKException("Not yet implemented");
    }

    @Override
    public Set<Entry<List<Object>, V>> entrySet() {
        throw new GATKException("Not yet implemented");
    }

    // -------------------------------------------------------------------------------------
    //
    // utilities
    //
    // -------------------------------------------------------------------------------------

    public static List<List<Object>> combineStates(final List<Object> first, final List<Object> second) {
        final List<List<Object>> combined = new ArrayList<List<Object>>(first.size());
        for ( int i = 0; i < first.size(); i++ ) {
            final Object firstI = first.get(i);
            final Object secondI = second.get(i);
            if ( firstI.equals(secondI) ) 
                combined.add(Collections.singletonList(firstI));
            else 
                combined.add(Arrays.asList(firstI, secondI));
        }
        return combined;
    }

    public interface Combiner<V> {
        /** take two values of type V and return a combined value of type V */
        public V combine(final V lhs, final V rhs);
    }

    /**
     * Remaps the stratifications from one stratification set to another, combining
     * the values in V according to the combiner function.
     *
     * stratifierToReplace defines a set of states S1, while newStratifier defines
     * a new set S2.  remappedStates is a map from all of S1 into at least some of
     * S2.  This function creates a new, fully initialized manager where all of the
     * data in this new manager is derived from the original data in this object
     * combined according to the mapping remappedStates.  When multiple
     * elements of S1 can map to the same value in S2, these are sequentially
     * combined by the function combiner.  Suppose for example at states s1, s2, and
     * s3 all map to N1.  Eventually the value associated with state N1 would be
     *
     *   value(N1) = combine(value(s1), combine(value(s2), value(s3))
     *
     * in some order for s1, s2, and s3, which is not defined.  Note that this function
     * only supports combining one stratification at a time, but in principle a loop over
     * stratifications and this function could do the multi-dimensional collapse.
     *
     * @param stratifierToReplace
     * @param newStratifier
     * @param combiner
     * @param remappedStates
     * @return
     */
    public StratificationManager<K, V> combineStrats(final K stratifierToReplace,
                                                     final K newStratifier,
                                                     final Combiner<V> combiner,
                                                     final Map<Object, Object> remappedStates) {
        // make sure the mapping is reasonable
        if ( ! newStratifier.getAllStates().containsAll(remappedStates.values()) )
            throw new GATKException("combineStrats: remapped states contains states not found in newStratifer state set");

        if ( ! remappedStates.keySet().containsAll(stratifierToReplace.getAllStates()) )
            throw new GATKException("combineStrats: remapped states missing mapping for some states");

        // the new strats are the old ones with the single replacement
        final List<K> newStrats = new ArrayList<K>(getStratifiers());
        final int stratOffset = newStrats.indexOf(stratifierToReplace);
        if ( stratOffset == -1 )
            throw new GATKException("Could not find strat to replace " + stratifierToReplace + " in existing strats " + newStrats);
        newStrats.set(stratOffset, newStratifier);

        // create an empty but fully initialized new manager
        final StratificationManager<K, V> combined = new StratificationManager<K, V>(newStrats);

        // for each key, get its state, update it according to the map, and update the combined manager
        for ( int key = 0; key < size(); key++ ) {
            // the new state is just the old one with the replacement
            final List<Object> newStates = new ArrayList<Object>(getStatesForKey(key));
            final Object oldState = newStates.get(stratOffset);
            final Object newState = remappedStates.get(oldState);
            newStates.set(stratOffset, newState);

            // look up the new key given the new state
            final int combinedKey = combined.getKey(newStates);
            if ( combinedKey == -1 ) throw new GATKException("Couldn't find key for states: " + Utils.join(",", newStates));

            // combine the old value with whatever new value is in combined already
            final V combinedValue = combiner.combine(combined.get(combinedKey), get(key));

            // update the value associated with combined key
            combined.set(combinedKey, combinedValue);
        }

        return combined;
    }
}