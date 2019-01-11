package org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications.manager;

import org.broadinstitute.hellbender.exceptions.GATKException;

import java.util.*;

/**
 * Helper class representing a tree of stratification splits, where leaf nodes
 * are given a unique integer key starting at 0 and incrementing up to the
 * number of leaves in the tree.  This allows you to use this tree to produce
 * a key to map into an array index mapped data structure.
 *
 * Suppose I have to strats, each with two values: A = 1, 2 and B = 3, 4
 *
 * This data structure creates a tree such as:
 *
 * root -> A -> 1 -> B -> 3  : 0
 *                |- B -> 4  : 1
 *      |- A -> 2 -> B -> 3  : 2
 *                |- B -> 4  : 3
 *
 * This code allows us to efficiently look up a state key (A=2, B=3) and map it
 * to a specific key (an integer) that's unique over the tree
 *
 * Note the structure of this tree is that the keys are -1 for all internal nodes, and
 * leafs are the only nodes with meaningful keys.  So for a tree with 2N nodes N of these
 * will be internal, with no keys, and meaningful maps from states -> subtrees.  The
 * other N nodes are leafs, with meaningful keys, empty maps, and null stratification objects
 *
 * @author Mark DePristo
 * @since 3/27/12
 */
class StratNode<T extends Stratifier<Object>> implements Iterable<StratNode<T>> {
    int key = -1;
    final T stratifier;
    final Map<Object, StratNode<T>> subnodes; // NOTE, because we don't iterate our best option is a HashMap

    protected StratNode() {
        this.subnodes = Collections.emptyMap();
        this.stratifier = null;
    }

    protected StratNode(final T stratifier, final Map<Object, StratNode<T>> subnodes) {
        this.stratifier = stratifier;
        // important to reallocate an unmodififable hashmap with this specific size for space and safety
        this.subnodes = Collections.unmodifiableMap(new HashMap<Object, StratNode<T>>(subnodes));
    }

    public void setKey(final int key) {
        if ( ! isLeaf() )
            throw new GATKException("Cannot set key of non-leaf node");
        this.key = key;
    }

    public int find(final List<Object> states, int offset) {
        if ( isLeaf() ) // we're here!
            return key;
        else {
            final Object state = states.get(offset);
            StratNode<T> subnode = subnodes.get(state);
            if ( subnode == null )
                return -1;
            else
                return subnode.find(states, offset+1);
        }
    }

    public void find(final List<List<Object>> multipleStates, final int offset, final HashSet<Integer> keys) {
        if ( isLeaf() ) // we're here!
            keys.add(key);
        else {
            for ( final Object state : multipleStates.get(offset) ) {
                // loop over all of the states at this offset
                final StratNode<T> subnode = subnodes.get(state);
                if ( subnode == null )
                    throw new GATKException("Couldn't find state for " + state + " at node " + this);
                else
                    subnode.find(multipleStates, offset+1, keys);
            }
        }
    }

    public int getKey() {
        if ( ! isLeaf() )
            throw new GATKException("Cannot get key of non-leaf node");
        else
            return key;
    }

    protected Map<Object, StratNode<T>> getSubnodes() {
        return subnodes;
    }

    public int size() {
        if ( isLeaf() )
            return 1;
        else {
            return subnodes.values().iterator().next().size() * subnodes.size();
        }
    }

    public T getSetOfStates() {
        return stratifier;
    }

    /**
     * @return true if this node is a leaf
     */
    public boolean isLeaf() {
        return stratifier == null;
    }

    /**
     * Returns an iterator over this node and all subnodes including internal and leaf nodes
     * @return
     */
    @Override
    public Iterator<StratNode<T>> iterator() {
        return new StratNodeIterator<T>(this);
    }
}
