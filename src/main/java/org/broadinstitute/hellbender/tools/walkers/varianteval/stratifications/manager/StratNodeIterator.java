package org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications.manager;

import org.broadinstitute.hellbender.exceptions.GATKException;

import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Queue;

/**
 * Helper class for creating iterators over all nodes in the stratification tree
 *
 * @author Mark DePristo
 * @since 3/27/12
 */
class StratNodeIterator<T extends Stratifier<Object>> implements Iterator<StratNode<T>> {
    Queue<Iterator<StratNode<T>>> iterators = new LinkedList<Iterator<StratNode<T>>>();
    Iterator<StratNode<T>> currentIterator;

    StratNodeIterator(final StratNode<T> root) {
        currentIterator = Collections.singleton(root).iterator();
        for ( final StratNode<T> subNode : root.subnodes.values() )
            iterators.add(new StratNodeIterator<T>(subNode));
    }

    @Override
    public boolean hasNext() {
        return currentIterator.hasNext() || ! iterators.isEmpty();
    }

    @Override
    public StratNode<T> next() {
        if ( currentIterator.hasNext() )
            return currentIterator.next();
        else if ( ! iterators.isEmpty() ) {
            currentIterator = iterators.poll();
            return next();
        } else {
            throw new IllegalStateException("Next called on empty iterator");
        }
    }

    @Override
    public void remove() {
        throw new GATKException("Cannot remove from StratNode iterator");
    }
}
