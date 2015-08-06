package org.broadinstitute.hellbender.utils;

import htsjdk.samtools.util.Locatable;

import java.util.Iterator;
import java.util.Set;

/**
 * Created by edwardk on 8/6/15.
 */
public class IntervalTree<V> implements Iterable<IntervalTree.Node<V>>
{

    @Override
    public Iterator<Node<V>> iterator() {
        return null;
    }

    public static class Node<V> {

    }


    //TODO: methods to be implemented

    public boolean overlaps(Locatable locatable) {
        return false;
    }

    Set<Locatable> getOverlapping(Locatable locatable) {
        return null;
    }

}
