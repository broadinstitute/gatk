package org.broadinstitute.hellbender.utils.pileup;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

abstract class PileupElementTracker<PE extends PileupElement> implements Iterable<PE> {
    public abstract int size();

    /**
     * Iterate through the PEs here, but in any order, which may improve performance
     * if you don't care about the underlying order the reads are coming to you in.
     * @return an iteratable over all pileup elements in this tracker
     */
    public abstract Iterable<PE> unorderedIterable();

    /**
     * Same as @see #unorderedIterable but the actual iterator itself
     * @return
     */
    public Iterator<PE> unorderedIterator() { return unorderedIterable().iterator(); }

    public abstract PileupElementTracker<PE> copy();
}

