package org.broadinstitute.hellbender.tools.spark.sv.utils;

import java.util.ConcurrentModificationException;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.Objects;

/**
 * A Red-Black tree with intervals for keys.
 * Intervals are kept in sorted order first by start value and then by end value. This cannot be overridden.
 * Not thread-safe, and cannot be made so.  You must synchronize externally.
 * <p>
 * There's some weird stuff about sentinel values for the put and remove methods.  Here's what's up with that:
 * When you update the value associated with some interval by doing a put with an interval that's already in the
 * tree, the old value is returned to you.  But maybe you've put some nulls into the tree as values.  (That's legal.)
 * In that case, when you get a null value returned by put you can't tell whether the interval was inserted into the tree
 * and there was no old value to return to you or whether you just updated an existing interval that had a null value
 * associated with it.  (Both situations return null.)  IF you're inserting nulls as values, and IF you need to be able
 * to tell whether the put operation did an insert or an update, you can do a special thing so that you can distinguish
 * these cases:  set the sentinel value for the tree to some singleton object that you never ever use as a legitimate
 * value.  Then when you call put you'll get your sentinel value back for an insert, but you'll get null back for an
 * update of a formerly-null value.  Same thing happens for remove:  set the sentinel IF you've used nulls for values,
 * and IF you need to be able to tell the difference between remove not finding the interval and remove removing an
 * interval that has a null value associated with it.
 * If you're not using nulls as values, or if you don't care to disambiguate these cases, then just forget about
 * all this weirdness.  The sentinel value is null by default, so put and remove will behave like you might expect them
 * to if you're not worrying about this stuff:  they'll return null for novel insertions and failed deletions.
 */
public final class SVIntervalTree<V> implements Iterable<SVIntervalTree.Entry<V>> {
    private Node<V> root;
    private V sentinel;

    /**
     * Return the number of intervals in the tree.
     *
     * @return The number of intervals.
     */
    public int size() {
        return root == null ? 0 : root.getSize();
    }

    /**
     * Remove all entries.
     */
    public void clear() {
        root = null;
    }

    /**
     * Put a new interval into the tree (or update the value associated with an existing interval).
     * If the interval is novel, the special sentinel value (which is null by default) is returned.
     *
     * @param interval The interval.
     * @param value    The associated value.
     * @return The old value associated with that interval, or the sentinel value.
     */
    public V put( final SVInterval interval, final V value ) {
        V result = sentinel;

        if ( root == null ) {
            root = new Node<>(interval, value);
        } else {
            Node<V> parent = null;
            Node<V> node = root;
            int cmpVal = 0;

            while ( node != null ) {
                parent = node; // last non-null node
                cmpVal = interval.compareTo(node.getInterval());
                if ( cmpVal == 0 ) {
                    break;
                }
                node = cmpVal < 0 ? node.getLeft() : node.getRight();
            }

            if ( cmpVal == 0 ) {
                result = parent.setValue(value);
            } else if ( cmpVal < 0 ) {
                root = parent.insertLeft(interval, value, root);
            } else {
                root = parent.insertRight(interval, value, root);
            }
        }

        return result;
    }

    /**
     * Remove an interval from the tree.
     * If the interval is not found, the special sentinel value (which is null by default) is returned.
     *
     * @param interval The interval to remove.
     * @return The value associated with the deleted interval, or the sentinel value.
     */
    public V remove( final SVInterval interval ) {
        V result = sentinel;
        Node<V> node = root;

        while ( node != null ) {
            final int cmpVal = interval.compareTo(node.getInterval());
            if ( cmpVal == 0 ) {
                result = node.getValue();
                root = node.remove(root);
                break;
            }
            node = cmpVal < 0 ? node.getLeft() : node.getRight();
        }

        return result;
    }

    /**
     * Find an interval.
     *
     * @param interval The interval sought.
     * @return The Node that represents that interval, or null.
     */
    public Entry<V> find( final SVInterval interval ) {
        Node<V> node = root;

        while ( node != null ) {
            final int cmpVal = interval.compareTo(node.getInterval());
            if ( cmpVal == 0 ) {
                break;
            }
            node = cmpVal < 0 ? node.getLeft() : node.getRight();
        }

        return node;
    }

    /**
     * Find the nth interval in the tree.
     *
     * @param idx The rank of the interval sought (from 0 to size()-1).
     * @return The Node that represents the nth interval.
     */
    public Entry<V> findByIndex( final int idx ) {
        return Node.findByRank(root, idx + 1);
    }

    /**
     * Find the rank of the specified interval.  If the specified interval is not in the
     * tree, then -1 is returned.
     *
     * @param interval The interval for which the index is sought.
     * @return The rank of that interval, or -1.
     */
    public int getIndex( final SVInterval interval ) { return Node.getRank(root, interval) - 1; }

    /**
     * Find the least interval in the tree.
     *
     * @return The earliest interval, or null if the tree is empty.
     */
    public Entry<V> min() {
        Node<V> result = null;
        Node<V> node = root;

        while ( node != null ) {
            result = node;
            node = node.getLeft();
        }

        return result;
    }

    /**
     * Find the earliest interval in the tree greater than or equal to the specified interval.
     *
     * @param interval The interval sought.
     * @return The earliest >= interval, or null if there is none.
     */
    public Entry<V> min( final SVInterval interval ) {
        Node<V> result = null;
        Node<V> node = root;
        int cmpVal = 0;

        while ( node != null ) {
            result = node;
            cmpVal = interval.compareTo(node.getInterval());
            if ( cmpVal == 0 ) {
                break;
            }
            node = cmpVal < 0 ? node.getLeft() : node.getRight();
        }

        if ( cmpVal > 0 ) {
            result = result.getNext();
        }

        return result;
    }

    /**
     * Find the earliest interval in the tree that overlaps the specified interval.
     *
     * @param interval The interval sought.
     * @return The earliest overlapping interval, or null if there is none.
     */
    public Entry<V> minOverlapper( final SVInterval interval ) {
        Node<V> result = null;
        Node<V> node = root;

        if ( node != null && !node.getMaxEndInterval().isUpstreamOf(interval) ) {
            while ( true ) {
                if ( node.getInterval().overlaps(interval) ) {
                    // this node overlaps.  however, there might be an earlier overlapper down the left sub-tree.
                    // no need to consider the right sub-tree:  even if there's an overlapper, if won't be minimal
                    result = node;
                    node = node.getLeft();
                    // if no left sub-tree or all nodes end too early, then break
                    if ( node == null || node.getMaxEndInterval().isUpstreamOf(interval) ) {
                        break;
                    }
                } else { // no overlap.  if there might be a left sub-tree overlapper, consider the left sub-tree.
                    final Node<V> left = node.getLeft();
                    if ( left != null && !left.getMaxEndInterval().isUpstreamOf(interval) ) {
                        node = left;
                    } else { // left sub-tree cannot contain an overlapper.  consider the right sub-tree.
                        // if everything in the right sub-tree is past the end of the query interval, then break
                        if ( interval.isUpstreamOf(node.getInterval()) ) {
                            break;
                        }

                        node = node.getRight();
                        // if no right sub-tree or all nodes end too early, then break
                        if ( node == null || node.getMaxEndInterval().isUpstreamOf(interval) ) {
                            break;
                        }
                    }
                }
            }
        }

        return result;
    }

    /**
     * Find the greatest interval in the tree.
     *
     * @return The latest interval, or null if the tree is empty.
     */
    public Entry<V> max() {
        Node<V> result = null;
        Node<V> node = root;

        while ( node != null ) {
            result = node;
            node = node.getRight();
        }

        return result;
    }

    /**
     * Find the latest interval in the tree less than or equal to the specified interval.
     *
     * @param interval The interval sought.
     * @return The latest <= interval, or null if there is none.
     */
    public Entry<V> max( final SVInterval interval ) {
        Node<V> result = null;
        Node<V> node = root;
        int cmpVal = 0;

        while ( node != null ) {
            result = node;
            cmpVal = interval.compareTo(node.getInterval());
            if ( cmpVal == 0 ) {
                break;
            }

            node = cmpVal < 0 ? node.getLeft() : node.getRight();
        }

        if ( cmpVal < 0 ) {
            result = result.getPrev();
        }

        return result;
    }

    /**
     * Return the interval having the largest ending value.
     * This will be null if the tree is empty.
     */
    public SVInterval maxEnd() {
        return root == null ? null : root.getMaxEndInterval();
    }

    /**
     * Return an iterator over the entire tree.
     *
     * @return An iterator.
     */
    public Iterator<Entry<V>> iterator() { return new FwdIterator((Node<V>)min()); }

    /**
     * Return an iterator over all intervals greater than or equal to the specified interval.
     *
     * @param interval The minimum interval.
     * @return An iterator.
     */
    public Iterator<Entry<V>> iterator( final SVInterval interval ) { return new FwdIterator((Node<V>)min(interval)); }

    /**
     * Return an iterator over all intervals overlapping the specified interval.
     *
     * @param interval Interval to overlap.
     * @return An iterator.
     */
    public Iterator<Entry<V>> overlappers( final SVInterval interval ) { return new OverlapIterator(interval); }

    /**
     * Return an iterator over the entire tree that returns intervals in reverse order.
     *
     * @return An iterator.
     */
    public Iterator<Entry<V>> reverseIterator() { return new RevIterator((Node<V>)max()); }

    /**
     * Return an iterator over all intervals less than or equal to the specified interval, in reverse order.
     *
     * @param interval The maximum interval.
     * @return An iterator.
     */
    public Iterator<Entry<V>> reverseIterator( final SVInterval interval ) {
        return new RevIterator((Node<V>)max(interval));
    }

    /**
     * Get the special sentinel value that will be used to signal novelty when putting a new interval
     * into the tree, or to signal "not found" when removing an interval.  This is null by default.
     *
     * @return The sentinel value.
     */
    public V getSentinel() {
        return sentinel;
    }

    /**
     * Set the special sentinel value that will be used to signal novelty when putting a new interval
     * into the tree, or to signal "not found" when removing an interval.
     *
     * @param sentinel The new sentinel value.
     * @return The old sentinel value.
     */
    public V setSentinel( final V sentinel ) {
        final V result = this.sentinel;
        this.sentinel = sentinel;
        return result;
    }

    void removeNode( final Node<V> node ) {
        root = node.remove(root);
    }

    public interface Entry<V1> {
        SVInterval getInterval();
        V1 getValue();
        V1 setValue( final V1 value );
    }

    static final class Node<V1> implements Entry<V1> {
        private final SVInterval interval;
        private V1 value;
        private Node<V1> parent;
        private Node<V1> left;
        private Node<V1> right;
        private int size;
        private SVInterval maxEndInterval; // interval in this sub-tree having the greatest endpoint
        private boolean isBlack;

        /** make a root node */
        Node( final SVInterval interval, final V1 value ) {
            this.interval = interval;
            this.value = value;
            size = 1;
            maxEndInterval = interval;
            isBlack = true;
        }

        /** make a leaf node */
        Node( final Node<V1> parent, final SVInterval interval, final V1 value ) {
            this.interval = interval;
            this.value = value;
            this.parent = parent;
            size = 1;
            maxEndInterval = interval;
        }

        @Override
        public SVInterval getInterval() { return interval; }

        @Override
        public V1 getValue() {
            return value;
        }

        @Override
        public V1 setValue( final V1 value ) {
            final V1 result = this.value;
            this.value = value;
            return result;
        }

        int getSize() {
            return size;
        }

        SVInterval getMaxEndInterval() {
            return maxEndInterval;
        }

        Node<V1> getLeft() {
            return left;
        }

        Node<V1> insertLeft( final SVInterval interval, final V1 value, final Node<V1> root ) {
            left = new Node<>(this, interval, value);
            return insertFixup(left, root);
        }

        Node<V1> getRight() {
            return right;
        }

        Node<V1> insertRight( final SVInterval interval, final V1 value, final Node<V1> root ) {
            right = new Node<>(this, interval, value);
            return insertFixup(right, root);
        }

        Node<V1> getNext() {
            Node<V1> result;

            if ( right != null ) {
                result = right;
                while ( result.left != null ) {
                    result = result.left;
                }
            } else {
                Node<V1> node = this;
                result = parent;
                while ( result != null && node == result.right ) {
                    node = result;
                    result = result.parent;
                }
            }

            return result;
        }

        Node<V1> getPrev() {
            Node<V1> result;

            if ( left != null ) {
                result = left;
                while ( result.right != null ) {
                    result = result.right;
                }
            } else {
                Node<V1> node = this;
                result = parent;
                while ( result != null && node == result.left ) {
                    node = result;
                    result = result.parent;
                }
            }

            return result;
        }

        boolean wasRemoved() {
            return size == 0;
        }

        Node<V1> remove( final Node<V1> initialRoot ) {
            Node<V1> root = initialRoot;
            if ( size == 0 ) {
                throw new IllegalStateException("Entry was already removed.");
            }

            if ( left == null ) {
                if ( right == null ) { // no children
                    if ( parent == null ) {
                        root = null;
                    } else if ( parent.left == this ) {
                        parent.left = null;
                        fixup(parent);

                        if ( isBlack ) {
                            root = removeFixup(parent, null, root);
                        }
                    } else {
                        parent.right = null;
                        fixup(parent);

                        if ( isBlack ) {
                            root = removeFixup(parent, null, root);
                        }
                    }
                } else { // single child on right
                    root = spliceOut(right, root);
                }
            } else if ( right == null ) { // single child on left
                root = spliceOut(left, root);
            } else { // two children
                final Node<V1> next = getNext();
                root = next.remove(root);

                // put next into tree in same position as this, effectively removing this
                if ( (next.parent = parent) == null ) {
                    root = next;
                } else if ( parent.left == this ) {
                    parent.left = next;
                } else {
                    parent.right = next;
                }

                if ( (next.left = left) != null ) {
                    left.parent = next;
                }

                if ( (next.right = right) != null ) {
                    right.parent = next;
                }

                next.isBlack = isBlack;
                next.size = size;
                fixup(next);
            }

            size = 0;
            return root;
        }

        static <V1> Node<V1> getNextOverlapper( final Node<V1> startingNode, final SVInterval interval ) {
            Node<V1> node = startingNode;
            do {
                Node<V1> nextNode = node.right;
                if ( nextNode != null && !nextNode.maxEndInterval.isUpstreamOf(interval) ) {
                    node = nextNode;
                    while ( (nextNode = node.left) != null && !nextNode.maxEndInterval.isUpstreamOf(interval) )
                        node = nextNode;
                } else {
                    nextNode = node;
                    while ( (node = nextNode.parent) != null && node.right == nextNode )
                        nextNode = node;
                }

                if ( node != null && interval.isUpstreamOf(node.interval) ) {
                    node = null;
                }
            }
            while ( node != null && interval.isDisjointFrom(node.interval) );

            return node;
        }

        static <V1> Node<V1> findByRank( final Node<V1> startingNode, final int initialRank ) {
            Node<V1> node = startingNode;
            int rank = initialRank;
            while ( node != null ) {
                final int nodeRank = node.getRank();
                if ( rank == nodeRank ) {
                    break;
                }

                if ( rank < nodeRank ) {
                    node = node.left;
                } else {
                    node = node.right;
                    rank -= nodeRank;
                }
            }

            return node;
        }

        static <V1> int getRank( final Node<V1> startingNode, final SVInterval interval ) {
            Node<V1> node = startingNode;
            int rank = 0;

            while ( node != null ) {
                final int cmpVal = interval.compareTo(node.getInterval());
                if ( cmpVal < 0 ) {
                    node = node.left;
                } else {
                    rank += node.getRank();
                    if ( cmpVal == 0 ) {
                        return rank; // EARLY RETURN!!!
                    }

                    node = node.right;
                }
            }

            return 0;
        }

        private int getRank() {
            int result = 1;
            if ( left != null ) {
                result = left.size + 1;
            }
            return result;
        }

        private Node<V1> spliceOut( final Node<V1> child, final Node<V1> initialRoot ) {
            Node<V1> root = initialRoot;
            if ( (child.parent = parent) == null ) {
                root = child;
                child.isBlack = true;
            } else {
                if ( parent.left == this ) {
                    parent.left = child;
                } else {
                    parent.right = child;
                }
                fixup(parent);

                if ( isBlack ) {
                    root = removeFixup(parent, child, root);
                }
            }

            return root;
        }

        private Node<V1> rotateLeft( final Node<V1> initialRoot ) {
            Node<V1> root = initialRoot;
            final Node<V1> child = right;

            final int childSize = child.size;
            child.size = size;
            size -= childSize;

            if ( (right = child.left) != null ) {
                right.parent = this;
                size += right.size;
            }

            if ( (child.parent = parent) == null ) {
                root = child;
            } else if ( this == parent.left ) {
                parent.left = child;
            } else {
                parent.right = child;
            }

            child.left = this;
            parent = child;

            setMaxEnd();
            child.setMaxEnd();

            return root;
        }

        private Node<V1> rotateRight( final Node<V1> initialRoot ) {
            Node<V1> root = initialRoot;
            final Node<V1> child = left;

            final int childSize = child.size;
            child.size = size;
            size -= childSize;

            if ( (left = child.right) != null ) {
                left.parent = this;
                size += left.size;
            }

            if ( (child.parent = parent) == null ) {
                root = child;
            } else if ( this == parent.left ) {
                parent.left = child;
            } else {
                parent.right = child;
            }

            child.right = this;
            parent = child;

            setMaxEnd();
            child.setMaxEnd();

            return root;
        }

        private void setMaxEnd() {
            maxEndInterval = interval;
            if ( left != null ) {
                maxEndInterval = laterEndingInterval(maxEndInterval, left.maxEndInterval);
            }
            if ( right != null ) {
                maxEndInterval = laterEndingInterval(maxEndInterval, right.maxEndInterval);
            }
        }

        private static SVInterval laterEndingInterval( final SVInterval interval1, final SVInterval interval2 ) {
            final int contig1 = interval1.getContig();
            final int contig2 = interval2.getContig();
            if ( contig1 > contig2 ) {
                return interval1;
            }
            if ( contig2 > contig1 ) {
                return interval2;
            }
            if ( interval2.getEnd() > interval1.getEnd() ) {
                return interval2;
            }
            return interval1;
        }

        private static <V1> void fixup( final Node<V1> initialNode ) {
            Node<V1> node = initialNode;
            do {
                node.size = 1;
                if ( node.left != null ) {
                    node.size += node.left.size;
                }
                if ( node.right != null ) {
                    node.size += node.right.size;
                }
                node.setMaxEnd();
            }
            while ( (node = node.parent) != null );
        }

        private static <V1> Node<V1> insertFixup( final Node<V1> initialDaughter, final Node<V1> initialRoot ) {
            Node<V1> daughter = initialDaughter;
            Node<V1> root = initialRoot;
            Node<V1> mom = daughter.parent;
            fixup(mom);

            while ( mom != null && !mom.isBlack ) {
                final Node<V1> gramma = mom.parent;
                Node<V1> auntie = gramma.left;
                if ( auntie == mom ) {
                    auntie = gramma.right;
                    if ( auntie != null && !auntie.isBlack ) {
                        mom.isBlack = true;
                        auntie.isBlack = true;
                        gramma.isBlack = false;
                        daughter = gramma;
                    } else {
                        if ( daughter == mom.right ) {
                            root = mom.rotateLeft(root);
                            mom = daughter;
                        }
                        mom.isBlack = true;
                        gramma.isBlack = false;
                        root = gramma.rotateRight(root);
                        break;
                    }
                } else {
                    if ( auntie != null && !auntie.isBlack ) {
                        mom.isBlack = true;
                        auntie.isBlack = true;
                        gramma.isBlack = false;
                        daughter = gramma;
                    } else {
                        if ( daughter == mom.left ) {
                            root = mom.rotateRight(root);
                            mom = daughter;
                        }
                        mom.isBlack = true;
                        gramma.isBlack = false;
                        root = gramma.rotateLeft(root);
                        break;
                    }
                }
                mom = daughter.parent;
            }
            root.isBlack = true;
            return root;
        }

        private static <V1> Node<V1> removeFixup( final Node<V1> initialParent,
                                                  final Node<V1> initialNode, final Node<V1> initialRoot ) {
            Node<V1> parent = initialParent;
            Node<V1> node = initialNode;
            Node<V1> root = initialRoot;
            do {
                if ( node == parent.left ) {
                    Node<V1> sister = parent.right;
                    if ( !sister.isBlack ) {
                        sister.isBlack = true;
                        parent.isBlack = false;
                        root = parent.rotateLeft(root);
                        sister = parent.right;
                    }
                    if ( (sister.left == null || sister.left.isBlack) && (sister.right == null || sister.right.isBlack) ) {
                        sister.isBlack = false;
                        node = parent;
                    } else {
                        if ( sister.right == null || sister.right.isBlack ) {
                            sister.left.isBlack = true;
                            sister.isBlack = false;
                            root = sister.rotateRight(root);
                            sister = parent.right;
                        }
                        sister.isBlack = parent.isBlack;
                        parent.isBlack = true;
                        sister.right.isBlack = true;
                        root = parent.rotateLeft(root);
                        node = root;
                    }
                } else {
                    Node<V1> sister = parent.left;
                    if ( !sister.isBlack ) {
                        sister.isBlack = true;
                        parent.isBlack = false;
                        root = parent.rotateRight(root);
                        sister = parent.left;
                    }
                    if ( (sister.left == null || sister.left.isBlack) && (sister.right == null || sister.right.isBlack) ) {
                        sister.isBlack = false;
                        node = parent;
                    } else {
                        if ( sister.left == null || sister.left.isBlack ) {
                            sister.right.isBlack = true;
                            sister.isBlack = false;
                            root = sister.rotateLeft(root);
                            sister = parent.left;
                        }
                        sister.isBlack = parent.isBlack;
                        parent.isBlack = true;
                        sister.left.isBlack = true;
                        root = parent.rotateRight(root);
                        node = root;
                    }
                }
                parent = node.parent;
            }
            while ( parent != null && node.isBlack );

            node.isBlack = true;
            return root;
        }
    }

    abstract class IteratorBase implements Iterator<Entry<V>> {
        protected Node<V> next;
        protected Node<V> last;

        protected IteratorBase( Node<V> node ) { next = node; }

        @Override
        public boolean hasNext() {
            return next != null;
        }

        @Override
        public void remove() {
            if ( last == null ) {
                throw new IllegalStateException("No entry to remove.");
            }

            removeNode(last);
            last = null;
        }

        /** equality of iterators is defined as having the same current position in the tree */
        @Override
        public boolean equals( final Object obj ) {
            return obj instanceof SVIntervalTree<?>.IteratorBase && ((SVIntervalTree<?>.IteratorBase)obj).next == next;
        }

        @Override
        public int hashCode() {
            return Objects.hashCode(next);
        }
    }

    public final class FwdIterator extends IteratorBase {

        public FwdIterator( final Node<V> node ) {
            super(node);
        }

        @Override
        public Node<V> next() {
            if ( next == null ) {
                throw new NoSuchElementException("No next element.");
            }

            if ( next.wasRemoved() ) {
                next = (Node<V>)min(next.getInterval());
                if ( next == null ) {
                    throw new ConcurrentModificationException("Current element was removed, and there are no more elements.");
                }
            }
            last = next;
            next = next.getNext();
            return last;
        }
    }

    public final class RevIterator extends IteratorBase {

        public RevIterator( final Node<V> node ) {
            super(node);
        }

        @Override
        public Node<V> next() {
            if ( next == null ) {
                throw new NoSuchElementException("No next element.");
            }

            if ( next.wasRemoved() ) {
                next = (Node<V>)max(next.getInterval());
                if ( next == null ) {
                    throw new ConcurrentModificationException("Current element was removed, and there are no more elements.");
                }
            }
            last = next;
            next = next.getPrev();
            return last;
        }
    }

    public final class OverlapIterator extends IteratorBase {
        private final SVInterval interval;

        public OverlapIterator( final SVInterval interval ) {
            super((Node<V>)minOverlapper(interval));
            this.interval = interval;
        }

        @Override
        public Node<V> next() {
            if ( next == null ) {
                throw new NoSuchElementException("No next element.");
            }

            if ( next.wasRemoved() ) {
                throw new ConcurrentModificationException("Current element was removed.");
            }

            last = next;
            next = Node.getNextOverlapper(next, interval);
            return last;
        }
    }

    public final static class ValuesIterator<V1> implements Iterator<V1> {
        private final Iterator<Entry<V1>> itr;

        public ValuesIterator( final Iterator<Entry<V1>> itr ) {
            this.itr = itr;
        }

        @Override
        public boolean hasNext() {
            return itr.hasNext();
        }

        @Override
        public V1 next() {
            return itr.next().getValue();
        }

        @Override
        public void remove() {
            itr.remove();
        }
    }
}
