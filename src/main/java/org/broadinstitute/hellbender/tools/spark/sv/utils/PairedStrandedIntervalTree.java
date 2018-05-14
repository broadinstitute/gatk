package org.broadinstitute.hellbender.tools.spark.sv.utils;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import org.broadinstitute.hellbender.exceptions.GATKException;
import scala.Tuple2;

import java.util.Collections;
import java.util.Iterator;
import java.util.NoSuchElementException;

@DefaultSerializer(PairedStrandedIntervalTree.Serializer.class)
public class PairedStrandedIntervalTree<V> implements Iterable<Tuple2<PairedStrandedIntervals, V>> {

    private SVIntervalTree<LeftEndEntry<V>> leftEnds = new SVIntervalTree<>();
    private int size = 0;

    public PairedStrandedIntervalTree() {}

    @SuppressWarnings("unchecked")
    public PairedStrandedIntervalTree(final Kryo kryo, final Input input) {
        leftEnds = (SVIntervalTree<LeftEndEntry<V>>) kryo.readClassAndObject(input);
    }

    public V put(PairedStrandedIntervals pair, V value) {
        final SVIntervalTree.Entry<LeftEndEntry<V>> leftEndEntry = leftEnds.find(pair.getLeft().getInterval());
        final LeftEndEntry<V> leftEntryValue;
        if (leftEndEntry == null) {
            leftEntryValue = new LeftEndEntry<>();
            leftEnds.put(pair.getLeft().getInterval(), leftEntryValue);
        } else {
            leftEntryValue = leftEndEntry.getValue();
        }

        RightEndEntry<V> rightValue;
        final SVIntervalTree<RightEndEntry<V>> rightEntryTree =
                pair.getLeft().getStrand() ? leftEntryValue.getPosStrandEntries() : leftEntryValue.getNegStrandEntries();
        final SVIntervalTree.Entry<RightEndEntry<V>> rightEndEntry = rightEntryTree.find(pair.getRight().getInterval());

        if (rightEndEntry == null) {
            rightValue = new RightEndEntry<>();
            rightEntryTree.put(pair.getRight().getInterval(), rightValue);
        } else {
            rightValue = rightEndEntry.getValue();
        }

        size = size + 1;
        if (pair.getRight().getStrand()) {
            final V existingValue = rightValue.getPosStrandValue();
            rightValue.setPosStrandValue(value);
            return existingValue;
        } else {
            final V existingValue = rightValue.getNegStrandValue();
            rightValue.setNegStrandValue(value);
            return existingValue;
        }

    }

    public int size() {
        return size;
    }

    public final class PairedStrandedIntervalTreeOverlapperIterator implements Iterator<Tuple2<PairedStrandedIntervals, V>> {

        private Iterator<SVIntervalTree.Entry<LeftEndEntry<V>>> leftOverlappers;
        private Iterator<SVIntervalTree.Entry<RightEndEntry<V>>> rightOverlappers;

        private SVIntervalTree.Entry<LeftEndEntry<V>> leftEntry;
        private SVIntervalTree.Entry<RightEndEntry<V>> rightEntry;

        private final PairedStrandedIntervals query;

        boolean emittedThisValue = false;
        boolean canRemove = false;

        public PairedStrandedIntervalTreeOverlapperIterator(final PairedStrandedIntervalTree<V> tree,
                                                            final PairedStrandedIntervals query) {
            this.query = query;
            leftOverlappers = tree.leftEnds.overlappers(query.getLeft().getInterval());
            if (leftOverlappers.hasNext()) {
                leftEntry = leftOverlappers.next();

                final SVIntervalTree<RightEndEntry<V>> rightEntries =
                        query.getLeft().getStrand() ? leftEntry.getValue().getPosStrandEntries() : leftEntry.getValue().getNegStrandEntries();
                rightOverlappers = rightEntries.overlappers(query.getRight().getInterval());
            } else {
                rightOverlappers = Collections.emptyIterator();
            }

            advance(query);
        }

        private void advance(PairedStrandedIntervals query) {
            emittedThisValue = false;
            while (rightOverlappers.hasNext()) {
                rightEntry = this.rightOverlappers.next();
                final V value = query.getRight().getStrand() ? rightEntry.getValue().getPosStrandValue() : rightEntry.getValue().getNegStrandValue();
                if (value != null) {
                    return;
                }
            }

            while (leftOverlappers.hasNext()) {
                leftEntry = leftOverlappers.next();
                final SVIntervalTree<RightEndEntry<V>> rightEntries =
                        query.getLeft().getStrand() ? leftEntry.getValue().getPosStrandEntries() : leftEntry.getValue().getNegStrandEntries();

                rightOverlappers = rightEntries.overlappers(query.getRight().getInterval());

                while (rightOverlappers.hasNext()) {
                    rightEntry = this.rightOverlappers.next();
                    final V value = query.getRight().getStrand() ? rightEntry.getValue().getPosStrandValue() : rightEntry.getValue().getNegStrandValue();
                    if (value != null) {
                        return;
                    }
                }
            }
            rightEntry = null;
        }

        @Override
        public boolean hasNext() {
            canRemove = false;
            if (emittedThisValue) {
                advance(query);
            }
            return rightEntry != null;
        }

        @Override
        public Tuple2<PairedStrandedIntervals, V> next() {
            if (rightEntry == null) {
                throw new NoSuchElementException();
            }
            Tuple2<PairedStrandedIntervals, V> nextVal = new Tuple2<>(
                    new PairedStrandedIntervals(
                            new StrandedInterval(leftEntry.getInterval(),
                                    query.getLeft().getStrand()),
                            new StrandedInterval(rightEntry.getInterval(),
                                    query.getRight().getStrand())),
                    query.getRight().getStrand() ? rightEntry.getValue().getPosStrandValue() : rightEntry.getValue().getNegStrandValue());
            emittedThisValue = true;
            canRemove = true;
            return nextVal;
        }

        @Override
        public void remove() {
            if (! canRemove) {
                throw new UnsupportedOperationException("Remove can only be called on this iterator immediately after a call to next. Calling remove after hasNext is unsupported.");
            }
            if (query.getRight().getStrand()) {
                rightEntry.getValue().setPosStrandValue(null);
            } else {
                rightEntry.getValue().setNegStrandValue(null);
            }
            if (rightEntry.getValue().isEmpty()) {
                rightOverlappers.remove();
            }
            if (leftEntry.getValue().getPosStrandEntries().size() == 0 &&
                    leftEntry.getValue().getNegStrandEntries().size() == 0) {
                leftOverlappers.remove();
            }
            canRemove = false;
            size = size - 1;
        }
    }

    public Iterator<Tuple2<PairedStrandedIntervals,V>> overlappers(PairedStrandedIntervals pair) {
        return new PairedStrandedIntervalTreeOverlapperIterator(this, pair);
    }

    public final class PairedStrandedIntervalTreeIterator implements Iterator<Tuple2<PairedStrandedIntervals, V>> {

        private final Iterator<SVIntervalTree.Entry<LeftEndEntry<V>>> leftEndIterator;
        private Iterator<SVIntervalTree.Entry<RightEndEntry<V>>> rightEndIterator;
        private SVIntervalTree.Entry<RightEndEntry<V>> rightEndEntry;
        private SVIntervalTree.Entry<LeftEndEntry<V>> leftEntry;
        private boolean posLeftStrand;
        private boolean posRightStrand;

        PairedStrandedIntervalTreeIterator(PairedStrandedIntervalTree<V> tree) {
            leftEndIterator = tree.leftEnds.iterator();
            if (leftEndIterator.hasNext()) {
                leftEntry = leftEndIterator.next();
                posLeftStrand = true;
                rightEndIterator = leftEntry.getValue().getPosStrandEntries().iterator();
            } else {
                leftEntry = null;
                rightEndIterator = Collections.emptyIterator();
            }
        }

        @Override
        public boolean hasNext() {
            if (rightEndEntry != null && posRightStrand && rightEndEntry.getValue().getNegStrandValue() != null) return true;
            final boolean rightIteratorHasNext = rightEndIterator.hasNext();
            if (rightIteratorHasNext) return true;
            if (posLeftStrand && leftEntry.getValue().getNegStrandEntries().iterator().hasNext()) return true;
            return leftEndIterator.hasNext();
        }

        private void advance() {
            // check if there's a neg strand value on the right entry
            if (rightEndEntry != null && posRightStrand && rightEndEntry.getValue().getNegStrandValue() != null) {
                posRightStrand = false;
                return;
            }

            // advance to the next right entry
            if (rightEndIterator.hasNext()) {
                rightEndEntry = rightEndIterator.next();
                posRightStrand = rightEndEntry.getValue().getPosStrandValue() != null;
                return;
            }

            // out of stuff from right end iterator
            if (posLeftStrand) {
                // check the neg strand entries on the left entry
                rightEndIterator = leftEntry.getValue().getNegStrandEntries().iterator();
                posLeftStrand = false;
                if (rightEndIterator.hasNext()) {
                    rightEndEntry = rightEndIterator.next();
                    posRightStrand = rightEndEntry.getValue().getPosStrandValue() != null;
                    return;
                }
            }

            // otherwise move to the next left entry
            if (leftEndIterator.hasNext()) {
                leftEntry = leftEndIterator.next();
                rightEndIterator = leftEntry.getValue().getPosStrandEntries().iterator();
                posLeftStrand = true;
                if (! rightEndIterator.hasNext()) {
                    rightEndIterator = leftEntry.getValue().getNegStrandEntries().iterator();
                    posLeftStrand = false;
                }
                if (rightEndIterator.hasNext()) {
                    rightEndEntry = rightEndIterator.next();
                    posRightStrand = rightEndEntry.getValue().getPosStrandValue() != null;
                    return;
                }
            }

            throw new GATKException.ShouldNeverReachHereException("Couldn't find a next element to advance too");
        }

        @Override
        public Tuple2<PairedStrandedIntervals, V> next() {
            advance();
            if (rightEndEntry == null) {
                throw new NoSuchElementException();
            }
            final V value = posRightStrand ? rightEndEntry.getValue().getPosStrandValue() : rightEndEntry.getValue().getNegStrandValue();
            final Tuple2<PairedStrandedIntervals, V> next = new Tuple2<>(new PairedStrandedIntervals(
                    new StrandedInterval(leftEntry.getInterval(),
                            posLeftStrand),
                    new StrandedInterval(rightEndEntry.getInterval(),
                            posRightStrand)),
                    value);
            return next;
        }

        @Override
        public void remove() {
            if (posRightStrand) {
                rightEndEntry.getValue().setPosStrandValue(null);
            } else {
                rightEndEntry.getValue().setNegStrandValue(null);
            }
            if (rightEndEntry.getValue().isEmpty()) {
                rightEndIterator.remove();
            }
            if (leftEntry.getValue().getPosStrandEntries().size() == 0 &&
                    leftEntry.getValue().getNegStrandEntries().size() == 0) {
                leftEndIterator.remove();
            }
            size = size - 1;
        }

    }

    public Iterator<Tuple2<PairedStrandedIntervals, V>> iterator() {
        return new PairedStrandedIntervalTreeIterator(this);
    }

    public boolean contains(PairedStrandedIntervals pair) {
        final int leftEndIndex = leftEnds.getIndex(pair.getLeft().getInterval());
        if (leftEndIndex == -1) return false;
        final SVIntervalTree.Entry<LeftEndEntry<V>> leftEndEntry = leftEnds.findByIndex(leftEndIndex);
        final LeftEndEntry<V> leftEndEntryValue = leftEndEntry.getValue();

        final SVIntervalTree<RightEndEntry<V>> rightEntries = pair.getLeft().getStrand() ? leftEndEntryValue.getPosStrandEntries() : leftEndEntryValue.getNegStrandEntries();
        final int rightIndex = rightEntries.getIndex(pair.getRight().getInterval());
        final SVIntervalTree.Entry<RightEndEntry<V>> rightEndEntry = rightEntries.findByIndex(rightIndex);
        return rightIndex != -1 && (pair.getRight().getStrand() ? rightEndEntry.getValue().getPosStrandValue() != null : rightEndEntry.getValue().getNegStrandValue() != null);
    }

    public static final class Serializer<T> extends com.esotericsoftware.kryo.Serializer<PairedStrandedIntervalTree<T>> {
        @Override
        public void write(final Kryo kryo, final Output output, final PairedStrandedIntervalTree<T> tree ) {
            tree.serialize(kryo, output);
        }

        @Override
        public PairedStrandedIntervalTree<T> read(final Kryo kryo, final Input input, final Class<PairedStrandedIntervalTree<T>> klass ) {
            return new PairedStrandedIntervalTree<>(kryo, input);
        }
    }

    static final class LeftEndEntry<V> {
        final SVIntervalTree<RightEndEntry<V>> posStrandEntries = new SVIntervalTree<>();
        final SVIntervalTree<RightEndEntry<V>> negStrandEntries = new SVIntervalTree<>();

        public SVIntervalTree<RightEndEntry<V>> getPosStrandEntries() {
            return posStrandEntries;
        }

        public SVIntervalTree<RightEndEntry<V>> getNegStrandEntries() {
            return negStrandEntries;
        }
    }

    static final class RightEndEntry<V> {
        V posStrandValue;
        V negStrandValue;

        public V getPosStrandValue() {
            return posStrandValue;
        }

        public void setPosStrandValue(final V posStrandValue) {
            this.posStrandValue = posStrandValue;
        }

        public V getNegStrandValue() {
            return negStrandValue;
        }

        public void setNegStrandValue(final V negStrandValue) {
            this.negStrandValue = negStrandValue;
        }

        public boolean isEmpty() {
            return posStrandValue == null && negStrandValue == null;
        }
    }

    private void serialize(final Kryo kryo, final Output output) {
        kryo.writeClassAndObject(output, this);
    }
}
