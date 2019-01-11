package org.broadinstitute.hellbender.tools.spark.sv.utils;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;

import java.util.Collections;
import java.util.Iterator;

@DefaultSerializer(PairedStrandedIntervalTree.Serializer.class)
public class PairedStrandedIntervalTree<V> implements Iterable<Tuple2<PairedStrandedIntervals, V>> {

    private SVIntervalTree<Tuple2<Boolean, SVIntervalTree<Tuple2<Boolean,V>>>> leftEnds = new SVIntervalTree<>();

    public PairedStrandedIntervalTree() {}

    @SuppressWarnings("unchecked")
    public PairedStrandedIntervalTree(final Kryo kryo, final Input input) {
        leftEnds = (SVIntervalTree<Tuple2<Boolean, SVIntervalTree<Tuple2<Boolean, V>>>>) kryo.readClassAndObject(input);
    }

    public boolean put(PairedStrandedIntervals pair, V value) {
        if (contains(pair)) return false;

        final SVIntervalTree.Entry<Tuple2<Boolean, SVIntervalTree<Tuple2<Boolean,V>>>> leftEntry = leftEnds.find(pair.getLeft().getInterval());
        if (leftEntry != null) {
            leftEntry.getValue()._2().put(pair.getRight().getInterval(), new Tuple2<>(pair.getRight().getStrand(), value));
        } else {
            final SVIntervalTree<Tuple2<Boolean, V>> rightEnds = new SVIntervalTree<>();
            rightEnds.put(pair.getRight().getInterval(), new Tuple2<>(pair.getRight().getStrand(), value));
            leftEnds.put(pair.getLeft().getInterval(), new Tuple2<>(pair.getLeft().getStrand(), rightEnds));
        }

        return true;
    }

    public int size() {
        return Utils.stream(leftEnds).mapToInt(e -> e.getValue()._2.size()).sum();
    }

    public final class PairedStrandedIntervalTreeOverlapperIterator implements Iterator<Tuple2<PairedStrandedIntervals, V>> {

        private Iterator<SVIntervalTree.Entry<Tuple2<Boolean, SVIntervalTree<Tuple2<Boolean, V>>>>> leftOverlappers;
        private Iterator<SVIntervalTree.Entry<Tuple2<Boolean, V>>> rightOverlappers;

        private SVIntervalTree.Entry<Tuple2<Boolean, SVIntervalTree<Tuple2<Boolean, V>>>> leftEntry;
        private SVIntervalTree.Entry<Tuple2<Boolean, V>> rightEntry;

        private final PairedStrandedIntervals query;

        public PairedStrandedIntervalTreeOverlapperIterator(PairedStrandedIntervalTree<V> tree, PairedStrandedIntervals query) {
            this.query = query;
            leftOverlappers = tree.leftEnds.overlappers(query.getLeft().getInterval());
        }

        private void advance(PairedStrandedIntervals query) {
            leftEntry = null;
            rightEntry = null;

            while (leftOverlappers.hasNext()) {
                leftEntry = leftOverlappers.next();
                if (! leftEntry.getValue()._1() == query.getLeft().getStrand()) {
                    continue;
                }
                rightOverlappers = leftEntry.getValue()._2().overlappers(query.getRight().getInterval());
                while (rightOverlappers.hasNext()) {
                    rightEntry = rightOverlappers.next();
                    if (! rightEntry.getValue()._1() == query.getRight().getStrand()) {
                        continue;
                    }
                    return;
                }
            }
            rightEntry = null;
        }

        @Override
        public boolean hasNext() {
            advance(query);
            return leftEntry != null && rightEntry != null;
        }

        @Override
        public Tuple2<PairedStrandedIntervals, V> next() {
            Tuple2<PairedStrandedIntervals, V> nextVal = new Tuple2<>(
                    new PairedStrandedIntervals(
                            new StrandedInterval(leftEntry.getInterval(),
                                    leftEntry.getValue()._1()),
                            new StrandedInterval(rightEntry.getInterval(),
                                    rightEntry.getValue()._1())),
                    rightEntry.getValue()._2);
            return nextVal;
        }

        @Override
        public void remove() {
            rightOverlappers.remove();
            if (leftEntry.getValue()._2().size() == 0) {
                leftOverlappers.remove();
            }
        }
    }

    public Iterator<Tuple2<PairedStrandedIntervals,V>> overlappers(PairedStrandedIntervals pair) {
        return new PairedStrandedIntervalTreeOverlapperIterator(this, pair);
    }

    public final class PairedStrandedIntervalTreeIterator implements Iterator<Tuple2<PairedStrandedIntervals, V>> {

        private final Iterator<SVIntervalTree.Entry<Tuple2<Boolean, SVIntervalTree<Tuple2<Boolean, V>>>>> leftEndIterator;
        private Iterator<SVIntervalTree.Entry<Tuple2<Boolean, V>>> rightEndIterator;
        private SVIntervalTree.Entry<Tuple2<Boolean, SVIntervalTree<Tuple2<Boolean, V>>>> leftEntry;

        PairedStrandedIntervalTreeIterator(PairedStrandedIntervalTree<V> tree) {
            leftEndIterator = tree.leftEnds.iterator();
            if (leftEndIterator.hasNext()) {
                leftEntry = leftEndIterator.next();
                rightEndIterator = leftEntry.getValue()._2().iterator();
            } else {
                leftEntry = null;
                rightEndIterator = Collections.emptyIterator();
            }
        }

        @Override
        public boolean hasNext() {
            return rightEndIterator.hasNext() || leftEndIterator.hasNext();
        }

        @Override
        public Tuple2<PairedStrandedIntervals, V> next() {
            SVIntervalTree.Entry<Tuple2<Boolean, V>> rightEntry;
            if (!rightEndIterator.hasNext()) {
                leftEntry = leftEndIterator.next();
                rightEndIterator = leftEntry.getValue()._2().iterator();
            }
            rightEntry = rightEndIterator.next();
            return new Tuple2<>(new PairedStrandedIntervals(
                    new StrandedInterval(leftEntry.getInterval(),
                            leftEntry.getValue()._1()),
                    new StrandedInterval(rightEntry.getInterval(),
                            rightEntry.getValue()._1())),
                    rightEntry.getValue()._2());
        }

        @Override
        public void remove() {
            rightEndIterator.remove();
            if (leftEntry.getValue()._2().size() == 0) {
                leftEndIterator.remove();
            }
        }

    }

    public Iterator<Tuple2<PairedStrandedIntervals, V>> iterator() {
        return new PairedStrandedIntervalTreeIterator(this);
    }

    public boolean contains(PairedStrandedIntervals pair) {
        final int leftEndIndex = leftEnds.getIndex(pair.getLeft().getInterval());
        if (leftEndIndex == -1) return false;
        final SVIntervalTree.Entry<Tuple2<Boolean, SVIntervalTree<Tuple2<Boolean, V>>>> leftEndEntry = leftEnds.findByIndex(leftEndIndex);
        final Tuple2<Boolean, SVIntervalTree<Tuple2<Boolean, V>>> storedValue = leftEndEntry.getValue();

        if (pair.getLeft().getStrand() != storedValue._1()) return false;

        final SVIntervalTree<Tuple2<Boolean, V>> rightEnds = storedValue._2();
        final int rightIndex = rightEnds.getIndex(pair.getRight().getInterval());
        return rightIndex != -1 && (pair.getRight().getStrand() == rightEnds.findByIndex(rightIndex).getValue()._1());
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

    private void serialize(final Kryo kryo, final Output output) {
        kryo.writeClassAndObject(output, this);
    }
}
