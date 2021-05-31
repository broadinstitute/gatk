package org.broadinstitute.hellbender.tools.dragstr;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.function.Consumer;

/**
 * Spliterator that splits by segragating even and odd elements as opposed to
 * first and second half which the standard (e.g {@link Collection#stream} delgates onto {@link Spliterators.ArraySpliterator}).
*
* This approach makes the load more balanced between threads and the excecution order follows more closelly
* the input source list order.
*
* For its single use in this tool, this makes it more likely that there won't be a need to hold
* the full genome in memory by the underlaying ReferenceSource and that the
* actual progress corresponds more closely to the one reported by the progress-meter (always following
* sreference coordinate order).
*/
class InterleavingListSpliterator<T> implements Spliterator<T> {

    private final List<T> source;
    private int next;
    private int increment;
    private int remaining;

    InterleavingListSpliterator(final List<T> source) {
        this(source, 0, 1);
    }

    private InterleavingListSpliterator(final List<T> source, final int start, final int increment) {
        Utils.nonNull(this.source = source);
        final int sourceSize = source.size();
        Utils.validate(start >= 0 && start <= sourceSize, "start is out of range");
        this.next = start;
        Utils.validate((this.increment = increment) > 0, "increment must be extrictly positive");
        // we add increment -1 to numerator to get the ceiling of the int div.
        this.remaining = (sourceSize - start + increment - 1) / increment;
    }

    @Override
    public boolean tryAdvance(Consumer<? super T> action) {
        // per parent spec, a NPE is required in this case:
        Objects.requireNonNull(action);
        if (remaining > 0) {
            action.accept(source.get(next));
            next += increment;
            remaining--;
            return true;
        } else {
            return false;
        }
    }

    @Override
    public Spliterator<T> trySplit() {
        if (remaining >= 2) {
            final int splitStart = next + increment;
            increment <<= 1;
            // numerator + 1 to calculate the ceiling of the x/2 int div:
            remaining = (remaining + 1) >> 1;
            return new InterleavingListSpliterator<>(source, splitStart, increment);
        } else {
            return null;
        }
    }

    @Override
    public long estimateSize() {
        return remaining;
    }

    @Override
    public long getExactSizeIfKnown() {
        return remaining;
    }

    @Override
    public int characteristics() {
        return Spliterator.SUBSIZED | Spliterator.IMMUTABLE | Spliterator.NONNULL | Spliterator.SIZED;
    }

    @Override
    public void forEachRemaining(final Consumer<? super T> action) {
        Objects.requireNonNull(action);
        while (remaining > 0) {
            action.accept(source.get(next));
            next += increment;
            remaining--;
        }
    }
}
