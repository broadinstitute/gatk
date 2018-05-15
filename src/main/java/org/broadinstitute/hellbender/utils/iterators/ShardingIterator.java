package org.broadinstitute.hellbender.utils.iterators;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.PeekableIterator;
import org.broadinstitute.hellbender.engine.Shard;
import org.broadinstitute.hellbender.engine.ShardBoundary;
import org.broadinstitute.hellbender.engine.ShardBoundaryShard;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Queue;

/**
 * Iterator over shards/windows (overlapping or not) of sorted {@link Locatable}.
 *
 * <p>Returned {@link Shard<T>} would be filled with {@link Locatable} overlapping the provided
 * {@link ShardBoundary}. If no {@link Locatable} overlaps, an empty {@link Shard<T>} would be returned.
 *
 * <p>Warning: for overlapping shards, the objects stored in the {@link Shard<T>} returned by
 * {@link #next()} are not clones, and thus shards are not independent. If the objects are modified
 * in each shard, the caller should clone to safely operate on the objects.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public final class ShardingIterator<T extends Locatable> implements Iterator<Shard<T>> {

    // underlying iterator - should be sorted by coordinate
    private final PeekableIterator<T> it;
    // queuing the shards - should be sorted by coordinate
    private final Queue<ShardBoundary> shards;
    // this wil be filled with the data for the shard in the first position of the shard queue
    private final Queue<T> nextData;
    // required to check sort order
    private final SAMSequenceDictionary dictionary;

    /**
     * Wraps an iterator to perform sharding.
     *
     * @param iterator   sorted iterator to iterate over windows.
     * @param boundaries sorted list of shard boundaries to iterate over.
     * @param dictionary dictionary for position comparison.
     */
    public ShardingIterator(final Iterator<T> iterator,
            final List<ShardBoundary> boundaries,
            final SAMSequenceDictionary dictionary) {
        // check arguments
        Utils.nonNull(iterator, "null iterator");
        Utils.nonEmpty(boundaries, "Empty shard boundaries");
        Utils.nonNull(dictionary, "null dictionary");

        // store the iterator and allow peek()
        this.it = new PeekableIterator<>(iterator);
        // convert the boundaries into a queue
        this.shards = new ArrayDeque<>(boundaries);
        // store the dictionary
        this.dictionary = dictionary;
        // create a new queue for the next data
        this.nextData = new ArrayDeque<>();

        // start the iteration
        advance();
    }

    @Override
    public boolean hasNext() {
        return !shards.isEmpty();
    }

    /**
     * {@inheritDoc}
     *
     * <p>Note: for overlapping shards, the objects stored in two (or more) consecutive shards are
     * references and not clones. If modified in the previous shard, the state would be kept in the
     * next returned value. It is up to the caller to make copies if necessary to avoid state in the
     * next iteration.
     */
    @Override
    public Shard<T> next() {
        if (!hasNext()) {
            throw new NoSuchElementException();
        }
        // copying the deque to being able to remove unnecesary data and fill in with the next
        final Shard<T> next = new ShardBoundaryShard<T>(shards.poll(),
                nextData.isEmpty() ? Collections.emptyList() : new ArrayList<>(nextData));
        advance();
        return next;
    }

    // advance to the next data
    private void advance() {
        // if we already finished the iteration, we clean up the data
        if (shards.isEmpty()) {
            nextData.clear();
        } else {
            // first we empty the queue if it has non-overlapping data with the next shard
            removeNonOverlapping();
            // then, fill in the data
            fillQueueWithOverlapping();
        }
    }

    // assumes that shards are not empty
    private void removeNonOverlapping() {
        final SimpleInterval nextInterval = shards.peek().getPaddedInterval();
        // first remove all the data that does not overlap
        // (e.g., different contig or before the start of the next interval)
        while (!nextData.isEmpty() && !nextData.peek().overlaps(nextInterval)) {
            nextData.poll();
        }


        // we should also check all the data coming from the iterator that is before the start
        // (and in the same contig) to ensure that it overlaps
        // we keep that data in the overlapping list, discarding the rest
        // this follows a similar approach as the FeatureCache.trimToNewStartPosition
        final List<T> overlapping = new ArrayList<>(128);
        while (!nextData.isEmpty() && nextData.peek().contigsMatch(nextInterval) && nextData.peek().getStart() < nextInterval.getStart()) {
            T beforeStart = nextData.poll();
            if (beforeStart.overlaps(nextInterval)) {
                overlapping.add(beforeStart);
            }
        }

        // now we add back the data which overlaps, in the same order
        for (int i = overlapping.size() - 1; i >= 0; --i) {
            nextData.add(overlapping.get(i));
        }
    }

    // assumes that shards are not empty
    private void fillQueueWithOverlapping() {
        while (it.hasNext()) {
            // if this one overlaps
            if (it.peek().overlaps(shards.peek().getPaddedInterval())) {
                nextData.offer(it.next());
            } else if (IntervalUtils
                    .isBefore(it.peek(), shards.peek().getPaddedInterval(), dictionary)) {
                // just ignore this, because is before the current interval
                it.next();
            } else {
                // break the while loop because it is after this window
                break;
            }
        }
    }
}
