package org.broadinstitute.hellbender.utils.downsampling;

import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;

/**
 * Leveling Downsampler: Given a set of Lists of arbitrary items and a target size, removes items from
 * the Lists in an even fashion until the total size of all Lists is <= the target size. Leveling
 * does not occur until all Lists have been submitted and signalEndOfInput() is called.
 *
 * The Lists should be LinkedLists for maximum efficiency during item removal, however other
 * kinds of Lists are also accepted (albeit at a slight performance penalty).
 *
 * Since this downsampler extends the Downsampler interface rather than the ReadsDownsampler interface,
 * the Lists need not contain reads.
 *
 * @param <T> the List type representing the stacks to be leveled
 * @param <E> the type of the elements of each List
 *
 * @author David Roazen
 */
public final class LevelingDownsampler<T extends List<E>, E> extends Downsampler<T> {
    private final int minElementsPerStack;

    private final long targetSize;

    private List<T> groups;

    private boolean groupsAreFinalized;

    /**
     * Construct a LevelingDownsampler
     *
     * Uses the default minElementsPerStack of 1
     *
     * @param targetSize the sum of the sizes of all individual Lists this downsampler is fed may not exceed
     *                   this value -- if it does, items are removed from Lists evenly until the total size
     *                   is <= this value
     */
    public LevelingDownsampler(final long targetSize ) {
        this(targetSize, 1);
    }

    /**
     * Construct a LevelingDownsampler
     *
     * @param targetSize the sum of the sizes of all individual Lists this downsampler is fed may not exceed
     *                   this value -- if it does, items are removed from Lists evenly until the total size
     *                   is <= this value
     * @param minElementsPerStack no stack will be reduced below this size during downsampling.  That is,
     *                            if a stack has only 3 elements and minElementsPerStack is 3, no matter what
     *                            we'll not reduce this stack below 3.
     */
    public LevelingDownsampler(final long targetSize, final int minElementsPerStack ) {
        Utils.validateArg( targetSize >= 0, "targetSize must be >= 0 but got " + targetSize);
        Utils.validateArg( minElementsPerStack >= 0, "minElementsPerStack must be >= 0 but got " + minElementsPerStack);

        this.targetSize = targetSize;
        this.minElementsPerStack = minElementsPerStack;
        clearItems();
        resetStats();
    }

    @Override
    public void submit( final T item ) {
        Utils.nonNull(item, "item");
        groups.add(item);
    }

    @Override
    public void submit( final Collection<T> items ) {
        Utils.nonNull(items, "items");
        Utils.validateArg(!items.contains(null), "null item");
        groups.addAll(items);
    }

    @Override
    public boolean hasFinalizedItems() {
        return groupsAreFinalized && !groups.isEmpty();
    }

    @Override
    public List<T> consumeFinalizedItems() {
        if ( ! hasFinalizedItems() ) {
            return new ArrayList<>();
        }

        // pass by reference rather than make a copy, for speed
        final List<T> toReturn = groups;
        clearItems();
        return toReturn;
    }

    @Override
    public boolean hasPendingItems() {
        return ! groupsAreFinalized && !groups.isEmpty();
    }

    @Override
    public T peekFinalized() {
        return hasFinalizedItems() ? groups.get(0) : null;
    }

    @Override
    public T peekPending() {
        return hasPendingItems() ? groups.get(0) : null;
    }

    @Override
    public int size() {
        return groups.stream().mapToInt(g -> g.size()).sum();
    }

    @Override
    public void signalEndOfInput() {
        levelGroups();
        groupsAreFinalized = true;
    }

    @Override
    public void clearItems() {
        groups = new ArrayList<>();
        groupsAreFinalized = false;
    }

    private void levelGroups() {
        final int[] groupSizes = new int[groups.size()];
        int totalSize = 0;
        int currentGroupIndex = 0;

        for ( final T group : groups ) {
            groupSizes[currentGroupIndex] = group.size();
            totalSize += groupSizes[currentGroupIndex];
            currentGroupIndex++;
        }

        if ( totalSize <= targetSize ) {
            return;    // no need to eliminate any items
        }

        // We will try to remove exactly this many items, however we will refuse to allow any
        // one group to fall below size 1, and so might end up removing fewer items than this
        long numItemsToRemove = totalSize - targetSize;

        currentGroupIndex = 0;
        int numConsecutiveUmodifiableGroups = 0;

        // Continue until we've either removed all the items we wanted to, or we can't
        // remove any more items without violating the constraint that all groups must
        // be left with at least one item
        while ( numItemsToRemove > 0 && numConsecutiveUmodifiableGroups < groupSizes.length ) {
            if ( groupSizes[currentGroupIndex] > minElementsPerStack ) {
                groupSizes[currentGroupIndex]--;
                numItemsToRemove--;
                numConsecutiveUmodifiableGroups = 0;
            }
            else {
                numConsecutiveUmodifiableGroups++;
            }

            currentGroupIndex = (currentGroupIndex + 1) % groupSizes.length;
        }

        // Now we actually go through and reduce each group to its new count as specified in groupSizes
        currentGroupIndex = 0;
        for ( final T group : groups ) {
            downsampleOneGroup(group, groupSizes[currentGroupIndex]);
            currentGroupIndex++;
        }
    }

    private void downsampleOneGroup( final T group, final int numItemsToKeep ) {
        if ( numItemsToKeep >= group.size() ) {
            return;
        }

        final BitSet itemsToKeep = new BitSet(group.size());
        for ( final Integer selectedIndex : MathUtils.sampleIndicesWithoutReplacement(group.size(), numItemsToKeep) ) {
            itemsToKeep.set(selectedIndex);
        }

        int currentIndex = 0;

        // If our group is a linked list, we can remove the desired items in a single O(n) pass with an iterator
        if ( group instanceof LinkedList) {
            final Iterator<E> iter = group.iterator();
            while ( iter.hasNext() ) {
                final E item = iter.next();

                if ( ! itemsToKeep.get(currentIndex) ) {
                    iter.remove();
                    incrementNumberOfDiscardedItems(1);
                }

                currentIndex++;
            }
        }
        // If it's not a linked list, it's more efficient to copy the desired items into a new list and back rather
        // than suffer O(n^2) of item shifting
        else {
            final List<E> keptItems = new ArrayList<>(group.size());

            for ( final E item : group ) {
                if ( itemsToKeep.get(currentIndex) ) {
                    keptItems.add(item);
                }
                currentIndex++;
            }
            incrementNumberOfDiscardedItems(group.size() - keptItems.size());
            group.clear();
            group.addAll(keptItems);
        }
    }
}
