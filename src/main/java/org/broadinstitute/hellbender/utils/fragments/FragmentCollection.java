package org.broadinstitute.hellbender.utils.fragments;

import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.*;
import java.util.function.Function;

/**
 * Represents the results of the reads -> fragment calculation.
 *
 * Contains singleton -- objects whose underlying reads do not overlap their mate pair
 * Contains overlappingPairs -- objects whose underlying reads do overlap their mate pair
 */
public final class FragmentCollection<T> {

    private final Collection<T> singletons;
    private final Collection<List<T>> overlappingPairs;

    /**
     * Makes a new collection.
     * Note: this collection stores live pointers to the argument collections.
     * The callers must not modify those arguments after handing them off to this collection.
     *
     * The constructor is private - use the factory method if you need an object.
     */
    private FragmentCollection(final Collection<T> singletons, final Collection<List<T>> overlappingPairs) {
        this.singletons = singletons == null ? Collections.emptyList() : singletons;
        this.overlappingPairs = overlappingPairs == null ? Collections.emptyList() : overlappingPairs;
    }

    /**
     * Gets the T elements not containing overlapping elements, in no particular order.
     * The returned collection is unmodifiable.
     */
    public Collection<T> getSingletonReads() {
        return Collections.unmodifiableCollection(singletons);
    }

    /**
     * Gets the T elements containing overlapping elements, in no particular order
     * The returned collection is unmodifiable.
     */
    public Collection<List<T>> getOverlappingPairs() {
        return Collections.unmodifiableCollection(overlappingPairs);
    }

    /**
     * Generic algorithm that takes an iterable over T objects, a getter routine to extract the reads in T,
     * and returns a FragmentCollection that contains the T objects whose underlying reads either overlap (or
     * not) with their mate pairs.
     *
     * @param readContainingObjects An iterator of objects that contain SAMRecords
     * @param nElements the number of elements to be provided by the iterator, which is usually known upfront and
     *                  greatly improves the efficiency of the fragment calculation
     * @param getter a helper function that takes an object of type T and returns is associated SAMRecord
     * @param <T>
     * @return a fragment collection
     */
    private static <T> FragmentCollection<T> create(final Iterable<T> readContainingObjects, final int nElements, final Function<T, GATKRead> getter) {
        Collection<T> singletons = null;
        Collection<List<T>> overlapping = null;
        Map<String, T> nameMap = null;

        int lastStart = -1;

        // build an initial map, grabbing all of the multi-read fragments
        for ( final T p : readContainingObjects ) {
            final GATKRead read = getter.apply(p);

            if ( read.getStart() < lastStart ) {
                throw new IllegalArgumentException(String.format(
                        "FragmentUtils.create assumes that the incoming objects are ordered by " +
                                "SAMRecord alignment start, but saw a read %s with alignment start " +
                                "%d before the previous start %d", read.getName(), read.getStart(), lastStart));
            }
            lastStart = read.getStart();


            if ( ! read.isPaired() || read.mateIsUnmapped() || read.getMateStart() == 0 || read.getMateStart() > read.getEnd() ) {
                // if we know that this read won't overlap its mate, or doesn't have one, jump out early
                if ( singletons == null ) {
                    singletons = new ArrayList<>(nElements); // lazy init
                }
                singletons.add(p);
            } else {
                // the read might overlap it's mate, or is the rightmost read of a pair
                final String readName = read.getName();
                final T pe1 = nameMap == null ? null : nameMap.get(readName);
                if ( pe1 != null ) {
                    // assumes we have at most 2 reads per fragment
                    if ( overlapping == null ) {
                        overlapping = new ArrayList<>(); // lazy init
                    }
                    overlapping.add(Arrays.asList(pe1, p));
                    nameMap.remove(readName);
                } else {
                    if ( nameMap == null ) {
                        nameMap = new LinkedHashMap<>(nElements); // lazy init
                    }
                    nameMap.put(readName, p);
                }
            }
        }

        // add all of the reads that are potentially overlapping but whose mate never showed
        // up to the oneReadPile
        if ( nameMap != null && ! nameMap.isEmpty() ) {
            if ( singletons == null ) {
                singletons = nameMap.values();
            } else {
                singletons.addAll(nameMap.values());
            }
        }

        return new FragmentCollection<>(singletons, overlapping);
    }

    /**
     * Create a FragmentCollection containing PileupElements from the ReadBackedPileup rbp
     * @param rbp a non-null read-backed pileup.  The elements in this ReadBackedPileup must be ordered
     * @return a non-null FragmentCollection
     */
    public static FragmentCollection<PileupElement> create(final ReadPileup rbp) {
        if ( rbp == null ) {
            throw new IllegalArgumentException("Pileup cannot be null");
        }
        return create(rbp::sortedIterator, rbp.size(), pileup -> pileup.getRead());
    }

    /**
     * Create a FragmentCollection containing SAMRecords from a list of reads
     *
     * @param reads a non-null list of reads, ordered by their start location
     * @return a non-null FragmentCollection
     */
    public static FragmentCollection<GATKRead> create(final List<GATKRead> reads) {
        if ( reads == null ) {
            throw new IllegalArgumentException("Pileup cannot be null");
        }
        return create(reads, reads.size(), read -> read);
    }
}

