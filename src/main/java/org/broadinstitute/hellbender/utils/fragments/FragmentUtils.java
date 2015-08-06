package org.broadinstitute.hellbender.utils.fragments;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.*;
import java.util.function.Function;
import java.util.function.Supplier;

/**
 * An easy to access fragment-based pileup, which contains two separate pileups.  The first
 * is a regular collection of PileupElements containing all of the reads in the original RBP
 * that uniquely info about a fragment.  The second are TwoReadPileupElements that, as the
 * name suggests, contain two reads that are sequenced from the same underlying fragment.
 *
 * Note: the order of the oneReadPileup and twoReadPileups are not
 * defined.  The algorithms that produce these lists are in fact producing
 * lists of Pileup elements *NOT* sorted by alignment start position of the underlying
 * reads.
 */
public final class FragmentUtils {

    public static final double DEFAULT_PCR_ERROR_RATE = 1e-4;
    private FragmentUtils() {} // private constructor

    public static FragmentCollection<PileupElement> create(final ReadPileup rbp) {
        if ( rbp == null ) throw new IllegalArgumentException("Pileup cannot be null");
        return create(rbp, rbp.size(), pe -> pe.getRead());
    }

    /**
     * Generic algorithm that takes an iterable over T objects, a getter routine to extract the reads in T,
     * and returns a FragmentCollection that contains the T objects whose underlying reads either overlap (or
     * not) with their mate pairs.
     *
     * @param readContainingObjects An iterator of objects that contain GATKSAMRecords
     * @param nElements the number of elements to be provided by the iterator, which is usually known upfront and
     *                  greatly improves the efficiency of the fragment calculation
     * @param getter a helper function that takes an object of type T and returns is associated GATKSAMRecord
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
                                "%d before the previous start %d", read.toString(), read.getStart(), lastStart));
            }
            lastStart = read.getStart();

            final int mateStart = read.getMateStart();
            if ( mateStart == 0 || mateStart > read.getEnd()) {
                // if we know that this read won't overlap its mate, or doesn't have one, jump out early
                if ( singletons == null ) singletons = new ArrayList<>(nElements); // lazy init
                singletons.add(p);
            } else {
                // the read might overlap it's mate, or is the rightmost read of a pair
                final String readName = read.getName();
                final T pe1 = nameMap == null ? null : nameMap.get(readName);
                if ( pe1 != null ) {
                    // assumes we have at most 2 reads per fragment
                    if ( overlapping == null ) overlapping = new ArrayList<>(); // lazy init
                    overlapping.add(Arrays.asList(pe1, p));
                    nameMap.remove(readName);
                } else {
                    if ( nameMap == null ) nameMap = new HashMap<>(nElements); // lazy init
                    nameMap.put(readName, p);
                }
            }
        }

        // add all of the reads that are potentially overlapping but whose mate never showed
        // up to the oneReadPile
        if ( nameMap != null && ! nameMap.isEmpty() ) {
            if ( singletons == null )
                singletons = nameMap.values();
            else
                singletons.addAll(nameMap.values());
        }

        return new FragmentCollection<>(singletons, overlapping);
    }


}
