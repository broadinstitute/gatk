package org.broadinstitute.hellbender.testutils;

import com.google.common.collect.Lists;
import htsjdk.samtools.*;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.CigarBuilder;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.*;
import java.util.stream.Collectors;

public final class ReadClipperTestUtils {
    //Should contain all the utils needed for tests to mass produce
    //reads, cigars, and other needed classes

    static final byte[] BASES = {'A', 'C', 'T', 'G'};
    static final byte[] QUALS = {2, 15, 25, 30};

    static final List<List<CigarElement>> LEADING_CLIPS = Arrays.asList(
            Collections.emptyList(),
            Collections.singletonList(new CigarElement(1, CigarOperator.HARD_CLIP)),
            Collections.singletonList(new CigarElement(1, CigarOperator.SOFT_CLIP)),
            Arrays.asList(new CigarElement(1, CigarOperator.HARD_CLIP), new CigarElement(1, CigarOperator.SOFT_CLIP)));

    static final List<List<CigarElement>> TRAILING_CLIPS = LEADING_CLIPS.stream().map(Lists::reverse).collect(Collectors.toList());

    static final CigarElement[] CORE_CIGAR_ELEMENTS = {
            new CigarElement(1, CigarOperator.INSERTION),
            new CigarElement(1, CigarOperator.DELETION),
            new CigarElement(1, CigarOperator.MATCH_OR_MISMATCH)};

    static final CigarElement[] CORE_CIGAR_ELEMENTS_INCLUDING_SKIPS = {
            new CigarElement(1, CigarOperator.INSERTION),
            new CigarElement(1, CigarOperator.DELETION),
            new CigarElement(1, CigarOperator.MATCH_OR_MISMATCH),
            new CigarElement(1, CigarOperator.SKIPPED_REGION)};

    /**
     * Make a read from the CIGAR string
     *
     * @param cigarString string used to create a CIGAR
     * @return artificial read
     */
    public static GATKRead makeReadFromCigar(String cigarString) {
        return makeReadFromCigar(TextCigarCodec.decode(cigarString));
    }

    /**
     * Make a read from the CIGAR.
     *
     * @param cigar
     * @return artificial read
     */
    public static GATKRead makeReadFromCigar(Cigar cigar) {
        return makeReadFromCigar(cigar, 0);
    }

    private static GATKRead makeReadFromCigar(Cigar cigar, int lengthChange) {
        int readLength = cigar.getReadLength();
        if (readLength >= -lengthChange) {
            readLength += lengthChange;
        }
        return ArtificialReadUtils.createArtificialRead(arrayFromArrayWithLength(BASES, readLength), arrayFromArrayWithLength(QUALS, readLength), cigar.toString());
    }

    private static byte [] arrayFromArrayWithLength(final byte[] array, final int length) {
        final byte [] output = new byte[length];
        for (int j = 0; j < length; j++) {
            output[j] = array[(j % array.length)];
        }
        return output;
    }

    /**
     * Make a read from the CIGAR
     *
     * @param cigarString  string used to create a CIGAR
     * @param lengthChange change in read length relative the CIGAR length
     * @return artificial read
     */
    public static GATKRead makeReadFromCigar(String cigarString, int lengthChange) {
        return makeReadFromCigar(TextCigarCodec.decode(cigarString), lengthChange);
    }

    /**
     * This function generates every valid permutation of cigar strings (with a given set of cigarElement) with a given length.
     * <p>
     * A valid cigar object obeys the following rules:
     * - No Hard/Soft clips in the middle of the read
     * - No deletions in the beginning / end of the read
     * - No repeated adjacent element (e.g. 1M2M, this should be 3M)
     * - No consecutive I/D elements
     *
     * @param maximumCigarElements the maximum number of elements in the cigar
     * @return a list with all valid Cigar objects
     */
    public static List<Cigar> generateCigarList(int maximumCigarElements, final boolean includeSkips) {
        final CigarElement[] coreElements = includeSkips ? CORE_CIGAR_ELEMENTS_INCLUDING_SKIPS : CORE_CIGAR_ELEMENTS;
        List<Cigar> cigarList = new ArrayList<>();

        for (final List<CigarElement> leadingClipElements : LEADING_CLIPS) {
            for (final List<CigarElement> trailingClipElements : TRAILING_CLIPS) {

                final int maximumCoreElements = maximumCigarElements - leadingClipElements.size() - trailingClipElements.size();
                if (maximumCoreElements < 1) {
                    continue;
                }

                int numCoreElements = coreElements.length;

                int[] cigarCombination = new int[maximumCoreElements];

                Arrays.fill(cigarCombination, 0);
                int currentIndex = 0;
                while (true) {
                    try {
                        final Cigar coreCigar = createCigarFromCombination(cigarCombination, coreElements);    // create the cigar
                        if (CigarUtils.isGood(coreCigar)) {// check if it's valid
                            final CigarBuilder builder = new CigarBuilder();
                            builder.addAll(leadingClipElements);
                            builder.addAll(coreCigar.getCigarElements());
                            builder.addAll(trailingClipElements);
                            cigarList.add(builder.make());                                      // add it
                        }
                    } catch (final Exception ex) {
                    } // we are creating random cigars, many of which are invalid.

                    boolean currentIndexChanged = false;
                    while (currentIndex < maximumCoreElements && cigarCombination[currentIndex] == numCoreElements - 1) {
                        currentIndex++;                                            // find the next index to increment
                        currentIndexChanged = true;                                // keep track of the fact that we have changed indices!
                    }

                    if (currentIndex == maximumCoreElements)                             // if we hit the end of the array, we're done.
                        break;

                    cigarCombination[currentIndex]++;                              // otherwise advance the current index

                    if (currentIndexChanged) {                                     // if we have changed index, then...
                        for (int i = 0; i < currentIndex; i++)
                            cigarCombination[i] = 0;                               // reset everything from 0->currentIndex
                        currentIndex = 0;                                          // go back to the first index
                    }
                }
            }
        }

        return cigarList;
    }

    private static Cigar createCigarFromCombination(int[] cigarCombination, CigarElement[] cigarElements) {
        return new Cigar(Arrays.stream(cigarCombination).mapToObj(n -> cigarElements[n]).collect(Collectors.toList()));
    }
}
