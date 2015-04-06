package org.broadinstitute.hellbender.utils;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.utils.read.ArtificialSAMUtils;
import org.broadinstitute.hellbender.utils.read.CigarUtils;

import java.util.LinkedList;
import java.util.List;

public class ReadClipperTestUtils {
    //Should contain all the utils needed for tests to mass produce
    //reads, cigars, and other needed classes

    final static byte[] BASES = {'A', 'C', 'T', 'G'};
    final static byte[] QUALS = {2, 15, 25, 30};
    final static String CIGAR = "4M";
    final static CigarElement[] cigarElements = {new CigarElement(1, CigarOperator.HARD_CLIP),
            new CigarElement(1, CigarOperator.SOFT_CLIP),
            new CigarElement(1, CigarOperator.INSERTION),
            new CigarElement(1, CigarOperator.DELETION),
            new CigarElement(1, CigarOperator.MATCH_OR_MISMATCH)};


    public static SAMRecord makeReadFromCigar(Cigar cigar) {
        return ArtificialSAMUtils.createArtificialRead(Utils.arrayFromArrayWithLength(BASES, cigar.getReadLength()), Utils.arrayFromArrayWithLength(QUALS, cigar.getReadLength()), cigar.toString());
    }

    /**
     * This function generates every valid permutation of cigar strings (with a given set of cigarElement) with a given length.
     * <p>
     * A valid cigar object obeys the following rules:
     * - No Hard/Soft clips in the middle of the read
     * - No deletions in the beginning / end of the read
     * - No repeated adjacent element (e.g. 1M2M -> this should be 3M)
     * - No consecutive I/D elements
     *
     * @param maximumLength the maximum number of elements in the cigar
     * @return a list with all valid Cigar objects
     */
    public static List<Cigar> generateCigarList(int maximumLength, CigarElement[] cigarElements) {
        int numCigarElements = cigarElements.length;
        LinkedList<Cigar> cigarList = new LinkedList<Cigar>();
        byte[] cigarCombination = new byte[maximumLength];

        Utils.fillArrayWithByte(cigarCombination, (byte) 0);               // we start off with all 0's in the combination array.
        int currentIndex = 0;
        while (true) {
            Cigar cigar = createCigarFromCombination(cigarCombination, cigarElements);    // create the cigar
            cigar = CigarUtils.combineAdjacentCigarElements(cigar);                   // combine adjacent elements
            if (CigarUtils.isCigarValid(cigar)) {                                     // check if it's valid
                cigarList.add(cigar);                                      // add it
            }

            boolean currentIndexChanged = false;
            while (currentIndex < maximumLength && cigarCombination[currentIndex] == numCigarElements - 1) {
                currentIndex++;                                            // find the next index to increment
                currentIndexChanged = true;                                // keep track of the fact that we have changed indices!
            }

            if (currentIndex == maximumLength)                             // if we hit the end of the array, we're done.
                break;

            cigarCombination[currentIndex]++;                              // otherwise advance the current index

            if (currentIndexChanged) {                                     // if we have changed index, then...
                for (int i = 0; i < currentIndex; i++)
                    cigarCombination[i] = 0;                               // reset everything from 0->currentIndex
                currentIndex = 0;                                          // go back to the first index
            }
        }

        return cigarList;
    }

    private static Cigar createCigarFromCombination(byte[] cigarCombination, CigarElement[] cigarElements) {
        Cigar cigar = new Cigar();
        for (byte i : cigarCombination) {
            cigar.add(cigarElements[i]);
        }
        return cigar;
    }

}
