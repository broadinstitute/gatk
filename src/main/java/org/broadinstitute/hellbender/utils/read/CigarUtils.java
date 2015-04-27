package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Stack;

public class CigarUtils {

    /**
     * Combines equal adjacent elements of a Cigar object
     *
     * @param rawCigar the cigar object
     * @return a combined cigar object
     */
    public static Cigar combineAdjacentCigarElements(Cigar rawCigar) {
        Cigar combinedCigar = new Cigar();
        CigarElement lastElement = null;
        int lastElementLength = 0;
        for (CigarElement cigarElement : rawCigar.getCigarElements()) {
            if (lastElement != null &&
                    ((lastElement.getOperator() == cigarElement.getOperator()) ||
                            (lastElement.getOperator() == CigarOperator.I && cigarElement.getOperator() == CigarOperator.D) ||
                            (lastElement.getOperator() == CigarOperator.D && cigarElement.getOperator() == CigarOperator.I)))
                lastElementLength += cigarElement.getLength();
            else
            {
                if (lastElement != null)
                    combinedCigar.add(new CigarElement(lastElementLength, lastElement.getOperator()));

                lastElement = cigarElement;
                lastElementLength = cigarElement.getLength();
            }
        }
        if (lastElement != null)
            combinedCigar.add(new CigarElement(lastElementLength, lastElement.getOperator()));

        return combinedCigar;
    }

    /**
     * Checks whether or not the read has any cigar element that is not H or S
     *
     * @param read the read
     * @return true if it has any M, I or D, false otherwise
     */
    public static boolean readHasNonClippedBases(SAMRecord read) {
        for (CigarElement cigarElement : read.getCigar().getCigarElements())
            if (cigarElement.getOperator() != CigarOperator.SOFT_CLIP && cigarElement.getOperator() != CigarOperator.HARD_CLIP)
                return true;
        return false;
    }

    public static Cigar invertCigar (Cigar cigar) {
        Stack<CigarElement> cigarStack = new Stack<>();
        for (CigarElement cigarElement : cigar.getCigarElements()) {
            cigarStack.push(cigarElement);
        }

        Cigar invertedCigar = new Cigar();
        while (!cigarStack.isEmpty()) {
            invertedCigar.add(cigarStack.pop());
        }

        return invertedCigar;
    }

    /**
    * A valid cigar object obeys the following rules:
    *  - No Hard/Soft clips in the middle of the read
    *  - No deletions in the beginning / end of the read
    *  - No repeated adjacent element (e.g. 1M2M -> this should be 3M)
    *  - No consecutive I/D elements
    **/
    public static boolean isCigarValid(Cigar cigar) {
        if (cigar.isValid(null, -1) == null) {                                                                          // This should take care of most invalid Cigar Strings (picard's "exhaustive" implementation)

            Stack<CigarElement> cigarElementStack = new Stack<>();                                          // Stack to invert cigar string to find ending operator
            CigarOperator startingOp = null;
            CigarOperator endingOp = null;

            // check if it doesn't start with deletions
            boolean readHasStarted = false;                                                                             // search the list of elements for the starting operator
            for (CigarElement cigarElement : cigar.getCigarElements()) {
                if (!readHasStarted) {
                    if (cigarElement.getOperator() != CigarOperator.SOFT_CLIP && cigarElement.getOperator() != CigarOperator.HARD_CLIP) {
                        readHasStarted = true;
                        startingOp = cigarElement.getOperator();
                    }
                }
                cigarElementStack.push(cigarElement);
            }

            while (!cigarElementStack.empty()) {
                CigarElement cigarElement = cigarElementStack.pop();
                if (cigarElement.getOperator() != CigarOperator.SOFT_CLIP && cigarElement.getOperator() != CigarOperator.HARD_CLIP) {
                    endingOp = cigarElement.getOperator();
                    break;
                }
            }

            if (startingOp != CigarOperator.DELETION && endingOp != CigarOperator.DELETION && startingOp != CigarOperator.SKIPPED_REGION && endingOp != CigarOperator.SKIPPED_REGION)
                return true;                                                                                          // we don't accept reads starting or ending in deletions (add any other constraint here)
        }

        return false;
    }

    public static int countRefBasesBasedOnCigar(final SAMRecord read, final int cigarStartIndex, final int cigarEndIndex){
        int result = 0;
        for(int i = cigarStartIndex; i<cigarEndIndex;i++){
            final CigarElement cigarElement = read.getCigar().getCigarElement(i);
            switch (cigarElement.getOperator()) {
                case M:
                case S:
                case D:
                case N:
                case H:
                    result += cigarElement.getLength();
                    break;
                case I:
                    break;
                default:
                    throw new GATKException("Unsupported cigar operator: " + cigarElement.getOperator());
            }
        }
        return result;
    }

    public static Cigar unclipCigar(Cigar cigar) {
        ArrayList<CigarElement> elements = new ArrayList<>(cigar.numCigarElements());
        for ( CigarElement ce : cigar.getCigarElements() ) {
            if ( !isClipOperator(ce.getOperator()) )
                elements.add(ce);
        }
        return new Cigar(elements);
    }

    private static boolean isClipOperator(CigarOperator op) {
        return op == CigarOperator.S || op == CigarOperator.H || op == CigarOperator.P;
    }

    public static Cigar reclipCigar(Cigar cigar, SAMRecord read) {
        ArrayList<CigarElement> elements = new ArrayList<>();

        int i = 0;
        int n = read.getCigar().numCigarElements();
        while ( i < n && isClipOperator(read.getCigar().getCigarElement(i).getOperator()) )
            elements.add(read.getCigar().getCigarElement(i++));

        elements.addAll(cigar.getCigarElements());

        i++;
        while ( i < n && !isClipOperator(read.getCigar().getCigarElement(i).getOperator()) )
            i++;

        while ( i < n && isClipOperator(read.getCigar().getCigarElement(i).getOperator()) )
            elements.add(read.getCigar().getCigarElement(i++));

        return new Cigar(elements);
    }

    public static boolean containsNOperator(final Cigar cigar) {
        for (final CigarElement ce : cigar.getCigarElements()) {
            if (ce.getOperator() == CigarOperator.N) {
                return true;
            }
        }
        return false;
    }

    /**
     * A good Cigar object obeys the following rules:
     *  - is valid as per SAM spec {@link Cigar#isValid(String, long)}.
     *  - has no consecutive I/D elements
     *  - does not start or end with deletions (with or without preceding clips).
     */
    public static boolean isGood(final Cigar c) {
        if (c == null){
            throw new IllegalArgumentException("null cigar");
        }
        //Note: the string comes from the SAMRecord so it must be a wellformed CIGAR (that is, in "\*|([0-9]+[MIDNSHPX=])+" as per SAM spec).
        //We don't have to check that
        if (c.isValid(null, -1) != null){  //if it's invalid, then it's not good
            return false;
        }
        final List<CigarElement> elems = c.getCigarElements();
        if (hasConsecutiveIndels(elems)){
            return false;
        }
        if (startsWithDeletionIgnoringClips(elems)){
            return false;
        }
        //revert the list and check deletions at the end
        final List<CigarElement> elemsRev = new ArrayList<>(elems);
        Collections.reverse(elemsRev);
        return !startsWithDeletionIgnoringClips(elemsRev);
    }

    /**
     * Checks if cigar has consecutive I/D elements.
     */
    private static boolean hasConsecutiveIndels(final List<CigarElement> elems) {
        boolean prevIndel = false;
        for (final CigarElement elem : elems) {
            final CigarOperator op = elem.getOperator();
            final boolean isIndel = (op == CigarOperator.INSERTION || op == CigarOperator.DELETION);
            if (prevIndel && isIndel) {
                return true;
            }
            prevIndel = isIndel;
        }
        return false;
    }

    /**
     * Checks if cigar starts with a deletion (ignoring any clips at the beginning).
     */
    private static boolean startsWithDeletionIgnoringClips(final List<CigarElement> elems) {
        final Iterator<CigarElement> iter = elems.iterator();
        boolean isClip = true;
        CigarOperator op = null;
        while(iter.hasNext() && isClip) { //consume clips at the beginning
            final CigarElement elem = iter.next();
            op = elem.getOperator();
            isClip = (op == CigarOperator.HARD_CLIP || op == CigarOperator.SOFT_CLIP);
        }
        //once all clips are consumed, is it a deletion or not?
        return op == CigarOperator.DELETION;
    }
}
