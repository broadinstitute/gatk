package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Stack;

public final class CigarUtils {
    private CigarUtils(){}

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
     * Checks whether the cigar has any element that is not H or S
     * @return true the cigar has elements other than S or H, false otherwise.
     */
    public static boolean hasNonClippedBases(final Cigar cigar) {
        if (cigar == null){
            throw new IllegalArgumentException("null cigar");
        }
        return cigar.getCigarElements().stream().anyMatch(el -> el.getOperator() != CigarOperator.SOFT_CLIP && el.getOperator() != CigarOperator.HARD_CLIP);
    }

    /**
     * Inverts the order of the operators in the cigar.
     * Eg 10M1D20M -> 20M1D10M
     */
    public static Cigar invertCigar (final Cigar cigar) {
        if (cigar == null){
            throw new IllegalArgumentException("null cigar");
        }
        final List<CigarElement>  els = new ArrayList<>(cigar.getCigarElements());
        Collections.reverse(els);
        return new Cigar(els);
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

    /**
     * Compute the number of reference bases between the start (inclusive) and end (exclusive) cigar elements.
     * Note: The method does NOT use CigarOperator.consumesReferenceBases, since it checks something different.
     * The idea is you remove some elements from the beginning of the cigar string,
     * and want to recalculate what if the new starting reference position,
     * you want to count all the elements that indicate existing bases in the reference
     * (everything but I and P).
     * For example original position = 10. cigar: 2M3I2D1M. If you remove the 2M the new starting position is 12.
     * If you remove the 2M3I it is still 12. If you remove 2M3I2D (not reasonable cigar), you will get position 14.
     */
    @SuppressWarnings("fallthru")
    public static int countRefBasesBasedOnCigar(final GATKRead read, final int cigarStartIndex, final int cigarEndIndex){
        if (read == null){
            throw new IllegalArgumentException("null read");
        }
        final List<CigarElement> elems = read.getCigar().getCigarElements();
        if (cigarStartIndex < 0 || cigarEndIndex > elems.size() || cigarStartIndex > cigarEndIndex){
            throw new IllegalArgumentException("invalid index:" + 0 + " -" + elems.size());
        }
        int result = 0;
        for(int i = cigarStartIndex; i < cigarEndIndex; i++){
            final CigarElement cigarElement = elems.get(i);
            switch (cigarElement.getOperator()) {
                case M:
                case D:
                case N:
                case EQ:
                case X:
                case S:
                case H:
                    result += cigarElement.getLength();
                    break;
                case I:
                case P:        //for these two, nothing happens.
                    break;
                default:
                    throw new GATKException("Unsupported cigar operator: " + cigarElement.getOperator());
            }
        }
        return result;
    }

    /**
     * Removes all clipping operators from the cigar.
     */
    public static Cigar trimReadToUnclippedBases(final Cigar cigar) {
        if (cigar == null){
            throw new IllegalArgumentException("null cigar");
        }
        final List<CigarElement> elements = new ArrayList<>(cigar.numCigarElements());
        for ( CigarElement ce : cigar.getCigarElements() ) {
            if ( !isClipOperator(ce.getOperator()) )
                elements.add(ce);
        }
        return new Cigar(elements);
    }

    private static boolean isClipOperator(final CigarOperator op) {
        return op == CigarOperator.S || op == CigarOperator.H || op == CigarOperator.P;
    }

    private static boolean isClipOperator(final CigarElement el) {
        return isClipOperator(el.getOperator());
    }

    /**
     * Given a cigar1 and a read with cigar2,
     * this method creates cigar3 such that it has flanking clip operators from cigar2
     * and it has all operators from cigar1 in the middle.
     *
     * In other words if:
     * cigar2 = leftClip2 + noclips2 + rightClip2
     *
     * then
     * cigar3 = leftClip2 + cigar1 + rightClip2
     */
    public static Cigar reclipCigar(final Cigar cigar, final GATKRead read) {
        if (cigar == null || read == null){
            throw new IllegalArgumentException("null argument");
        }

        final List<CigarElement> elements = new ArrayList<>();

        int i = 0;
        final Cigar readCigar = read.getCigar();
        final int n = readCigar.numCigarElements();
        final List<CigarElement> readEls = readCigar.getCigarElements();

        //copy head clips
        while ( i < n && isClipOperator(readEls.get(i)) ) {
            elements.add(readEls.get(i));
            i++;
        }

        elements.addAll(cigar.getCigarElements());

        //skip over non-clips
        i++;
        while ( i < n && !isClipOperator(readEls.get(i)) ) {
            i++;
        }

        //copy tail clips
        while ( i < n && isClipOperator(readEls.get(i)) ) {
            elements.add(readEls.get(i));
            i++;
        }

        return new Cigar(elements);
    }

    /**
     * Returns whether the cigar has any N operators.
     */
    public static boolean containsNOperator(final Cigar cigar) {
        if (cigar == null){
            throw new IllegalArgumentException("null cigar");
        }
        return cigar.getCigarElements().stream().anyMatch(el -> el.getOperator() == CigarOperator.N);
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
