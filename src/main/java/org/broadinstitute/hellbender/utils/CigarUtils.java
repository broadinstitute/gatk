/*
* Copyright (c) 2012 The Broad Institute
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.hellbender.utils;

import htsjdk.samtools.*;
import org.broadinstitute.hellbender.utils.sam.AlignmentUtils;
import org.broadinstitute.hellbender.utils.smithwaterman.Parameters;
import org.broadinstitute.hellbender.utils.smithwaterman.SWPairwiseAlignment;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWaterman;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.util.Arrays;
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

    public static Cigar invertCigar (Cigar cigar) {
        Stack<CigarElement> cigarStack = new Stack<CigarElement>();
        for (CigarElement cigarElement : cigar.getCigarElements())
            cigarStack.push(cigarElement);

        Cigar invertedCigar = new Cigar();
        while (!cigarStack.isEmpty())
            invertedCigar.add(cigarStack.pop());

        return invertedCigar;
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

    public static Cigar cigarFromString(String cigarString) {
        return TextCigarCodec.getSingleton().decode(cigarString);
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

            Stack<CigarElement> cigarElementStack = new Stack<CigarElement>();                                          // Stack to invert cigar string to find ending operator
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

    public static final int countRefBasesBasedOnCigar(final SAMRecord read, final int cigarStartIndex, final int cigarEndIndex){
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

    // used in the bubble state machine to apply Smith-Waterman to the bubble sequence
    // these values were chosen via optimization against the NA12878 knowledge base
    public static final Parameters NEW_SW_PARAMETERS = new Parameters(200, -150, -260, -11);

    private final static String SW_PAD = "NNNNNNNNNN";

    /**
     * Calculate the cigar elements for this path against the reference sequence
     *
     * @param refSeq the reference sequence that all of the bases in this path should align to
     * @return a Cigar mapping this path to refSeq, or null if no reasonable alignment could be found
     */
    public static Cigar calculateCigar(final byte[] refSeq, final byte[] altSeq) {
        if ( altSeq.length == 0 ) {
            // horrible edge case from the unit tests, where this path has no bases
            return new Cigar(Arrays.asList(new CigarElement(refSeq.length, CigarOperator.D)));
        }

        final Cigar nonStandard;

        final String paddedRef = SW_PAD + new String(refSeq) + SW_PAD;
        final String paddedPath = SW_PAD + new String(altSeq) + SW_PAD;
        final SmithWaterman alignment = new SWPairwiseAlignment( paddedRef.getBytes(), paddedPath.getBytes(), NEW_SW_PARAMETERS);

        if ( isSWFailure(alignment) ) {
            return null;
        }


        // cut off the padding bases
        final int baseStart = SW_PAD.length();
        final int baseEnd = paddedPath.length() - SW_PAD.length() - 1; // -1 because it's inclusive
        nonStandard = AlignmentUtils.trimCigarByBases(alignment.getCigar(), baseStart, baseEnd);

        if ( nonStandard.getReferenceLength() != refSeq.length ) {
            nonStandard.add(new CigarElement(refSeq.length - nonStandard.getReferenceLength(), CigarOperator.D));
        }

        // finally, return the cigar with all indels left aligned
        return leftAlignCigarSequentially(nonStandard, refSeq, altSeq, 0, 0);
    }

    /**
     * Make sure that the SW didn't fail in some terrible way, and throw exception if it did
     */
    private static boolean isSWFailure(final SmithWaterman alignment) {
        // check that the alignment starts at the first base, which it should given the padding
        if ( alignment.getAlignmentStart2wrt1() > 0 ) {
            return true;
//          throw new IllegalStateException("SW failure ref " + paddedRef + " vs. " + paddedPath + " should always start at 0, but got " + alignment.getAlignmentStart2wrt1() + " with cigar " + alignment.getCigar());
        }

        // check that we aren't getting any S operators (which would be very bad downstream)
        for ( final CigarElement ce : alignment.getCigar().getCigarElements() ) {
            if ( ce.getOperator() == CigarOperator.S )
                return true;
            // soft clips at the end of the alignment are really insertions
//                throw new IllegalStateException("SW failure ref " + paddedRef + " vs. " + paddedPath + " should never contain S operators but got cigar " + alignment.getCigar());
        }

        return false;
    }

    /**
     * Left align the given cigar sequentially. This is needed because AlignmentUtils doesn't accept cigars with more than one indel in them.
     * This is a target of future work to incorporate and generalize into AlignmentUtils for use by others.
     * @param cigar     the cigar to left align
     * @param refSeq    the reference byte array
     * @param readSeq   the read byte array
     * @param refIndex  0-based alignment start position on ref
     * @param readIndex 0-based alignment start position on read
     * @return          the left-aligned cigar
     */
    public static Cigar leftAlignCigarSequentially(final Cigar cigar, final byte[] refSeq, final byte[] readSeq, int refIndex, int readIndex) {
        final Cigar cigarToReturn = new Cigar();
        Cigar cigarToAlign = new Cigar();
        for (int i = 0; i < cigar.numCigarElements(); i++) {
            final CigarElement ce = cigar.getCigarElement(i);
            if (ce.getOperator() == CigarOperator.D || ce.getOperator() == CigarOperator.I) {
                cigarToAlign.add(ce);
                final Cigar leftAligned = AlignmentUtils.leftAlignSingleIndel(cigarToAlign, refSeq, readSeq, refIndex, readIndex, false);
                for ( final CigarElement toAdd : leftAligned.getCigarElements() ) { cigarToReturn.add(toAdd); }
                refIndex += cigarToAlign.getReferenceLength();
                readIndex += cigarToAlign.getReadLength();
                cigarToAlign = new Cigar();
            } else {
                cigarToAlign.add(ce);
            }
        }
        if( !cigarToAlign.isEmpty() ) {
            for( final CigarElement toAdd : cigarToAlign.getCigarElements() ) {
                cigarToReturn.add(toAdd);
            }
        }

        final Cigar result = AlignmentUtils.consolidateCigar(cigarToReturn);
        if( result.getReferenceLength() != cigar.getReferenceLength() )
            throw new IllegalStateException("leftAlignCigarSequentially failed to produce a valid CIGAR.  Reference lengths differ.  Initial cigar " + cigar + " left aligned into " + result);
        return result;
    }
}
