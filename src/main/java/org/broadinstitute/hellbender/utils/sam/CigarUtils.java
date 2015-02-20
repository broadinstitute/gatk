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

package org.broadinstitute.hellbender.utils.sam;

import htsjdk.samtools.*;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.util.ArrayList;
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
}
