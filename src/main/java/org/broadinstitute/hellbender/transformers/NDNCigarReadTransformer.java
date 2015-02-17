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

package org.broadinstitute.hellbender.transformers;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.tools.walkers.bqsr.ReadTransformer;

/**
 * A read transformer that refactor NDN cigar elements to one N element.
 *
 *  <p>
 *     This read transformer will refactor cigar strings that contain N-D-N elements to one N element (with total length of the three refactored elements).
 *     This is intended primarily for users of RNA-Seq data handling programs such as TopHat2.
 *     Currently we consider that the internal N-D-N motif is illegal and we error out when we encounter it. By refactoring the cigar string of
 *     those specific reads, users of TopHat and other tools can circumvent this problem without affecting the rest of their dataset.
 *
 *     NOTE: any walker that need that functionality should apply that read transformer in its map function, since it won't be activated by the GATK engine.
 *
 *     The engine parameter that activate this read transformer is --refactor_NDN_cigar_string or -fixNDN
 *  </p>
 *
 */

public final class NDNCigarReadTransformer implements ReadTransformer {

    @Override
    public SAMRecord apply(final SAMRecord read) {
        final Cigar originalCigar = read.getCigar();
        if (originalCigar.isValid(read.getReadName(),-1) != null)
            throw new IllegalArgumentException("try to transform a read with non-valid cigar string: readName: "+read.getReadName()+" Cigar String: "+originalCigar);
        read.setCigar(refactorNDNtoN(originalCigar));
        return read;
    }

    private Cigar refactorNDNtoN(final Cigar originalCigar) {
        final Cigar refactoredCigar = new Cigar();
        final int cigarLength = originalCigar.numCigarElements();
        for(int i = 0; i < cigarLength; i++){
            final CigarElement element = originalCigar.getCigarElement(i);
            if(element.getOperator() == CigarOperator.N && thereAreAtLeast2MoreElements(i,cigarLength)){
                final CigarElement nextElement = originalCigar.getCigarElement(i+1);
                final CigarElement nextNextElement = originalCigar.getCigarElement(i+2);

                // if it is N-D-N replace with N (with the total length) otherwise just add the first N.
                if(nextElement.getOperator() == CigarOperator.D && nextNextElement.getOperator() == CigarOperator.N){
                    final int threeElementsLength = element.getLength() + nextElement.getLength() + nextNextElement.getLength();
                    final CigarElement refactoredElement = new CigarElement(threeElementsLength,CigarOperator.N);
                    refactoredCigar.add(refactoredElement);
                    i += 2; //skip the elements that were refactored
                }
                else
                    refactoredCigar.add(element);  // add only the first N
            }
            else
                refactoredCigar.add(element);  // add any non-N element
        }
        return refactoredCigar;
    }

    private boolean thereAreAtLeast2MoreElements(final int index, final int cigarLength){
        return index < cigarLength - 2;
    }

}

