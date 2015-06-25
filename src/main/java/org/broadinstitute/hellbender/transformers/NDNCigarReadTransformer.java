package org.broadinstitute.hellbender.transformers;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

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
    private static final long serialVersionUID = 1L;

    @Override
    public SAMRecord apply(final SAMRecord read) {
        final Cigar originalCigar = read.getCigar();
        if (originalCigar.isValid(read.getReadName(),-1) != null)
            throw new IllegalArgumentException("try to transform a read with non-valid cigar string: readName: "+read.getReadName()+" Cigar String: "+originalCigar);
        read.setCigar(refactorNDNtoN(originalCigar));
        return read;
    }

    /**
     * Refactor cigar strings that contain N-D-N elements to one N element (with total length of the three refactored elements).
     */
    @VisibleForTesting
    protected static Cigar refactorNDNtoN(final Cigar originalCigar) {
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
                } else {
                    refactoredCigar.add(element);  // add only the first N
                }
            } else {
                refactoredCigar.add(element);  // add any non-N element
            }
        }
        return refactoredCigar;
    }

    private static boolean thereAreAtLeast2MoreElements(final int index, final int cigarLength){
        return index < cigarLength - 2;
    }

}

