package org.broadinstitute.hellbender.utils.locusiterator;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.broadinstitute.hellbender.utils.read.GATKRead;

public final class LIBS_position {

    GATKRead read;

    final int numOperators;
    int currentOperatorIndex = 0;
    int currentPositionOnOperator = 0;
    int currentReadOffset = 0;
    int currentGenomeOffset = 0;

    public boolean isBeforeDeletionStart = false;
    public boolean isBeforeDeletedBase = false;
    public boolean isAfterDeletionEnd = false;
    public boolean isAfterDeletedBase = false;
    public boolean isBeforeInsertion = false;
    public boolean isAfterInsertion = false;
    public boolean isNextToSoftClip = false;

    boolean sawMop = false;

    public LIBS_position(final GATKRead read) {
        this.read = read;
        numOperators = read.numCigarElements();
    }

    public int getCurrentReadOffset() {
        return Math.max(0, currentReadOffset - 1);
    }

    public int getCurrentPositionOnOperatorBase0() {
        return currentPositionOnOperator - 1;
    }

    public int getCurrentGenomeOffsetBase0() {
        return currentGenomeOffset - 1;
    }

    /**
     * Steps forward on the genome.  Returns false when done reading the read, true otherwise.
     */
    @SuppressWarnings("fallthrough")
    public boolean stepForwardOnGenome() {
        if ( currentOperatorIndex == numOperators )
            return false;

        CigarElement curElement = read.getCigar().getCigarElement(currentOperatorIndex);
        if ( currentPositionOnOperator >= curElement.getLength() ) {
            if ( ++currentOperatorIndex == numOperators )
                return false;

            curElement = read.getCigar().getCigarElement(currentOperatorIndex);
            currentPositionOnOperator = 0;
        }

        switch ( curElement.getOperator() ) {
            case I: // insertion w.r.t. the reference
//                if ( !sawMop )
//                    break;
            case S: // soft clip
                currentReadOffset += curElement.getLength();
            case H: // hard clip
            case P: // padding
                currentOperatorIndex++;
                return stepForwardOnGenome();

            case D: // deletion w.r.t. the reference
            case N: // reference skip (looks and gets processed just like a "deletion", just different logical meaning)
                currentPositionOnOperator++;
                currentGenomeOffset++;
                break;

            case M:
            case EQ:
            case X:
                sawMop = true;
                currentReadOffset++;
                currentPositionOnOperator++;
                currentGenomeOffset++;
                break;
            default:
                throw new IllegalStateException("No support for cigar op: " + curElement.getOperator());
        }

        final boolean isFirstOp = currentOperatorIndex == 0;
        final boolean isLastOp = currentOperatorIndex == numOperators - 1;
        final boolean isFirstBaseOfOp = currentPositionOnOperator == 1;
        final boolean isLastBaseOfOp = currentPositionOnOperator == curElement.getLength();

        isBeforeDeletionStart = isBeforeOp(read.getCigar(), currentOperatorIndex, CigarOperator.D, isLastOp, isLastBaseOfOp);
        isBeforeDeletedBase = isBeforeDeletionStart || (!isLastBaseOfOp && curElement.getOperator() == CigarOperator.D);
        isAfterDeletionEnd = isAfterOp(read.getCigar(), currentOperatorIndex, CigarOperator.D, isFirstOp, isFirstBaseOfOp);
        isAfterDeletedBase  = isAfterDeletionEnd || (!isFirstBaseOfOp && curElement.getOperator() == CigarOperator.D);
        isBeforeInsertion = isBeforeOp(read.getCigar(), currentOperatorIndex, CigarOperator.I, isLastOp, isLastBaseOfOp)
                || (!sawMop && curElement.getOperator() == CigarOperator.I);
        isAfterInsertion = isAfterOp(read.getCigar(), currentOperatorIndex, CigarOperator.I, isFirstOp, isFirstBaseOfOp);
        isNextToSoftClip = isBeforeOp(read.getCigar(), currentOperatorIndex, CigarOperator.S, isLastOp, isLastBaseOfOp)
                || isAfterOp(read.getCigar(), currentOperatorIndex, CigarOperator.S, isFirstOp, isFirstBaseOfOp);

        return true;
    }

    private static boolean isBeforeOp(final Cigar cigar,
                                      final int currentOperatorIndex,
                                      final CigarOperator op,
                                      final boolean isLastOp,
                                      final boolean isLastBaseOfOp) {
        return  !isLastOp && isLastBaseOfOp && cigar.getCigarElement(currentOperatorIndex+1).getOperator() == op;
    }

    private static boolean isAfterOp(final Cigar cigar,
                                     final int currentOperatorIndex,
                                     final CigarOperator op,
                                     final boolean isFirstOp,
                                     final boolean isFirstBaseOfOp) {
        return  !isFirstOp && isFirstBaseOfOp && cigar.getCigarElement(currentOperatorIndex-1).getOperator() == op;
    }
}
