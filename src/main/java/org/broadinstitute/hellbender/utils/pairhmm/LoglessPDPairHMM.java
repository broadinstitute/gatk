package org.broadinstitute.hellbender.utils.pairhmm;

import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.haplotype.PartiallyDeterminedHaplotype;

import java.util.Arrays;

import static org.broadinstitute.hellbender.utils.pairhmm.PairHMMModel.*;

public class LoglessPDPairHMM extends N2MemoryPDPairHMM {
    static final double INITIAL_CONDITION = Math.pow(2, 1020);
    static final double INITIAL_CONDITION_LOG10 = Math.log10(INITIAL_CONDITION);

    // we divide e by 3 because the observed base could have come from any of the non-observed alleles
    static final double TRISTATE_CORRECTION = 3.0;

    //Branch Stat Arrays
    protected double[][] branchMatchMatrix = null;
    protected double[][] branchInsertionMatrix = null;
    protected double[][] branchDeletionMatrix = null;

    enum HMMState {
        // The regular state
        NORMAL,

        // Idicating that we must be copying the array elements to the right
        INSIDE_DEL,

        // Indicating that we must handle the special case for merging events after the del
        AFTER_DEL,
    }

    public double subComputeReadLikelihoodGivenHaplotypeLog10( final byte[] haplotypeBases,
        final byte[] haplotypePDBases,
        final byte[] readBases,
        final byte[] readQuals,
        final byte[] insertionGOP,
        final byte[] deletionGOP,
        final byte[] overallGCP,
        final int hapStartIndex,
        final boolean recacheReadValues,
        final int nextHapStartIndex) {

        if (previousHaplotypeBases == null || previousHaplotypeBases.length != haplotypeBases.length) {
            final double initialValue = INITIAL_CONDITION / haplotypeBases.length;
            // set the initial value (free deletions in the beginning) for the first row in the deletion matrix
            for( int j = 0; j < paddedHaplotypeLength; j++ ) {
                deletionMatrix[0][j] = initialValue;
            }
        }

        if ( ! constantsAreInitialized || recacheReadValues ) {
            initializeProbabilities(transition, insertionGOP, deletionGOP, overallGCP);

            // note that we initialized the constants
            constantsAreInitialized = true;
        }

        initializePriors(haplotypeBases, haplotypePDBases, readBases, readQuals, hapStartIndex);

        HMMState currentState = HMMState.NORMAL;
        for (int i = 1; i < paddedReadLength; i++) {
            for (int j = hapStartIndex+1; j < paddedHaplotypeLength; j++) {
                switch (currentState) {
                    case NORMAL:
                        // Pre-emptively copy over from the left
                        branchMatchMatrix[i][j] = matchMatrix[i][j - 1];
                        branchDeletionMatrix[i][j] = deletionMatrix[i][j - 1];
                        branchInsertionMatrix[i][j] = insertionMatrix[i][j - 1];

                        //Inlined the code from updateCell - helps JIT to detect hotspots and produce good native code
                        matchMatrix[i][j] = prior[i][j] * (matchMatrix[i - 1][j - 1] * transition[i][matchToMatch] +
                                insertionMatrix[i - 1][j - 1] * transition[i][indelToMatch] +
                                deletionMatrix[i - 1][j - 1] * transition[i][indelToMatch]);
                        deletionMatrix[i][j] = matchMatrix[i][j - 1] * transition[i][matchToDeletion] + deletionMatrix[i][j - 1] * transition[i][deletionToDeletion];

                        // Special case since we MIGHT be at the deletion end and have to be handling jump states from above
                        if ((haplotypePDBases[j-1] & PartiallyDeterminedHaplotype.DEL_END) == PartiallyDeterminedHaplotype.DEL_END) {
                            insertionMatrix[i][j] = Math.max(branchMatchMatrix[i - 1][j], matchMatrix[i - 1][j]) * transition[i][matchToInsertion]
                                    + Math.max(branchInsertionMatrix[i - 1][j], insertionMatrix[i - 1][j]) * transition[i][insertionToInsertion];
                        } else {
                            insertionMatrix[i][j] = matchMatrix[i - 1][j] * transition[i][matchToInsertion] + insertionMatrix[i - 1][j] * transition[i][insertionToInsertion];
                        }
                        break;

                    case INSIDE_DEL:
                        // If we are in a deletion copy from the left
                        branchMatchMatrix[i][j] = branchMatchMatrix[i][j - 1];
                        branchDeletionMatrix[i][j] = branchDeletionMatrix[i][j - 1];
                        branchInsertionMatrix[i][j] = branchInsertionMatrix[i][j - 1];


                        //Inlined the code from updateCell - helps JIT to detect hotspots and produce good native code
                        matchMatrix[i][j] = prior[i][j] * (matchMatrix[i - 1][j - 1] * transition[i][matchToMatch] +
                                insertionMatrix[i - 1][j - 1] * transition[i][indelToMatch] +
                                deletionMatrix[i - 1][j - 1] * transition[i][indelToMatch]);
                        deletionMatrix[i][j] = matchMatrix[i][j - 1] * transition[i][matchToDeletion] + deletionMatrix[i][j - 1] * transition[i][deletionToDeletion];

                        // Special case since we MIGHT be at the deletion end and have to be handling jump states from above
                        // TODO these indexes should probably be -1
                        if ((haplotypePDBases[j-1] & PartiallyDeterminedHaplotype.DEL_END) == PartiallyDeterminedHaplotype.DEL_END) {
                            insertionMatrix[i][j] = Math.max(branchMatchMatrix[i - 1][j], matchMatrix[i - 1][j]) * transition[i][matchToInsertion]
                                    + Math.max(branchInsertionMatrix[i - 1][j], insertionMatrix[i - 1][j]) * transition[i][insertionToInsertion];
                        } else {
                            insertionMatrix[i][j] = matchMatrix[i - 1][j] * transition[i][matchToInsertion] + insertionMatrix[i - 1][j] * transition[i][insertionToInsertion];
                        }
                        break;

                    case AFTER_DEL:
                        // Pre-emptively copy over from the left
                        branchMatchMatrix[i][j] = Math.max(branchMatchMatrix[i][j - 1], matchMatrix[i][j - 1]);
                        branchDeletionMatrix[i][j] = Math.max(branchDeletionMatrix[i][j - 1], deletionMatrix[i][j - 1]);
                        branchInsertionMatrix[i][j] = Math.max(branchInsertionMatrix[i][j - 1],  insertionMatrix[i][j - 1]);

                        //Inlined the code from updateCell - helps JIT to detect hotspots and produce good native code
                        matchMatrix[i][j] = prior[i][j] * (Math.max(branchMatchMatrix[i - 1][j - 1], matchMatrix[i - 1][j - 1]) * transition[i][matchToMatch] +
                                Math.max(branchInsertionMatrix[i - 1][j - 1], insertionMatrix[i - 1][j - 1]) * transition[i][indelToMatch] +
                                Math.max(branchDeletionMatrix[i - 1][j - 1], deletionMatrix[i - 1][j - 1]) * transition[i][indelToMatch]);
                        deletionMatrix[i][j] = Math.max(branchMatchMatrix[i][j - 1], matchMatrix[i][j - 1]) * transition[i][matchToDeletion]
                                + Math.max(branchDeletionMatrix[i][j - 1], deletionMatrix[i][j - 1]) * transition[i][deletionToDeletion];

                        // Special case since we MIGHT be at the deletion end and have to be handling jump states from above
                        if ((haplotypePDBases[j-1] & PartiallyDeterminedHaplotype.DEL_END) == PartiallyDeterminedHaplotype.DEL_END) {
                            insertionMatrix[i][j] = Math.max(branchMatchMatrix[i - 1][j], matchMatrix[i - 1][j]) * transition[i][matchToInsertion]
                                    + Math.max(branchInsertionMatrix[i - 1][j], insertionMatrix[i - 1][j]) * transition[i][insertionToInsertion];
                        } else {
                            insertionMatrix[i][j] = matchMatrix[i - 1][j] * transition[i][matchToInsertion] + insertionMatrix[i - 1][j] * transition[i][insertionToInsertion];
                        }
                        currentState = HMMState.NORMAL;
                        break;
                }
                // If we are at a deletion start base, start copying the branch states
                if ((haplotypePDBases[j-1] & PartiallyDeterminedHaplotype.DEL_START) == PartiallyDeterminedHaplotype.DEL_START) {
                    currentState = HMMState.INSIDE_DEL;
                }
                // Being after a deltion overrides being inside a deletion by virtue of the fact that we allow single element deletions
                if ((haplotypePDBases[j-1] & PartiallyDeterminedHaplotype.DEL_END) == PartiallyDeterminedHaplotype.DEL_END) {
                    currentState = HMMState.AFTER_DEL;
                }
            }
        }

        // final log probability is the log10 sum of the last element in the Match and Insertion state arrays
        // this way we ignore all paths that ended in deletions! (huge)
        // but we have to sum all the paths ending in the M and I matrices, because they're no longer extended.
        final int endI = paddedReadLength - 1;
        double finalSumProbabilities = 0.0;
        for (int j = 1; j < paddedHaplotypeLength; j++) {
            finalSumProbabilities += matchMatrix[endI][j] + insertionMatrix[endI][j];
        }
        return Math.log10(finalSumProbabilities) - INITIAL_CONDITION_LOG10;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void initialize(final int readMaxLength, final int haplotypeMaxLength ) {
        super.initialize(readMaxLength, haplotypeMaxLength);

        branchMatchMatrix = new double[paddedMaxReadLength][paddedMaxHaplotypeLength];
        branchInsertionMatrix = new double[paddedMaxReadLength][paddedMaxHaplotypeLength];
        branchDeletionMatrix = new double[paddedMaxReadLength][paddedMaxHaplotypeLength];
    }


    /**
     * Initializes the matrix that holds all the constants related to the editing
     * distance between the read and the haplotype.
     *
     * @param haplotypeBases the bases of the haplotype
     * @param readBases      the bases of the read
     * @param readQuals      the base quality scores of the read
     * @param startIndex     where to start updating the distanceMatrix (in case this read is similar to the previous read)
     */
    void initializePriors(final byte[] haplotypeBases, final byte[] haplotypePDBases, final byte[] readBases, final byte[] readQuals, final int startIndex) {

        // initialize the prior matrix for all combinations of read x haplotype bases
        // Abusing the fact that java initializes arrays with 0.0, so no need to fill in rows and columns below 2.

        for (int i = 0; i < readBases.length; i++) {
            final byte x = readBases[i];
            final byte qual = readQuals[i];
            for (int j = startIndex; j < haplotypeBases.length; j++) {
                final byte y = haplotypeBases[j];
                final byte hapPDBase = haplotypePDBases[j];
                prior[i+1][j+1] = ( x == y || x == (byte) 'N' || y == (byte) 'N' || isBasePDMatching(x,hapPDBase) ?
                        QualityUtils.qualToProb(qual) : (QualityUtils.qualToErrorProb(qual) / (doNotUseTristateCorrection ? 1.0 : TRISTATE_CORRECTION)) );
            }
        }
    }
    //Helper method for handling if a pd base counts as matching
    static boolean isBasePDMatching(byte x, byte hapPDBases) {
        if ((hapPDBases & PartiallyDeterminedHaplotype.SNP) != 0) {
            switch (x) {
                case 'A':
                case 'a':
                    return ((hapPDBases & PartiallyDeterminedHaplotype.A) != 0);
                case 'C':
                case 'c':
                    return ((hapPDBases & PartiallyDeterminedHaplotype.C) != 0);
                case 'T':
                case 't':
                    return ((hapPDBases & PartiallyDeterminedHaplotype.T) != 0);
                case 'G':
                case 'g':
                    return ((hapPDBases & PartiallyDeterminedHaplotype.G) != 0);
                default:
                    throw new RuntimeException("Found unexpected base in alt alleles");
            }
        }
        return false;
    }


    /**
     * Initializes the matrix that holds all the constants related to quality scores.
     *
     * @param insertionGOP   insertion quality scores of the read
     * @param deletionGOP    deletion quality scores of the read
     * @param overallGCP     overall gap continuation penalty
     */
    static void initializeProbabilities(final double[][] transition, final byte[] insertionGOP, final byte[] deletionGOP, final byte[] overallGCP) {
        PairHMMModel.qualToTransProbs(transition,insertionGOP,deletionGOP,overallGCP);
    }
}
