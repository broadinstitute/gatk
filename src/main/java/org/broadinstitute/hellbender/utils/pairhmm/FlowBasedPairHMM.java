package org.broadinstitute.hellbender.utils.pairhmm;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.FlowBasedAlignmentLikelihoodEngine;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.haplotype.FlowBasedHaplotype;
import org.broadinstitute.hellbender.utils.read.FlowBasedKeyCodec;
import org.broadinstitute.hellbender.utils.read.FlowBasedRead;

import java.util.List;

import static org.broadinstitute.hellbender.utils.pairhmm.PairHMMModel.*;
import static org.broadinstitute.hellbender.utils.pairhmm.PairHMMModel.deletionToDeletion;
import static org.broadinstitute.hellbender.utils.pairhmm.PairHMMModel.matchToDeletion;

/**
 * Class for performing the pair HMM for global alignment in FlowSpace. See {@link FlowBasedAlignmentLikelihoodEngine} and {@link LoglessPairHMM}.
 */
public class FlowBasedPairHMM extends PairHMM {
    protected static final Logger logger = LogManager.getLogger(FlowBasedPairHMM.class);
    private static final int DEFAULT_INDEL_NO_DATA_FILL = 40;
    private static final int DEFAULT_MATCH_NO_DATA_FILL = 10;

    protected boolean constantsAreInitialized = false;

    protected int[] previousHaplotypeBases;

    @Override
    protected double subComputeReadLikelihoodGivenHaplotypeLog10(byte[] haplotypeBases, byte[] readBases, byte[] readQuals, byte[] insertionGOP, byte[] deletionGOP, byte[] overallGCP, int hapStartIndex, boolean recacheReadValues, int nextHapStartIndex) {
        throw new UnsupportedOperationException("this should not be called");
    }

    protected double[][] transition = null; // The transition probabilities cache
    protected double[][] prior = null;      // The prior probabilities cache
    protected double[][] matchMatrix = null;
    protected double[][] insertionMatrix = null;
    protected double[][] deletionMatrix = null;

    //NOTE: if you ever find yourself wanting to change this parameter I would seriously reconsider if it will be possible to adapt this
    //      model to work in that case. The alignment and model math is rooted in the assumption that frame-shifts necessarily happen in
    //      discrete units of 4 and if that assumption is broken then this HMM is unlikely to be representative without rewriting.
    final static public int FLOW_SIZE = 4;

    /**
     * Initializes the matrix that holds all the constants related to the editing
     * distance between the read and the haplotype.
     *
     * @param haplotypeFlows the bases of the haplotype
     * @param readFlows      the bases of the read
     * @param read           the read with the ability to compute qual
     */
    void initializePriors(final int[] haplotypeFlows, final byte[] haplotypeFlowOrder, final int[] readFlows, final byte[] readFlowOrder, final FlowBasedRead read) {

        // initialize the prior matrix for all combinations of read x haplotype bases
        // Abusing the fact that java initializes arrays with 0.0, so no need to fill in rows and columns below 2.

        for (int i = 0; i < readFlows.length; i++) {
            for (int j = 0; j < haplotypeFlows.length; j++) {
                if (haplotypeFlowOrder[j] == readFlowOrder[i]) {
                    prior[i + 1 + FLOW_SIZE][j + 1 + FLOW_SIZE] = read.getProb(i, haplotypeFlows[j]);
                } else {
                    prior[i + 1 + FLOW_SIZE][j + 1 + FLOW_SIZE] = 0;
                }
            }
        }
    }

    static void initializeProbabilities(final double[][] transition, final byte[] insertionGOP, final byte[] deletionGOP, final byte[] overallGCP) {
        qualToTransProbs(transition,insertionGOP,deletionGOP,overallGCP);
    }

    static int findMaxReadLengthFlow(final List<FlowBasedRead> reads) {

        return reads.stream()
                .map(read -> read.getKey().length)
                .mapToInt(value -> value)
                .max().orElse(0);
    }

    private static int findMaxAlleleLengthFlow(final List<FlowBasedHaplotype> alleles) {

        return alleles.stream()
                .map(allele -> allele.getKeyLength())
                .mapToInt(value -> value)
                .max().orElse(0);
    }

    /**
     *  Given a list of reads and haplotypes, for every read compute the total probability of said read arising from
     *  each haplotype given base substitution, insertion, and deletion probabilities.
     *
     * @param processedReads reads to analyze instead of the ones present in the destination read-likelihoods.
     * @param logLikelihoods where to store the log likelihoods where position [a][r] is reserved for the log likelihood of {@code reads[r]}
     *             conditional to {@code alleles[a]}.
     */
    public void computeLog10LikelihoodsFlowBased(final LikelihoodMatrix<GATKRead, Haplotype> logLikelihoods,
                                                 final List<FlowBasedRead> processedReads,
                                                 final List<FlowBasedHaplotype> processedHaplotypes) {
        if (processedReads.isEmpty()) {
            return;
        }
        if(doProfiling) {
            startTime = System.nanoTime();
        }
        // (re)initialize the pairHMM only if necessary
        final int readMaxLength = findMaxReadLengthFlow(processedReads);
        final int haplotypeMaxLength = findMaxAlleleLengthFlow(processedHaplotypes);
        if (!initialized || readMaxLength > maxReadLength || haplotypeMaxLength > maxHaplotypeLength) {
            initialize(readMaxLength, haplotypeMaxLength);
        }

        final int readCount = processedReads.size();
        final List<Haplotype> alleles = logLikelihoods.alleles();
        final int alleleCount = alleles.size();
        mLogLikelihoodArray = new double[readCount * alleleCount];
        int idx = 0;
        int readIndex = 0;
        for(final FlowBasedRead read : processedReads){

            final int[] readKey = read.getKey();

            final byte[] readInsQuals = read.getReadInsQuals();
            final byte[] readDelQuals = read.getReadDelQuals();
            final byte[] overallGCP = read.getOverallGCP();

            final byte[] flowReadInsQuals = FlowBasedKeyCodec.baseArrayToKeySpace(read.getBases(), read.getKeyLength(), readInsQuals, (byte) DEFAULT_INDEL_NO_DATA_FILL, read.getFlowOrder());
            final byte[] flowReadDelQuals = FlowBasedKeyCodec.baseArrayToKeySpace(read.getBases(), read.getKeyLength(), readDelQuals, (byte) DEFAULT_INDEL_NO_DATA_FILL, read.getFlowOrder());
            final byte[] flowReadGCPQuals = FlowBasedKeyCodec.baseArrayToKeySpace(read.getBases(), read.getKeyLength(), overallGCP, (byte) DEFAULT_MATCH_NO_DATA_FILL, read.getFlowOrder());

            // peek at the next haplotype in the list (necessary to get nextHaplotypeBases, which is required for caching in the array implementation)
            final boolean isFirstHaplotype = true;
            for (int a = 0; a < alleleCount; a++) {
                final FlowBasedHaplotype allele = processedHaplotypes.get(a);
                final int[] alleleKey = allele.getKey();

                final byte[] flowOrder = allele.getFlowOrderArray();
                int startingPoint = 0;
                for (int i = 0; i < flowOrder.length; i++) {
                    if (flowOrder[i] == read.getFlowOrderArray()[0]) {
                        startingPoint = i;
                        break;
                    }
                }

                final double lk =  subComputekeayLikelihoodGivenHaplotypeKeysLog10(startingPoint, alleleKey, flowOrder,
                        readKey, read.getFlowOrderArray(), read, flowReadInsQuals, flowReadDelQuals, flowReadGCPQuals, isFirstHaplotype);
                logLikelihoods.set(a, readIndex, lk);
                mLogLikelihoodArray[idx++] = lk;
            }
            readIndex++;
        }
        if(doProfiling) {
            threadLocalPairHMMComputeTimeDiff = (System.nanoTime() - startTime);
            {
                pairHMMComputeTime += threadLocalPairHMMComputeTimeDiff;
            }
        }
    }




    protected double subComputekeayLikelihoodGivenHaplotypeKeysLog10(int hapStartIndex, int[] haplotypeKey, final byte[] haplotypeFlowOrder, int[] readKey, final byte[] readFlowOrder,
                                                                     final FlowBasedRead read, byte[] insertionGOP, byte[] deletionGOP, byte[] overallGCP, boolean recacheReadValues) {

        int thisReadPaddedLength = readKey.length + 1 + FLOW_SIZE;
        int thisHaplotypePaddedLenth = haplotypeKey.length + 1 + FLOW_SIZE;

//        //TODO THIS VERSION MIGHT BE NECESSARY FOR DEBUGGING IN THE FUTURE. NOTE THIS, IT MIGHT SAVE YOU
        // In a previous iteraiton of this code there was the risk that likelihoods scores from previous runs of this tool
        // were accidentally making their way into the math for the reads. This is a new risk compared to the regular PairHMM
        // since for any read/haplotype combination only 1/4th of the array fields actually get filled. This means there is
        // an outside risk that if the wrong array elements are pulled upon in the following code noise could be introduced
        // for reads/haplotypes based on previous executions of the hmm. If its happening now it has proven very difficult to track.
//        for( int i = FLOW_SIZE + 1; i < thisReadPaddedLength; i++ ) {
//            for( int j = FLOW_SIZE + 1; j < thisHaplotypePaddedLenth; j++ ) {
//                deletionMatrix[i][j] = 0.0;
//                insertionMatrix[i][j] = 0.0;
//                matchMatrix[i][j] = 0.0;
//            }
//
//        }
        // Zero out the array elements that will later be summed summed across likelihoods to eliminate the risk that
        // scores from previous reads/haplotypes might affect these results.
        final int endI = thisReadPaddedLength - 1;
        for (int j = 1; j < thisHaplotypePaddedLenth; j++) {
            matchMatrix[endI][j] = 0;
            insertionMatrix[endI][j] = 0;
        }


        // TODO maybe handle previous haplotype computation...
        if (previousHaplotypeBases == null || previousHaplotypeBases.length != haplotypeKey.length) {
            final double initialValue = LoglessPairHMM.INITIAL_CONDITION / haplotypeKey.length;
            // set the initial value (free deletions in the beginning) for the first row in the deletion matrix
            for( int j = 0; j < thisHaplotypePaddedLenth; j++ ) {
                deletionMatrix[0][j] = initialValue;
                deletionMatrix[1][j] = initialValue;
                deletionMatrix[2][j] = initialValue;
                deletionMatrix[3][j] = initialValue;
                deletionMatrix[4][j] = initialValue;
            }
        }
        previousHaplotypeBases = haplotypeKey;

        if ( ! constantsAreInitialized || recacheReadValues ) {
            initializeProbabilities(transition, insertionGOP, deletionGOP, overallGCP);

            // note that we initialized the constants
            constantsAreInitialized = true;
        }

        initializePriors(haplotypeKey, haplotypeFlowOrder, readKey, readFlowOrder, read);

        for (int i = 1+FLOW_SIZE; i < thisReadPaddedLength; i++) {
            // +1 here is because hapStartIndex is 0-based, but our matrices are 1 based
            // Only loop over array elements that correspond to valid flow-flow base comparisons
            for (int j = ((hapStartIndex + i)%FLOW_SIZE) + FLOW_SIZE; j < thisHaplotypePaddedLenth; j+=FLOW_SIZE) {
                //Inlined the code from updateCell - helps JIT to detect hotspots and produce good native code -- this reflects the LoglessPairHMM optimizations (which this code is based off of).
                matchMatrix[i][j] =  prior[i][j] * ( matchMatrix[i - 1][j - 1] * transition[i - FLOW_SIZE][matchToMatch] +
                        insertionMatrix[i - 1][j - 1] * transition[i - FLOW_SIZE][indelToMatch] +
                        deletionMatrix[i - 1][j - 1] * transition[i - FLOW_SIZE][indelToMatch] );
                insertionMatrix[i][j] = matchMatrix[i - FLOW_SIZE][j] * transition[i - FLOW_SIZE][matchToInsertion] + insertionMatrix[i - FLOW_SIZE][j] * transition[i - FLOW_SIZE][insertionToInsertion];
                deletionMatrix[i][j] = matchMatrix[i][j - FLOW_SIZE] * transition[i - FLOW_SIZE][matchToDeletion] + deletionMatrix[i][j - FLOW_SIZE] * transition[i - FLOW_SIZE][deletionToDeletion];
            }
        }

        // final log probability is the log10 sum of the last element in the Match and Insertion state arrays
        // this way we ignore all paths that ended in deletions! (huge)
        // but we have to sum all the paths ending in the M and I matrices, because they're no longer extended.
        double finalSumProbabilities = 0.0;
        for (int j = 1; j < thisHaplotypePaddedLenth; j++) {
            finalSumProbabilities += matchMatrix[endI][j] + insertionMatrix[endI][j];
        }

        return Math.log10(finalSumProbabilities) - LoglessPairHMM.INITIAL_CONDITION_LOG10;
    }


    public void initialize( final int readMaxLength, final int haplotypeMaxLength ) throws IllegalArgumentException {
        Utils.validateArg(haplotypeMaxLength > 0, () -> "haplotypeMaxLength must be > 0 but got " + haplotypeMaxLength);

        maxHaplotypeLength = haplotypeMaxLength;
        maxReadLength = readMaxLength;

        // M, X, and Y arrays are of size read and haplotype + 1 because of an extra column for initial conditions and + 1 to consider the final base in a non-global alignment
        // NOTE we add FLOW_SIZE here to account for backtracing indels
        paddedMaxReadLength = readMaxLength + 1 + FLOW_SIZE;
        paddedMaxHaplotypeLength = haplotypeMaxLength + 1 + FLOW_SIZE;

        previousHaplotypeBases = null;

        constantsAreInitialized = false;
        initialized = true;

        matchMatrix = new double[paddedMaxReadLength][paddedMaxHaplotypeLength];
        insertionMatrix = new double[paddedMaxReadLength][paddedMaxHaplotypeLength];
        deletionMatrix = new double[paddedMaxReadLength][paddedMaxHaplotypeLength];

        transition = PairHMMModel.createTransitionMatrix(paddedMaxReadLength);
        prior = new double[paddedMaxReadLength][paddedMaxHaplotypeLength];
    }
}
