package org.broadinstitute.hellbender.utils.pairhmm.pseudo;

import org.broadinstitute.gatk.nativebindings.pairhmm.HaplotypeDataHolder;
import org.broadinstitute.gatk.nativebindings.pairhmm.PairHMMNativeArguments;
import org.broadinstitute.gatk.nativebindings.pairhmm.PairHMMNativeBinding;
import org.broadinstitute.gatk.nativebindings.pairhmm.ReadDataHolder;
import org.broadinstitute.hellbender.utils.MathUtils;

import java.io.File;

/**
 * This class is a naive implementation of the native pairhmm bindings. It is NOT actually native and
 * it should NOT be used in production because it is super slow.
 * The code is extracted from GATK's LoglessPairHMM and made self-contained ie not depenenent on GATK.
 * It is only intended for help in unit testing pairhmm implementations.
 */
public final class PseudoNativePairHMM implements PairHMMNativeBinding {

    static final double INITIAL_CONDITION = Math.pow(2, 1020);
    static final double INITIAL_CONDITION_LOG10 = Math.log10(INITIAL_CONDITION);
    // we divide e by 3 because the observed base could have come from any of the non-observed alleles
    static final double TRISTATE_CORRECTION = 3.0;

    /**
     * Length of the standard transition probability array.
     */
    public static final int TRANS_PROB_ARRAY_LENGTH = 6;

    /**
     * Position in the transition probability array for the Match-to-Match transition.
     */
    public static final int matchToMatch = 0;

    /**
     * Position in the transition probability array for the Indel-to-Match transition.
     */
    public static final int indelToMatch = 1;

    /**
     * Position in the transition probability array for the Match-to-Insertion transition.
     */
    public static final int matchToInsertion = 2;

    /**
     * Position in the transition probability array for the Insertion-to-Insertion transition.
     */
    public static final int insertionToInsertion = 3;

    /**
     * Position in the transition probability array for the Match-to-Deletion transition.
     */
    public static final int matchToDeletion = 4;

    /**
     * Position in the transition probability array for the Deletion-to-Deletion transition.
     */
    public static final int deletionToDeletion = 5;

    /**
     * Convenient ln10 constant.
     */
    private static final double LN10 = Math.log(10);

    /**
     * Convenient (ln10)^-1 constant.
     */
    private static final double INV_LN10 = 1.0 / LN10;

    @Override
    public void initialize(final PairHMMNativeArguments args) {
        //ignore them
    }

    @Override
    public void computeLikelihoods(final ReadDataHolder[] readDataArray,
                                   final HaplotypeDataHolder[] haplotypeDataArray,
                                   final double[] likelihoodArray) {
        for (int i = 0; i < readDataArray.length; i++) {
            for (int j = 0; j < haplotypeDataArray.length; j++) {
                likelihoodArray[i + j] = computeLikelihood(readDataArray[i], haplotypeDataArray[j]);
            }
        }
    }

    private double computeLikelihood(final ReadDataHolder read, final HaplotypeDataHolder hap) {
        return computeLikelihood(read.readBases, read.readQuals, read.insertionGOP, read.deletionGOP, read.overallGCP, hap.haplotypeBases);
    }

    /**
     * Compute 1 likelihood of read vs haplotype.
     */
    private double computeLikelihood(final byte[] readBases,
                                     final byte[] readQuals,
                                     final byte[] insertionGOP,
                                     final byte[] deletionGOP,
                                     final byte[] overallGCP,
                                     final byte[] haplotypeBases) {

        double[][] transition = null; // The transition probabilities cache
        double[][] prior = null;      // The prior probabilities cache
        double[][] matchMatrix = null;
        double[][] insertionMatrix = null;
        double[][] deletionMatrix = null;

        final int paddedReadLength = readBases.length + 1;
        final int paddedHaplotypeLength = haplotypeBases.length + 1;

        matchMatrix = new double[paddedReadLength][paddedHaplotypeLength];
        insertionMatrix = new double[paddedReadLength][paddedHaplotypeLength];
        deletionMatrix = new double[paddedReadLength][paddedHaplotypeLength];

        transition = createTransitionMatrix(paddedReadLength);
        prior = new double[paddedReadLength][paddedHaplotypeLength];

        final double initialValue = INITIAL_CONDITION / haplotypeBases.length;
        // set the initial value (free deletions in the beginning) for the first row in the deletion matrix
        for( int j = 0; j < paddedHaplotypeLength; j++ ) {
            deletionMatrix[0][j] = initialValue;
        }

        //initialize probabilities
        PairHMMModel.qualToTransProbs(transition,insertionGOP,deletionGOP,overallGCP);
        final int hapStartIndex = 0;
        initializePriors(readBases, readQuals, haplotypeBases, prior, hapStartIndex);

        for (int i = 1; i < paddedReadLength; i++) {
            // +1 here is because hapStartIndex is 0-based, but our matrices are 1 based
            for (int j = hapStartIndex+1; j < paddedHaplotypeLength; j++) {
                //Inlined the code from updateCell - helps JIT to detect hotspots and produce good native code
                matchMatrix[i][j] = prior[i][j] * ( matchMatrix[i - 1][j - 1] * transition[i][matchToMatch] +
                        insertionMatrix[i - 1][j - 1] * transition[i][indelToMatch] +
                        deletionMatrix[i - 1][j - 1] * transition[i][indelToMatch] );
                insertionMatrix[i][j] = matchMatrix[i - 1][j] * transition[i][matchToInsertion] + insertionMatrix[i - 1][j] * transition[i][insertionToInsertion];
                deletionMatrix[i][j] = matchMatrix[i][j - 1] * transition[i][matchToDeletion] + deletionMatrix[i][j - 1] * transition[i][deletionToDeletion];
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

    private static double[][] createTransitionMatrix(final int maxReadLength) {
        return new double[maxReadLength + 1][TRANS_PROB_ARRAY_LENGTH];
    }

    private void initializePriors(byte[] readBases, byte[] readQuals, byte[] haplotypeBases, double[][] prior, int hapStartIndex) {

        // initialize the prior matrix for all combinations of read x haplotype bases
        // Abusing the fact that java initializes arrays with 0.0, so no need to fill in rows and columns below 2.

        for (int i = 0; i < readBases.length; i++) {
            final byte x = readBases[i];
            final byte qual = readQuals[i];
            for (int j = hapStartIndex; j < haplotypeBases.length; j++) {
                final byte y = haplotypeBases[j];
                prior[i+1][j+1] = ( x == y || x == (byte) 'N' || y == (byte) 'N' ?
                        qualToProb(qual) : (qualToErrorProb(qual) / TRISTATE_CORRECTION) );
            }
        }
    }

    @Override
    public void done() {
        //do nothing
    }

    @Override
    public boolean load(File tmpDir) {
        return true;//yes, we are always 'supported'
    }

    /**
     * Maximum sense quality value.
     */
    public static final int MAX_QUAL = 254;

    /**
     * Cached values for qual as byte calculations so they are very fast
     */
    private static final double[] qualToErrorProbCache = new double[MAX_QUAL + 1];
    private static final double[] qualToProbLog10Cache = new double[MAX_QUAL + 1];


    static {
        for (int i = 0; i <= MAX_QUAL; i++) {
            qualToErrorProbCache[i] = qualToErrorProb((double) i);
            qualToProbLog10Cache[i] = Math.log10(1.0 - qualToErrorProbCache[i]);
        }
    }

    /**
     * COPIED FROM QualityUtils in GATK
     *
     * Convert a phred-scaled quality score to its probability of being true (Q30 => 0.999)
     *
     * This is the Phred-style conversion, *not* the Illumina-style conversion.
     *
     * Because the input is a discretized byte value, this function uses a cache so is very efficient
     *
     * WARNING -- because this function takes a byte for maxQual, you must be careful in converting
     * integers to byte.  The appropriate way to do this is ((byte)(myInt & 0xFF))
     *
     * @param qual a quality score (0-255)
     * @return a probability (0.0-1.0)
     */
    public static double qualToProb(final byte qual) {
        return 1.0 - qualToErrorProb(qual);
    }

    /**
     * COPIED FROM QualityUtils in GATK
     * Convert a phred-scaled quality score to its probability of being wrong (Q30 => 0.001)
     *
     * This is the Phred-style conversion, *not* the Illumina-style conversion.
     *
     * Because the input is a double value, this function must call Math.pow so can be quite expensive
     *
     * @param qual a phred-scaled quality score encoded as a double.  Can be non-integer values (30.5)
     * @return a probability (0.0-1.0)
     */
    public static double qualToErrorProb(final double qual) {
        return Math.pow(10.0, qual / -10.0);
    }

    /**
     * COPIED FROM QualityUtils in GATK
     *
     * Convert a phred-scaled quality score to its probability of being wrong (Q30 => 0.001)
     *
     * This is the Phred-style conversion, *not* the Illumina-style conversion.
     *
     * Because the input is a byte value, this function uses a cache so is very efficient
     *
     * WARNING -- because this function takes a byte for maxQual, you must be careful in converting
     * integers to byte.  The appropriate way to do this is ((byte)(myInt & 0xFF))
     *
     * @param qual a phred-scaled quality score encoded as a byte
     * @return a probability (0.0-1.0)
     */
    public static double qualToErrorProb(final byte qual) {
        return qualToErrorProbCache[(int)qual & 0xff]; // Map: 127 -> 127; -128 -> 128; -1 -> 255; etc.
    }

    public static void qualToTransProbs(final double[][] dest, final byte[] insQuals, final byte[] delQuals, final byte[] gcps) {
        final int readLength = insQuals.length;
        if (delQuals.length != readLength) throw new IllegalArgumentException("deletion quality array length does not match insert quality array length: " + readLength + " != " + delQuals.length);
        if (gcps.length != readLength) throw new IllegalArgumentException("deletion quality array length does not match insert quality array length: " + readLength + " != " + gcps.length);
        if (dest.length < readLength + 1) throw new IllegalArgumentException("destination length is not enough for the read length: " + dest.length + " < " + readLength + " + 1");

        for (int i = 0; i < readLength; i++) {
            qualToTransProbs(dest[i + 1], insQuals[i], delQuals[i], gcps[i]);
        }
    }

    public static void qualToTransProbs(final double[] dest, final byte insQual, final byte delQual, final byte gcp) {
        if (insQual < 0) throw new IllegalArgumentException("insert quality cannot less than 0: " + insQual);
        if (delQual < 0) throw new IllegalArgumentException("deletion quality cannot be less than 0: " + delQual);
        if (gcp < 0) throw new IllegalArgumentException("gcp cannot be less than 0: " + gcp);
        dest[matchToMatch] = matchToMatchProb(insQual, delQual);
        dest[matchToInsertion] = qualToErrorProb(insQual);
        dest[matchToDeletion] = qualToErrorProb(delQual);
        dest[indelToMatch] = qualToProb(gcp);
        dest[insertionToInsertion] = dest[deletionToDeletion] = qualToErrorProb(gcp);
    }

    private static double matchToMatchProb(final byte insQual, final byte delQual) {
        return matchToMatchProb((insQual & 0xFF), (delQual & 0xFF));
    }

    private static double matchToMatchProb(final int insQual, final int delQual) {
        final int minQual;
        final int maxQual;
        if (insQual <= delQual) {
            minQual = insQual;
            maxQual = delQual;
        } else {
            minQual = delQual;
            maxQual = insQual;
        }

        if (minQual < 0) throw new IllegalArgumentException("quality cannot be negative: " + minQual + " and " + maxQual);

        return (QualityUtils.MAX_QUAL < maxQual) ?  1.0 - Math.pow(10, MathUtils.approximateLog10SumLog10(-0.1 * minQual, -0.1 * maxQual)) :
                matchToMatchProb[((maxQual * (maxQual + 1)) >> 1) + minQual];
    }

    /**
     * Holds pre-calculated the matchToMatch log10-probabilities.
     *
     * <p/>
     * This is a triangular matrix stored in a unidimentional array like so:
     * <p/>
     * (0,0), (0,1), (1,1), (0,2), (1,2), (2,2), (0,3) ... ({@link QualityUtils#MAX_QUAL},{@link QualityUtils#MAX_QUAL})
     */
    private static final double[] matchToMatchLog10 = new double[((QualityUtils.MAX_QUAL + 1) * (QualityUtils.MAX_QUAL + 2)) >> 1];

    /**
     * Holds pre-calculated the matchToMath probability values in linear scale.
     *
     * <p/>
     * This is a triangular matrix stored in a unidimentional array like so:
     * <p/>
     * (0,0), (0,1), (1,1), (0,2), (1,2), (2,2), (0,3) ... ({@link QualityUtils#MAX_QUAL},{@link QualityUtils#MAX_QUAL})
     */
    private static final double[] matchToMatchProb = new double[((QualityUtils.MAX_QUAL + 1) * (QualityUtils.MAX_QUAL + 2)) >> 1];

    /**
     * Initialize matchToMatch cache tables {@link #matchToMatch} and {@link #matchToMatchLog}
     */
    static {
        for (int i = 0, offset = 0; i <= QualityUtils.MAX_QUAL; offset += ++i)
            for (int j = 0; j <= i; j++) {
                final double log10Sum = MathUtils.approximateLog10SumLog10(-0.1 * i, -0.1 * j);
                matchToMatchLog10[offset + j] =
                        Math.log1p(-Math.min(1, Math.pow(10, log10Sum))) * INV_LN10;
                matchToMatchProb[offset + j] = Math.pow(10, matchToMatchLog10[offset + j]);
            }
    }
}
