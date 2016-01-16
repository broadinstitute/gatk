package org.broadinstitute.hellbender.utils.pairhmm;

import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;

/**
 * Helper class that implement calculations required to implement the PairHMM Finite State Automation (FSA) model.
 */
public final class PairHMMModel {

    /**
     * Prevents instantiation of this class
     */
    private PairHMMModel() {
    }

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
     * Holds pre-calculated the matchToMatch log10-probabilities.
     *
     * <p/>
     * This is a triangular matrix stored in a unidimentional array like so:
     * <p/>
     * (0,0), (0,1), (1,1), (0,2), (1,2), (2,2), (0,3) ... ({@link QualityUtils#MAX_QUAL},{@link QualityUtils#MAX_QUAL})
     */
    private static final double[] matchToMatchLog10 = new double[((QualityUtils.MAX_QUAL + 1) * (QualityUtils.MAX_QUAL + 2)) >> 1];

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

    /**
     * Fills a transition probability array given the different quality scores affecting a read site
     *
     * @param insQual the insertion quality score as a byte.
     * @param delQual the deletion quality score as a byte.
     * @param gcp the gap-continuation-penalty score as a byte.
     *
     * @throws NullPointerException if {@code dest} is {@code null}.
     * @throws ArrayIndexOutOfBoundsException if {@code dest} is not large enough.
     * @throws IllegalArgumentException if {@code insQual}, {@code delQual} or {@code gcp} is less than negative.
     */
    public static void qualToTransProbs(final double[] dest, final byte insQual, final byte delQual, final byte gcp) {
        Utils.nonNull(dest, "dest array null");
        if (insQual < 0) throw new IllegalArgumentException("insert quality cannot less than 0: " + insQual);
        if (delQual < 0) throw new IllegalArgumentException("deletion quality cannot be less than 0: " + delQual);
        if (gcp < 0) throw new IllegalArgumentException("gcp cannot be less than 0: " + gcp);
        dest[matchToMatch] = matchToMatchProb(insQual, delQual);
        dest[matchToInsertion] = QualityUtils.qualToErrorProb(insQual);
        dest[matchToDeletion] = QualityUtils.qualToErrorProb(delQual);
        dest[indelToMatch] = QualityUtils.qualToProb(gcp);
        dest[insertionToInsertion] = dest[deletionToDeletion] = QualityUtils.qualToErrorProb(gcp);
    }

    /**
     * Returns a transition probability array given the different quality scores affecting a read site.
     *
     * @param insQual the insertion quality score as a byte.
     * @param delQual the deletion quality score as a byte.
     * @param gcp the gap-continuation-penalty score as a byte.
     *
     * @throws NullPointerException if {@code dest} is {@code null}.
     * @throws ArrayIndexOutOfBoundsException if {@code dest} is not large enough.
     * @throws IllegalArgumentException if {@code insQual}, {@code delQual} or {@code gcp} is less than negative.
     *
     * @return never {@code null}. An array of length {@link #TRANS_PROB_ARRAY_LENGTH}.
     */
    public static double[] qualToTransProbs(final byte insQual, final byte delQual, final byte gcp) {
        final double[] dest = new double[TRANS_PROB_ARRAY_LENGTH];
        qualToTransProbs(dest,insQual,delQual,gcp);
        return dest;
    }

    /**
     * Fills ax matrix with the transition probabilities for a number of bases.
     *
     * <p/>
     * The first dimension of the matrix correspond to the different bases where the first one is stored in position 1.
     * Thus the position 0 is left empty and the length of the resulting matrix is actually {@code insQual.length + 1}.
     * <p/>
     * Each entry is the transition probability array for that base with a length of {@link #TRANS_PROB_ARRAY_LENGTH}.
     *
     * @param dest the matrix to update
     * @param insQuals insertion qualities.
     * @param delQuals deletion qualities.
     * @param gcps gap-continuation penalty qualities.
     *
     * @throws NullPointerException if any of the input arrays, matrices is {@code null} or any entry in {@code dest} is {@code null}.
     * @throws IllegalArgumentException if {@code IllegalArgumentException}
     *  if the input array don't have the same length.
     * @throws ArrayIndexOutOfBoundsException if {@code dest} or any of its elements is not large enough to contain the
     *  transition  matrix.
     */
    public static void qualToTransProbs(final double[][] dest, final byte[] insQuals, final byte[] delQuals, final byte[] gcps) {
        Utils.nonNull(dest,     "dest array null");
        Utils.nonNull(insQuals, "insQuals array null");
        Utils.nonNull(delQuals, "delQuals array null");
        Utils.nonNull(gcps,     "gcps array null");

        final int readLength = insQuals.length;
        if (delQuals.length != readLength) throw new IllegalArgumentException("deletion quality array length does not match insert quality array length: " + readLength + " != " + delQuals.length);
        if (gcps.length != readLength) throw new IllegalArgumentException("deletion quality array length does not match insert quality array length: " + readLength + " != " + gcps.length);
        if (dest.length < readLength + 1) throw new IllegalArgumentException("destination length is not enough for the read length: " + dest.length + " < " + readLength + " + 1");

        for (int i = 0; i < readLength; i++) {
            qualToTransProbs(dest[i + 1], insQuals[i], delQuals[i], gcps[i]);
        }
    }

    /**
     * Returns a matrix with the transition probabilities for a number of bases.
     *
     * <p/>
     * The first dimension of the matrix correspond to the different bases where the first one is stored in position 1.
     * Thus the position 0 is left empty and the length of the resulting matrix is actually {@code insQual.length + 1}.
     * <p/>
     * Each entry is the transition probability array for that base with a length of {@link #TRANS_PROB_ARRAY_LENGTH}.
     *
     * @param insQuals insertion qualities.
     * @param delQuals deletion qualities.
     * @param gcps gap-continuation penalty qualities.
     *
     * @throws NullPointerException if any of the input arrays is {@code null}.
     * @throws IllegalArgumentException if {@code IllegalArgumentException}
     *  if the input array don't have the same length.
     *
     * @return never {@code null}, an matrix of the dimensions explained above.
     */
    public static double[][] qualToTransProbs(final byte[] insQuals, final byte[] delQuals, final byte[] gcps) {
        Utils.nonNull(insQuals, "insQuals array null");
        Utils.nonNull(delQuals, "delQuals array null");
        Utils.nonNull(gcps,     "gcps array null");

        final double[][] dest = createTransitionMatrix(insQuals.length);
        qualToTransProbs(dest,insQuals,delQuals,gcps);
        return dest;
    }

    /**
     * Fills a transition log10-probability array given the different quality scores affecting a read site.
     *
     * @param insQual the insertion quality score as a byte.
     * @param delQual the deletion quality score as a byte.
     * @param gcp the gap-continuation-penalty score as a byte.
     *
     * @throws NullPointerException if {@code dest} is {@code null}.
     * @throws ArrayIndexOutOfBoundsException if {@code dest} is not large enough.
     * @throws IllegalArgumentException if {@code insQual}, {@code delQual} or {@code gcp} is less than negative.
     */
    public static void qualToTransProbsLog10(final double[] dest, final byte insQual, final byte delQual, final byte gcp) {
        Utils.nonNull(dest, "dest array null");
        if (insQual < 0) throw new IllegalArgumentException("insert quality cannot less than 0: " + insQual);
        if (delQual < 0) throw new IllegalArgumentException("deletion quality cannot be less than 0: " + delQual);
        if (gcp < 0) throw new IllegalArgumentException("gcp cannot be less than 0: " + gcp);
        dest[matchToMatch] = matchToMatchProbLog10(insQual, delQual);
        dest[matchToInsertion] = QualityUtils.qualToErrorProbLog10(insQual);
        dest[matchToDeletion] = QualityUtils.qualToErrorProbLog10(delQual);
        dest[indelToMatch] = QualityUtils.qualToProbLog10(gcp);
        dest[insertionToInsertion] = QualityUtils.qualToErrorProbLog10(gcp);
        dest[deletionToDeletion] = QualityUtils.qualToErrorProbLog10(gcp);
    }

    /**
     * Returns a transition log10 probability array given the different quality scores affecting a read site.
     *
     * @param insQual the insertion quality score as a byte.
     * @param delQual the deletion quality score as a byte.
     * @param gcp the gap-continuation-penalty score as a byte.
     *
     * @throws NullPointerException if {@code dest} is {@code null}.
     * @throws ArrayIndexOutOfBoundsException if {@code dest} is not large enough.
     * @throws IllegalArgumentException if {@code insQual}, {@code delQual} or {@code gcp} is less than negative.
     *
     * @return never {@code null}. An array of length {@link #TRANS_PROB_ARRAY_LENGTH}.
     */
    public static double[] qualToTransProbsLog10(final byte insQual, final byte delQual, final byte gcp) {
        final double[] dest = new double[TRANS_PROB_ARRAY_LENGTH];
        qualToTransProbsLog10(dest, insQual, delQual, gcp);
        return dest;
    }

    /**
     * Fills a matrix with the log10 transition probabilities for a number of bases.
     *
     * <p/>
     * The first dimension of the matrix correspond to the different bases where the first one is stored in position 1.
     * Thus the position 0 is left empty and the length of the resulting matrix is actually {@code insQual.length + 1}.
     * <p/>
     * Each entry is the transition probability array for that base with a length of {@link #TRANS_PROB_ARRAY_LENGTH}.
     *
     * @param insQuals insertion qualities.
     * @param delQuals deletion qualities.
     * @param gcps gap-continuation penalty qualities.
     *
     * @throws NullPointerException if any of the input arrays, matrices is {@code null} or any entry in {@code dest} is {@code null}.
     * @throws IllegalArgumentException if {@code IllegalArgumentException}
     *  if the input array don't have the same length.
     * @throws ArrayIndexOutOfBoundsException if {@code dest} or any of its elements is not large enough to contain the
     *  transition  matrix.
     */
    public static void qualToTransProbsLog10(final double[][] dest, final byte[] insQuals, final byte[] delQuals, final byte[] gcps) {
        Utils.nonNull(dest,     "dest array null");
        Utils.nonNull(insQuals, "insQuals array null");
        Utils.nonNull(delQuals, "delQuals array null");
        Utils.nonNull(gcps,     "gcps array null");

        final int readLength = insQuals.length;
        if (delQuals.length != readLength) throw new IllegalArgumentException("deletion quality array length does not match insert quality array length: " + readLength + " != " + delQuals.length);
        if (gcps.length != readLength) throw new IllegalArgumentException("deletion quality array length does not match insert quality array length: " + readLength + " != " + gcps.length);
        if (dest.length < readLength + 1) throw new IllegalArgumentException("destination length is not enough for the read length: " + dest.length + " < " + readLength + " + 1");

        for (int i = 0; i < readLength; i++) {
            qualToTransProbsLog10(dest[i + 1], insQuals[i], delQuals[i], gcps[i]);
        }
    }

    /**
     * Returns a matrix with the log10 transition probabilities for a number of bases.
     *
     * <p/>
     * The first dimension of the matrix correspond to the different bases where the first one is stored in position 1.
     * Thus the position 0 is left empty and the length of the resulting matrix is actually {@code insQual.length + 1}.
     * <p/>
     * Each entry is the transition probability array for that base with a length of {@link #TRANS_PROB_ARRAY_LENGTH}.
     *
     * @param insQuals insertion qualities.
     * @param delQuals deletion qualities.
     * @param gcps gap-continuation penalty qualities.
     *
     * @throws NullPointerException if any of the input arrays is {@code null}.
     * @throws IllegalArgumentException if {@code IllegalArgumentException}
     *  if the input array don't have the same length.
     *
     * @return never {@code null}, an matrix of the dimensions explained above.
     */
    public static double[][] qualToTransProbsLog10(final byte[] insQuals, final byte[] delQuals, final byte[] gcps) {
        Utils.nonNull(insQuals, "insQuals array null");
        Utils.nonNull(delQuals, "delQuals array null");
        Utils.nonNull(gcps,     "gcps array null");

        final double[][] dest = createTransitionMatrix(insQuals.length);
        qualToTransProbsLog10(dest,insQuals,delQuals,gcps);
        return dest;
    }

    /**
     * Creates a transition probability matrix large enough to work with sequences of a particular length.
     *
     * @param maxReadLength the maximum read length for the transition matrix.
     *
     * @return never {@code null}. A matrix of {@code maxReadLength + 1} by {@link #TRANS_PROB_ARRAY_LENGTH} positions.
     */
    public static double[][] createTransitionMatrix(final int maxReadLength) {
        return new double[maxReadLength + 1][TRANS_PROB_ARRAY_LENGTH];
    }

    /**
     * Returns the probability that neither of two event takes place.
     * <p/>
     *
     * We assume that both event never occur together and that delQual is the conditional probability
     * (qual. encoded) of the second event, given the first event didn't took place. So that the
     * probability of no event is: <br/>
     *
     * We assume that both event never occur together so that the probability of no event is: <br/>
     *
     * <code>1 - ProbErr(insQual) - ProbErr(delQual)</code> <br/>
     *
     * @param insQual PhRED scaled quality/probability of the first event.
     * @param delQual PhRED scaled quality/probability of the second event.
     *
     * @return a value between 0 and 1.
     */
    public static double matchToMatchProb(final byte insQual, final byte delQual) {
        return matchToMatchProb((insQual & 0xFF), (delQual & 0xFF));
    }

    /**
     * Returns the log10 probability that neither of two events, insertion and deletion, takes place.
     * <p/>
     *
     * We assume that both event never occur together so that the probability of no event is: <br/>
     *
     * <code>1 - ProbErr(insQual) - ProbErr(delQual)</code> <br/>
     *
     * @param insQual PhRED scaled quality/probability of an insertion.
     * @param delQual PhRED scaled quality/probability of a deletion.
     *
     * @return a value between 0 and -Inf.
     */
    public static double matchToMatchProbLog10(final byte insQual, final byte delQual) {
        return matchToMatchProbLog10((insQual & 0xFF), (delQual & 0xFF));
    }

    /**
     * Returns the probability that neither of two events, insertion and deletion, takes place.
     * <p/>
     *
     * We assume that both event never occur together and that delQual is the conditional probability
     * (qual. encoded) of the second event, given the first event didn't took place. So that the
     * probability of no event is: <br/>
     *
     * <code>1 - ProbErr(insQual) - ProbErr(delQual)</code> <br/>
     *
     * @param insQual PhRED scaled quality/probability of an insertion.
     * @param delQual PhRED scaled quality/probability of a deletion.
     * @return a value between 0 and 1.
     */
    public static double matchToMatchProb(final int insQual, final int delQual) {
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
     * Returns the log-probability that neither of two events takes place.
     * <p/>
     *
     * We assume that both events never occur together and that delQual is the conditional probability (qual. encoded)
     * of the second event, given the first event didn't took place. So that the probability of no event is: <br/>
     *
     * <code>1 - ProbErr(insQual) - ProbErr(delQual)</code> <br/>
     *
     * @param insQual PhRED scaled quality/probability of an insertion.
     * @param delQual PhRED scaled quality/probability of a deletion.
     *
     * @return a value between 0 and -Inf.
     */
    public static double matchToMatchProbLog10(final int insQual, final int delQual) {
        final int minQual;
        final int maxQual;
        if (insQual <= delQual) {
            minQual = insQual;
            maxQual = delQual;
        } else {
            minQual = delQual;
            maxQual = insQual;
        }
        return (QualityUtils.MAX_QUAL < maxQual) ? Math.log1p (
                -Math.min(1, Math.pow(10,
                        MathUtils.approximateLog10SumLog10(-.1 * minQual, -.1 * maxQual)))) * INV_LN10 :
                matchToMatchLog10[((maxQual * (maxQual + 1)) >> 1) + minQual];
    }
}
