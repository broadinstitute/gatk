package org.broadinstitute.hellbender.tools.exome.allelicbalancecaller;

import com.google.common.annotations.VisibleForTesting;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.BaseAbstractUnivariateIntegrator;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.SearchInterval;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;
import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.apache.commons.math3.util.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.tools.exome.ACNVModeledSegment;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.variant.HomoSapiensConstants;

import java.io.Serializable;
import java.util.Arrays;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;


/**
 * <p>This class makes calls for each segment on whether the segment is balanced (MAF=0.5) or CNLoH.  Please see the
 *  docs for detailed information.
 * </p>
 *
 * <p>Important note:  The CNLoH calls are not evaluated and preliminary performance indicated low sensitivity and low
 * precision.</p>
 *
 * <p>Throughout this class, rho is the CCF * purity.  M is the copy number of one allele and N is the copy number of the other allele.
 * There is no assumption about the relative relationship of M and N, though typically, we try to maintain M >= N,
 *  though this is not guaranteed.
 *</p>
 *
 *  <p>Many of the variable names assume that you have read the CNV methods pdf in the docs directory.</p>
 *
 * Variables:
 * <ul>
 *     <li>rho -- CCF * purity.  We use a set number of rhos that each segment can have.</li>
 *     <li>M -- absolute copy number for allele 1</li>
 *     <li>N -- absolute copy number for allele 2.  NOT assumed to have same or fewer copies than allele 1</li>
 *     <li>pi -- discrete pdf for values of M and N</li>
 *     <li>phi -- discrete pdf for each rho value</li>
 *     <li>maf -- minor allelic fraction</li>
 *     <li>fMaf -- likelihood of a minor allelic fraction, given other parameters</li>
 *     <li>cr -- copy ratio</li>
 *     <li>fCr -- likelihood of a copy ratio, given other parameters</li>
 *     <li>alpha -- concentration parameter of dirichlet prior on z_s</li>
 *     <li>z_s -- indicator variable (binary) whether segment s has a CNV on segment s.  This is not explicitly used in the code.</li>
 * </ul>
 *
 *  NOTES:
 *  <ul>
 *      <li>{@code lambda} (CCF * ploidy) is present, but the value is never changed.  In future versions, we may add this to the model.</li>
 *  </ul>
 */
public class CNLOHCaller implements Serializable {

    private static final Logger logger = LogManager.getLogger(CNLOHCaller.class);
    protected static double MIN_L = 1e-50;
    protected static double MIN_DIST_FROM_MODE = 1e-10;
    static final long serialVersionUID = 337337337123L;

    private double segmentMeanBiasInCR;
    private double segmentMeanVarianceInCR;

    private double rhoThreshold;

    public static double RHO_THRESHOLD_DEFAULT = 0.2;

    private static double NUM_STD_95_CONF_INTERVAL_GAUSSIAN = 1.96;

    private static double CLOSE_ENOUGH_TO_COPY_NEUTRAL_IN_CR = 0.1;

    /** Minimum value to integrate over in order to calculate E_alpha*/
    private static double MIN_E_ELPHA_INTEGRATION_RANGE = 0.0;

    /** Maximum value to integrate over in order to calculate E_alpha*/
    private static double MAX_E_ELPHA_INTEGRATION_RANGE = 10.0;

    /** How many copies expected in a normal.  For humans, this is 2 */
    private int normalNumCopies;

    /** Number of discrete rho values to use. A.k.a "K"*/
    protected static int NUM_RHOS = 25;

    public CNLOHCaller(){
        rhoThreshold = RHO_THRESHOLD_DEFAULT;
        normalNumCopies = HomoSapiensConstants.DEFAULT_PLOIDY;
    }

    public int getNormalNumCopies() {
        return normalNumCopies;
    }

    public void setNormalNumCopies(final int normalNumCopies) {
        this.normalNumCopies = normalNumCopies;
    }

    public double getRhoThreshold() {
        return rhoThreshold;
    }

    public void setRhoThreshold(final double rhoThreshold) {
        this.rhoThreshold = rhoThreshold;
    }

    /**
     * Calculate the likelihood of a MAF.
     *
     * This is done by assuming with two half-Gaussians each modeling below the mode and above the mode.  While this is not a
     *  particularly good form for a distribution, it is quick to calculate and generates reasonable results.  This is
     *  also done because we do not know the exact form of the distribution when coming from ACNV.
     *
     */
    @VisibleForTesting
    static double calculateFmaf(final double rho, final int m, final int n, final double credibleMode,
                                   final double credibleLow, final double credibleHigh, final double normalPloidy) {
        ParamUtils.inRange(rho, 0.0, 1.0, "Invalid rho value: " + rho + ".  Must be [0,1]");
        ParamUtils.isPositiveOrZero(m, "M must be positive.");
        ParamUtils.isPositiveOrZero(n, "N must be positive.");

        if (Double.isNaN(credibleMode) || Double.isNaN(credibleLow) || Double.isNaN(credibleHigh)){
            return 1;
        }

        final double maf = calculateMaf(rho, m, n, normalPloidy);
        return Math.max(calculateDoubleGaussian(maf, credibleMode, credibleLow, credibleHigh), MIN_L);

    }

    /**
     * Calculates the minor allelic fraction for given rho, M, and N.
     *
     *  Output is truncated at MIN_L
     *
     *  @param rho CCF*purity
     * @param m absolute copy number of allele1
     * @param n absolute copy number of allele2
     * @return calculated minor allele fraction
     */
    @VisibleForTesting
    static double calculateMaf(final double rho, final int m, final int n, final double normalPloidy) {
        ParamUtils.inRange(rho, 0.0, 1.0, "Invalid rho value: " + rho + ".  Must be [0,1]");
        ParamUtils.isPositiveOrZero(m, "M must be >= 0");
        ParamUtils.isPositiveOrZero(n, "N must be >= 0");

        double result = ((1-rho) + rho * (Math.min(n,m))) / ( (1-rho)*normalPloidy + rho * (n+m) );

        // If result is NaN, that means that we had rho = 1 and M,N = 0
        return (result < MIN_L) || (Double.isNaN(result)) ? MIN_L : result;
    }

    /**
     * This function takes two half-Gaussian distributions and uses one for the values below the mode and the other for values above the
     *  mode.
     *
     *  If the mode is outside the credible interval, the calculation is done as if the interval boundary is
     *   {@link CNLOHCaller#MIN_DIST_FROM_MODE} above/below the mode.
     *
     *
     * @param val value to calculate the "pdf"
     * @param credibleMode mode of the presumed distribution
     * @param credibleLow 95% confidence interval on the low tail
     * @param credibleHigh 95% confidence interval on the high tail
     * @return pdf with a minimum value of {@link CNLOHCaller#MIN_L}
     */
    @VisibleForTesting
    static double calculateDoubleGaussian(final double val, final double credibleMode, final double credibleLow,
                                             final double credibleHigh) {

        final double hiDist = Math.max(credibleHigh - credibleMode, MIN_DIST_FROM_MODE);
        final double loDist = Math.max(credibleMode - credibleLow, MIN_DIST_FROM_MODE);
        return new NormalDistribution(null, credibleMode, ( val >= credibleMode ? hiDist : loDist)/NUM_STD_95_CONF_INTERVAL_GAUSSIAN).density(val);
    }

    /**
     * Calculate a copy ratio.
     *
     * @param rho CCF * purity
     * @param m Absolute copy number of allele 1
     * @param n Absolute copy number of allele 2
     * @param lambda CCF * ploidy
     * @param normalPloidy ploidy for the normal cells.  For humans, this is 2.
     * @return copy ratio in CR space
     */
    @VisibleForTesting
    static double calculateCopyRatio(final double rho, final int m, final int n, final double lambda, final double normalPloidy) {
        return ((1-rho) * normalPloidy + rho*(n+m)) / lambda;
    }

    /**
     * Calculate the likelihood for the copy ratio.  Models the credible interval as a Gaussian.
     *
     * @param rho -- see class docs
     * @param m -- see class docs
     * @param n -- see class docs
     * @param lambda -- see class docs
     * @param credibleMode mode of the credible interval from for the copy ratio.
     * @param credibleLow  low point of the credible interval from for the copy ratio.
     * @param credibleHigh  high point of the credible interval from for the copy ratio.
     * @param additionalVariance Additional variance that should be included in the segment mean pdf
     * @param normalPloidy ploidy for the normal cells.  For humans, this is 2.
     * @return likelihood, with a minimum value of {@link CNLOHCaller#MIN_L}
     */
    @VisibleForTesting
    static double calculateFcr(final double rho, final int m, final int n, final double lambda, final double credibleMode,
                         final double credibleLow, final double credibleHigh, final double additionalVariance, final double normalPloidy) {
        ParamUtils.inRange(rho, 0.0, 1.0, "Invalid rho value: " + rho + ".  Must be [0,1]");
        ParamUtils.isPositiveOrZero(m, "M must be positive.");
        ParamUtils.isPositiveOrZero(n, "N must be positive.");
        ParamUtils.isPositiveOrZero(additionalVariance, "additional variance must be positive.");
        if (Double.isNaN(credibleMode) || Double.isNaN(credibleLow) || Double.isNaN(credibleHigh)){
            return 1;
        }
        final double cr = calculateCopyRatio(rho, m, n, lambda, normalPloidy);
        final double avgDist = ((credibleMode - credibleLow) + (credibleHigh - credibleMode))/2.0;
        return Math.max(new NormalDistribution(null, credibleMode,
                (avgDist / NUM_STD_95_CONF_INTERVAL_GAUSSIAN) + Math.sqrt(additionalVariance)).density(cr), MIN_L);
    }

    /**
     * Attempt to get an idea of segment mean variance near copy neutral.
     *
     * @param segments Never {@code null}
     * @return variance of segment mean (in CR space) of segments that are "close enough" to copy neutral.
     *   Zero if no segments are "close enough"
     */
    private double calculateVarianceOfCopyNeutralSegmentMeans(final List<ACNVModeledSegment> segments, final double meanBiasInCR) {
        Utils.nonNull(segments);

        // Only consider values "close enough" to copy neutral (CR == 1).
        final double neutralCR = 1 + meanBiasInCR;
        final double[] neutralSegmentMeans = segments.stream()
                .mapToDouble(ACNVModeledSegment::getSegmentMeanInCRSpace)
                .filter(m -> Math.abs(m - neutralCR) < CLOSE_ENOUGH_TO_COPY_NEUTRAL_IN_CR).toArray();
        return new Variance().evaluate(neutralSegmentMeans);
    }

    private double calculateSegmentMeanBiasInCRSpace(final List<ACNVModeledSegment> segments) {
        Utils.nonNull(segments);

        final double neutralCRApprox = 1;

        // Only consider values "close enough" to copy neutral (CR == 1).
        final double[] neutralSegmentMeans = segments.stream().mapToDouble(ACNVModeledSegment::getSegmentMeanInCRSpace)
                .filter(x -> Math.abs(x - neutralCRApprox) < CLOSE_ENOUGH_TO_COPY_NEUTRAL_IN_CR)
                .toArray();

        return new Percentile().evaluate(neutralSegmentMeans) - 1;
    }

    /**
     * Make both CNLoH and balanced segment calls.
     *
     * This method is not thread-safe, since it initializes member attributes.
     *
     * @param segments The segments with segment mean and minor allelic fraction estimates (typically from ACNV).
     * @param numIterations Number of EM iteration to run
     * @param ctx This call requires Spark, so not {@code null}
     * @return Never {@code null}
     */
    public List<AllelicSplitCall> makeCalls(final List<ACNVModeledSegment> segments, final int numIterations, final JavaSparkContext ctx) {
        ParamUtils.isPositive(numIterations, "Must be more than zero iterations.");
        Utils.nonNull(segments);
        Utils.nonNull(ctx, "Java SparkContext can't be null when attempting to make CNLOH and balance calls.");

        segmentMeanBiasInCR = calculateSegmentMeanBiasInCRSpace(segments);
        segmentMeanVarianceInCR = calculateVarianceOfCopyNeutralSegmentMeans(segments, segmentMeanBiasInCR);

        final AllelicBalanceCallerModelState state = AllelicBalanceCallerModelState.createInitialCNLOHCallerModelState(rhoThreshold, segments,
                normalNumCopies, NUM_RHOS);

        // Create a Spark RDD for the segments.
        final JavaRDD<ACNVModeledSegment> segs = ctx.parallelize(segments);

        List<AllelicSplitCall> allelicSplitCalls = null;

        // for each iteration
        for (int i = 1; i <= numIterations; i++) {

            logger.info("Starting iteration " + i + " of " + numIterations + " ... ");
            logger.info("Calculating responsibilities (E_zsk_vsm_wsn) for each segment ... ");

            // We are not broadcasting state, since it will only be used once before values change.  In that case,
            //  there is no gain for broadcasting

            // Effectively, the responsibilities are a 4D array:  S x K x M x N.  We do the responsibility calculation
            //  over the segments.  So the JavaRDD will be of length S.
            // Update E_zsk_vsm_wsn (responsibilities) for each segment (Steps 1 & 2)  {K x M x N} in S entries.  The RDD is of length S.
            final JavaRDD<double[][][]> responsibilitiesForSegs = segs.map(s -> calculateResponsibilities(state.getEffectivePhis(), state.getEffectivePis(), state.getRhos(),
                    s.getMinorAlleleFractionPosteriorSummary().getCenter(),
                    s.getMinorAlleleFractionPosteriorSummary().getLower(),
                    s.getMinorAlleleFractionPosteriorSummary().getUpper(),
                    Math.pow(2, s.getSegmentMeanPosteriorSummary().getCenter()) - segmentMeanBiasInCR,
                    Math.pow(2, s.getSegmentMeanPosteriorSummary().getLower()) - segmentMeanBiasInCR,
                    Math.pow(2, s.getSegmentMeanPosteriorSummary().getUpper()) - segmentMeanBiasInCR,
                    state.getLambda(), state.getmVals(), state.getnVals()));

            // Create an array (for each segment) that is a sum of responsibilities for each rho.  (length:K for S entries)
            final JavaRDD<double[]> responsibilitiesByRhoForSegs = responsibilitiesForSegs.map(da -> sumOverSecondAndThirdDimension(da));
            final JavaRDD<double[]> responsibilitiesByAllele1ForSegs = responsibilitiesForSegs.map(da -> sumOverFirstAndThirdDimension(da));
            final JavaRDD<double[]> responsibilitiesByAllele2ForSegs = responsibilitiesForSegs.map(da -> sumOverFirstAndSecondDimension(da));

            // Create arrays summed in one dimension (incl. over segs).
            final double[] responsibilitiesByRho = sumOverSegments(responsibilitiesByRhoForSegs);
            final double[] responsibilitiesByAllele1 = sumOverSegments(responsibilitiesByAllele1ForSegs);
            final double[] responsibilitiesByAllele2 = sumOverSegments(responsibilitiesByAllele2ForSegs);

            final List<double[][][]> responsibilitiesForSegsAsList = responsibilitiesForSegs.collect();

            // Update EffectivePhis (Step 3)
            logger.info("Updating effective phis (iteration " + i + " of " + numIterations +  ") ... ");
            state.setEffectivePhis(calcEffectivePhis(state.getEffectiveAlpha(), responsibilitiesByRho));

            // Update pis (Step 4)
            logger.info("Updating pis (iteration " + i + " of " + numIterations + ") ... ");
            state.setEffectivePis(calcPis(responsibilitiesByAllele1, responsibilitiesByAllele2));

            // Update E_alpha (Step 5)
            logger.info("Updating E_alpha (iteration " + i + " of " + numIterations + ") ... ");
            state.setEffectiveAlpha(calcEffectiveAlpha(state.getEffectivePhis()));

            // Update rhos (Step 6)
            logger.info("Updating rhos (iteration " + i + " of " + numIterations +  ") ... ");
            state.setRhos(calcNewRhos(segments, responsibilitiesForSegsAsList, state.getLambda(),
                    state.getRhos(), state.getmVals(), state.getnVals(), ctx));

            // This only really needs to  be done once, but is necessary here due to scoping issues that will cause
            //  compiler issues.
            logger.info("Updating calls (iteration " + i + " of " + numIterations +  ") ... ");
            allelicSplitCalls = createCNLoHCalls(segments, state, responsibilitiesForSegsAsList);
        }
        logger.info("ACNV segments bias, variance (around copy neutral only): " + segmentMeanBiasInCR + ", " + segmentMeanVarianceInCR);
        logger.info("Output file is not adjusted for the bias and variance.  No user action required.");
        return allelicSplitCalls;
    }

    private List<AllelicSplitCall> createCNLoHCalls(final List<ACNVModeledSegment> segments, final AllelicBalanceCallerModelState state, final List<double[][][]> responsibilitiesForSegs) {

        final List<AllelicSplitCall> allelicSplitCalls = segments.stream().map(AllelicSplitCall::new).collect(Collectors.toList());

        final List<Pair<Integer, Integer>> maxMIdxNIdxPerSeg = IntStream.range(0, segments.size()).boxed()
                .map(s -> max2dIndices(sumOverFirstDimension(responsibilitiesForSegs.get(s))))
                .collect(Collectors.toList());

        final List<double[]> rhosForBestMN = IntStream.range(0, segments.size()).boxed()
                .map(s -> getFirstDimensionArray(responsibilitiesForSegs.get(s), maxMIdxNIdxPerSeg.get(s).getFirst(), maxMIdxNIdxPerSeg.get(s).getSecond()))
                .collect(Collectors.toList());

        final double[] calledRhoPerSeg = rhosForBestMN.stream().map(MathUtils::maxElementIndex).mapToDouble(ri -> state.getRhos()[ri]).toArray();
        final int[] calledMPerSeg = maxMIdxNIdxPerSeg.stream().mapToInt(p -> state.getmVals()[p.getFirst()]).toArray();
        final int[] calledNPerSeg = maxMIdxNIdxPerSeg.stream().mapToInt(p -> state.getnVals()[p.getSecond()]).toArray();

        // Set the rho, M, N, FCr, and FMaf for the call
        for (int s = 0; s < allelicSplitCalls.size(); s ++) {
            final AllelicSplitCall allelicSplitCallForSeg = allelicSplitCalls.get(s);
            allelicSplitCallForSeg.setRho(calledRhoPerSeg[s]);
            allelicSplitCallForSeg.setM(calledMPerSeg[s]);
            allelicSplitCallForSeg.setN(calledNPerSeg[s]);
            allelicSplitCallForSeg.setfCr(
                    calculateFcr(calledRhoPerSeg[s], calledMPerSeg[s], calledNPerSeg[s], state.getLambda(),
                            Math.pow(2, allelicSplitCallForSeg.getAcnvSegment().getSegmentMeanPosteriorSummary().getCenter()),
                            Math.pow(2, allelicSplitCallForSeg.getAcnvSegment().getSegmentMeanPosteriorSummary().getLower()),
                            Math.pow(2, allelicSplitCallForSeg.getAcnvSegment().getSegmentMeanPosteriorSummary().getUpper()),
                            segmentMeanVarianceInCR, normalNumCopies
                    ));
            allelicSplitCallForSeg.setfMaf(
                    calculateFmaf(calledRhoPerSeg[s], calledMPerSeg[s], calledNPerSeg[s],
                            allelicSplitCallForSeg.getAcnvSegment().getMinorAlleleFractionPosteriorSummary().getCenter(),
                            allelicSplitCallForSeg.getAcnvSegment().getMinorAlleleFractionPosteriorSummary().getLower(),
                            allelicSplitCallForSeg.getAcnvSegment().getMinorAlleleFractionPosteriorSummary().getUpper(),
                            normalNumCopies)
            );

            // Set the balanced call
            allelicSplitCallForSeg.setBalancedCall(isBalanced(calledRhoPerSeg[s], calledMPerSeg[s], calledNPerSeg[s]));

            // Set the CNLoH call
            allelicSplitCallForSeg.setCnlohCall(isCNLoH(calledRhoPerSeg[s], calledMPerSeg[s], calledNPerSeg[s], allelicSplitCallForSeg.getAcnvSegment()));
        }

        return allelicSplitCalls;
    }

    private AllelicBalanceCall isBalanced(final double rho, final int m, final int n) {
        return ( (m == n) || (rho <= getRhoThreshold()) ? AllelicBalanceCall.BALANCED : AllelicBalanceCall.NOT_BALANCED);
    }

    private CNLoHCall isCNLoH(final double rho, final int m, final int n, final ACNVModeledSegment segment) {
        if (Double.isNaN(segment.getMinorAlleleFractionPosteriorSummary().getCenter())) {
            return CNLoHCall.NO_CALL;
        } else if (((m == getNormalNumCopies()) && (n == 0)) || ((n == getNormalNumCopies()) && (m == 0))
            && (rho > getRhoThreshold())) {
            return CNLoHCall.CNLOH;
        } else {
            return CNLoHCall.NOT_CNLOH;
        }
    }

    private double[] sumOverSegments(final JavaRDD<double[]> eZskBySeg) {
        final double[][] eZskBySeg2D = eZskBySeg.collect().stream().toArray(double[][]::new); // S x K
        final RealMatrix eZskBySegMatrix = MatrixUtils.createRealMatrix(eZskBySeg2D);
        return GATKProtectedMathUtils.columnSums(eZskBySegMatrix);
    }

    private double[] sumOverSecondAndThirdDimension(final double[][][] data) {
        final double [] result = new double[data.length];
        for (int i = 0; i < data.length; i++) {
            result[i] = 0;
            for (int j = 0; j < data[0].length; j++) {
                for (int k = 0; k < data[0][0].length; k++) {
                    result[i] += data[i][j][k];
                }
            }
        }

        return result;
    }

    @VisibleForTesting
    static Pair<Integer, Integer> max2dIndices(final double[][] data) {
        double maxValue = Double.NEGATIVE_INFINITY;
        Pair<Integer, Integer> result = null;
        for (int i=0; i < data.length; i++) {
            for (int j=0; j < data[0].length; j++) {
                if (data[i][j] > maxValue) {
                    maxValue = data[i][j];
                    result = Pair.create(i, j);
                }
            }
        }
        return result;
    }

    private double[] getFirstDimensionArray(final double[][][] array, final int dim2Index, final int dim3Index) {
        return Arrays.stream(array).mapToDouble(array2D -> array2D[dim2Index][dim3Index]).toArray();
    }

    @VisibleForTesting
    static double[][] sumOverFirstDimension(final double[][][] data) {
        final double [][] result = new double[data[0].length][data[0][0].length];

        for (int j = 0; j < data[0].length; j++) {
            for (int k = 0; k < data[0][0].length; k++) {
                result[j][k] = 0;
                for (int i = 0; i < data.length; i++) {
                    result[j][k] += data[i][j][k];
                }
            }
        }

        return result;
    }

    private double[] sumOverFirstAndThirdDimension(final double[][][] data) {
        final double [] result = new double[data[0].length];
        for (int j = 0; j < data[0].length; j++) {
            result[j] = 0;
            for (int i = 0; i < data.length; i++) {
                for (int k = 0; k < data[0][0].length; k++) {
                    result[j] += data[i][j][k];
                }
            }
        }

        return result;
    }

    private double[] sumOverFirstAndSecondDimension(final double[][][] data) {
        final double [] result = new double[data[0][0].length];
        for (int k = 0; k < data[0][0].length; k++) {
            result[k] = 0;
            for (int i = 0; i < data.length; i++) {
                for (int j = 0; j < data[0].length; j++) {
                    result[k] += data[i][j][k];
                }
            }
        }

        return result;
    }

    // Assumes that the pis parameter is a SortedMap.
    protected double[][][] calculateResponsibilities(final double[] effectivePhis, final double[] pis,
                                                     final double[] rhos, final double mafCredibleMode,
                                                     final double mafCredibleLow, final double mafCredibleHigh,
                                                     final double crCredibleMode, final double crCredibleLow,
                                                     final double crCredibleHigh, final double lambda, final int[] mVals,
                                                     final int[] nVals) {

        final double[][][] result = new double[rhos.length][mVals.length][nVals.length];
        double total = 0.0;
        for (int rhoi = 0; rhoi < rhos.length; rhoi++) {
            for (int mi = 0; mi < mVals.length; mi++) {
                for (int ni = 0; ni < nVals.length; ni++) {

                    final double rho = rhos[rhoi];
                    final int m = mVals[mi];
                    final int n = nVals[ni];

                    // Do not allow rho to be less than zero or greater than one, since it is a proportion.
                    if ((rho < 0) || (rho > 1)) {
                        result[rhoi][mi][ni] = 0;

                        // Do not allow rho to be zero and M = N = 1, since this would imply a rho of zero on a variant.
                    } else if ((rho == 0) && ((m != 1) || (n != 1))) {
                        result[rhoi][mi][ni] = 0;

                        // Do not allow rho to take on a value less than the threshold, unless it is zero
                    } else if ((rho > 0) && (rho < rhoThreshold)) {
                        result[rhoi][mi][ni] = 0;

                        // Do not allow rho to be greater than zero and M = N = 1, since this would imply a rho of non-zero
                        //  on a normal ploidy.
                    } else if ((rho != 0) && ((m == 1) && (n == 1))) {
                        result[rhoi][mi][ni] = 0;
                    } else {
                        result[rhoi][mi][ni] = effectivePhis[rhoi] *
                                pis[mi] * pis[ni] *
                                calculateFmaf(rho, m, n, mafCredibleMode, mafCredibleLow, mafCredibleHigh, normalNumCopies) *
                                calculateFcr(rho, m, n, lambda, crCredibleMode, crCredibleLow, crCredibleHigh,
                                        segmentMeanVarianceInCR, normalNumCopies);
                    }

                    total += result[rhoi][mi][ni];
                }
            }
        }

        // Normalize by total
        for (int rhoi = 0; rhoi < rhos.length; rhoi++) {
            for (int mi = 0; mi < mVals.length; mi++) {
                for (int ni = 0; ni < nVals.length; ni++) {
                    result[rhoi][mi][ni] /= total;
                }
            }
        }

        return result;
    }

    protected double[] calcEffectivePhis(final double E_alpha, final double[] responsibilitiesByRho) {

        final double sumResponsibilities = MathUtils.sum(responsibilitiesByRho);

        final double[] result = new double[responsibilitiesByRho.length];
        final int k = responsibilitiesByRho.length;

        // Initialize all pseudocounts to 1, except for index 0, which is 20;
        //  This artificially increases the odds for a rho = 0.
        final RealVector pseudocounts = new ArrayRealVector(responsibilitiesByRho.length);
        pseudocounts.set(1.0);
        pseudocounts.setEntry(0, 20.0);

        final double sumPseudocounts = MathUtils.sum(pseudocounts.toArray());
        final double term2 = Gamma.digamma(E_alpha + sumPseudocounts + sumResponsibilities);
        for (int i=0; i < result.length; i++) {
            final double term1 = Gamma.digamma(E_alpha/k + pseudocounts.getEntry(i) + responsibilitiesByRho[i]);
            result[i] = Math.exp(term1 - term2);
        }

        return result;
    }

    // Sum (and normalize) the responsibilities for each possible value of the absolute copy numbers.
    protected double[] calcPis(final double[] responsibilitiesByAllele1, final double[] responsibilitiesByAllele2) {

        final double total = MathUtils.sum(responsibilitiesByAllele1) + MathUtils.sum(responsibilitiesByAllele2);

        final double[] result = new double[responsibilitiesByAllele1.length];
        for (int i = 0; i < result.length; i++) {
            result[i] = (responsibilitiesByAllele1[i] + responsibilitiesByAllele2[i]) / total;
        }
        return result;
    }

    protected double calcEffectiveAlpha(final double[] effectivePhis) {

        final QAlphaUnivariateFunction qAlpha = new QAlphaUnivariateFunction(effectivePhis);
        final AQAlphaUnivariateFunction aQAlpha = new AQAlphaUnivariateFunction(effectivePhis);
        final BaseAbstractUnivariateIntegrator integrator = new SimpsonIntegrator();

        final double numerator = integrator.integrate(10, aQAlpha, MIN_E_ELPHA_INTEGRATION_RANGE, MAX_E_ELPHA_INTEGRATION_RANGE);
        final double denominator = integrator.integrate(10, qAlpha, MIN_E_ELPHA_INTEGRATION_RANGE, MAX_E_ELPHA_INTEGRATION_RANGE);

        return numerator/denominator;
    }

    private double[] calcNewRhos(final List<ACNVModeledSegment> segments,
                                 final List<double[][][]> responsibilitiesBySeg,
                                 final double lambda, final double[] rhos, final int[] mVals, final int[] nVals,
                                 final JavaSparkContext ctx) {

        // Since, we pass in the entire responsibilities matrix, we need the correct index for each rho.  That, and the
        //  fact that this is a univariate objective function, means we need to create an instance for each rho.  And
        //  then we blast across Spark.

        final List<Pair<? extends Function<Double, Double>, SearchInterval>> objectives = IntStream.range(0, rhos.length)
                .mapToObj(i -> new Pair<>(
                        new Function<Double, Double>() {
                            @Override
                            public Double apply(Double rho) {
                                return calculateESmnObjective(rho, segments, responsibilitiesBySeg, mVals, nVals, lambda, i);
                            }
                        },
                        new SearchInterval(0.0, 1.0, rhos[i])))
                .collect(Collectors.toList());

        final JavaRDD<Pair<? extends Function<Double, Double>, SearchInterval>> objectivesRDD = ctx.parallelize(objectives);

        final List<Double> resultsAsDouble = objectivesRDD
                .map(objective -> optimizeIt(objective.getFirst(), objective.getSecond()))
                .collect();

        return resultsAsDouble.stream().mapToDouble(Double::doubleValue).toArray();
    }

    // Returns a single number that is the total likelihood.  Note that it does this for the given rho, even though
    //  the entire responsibility 4d array (list of 3d, where list is per segment, then each array is KxMxN)
    //  is passed in.
    //
    //  HACK: All rhos have to be passed in with the index of interest.  Under the hood, all other values of rho are ignored.
    private double calculateESmnObjective(final double rho, List<ACNVModeledSegment> segments,
                                          final List<double[][][]> responsibilitiesForSegsAsList,
                                          final int[] mVals, final int[] nVals, final double lambda, final int rhoIndex) {

        // We will want to sum an entire matrix that is S x M x N for the given rho.
        final double[][][] eSMN = new double[responsibilitiesForSegsAsList.size()][mVals.length][nVals.length];

        // Populate eSMN
        for (int s=0; s<responsibilitiesForSegsAsList.size(); s++) {
            final ACNVModeledSegment seg = segments.get(s);
            final double mafMode = seg.getMinorAlleleFractionPosteriorSummary().getCenter();
            final double mafLow = seg.getMinorAlleleFractionPosteriorSummary().getLower();
            final double mafHigh = seg.getMinorAlleleFractionPosteriorSummary().getUpper();
            final double crMode = Math.pow(2, seg.getSegmentMeanPosteriorSummary().getCenter()) - segmentMeanBiasInCR;
            final double crLow = Math.pow(2, seg.getSegmentMeanPosteriorSummary().getLower()) - segmentMeanBiasInCR;
            final double crHigh = Math.pow(2, seg.getSegmentMeanPosteriorSummary().getUpper()) - segmentMeanBiasInCR;

            for (int m=0; m<mVals.length; m++) {
                for (int n=0; n<nVals.length; n++) {
                    final double mafLikelihood = calculateFmaf(rho, mVals[m], nVals[n], mafMode, mafLow, mafHigh, normalNumCopies);
                    final double crLikelihood = calculateFcr(rho, mVals[m], nVals[n], lambda, crMode, crLow, crHigh,
                            segmentMeanVarianceInCR, normalNumCopies);
                    if ( ((rho > 1) || (rho < 0)) || ((rho > 0) && (rho < rhoThreshold)))  {
                        eSMN[s][m][n] = MIN_L;
                    } else {
                        eSMN[s][m][n] = responsibilitiesForSegsAsList.get(s)[rhoIndex][m][n] * Math.log(mafLikelihood * crLikelihood);
                    }
                }
            }
        }

        return GATKProtectedMathUtils.sum(eSMN);
    }

    private class QAlphaUnivariateFunction implements UnivariateFunction {
        private static final double GAMMA_SHAPE = 1.0;
        private static final double GAMMA_SCALE = 1.0;

        private int numEffectivePhis;
        private double sumLogEffectivePhis;
        private GammaDistribution gammaDistribution;

        public QAlphaUnivariateFunction(final double[] effectivePhis) {
            this.numEffectivePhis = effectivePhis.length;
            this.sumLogEffectivePhis = Arrays.stream(effectivePhis).map(Math::log).sum();
            this.gammaDistribution = new GammaDistribution(null, GAMMA_SHAPE, GAMMA_SCALE);
        }

        @Override
        public double value(double a) {

            // Cannot evaluate a gamma at zero, so set a to a minimum value..
            final double aPrime = Math.max(a, MIN_L);
            final double pAlpha = gammaDistribution.density(aPrime);
            final double gammaAlpha = Gamma.gamma(aPrime);

            final double numerator = pAlpha * gammaAlpha;
            final double denominator = Math.pow(Gamma.gamma(aPrime/numEffectivePhis), numEffectivePhis);
            final double coeff = numerator/denominator;

            return coeff * Math.exp((aPrime/numEffectivePhis) * sumLogEffectivePhis);
        }
    }

    private class AQAlphaUnivariateFunction implements UnivariateFunction {
        private QAlphaUnivariateFunction qAlpha;

        public AQAlphaUnivariateFunction(final double[] effectivePhis) {
            qAlpha = new QAlphaUnivariateFunction(effectivePhis);
        }

        @Override
        public double value(double a) {
            return a * qAlpha.value(a);
        }
    }

    private double optimizeIt(final Function<Double, Double> objectiveFxn, final SearchInterval searchInterval) {
        final MaxEval BRENT_MAX_EVAL = new MaxEval(1000);
        final double RELATIVE_TOLERANCE = 0.001;
        final double ABSOLUTE_TOLERANCE = 0.001;
        final BrentOptimizer OPTIMIZER = new BrentOptimizer(RELATIVE_TOLERANCE, ABSOLUTE_TOLERANCE);

        final UnivariateObjectiveFunction objective = new UnivariateObjectiveFunction(x -> objectiveFxn.apply(x));
        return OPTIMIZER.optimize(objective, GoalType.MAXIMIZE, searchInterval, BRENT_MAX_EVAL).getPoint();
    }
}
