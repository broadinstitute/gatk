package org.broadinstitute.hellbender.tools.exome.cnlohcaller;

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
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.serializer.KryoSerializer;
import org.apache.spark.serializer.SerializerInstance;
import org.broadinstitute.hellbender.tools.exome.ACNVModeledSegment;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import scala.reflect.ClassTag;
import scala.reflect.ClassTag$;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.SortedMap;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;


/**
 * <p>This class makes calls for each segment on whether the segment is CNLoH and/or balanced (MAF=0.5).  Please see the
 *  docs for detailed information.
 * </p>
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
 *     <li>phi -- likelihood of each rho value</li>
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

    /** How many copies expected in a normal.  For humans, this is 2 */
    private int normalNumCopies;

    public CNLOHCaller(){
        rhoThreshold = 0.2;
        normalNumCopies = 2;
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
     * This is done by assuming with two Gaussians each modeling below the mode and above the mode.  While this is not a
     *  particularly good form for a distribution, it is quick to calculate and generates reasonable results.  This is
     *  also done, because we do not know the exact form of the distribution when coming from ACNV.
     *
     */
    @VisibleForTesting
    static double calculateFmaf(final double rho, final int m, final int n, final double credibleMode,
                                   final double credibleLow, final double credibleHigh) {
        ParamUtils.inRange(rho, 0.0, 1.0, "Invalid rho value: " + rho + ".  Must be [0,1]");
        ParamUtils.isPositiveOrZero(m, "M must be positive.");
        ParamUtils.isPositiveOrZero(n, "N must be positive.");

        if (Double.isNaN(credibleMode) || Double.isNaN(credibleLow) || Double.isNaN(credibleHigh)){
            return 1;
        }

        final double maf = calculateMaf(rho, m, n);
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
    static double calculateMaf(final double rho, final int m, final int n) {

        ParamUtils.inRange(rho, 0.0, 1.0, "Invalid rho value: " + rho + ".  Must be [0,1]");
        ParamUtils.isPositiveOrZero(m, "M must be >= 0");
        ParamUtils.isPositiveOrZero(n, "N must be >= 0");

        double result = ((1-rho) + rho * (Math.min(n,m))) / ( (1-rho)*2 + rho * (n+m) );

        // If result is NaN, that means that we had rho = 1 and M,N = 0
        if ((result < MIN_L) || (Double.isNaN(result))) {
            result = MIN_L;
        }
        return result;
    }

    /**
     * This function takes two Gaussian distributions and uses one for the values below the mode and the other for values above the
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

        double hiDist = credibleHigh - credibleMode;
        if (credibleMode >= credibleHigh) {
            hiDist = MIN_DIST_FROM_MODE;
        }
        double loDist = credibleMode - credibleLow;
        if (credibleMode <= credibleLow) {
            loDist = MIN_DIST_FROM_MODE;
        }

        if (val >= credibleMode) {
            return new NormalDistribution(credibleMode, hiDist / 1.96).density(val);
        } else {
            return new NormalDistribution(credibleMode, loDist / 1.96).density(val);
        }
    }

    /**
     * Calculate a copy ratio.
     *
     * @param rho CCF * purity
     * @param m Absolute copy number of allele 1
     * @param n Absolute copy number of allele 2
     * @param lambda CCF * ploidy
     * @return copy ratio in CR space
     */
    @VisibleForTesting
    static double calculateCopyRatio(final double rho, final int m, final int n, final double lambda) {
        return ((1-rho) * 2 + rho*(n+m)) / lambda;
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
     * @return likelihood, with a minimum value of {@link CNLOHCaller#MIN_L}
     */
    @VisibleForTesting
    static double calculateFcr(final double rho, final int m, final int n, final double lambda, final double credibleMode,
                         final double credibleLow, final double credibleHigh, final double additionalVariance) {
        ParamUtils.inRange(rho, 0.0, 1.0, "Invalid rho value: " + rho + ".  Must be [0,1]");
        ParamUtils.isPositiveOrZero(m, "M must be positive.");
        ParamUtils.isPositiveOrZero(n, "N must be positive.");
        ParamUtils.isPositiveOrZero(additionalVariance, "additional variance must be positive.");
        if (Double.isNaN(credibleMode) || Double.isNaN(credibleLow) || Double.isNaN(credibleHigh)){
            return 1;
        }
        final double cr = calculateCopyRatio(rho, m, n, lambda);
        final double avgDist = ((credibleMode - credibleLow) + (credibleHigh - credibleMode))/2.0;
        return Math.max(new NormalDistribution(credibleMode, (avgDist / 1.96) + Math.sqrt(additionalVariance)).density(cr), MIN_L);
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

        final Variance varianceCalculator = new Variance();

        // Only consider values "close enough" to copy neutral (CR == 1).
        return varianceCalculator.evaluate(segments.stream()
                .filter(s -> Math.abs(s.getSegmentMeanInCRSpace() - 1 + meanBiasInCR) < 0.1)
                .mapToDouble(s -> s.getSegmentMeanInCRSpace()).toArray());
    }

    private double calculateSegmentMeanBiasInCRSpace(final List<ACNVModeledSegment> segments) {
        Utils.nonNull(segments);
        final Percentile percentile = new Percentile();

        // Only consider values "close enough" to copy neutral (CR == 1).
        final double modeInCRSpace =  percentile.evaluate(segments.stream()
                .filter(s -> Math.abs(s.getSegmentMeanInCRSpace() - 1) < 0.1)
                .mapToDouble(s -> s.getSegmentMeanInCRSpace()).toArray(), 50);

        return modeInCRSpace - 1;
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
    public List<CNLOHCall> makeCalls(final List<ACNVModeledSegment> segments, final int numIterations, final JavaSparkContext ctx) {

        ParamUtils.isPositive(numIterations, "Must be more than zero iterations.");
        Utils.nonNull(segments);
        Utils.nonNull(ctx, "Java SparkContext can't be null when attempting to make CNLOH and balance calls.");

        segmentMeanBiasInCR = calculateSegmentMeanBiasInCRSpace(segments);
        segmentMeanVarianceInCR = calculateVarianceOfCopyNeutralSegmentMeans(segments, segmentMeanBiasInCR);

        CNLOHCallerModelState state = CNLOHCallerModelState.createInitialCNLOHCallerModelState(rhoThreshold, segments);

        // Create a Spark RDD for the segments.
        final JavaRDD<ACNVModeledSegment> segs = ctx.parallelize(segments);

        List<CNLOHCall> cnlohCalls = null;

        // for each iteration
        for (int i = 1; i <= numIterations; i++) {

            logger.info("Starting iteration " + i + " of " + numIterations + " ... ");
            logger.info("Calculating big E_zsk_vsm_wsn for each segment (iteration " + i + " of " + numIterations + ") ... ");

            // We are not broadcasting state, since it will only be used once before values change.  In that case,
            //  there is no gain for broadcasting

            // Effectively, the responsibilities are a 4D array:  S x K x M x N.  We do the responsibility calculation
            //  over the segments.  So the JavaRDD will be of length S.
            // Update E_zsk_vsm_wsn (responsibilities) for each segment (Steps 1 & 2)  {K x M x N} in S entries.  The RDD is of length S.
            final JavaRDD<double[][][]> eZskVsmWsnBySeg = segs.map(s -> calcE_zsk_vsm_wsn(state.getELnPhiK(), state.getCnToPiMap(), state.getRhos(),
                    s.getMinorAlleleFractionPosteriorSummary().getCenter(),
                    s.getMinorAlleleFractionPosteriorSummary().getLower(),
                    s.getMinorAlleleFractionPosteriorSummary().getUpper(),
                    Math.pow(2, s.getSegmentMeanPosteriorSummary().getCenter()) - segmentMeanBiasInCR,
                    Math.pow(2, s.getSegmentMeanPosteriorSummary().getLower()) - segmentMeanBiasInCR,
                    Math.pow(2, s.getSegmentMeanPosteriorSummary().getUpper()) - segmentMeanBiasInCR,
                    state.getLambda()));

            // Create an array (for each segment) that is a sum of eZskVsmWsnBySeg on the rhos.  (length:K for S entries)
            final JavaRDD<double[]> eZskBySeg = eZskVsmWsnBySeg.map(da -> sumOverFirstDimension(da));
            final JavaRDD<double[]> eVsmBySeg = eZskVsmWsnBySeg.map(da -> sumOverSecondDimension(da));
            final JavaRDD<double[]> eWsnBySeg = eZskVsmWsnBySeg.map(da -> sumOverThirdDimension(da));

            // Create arrays summed in one dimension (incl. over segs).
            final double[] eZskSummedOverSegs = sumOverSegs(eZskBySeg);
            final double[] eVsmSummedOverSegs = sumOverSegs(eVsmBySeg);
            final double[] eWsnSummedOverSegs = sumOverSegs(eWsnBySeg);

            final List<double[][][]> eZskVsmWsnBySegAsList = eZskVsmWsnBySeg.collect();

            // Update E_ln_phi_k (Step 3)
            logger.info("Updating E_ln_phi_k (iteration " + i + " of " + numIterations +  ") ... ");
            state.setELnPhiK(calc_E_ln_phi_k(state.getEAlpha(), eZskSummedOverSegs));

            // Update pis (Step 4)
            logger.info("Updating pis (iteration " + i + " of " + numIterations +  ") ... ");
            final double[] pis = calcPis(eVsmSummedOverSegs, eWsnSummedOverSegs);
            final List<Integer> keys = new ArrayList<>(state.getCnToPiMap().keySet());
            for (int p = 0; p < pis.length; p++) {
                state.getCnToPiMap().put(keys.get(p), pis[p]);
            }

            // Update E_alpha (Step 5)
            logger.info("Updating E_alpha (iteration " + i + " of " + numIterations + ") ... ");
            final int minCN = state.getCnToPiMap().keySet().stream().min(Integer::compare).get();
            final int maxCN = state.getCnToPiMap().keySet().stream().max(Integer::compare).get();
            state.setEAlpha(calcEAlpha(state.getELnPhiK(), minCN, maxCN));

            // Update rhos (Step 6)
            logger.info("Updating rhos (iteration " + i + ") ... ");
            state.setRhos(calcNewRhos(segments, eZskVsmWsnBySegAsList, state.getCnToPiMap(), state.getLambda(),
                    state.getRhos(), ctx));

            // This only really needs to  be done once, but is necessary here due to scoping issues that will cause
            //  compiler issues.
            cnlohCalls = createCNLoHCalls(segments, state, eZskVsmWsnBySegAsList);
        }
        logger.info("ACNV segments bias, variance (around copy neutral only): " + segmentMeanBiasInCR + ", " + segmentMeanVarianceInCR);
        logger.info("Output file is not adjusted for the bias and variance.  No user action required.");
        return cnlohCalls;
    }

    private List<CNLOHCall> createCNLoHCalls(List<ACNVModeledSegment> segments, CNLOHCallerModelState state, List<double[][][]> eZskVsmWsnBySegAsList) {

        final List<CNLOHCall> cnlohCalls = segments.stream().map(s -> new CNLOHCall(s)).collect(Collectors.toList());

        final List<Pair<Integer, Integer>> maxMiNiPerSeg = IntStream.range(0, segments.size()).boxed()
                .map(s -> max2dIndices(sumOnlyFirstDimension(eZskVsmWsnBySegAsList.get(s))))
                .collect(Collectors.toList());

        final List<double[]> rhosForBestMN = IntStream.range(0, segments.size()).boxed()
                .map(s -> getFirstDimensionArray(eZskVsmWsnBySegAsList.get(s), maxMiNiPerSeg.get(s).getFirst(), maxMiNiPerSeg.get(s).getSecond()))
                .collect(Collectors.toList());

        final List<Integer> maxRhoi = rhosForBestMN.stream().map(da -> maxIndex(da)).collect(Collectors.toList());
        final int[] mVals = state.getCnToPiMap().keySet().stream().mapToInt(Integer::intValue).toArray();
        final int[] nVals = state.getCnToPiMap().keySet().stream().mapToInt(Integer::intValue).toArray();

        final int[] calledMPerSeg = maxMiNiPerSeg.stream().mapToInt(p -> mVals[p.getFirst()]).toArray();
        final int[] calledNPerSeg = maxMiNiPerSeg.stream().mapToInt(p -> nVals[p.getSecond()]).toArray();
        final double[] calledRhoPerSeg = maxRhoi.stream().mapToDouble(ri -> state.getRhos()[ri]).toArray();

        // Set the rho, M, N, FCr, and FMaf for the call
        IntStream.range(0, cnlohCalls.size()).forEach(s -> cnlohCalls.get(s).setRho(calledRhoPerSeg[s]));
        IntStream.range(0, cnlohCalls.size()).forEach(s -> cnlohCalls.get(s).setM(calledMPerSeg[s]));
        IntStream.range(0, cnlohCalls.size()).forEach(s -> cnlohCalls.get(s).setN(calledNPerSeg[s]));
        IntStream.range(0, cnlohCalls.size()).forEach(s -> cnlohCalls.get(s)
                        .setfCr(calculateFcr(calledRhoPerSeg[s], calledMPerSeg[s], calledNPerSeg[s], state.getLambda(),
                                Math.pow(2, cnlohCalls.get(s).getAcnvSegment().getSegmentMeanPosteriorSummary().getCenter()),
                                Math.pow(2, cnlohCalls.get(s).getAcnvSegment().getSegmentMeanPosteriorSummary().getLower()),
                                Math.pow(2, cnlohCalls.get(s).getAcnvSegment().getSegmentMeanPosteriorSummary().getUpper()),
                                segmentMeanVarianceInCR
                        ))
        );
        IntStream.range(0, cnlohCalls.size()).forEach(s -> cnlohCalls.get(s)
                        .setfMaf(calculateFmaf(calledRhoPerSeg[s], calledMPerSeg[s], calledNPerSeg[s],
                                        cnlohCalls.get(s).getAcnvSegment().getMinorAlleleFractionPosteriorSummary().getCenter(),
                                        cnlohCalls.get(s).getAcnvSegment().getMinorAlleleFractionPosteriorSummary().getLower(),
                                        cnlohCalls.get(s).getAcnvSegment().getMinorAlleleFractionPosteriorSummary().getUpper())
                        )
        );

        // Set the balanced calls
        IntStream.range(0, cnlohCalls.size()).forEach(s -> cnlohCalls.get(s).setBalancedCall(isBalanced(calledRhoPerSeg[s], calledMPerSeg[s], calledNPerSeg[s])));

        // Set the CNLoH calls
        IntStream.range(0, cnlohCalls.size()).forEach(s -> cnlohCalls.get(s).setCnlohCall(isCNLoH(calledRhoPerSeg[s], calledMPerSeg[s], calledNPerSeg[s], cnlohCalls.get(s).getAcnvSegment())));

        return cnlohCalls;
    }

    private CNLOHBalancedCallEnum isBalanced(final double rho, final int m, final int n) {
        if ((m == n) || (rho <= getRhoThreshold())) {
            return CNLOHBalancedCallEnum.BALANCED;
        } else {
            return CNLOHBalancedCallEnum.NOT_BALANCED;
        }
    }

    private CNLOHLoHCallEnum isCNLoH(final double rho, final int m, final int n, final ACNVModeledSegment segment) {
        if (Double.isNaN(segment.getMinorAlleleFractionPosteriorSummary().getCenter())) {
            return CNLOHLoHCallEnum.NO_CALL;
        } else if (((m == getNormalNumCopies()) && (n == 0)) || ((n == getNormalNumCopies()) && (m == 0))
            && (rho > getRhoThreshold())) {
            return CNLOHLoHCallEnum.CNLOH;
        } else {
            return CNLOHLoHCallEnum.NOT_CNLOH;
        }
    }

    private double[] sumOverSegs(JavaRDD<double[]> eZskBySeg) {
        final double[][] eZskBySeg2D = eZskBySeg.collect().stream().toArray(double[][]::new); // S x K
        final RealMatrix eZskBySegMatrix = MatrixUtils.createRealMatrix(eZskBySeg2D);
        return GATKProtectedMathUtils.columnSum(eZskBySegMatrix);
    }

    private double[] sumOverFirstDimension(final double[][][] data) {
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

    private int maxIndex(final double[] data) {
        double maxValue = Double.NEGATIVE_INFINITY;
        int result = -1;
        for (int i=0; i < data.length; i++) {
            if (data[i] > maxValue) {
                maxValue = data[i];
                result = i;
            }
        }
        return result;
    }

    private double[] getFirstDimensionArray(final double[][][] array, final int dim2Index, final int dim3Index) {
        final double[] result = new double[array.length];
        for (int i = 0; i < result.length; i++) {
            result[i] = array[i][dim2Index][dim3Index];
        }
        return result;
    }

    @VisibleForTesting
    static double[][] sumOnlyFirstDimension(final double[][][] data) {
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

    private double[] sumOverSecondDimension(final double[][][] data) {
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

    private double[] sumOverThirdDimension(final double[][][] data) {
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
    protected double[][][] calcE_zsk_vsm_wsn(final double[] E_ln_phi_k, final SortedMap<Integer, Double> pis,
                                             final double[] rhos, final double mafCredibleMode,
                                             final double mafCredibleLow, final double mafCredibleHigh,
                                             final double crCredibleMode, final double crCredibleLow,
                                             final double crCredibleHigh, final double lambda) {


        final int[] mVals = pis.keySet().stream().mapToInt(Integer :: intValue).toArray();
        final int[] nVals = pis.keySet().stream().mapToInt(Integer :: intValue).toArray();

        final double[][][] result = new double[rhos.length][mVals.length][nVals.length];
        double total = 0.0;
        for (int rhoi = 0; rhoi < rhos.length; rhoi++) {
            for (int mi = 0; mi < mVals.length; mi++) {
                for (int ni = 0; ni < nVals.length; ni++) {

                    final double rho = rhos[rhoi];
                    final int m = mVals[mi];
                    final int n = nVals[ni];

                    if ((rho < 0) || (rho > 1)) {
                        result[rhoi][mi][ni] = 0;
                    } else if ((rho == 0) && ((m != 1) || (n != 1))) {
                        result[rhoi][mi][ni] = 0;
                    } else if ((rho > 0) && (rho < rhoThreshold)) {
                        result[rhoi][mi][ni] = 0;
                    } else if ((rho != 0) && ((m == 1) && (n == 1))) {
                        result[rhoi][mi][ni] = 0;
                    } else {
                        result[rhoi][mi][ni] = Math.exp(E_ln_phi_k[rhoi]) *
                                pis.get(m) * pis.get(n) *
                                calculateFmaf(rho, m, n, mafCredibleMode, mafCredibleLow, mafCredibleHigh) *
                                calculateFcr(rho, m, n, lambda, crCredibleMode, crCredibleLow, crCredibleHigh, segmentMeanVarianceInCR);
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

    protected double[] calc_E_ln_phi_k(final double E_alpha, final double[] eZskSummedOverSegs) {

        /**
         *
         E_zsk_summed_over_segs = sum(E_zsk_for_each_seg, 2); % K x 1
         E_zsk_summed_over_k_and_segs = sum(E_zsk_summed_over_segs);

         % Slight modifications to stabilize values (+1 and +K) and add a bias towards rho = 0.
         pseudocounts = ones(size(E_zsk_summed_over_segs));
         pseudocounts(1) = 20;
         term1 = psi(E_alpha/K + pseudocounts + E_zsk_summed_over_segs);
         term2 = psi(E_alpha + sum(pseudocounts) + E_zsk_summed_over_k_and_segs);

         result = term1 - term2;

         */

        final double eZskSummedOverAll = GATKProtectedMathUtils.sum(eZskSummedOverSegs);

        final double[] result = new double[eZskSummedOverSegs.length];
        final int k = eZskSummedOverSegs.length;

        // Initialize all pseudocounts to 1, except for index 0, which is 20;
        //  This artificially increases the odds for a rho = 0.
        final RealVector pseudocounts = new ArrayRealVector(eZskSummedOverSegs.length);
        pseudocounts.set(1.0);
        pseudocounts.setEntry(0, 20.0);
        final double sumPseudocounts = GATKProtectedMathUtils.sum(pseudocounts.toArray());

        for (int i=0; i < result.length; i++) {
            final double term1 = Gamma.digamma(E_alpha/k + pseudocounts.getEntry(i) + eZskSummedOverSegs[i]);
            final double term2 = Gamma.digamma(E_alpha + sumPseudocounts + eZskSummedOverAll);
            result[i] = term1 - term2;
        }

        return result;
    }

    protected double[] calcPis(final double[] eVsm, final double[] eWsn) {

        // Sum (and normalize) the responsibilities for each possible value of the absolute copy numbers.

        double total = 0;
        for (int i = 0; i < eVsm.length; i++) {
            total += eVsm[i] + eWsn[i];
        }

        final double[] result = new double[eVsm.length];
        for (int i = 0; i < result.length; i++) {
            result[i] = (eVsm[i] + eWsn[i]) / total;
        }
        return result;
    }

    protected double calcEAlpha(final double[] E_ln_phi_k, final int minCN, final int maxCN) {
        /**
         *
         numerator = integral(@(a) (a.*q_alpha(a, K, E_ln_phi_k)), 0, 5);
         denominator = integral(@(a) (q_alpha(a, K, E_ln_phi_k)), 0, 5);
         if numerator == 0
         result = 0;
         else
         result = numerator/denominator;
         end
         */

        final QAlphaUnivariateFunction qAlpha = new QAlphaUnivariateFunction(E_ln_phi_k);
        final AQAlphaUnivariateFunction aQAlpha = new AQAlphaUnivariateFunction(E_ln_phi_k);
        final BaseAbstractUnivariateIntegrator integrator = new SimpsonIntegrator();

        final double numerator = integrator.integrate(10, aQAlpha, minCN, maxCN);
        final double denominator = integrator.integrate(10, qAlpha, minCN, maxCN);

        // If denominator and numerator are both zero, we want to return zero.  But that simplifies to the
        //  numerator being zero.
        if (numerator == 0) {
            return 0;
        } else {
            return numerator/denominator;
        }
    }

    private double[] calcNewRhos(final List<ACNVModeledSegment> segments,
                               final List<double[][][]> eZskVsmWsnBySeg, final SortedMap<Integer, Double> pis,
                               final double lambda, final double[] rhos, final JavaSparkContext ctx) {

        final int[] mVals = pis.keySet().stream().mapToInt(i->i).toArray();
        final int[] nVals = pis.keySet().stream().mapToInt(i->i).toArray();

        // Since, we pass in the entire responsibilities matrix, we need the correct index for each rho.  That, and the
        //  fact that this is a univariate objective function, means we need to create an instance for each rho.  And
        //  then we blast across Spark.

        final List<Pair<? extends Function<Double, Double>, SearchInterval>> objectives = IntStream.range(0, rhos.length)
                .mapToObj(i -> new Pair<>(
                        new Function<Double, Double>() {
                            @Override
                            public Double apply(Double rho) {
                                return calculateESmnObjective(rho, segments, eZskVsmWsnBySeg, mVals, nVals, lambda, i);
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
                                          final List<double[][][]> eZskVsmWsnBySeg,
                                          final int[] m_vals, final int[] n_vals, final double lambda, final int rhoIndex) {

        // We will want to sum an entire matrix that is S x M x N for the given rho.
        final double[][][] eSMN = new double[eZskVsmWsnBySeg.size()][m_vals.length][n_vals.length];

        // Populate eSMN
        for (int s=0; s<eZskVsmWsnBySeg.size(); s++) {
            final ACNVModeledSegment seg = segments.get(s);
            final double mafMode = seg.getMinorAlleleFractionPosteriorSummary().getCenter();
            final double mafLow = seg.getMinorAlleleFractionPosteriorSummary().getLower();
            final double mafHigh = seg.getMinorAlleleFractionPosteriorSummary().getUpper();
            final double crMode = Math.pow(2, seg.getSegmentMeanPosteriorSummary().getCenter()) - segmentMeanBiasInCR;
            final double crLow = Math.pow(2, seg.getSegmentMeanPosteriorSummary().getLower()) - segmentMeanBiasInCR;
            final double crHigh = Math.pow(2, seg.getSegmentMeanPosteriorSummary().getUpper()) - segmentMeanBiasInCR;

            for (int m=0; m<m_vals.length; m++) {
                for (int n=0; n<n_vals.length; n++) {
                    final double mafLikelihood = calculateFmaf(rho, m_vals[m], n_vals[n], mafMode, mafLow, mafHigh);
                    final double crLikelihood = calculateFcr(rho, m_vals[m], n_vals[n], lambda, crMode, crLow, crHigh, segmentMeanVarianceInCR);
                    if ( ((rho > 1) || (rho < 0)) || ((rho > 0) && (rho < rhoThreshold)))  {
                        eSMN[s][m][n] = MIN_L;
                    } else {
                        eSMN[s][m][n] = eZskVsmWsnBySeg.get(s)[rhoIndex][m][n] * Math.log(mafLikelihood * crLikelihood);
                    }
                }
            }
        }

        return GATKProtectedMathUtils.sum(eSMN);
    }

    private class QAlphaUnivariateFunction implements UnivariateFunction {
        private static final double GAMMA_SHAPE = 1.0;
        private static final double GAMMA_SCALE = 1.0;

        private double[] E_ln_phi_k;
        private double sumE_ln_phi_k;
        private GammaDistribution gammaDistribution;

        public QAlphaUnivariateFunction(final double[] E_ln_phi_k) {
            this.E_ln_phi_k = E_ln_phi_k;
            this.sumE_ln_phi_k = GATKProtectedMathUtils.sum(E_ln_phi_k);
            this.gammaDistribution = new GammaDistribution(GAMMA_SHAPE, GAMMA_SCALE);
        }

        @Override
        public double value(double a) {
            /**
             * P_alpha = gampdf(a, 1, 1);
             gamma_alpha = gamma(max(a, eps));
             K = length(k);

             numerator = P_alpha .* gamma_alpha;
             denominator = gamma(a/K).^K;

             coeff = numerator ./ denominator;

             result = coeff .* exp( (a/K) * sum(E_ln_phi_k) );
             */
            // Cannot evaluate a gamma at zero, so set a to a minimum value..
            final double aPrime = Math.max(a, MIN_L);
            final double pAlpha = gammaDistribution.probability(aPrime);
            final double gammaAlpha = Gamma.gamma(aPrime);
            final int k = E_ln_phi_k.length;

            final double numerator = pAlpha * gammaAlpha;
            final double denominator = Math.pow(Gamma.gamma(aPrime/k), k);
            final double coeff = numerator/denominator;

            return coeff * Math.exp((aPrime/k) * sumE_ln_phi_k);
        }
    }

    private class AQAlphaUnivariateFunction implements UnivariateFunction {
        private QAlphaUnivariateFunction qAlpha;

        public AQAlphaUnivariateFunction(final double[] E_ln_phi_k) {
            qAlpha = new QAlphaUnivariateFunction(E_ln_phi_k);
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

    public static <T> T roundTripInKryo(T input, Class<?> inputClazz, final SparkConf conf) {
        KryoSerializer kryoSerializer = new KryoSerializer(conf);
        SerializerInstance sparkSerializer = kryoSerializer.newInstance();
        final ClassTag<T> tag = ClassTag$.MODULE$.apply(inputClazz);
        return sparkSerializer.deserialize(sparkSerializer.serialize(input, tag), tag);
    }
}
