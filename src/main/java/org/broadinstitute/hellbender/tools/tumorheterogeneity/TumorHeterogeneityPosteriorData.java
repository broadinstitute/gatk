package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import com.google.common.primitives.Doubles;
import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.MultivariateOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.NelderMeadSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.apache.commons.math3.util.FastMath;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.exome.ACNVModeledSegment;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.mcmc.DataCollection;
import org.broadinstitute.hellbender.utils.mcmc.DecileCollection;

import java.util.Arrays;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * {@link DataCollection} for the tumor-heterogeneity model that allows the calculation of log posterior probabilities
 * for (copy ratio, minor-allele fraction) for each {@link ACNVModeledSegment}.  Fits a normal distribution to
 * the log_2 copy-ratio posterior deciles and a scaled beta distribution to the minor-allele fraction posterior deciles.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class TumorHeterogeneityPosteriorData implements DataCollection {
    private static final double LN2 = Math.log(2.);
    private static final double INV_LN2 = 1. / LN2;
    private static final double LN_LN2 = Math.log(LN2);
    private static final double VARIANCE_SAMPLE_TO_POPULATION_CONVERSION_FACTOR =
            (DecileCollection.NUM_DECILES - 2.) / (DecileCollection.NUM_DECILES - 3.);  //n / (n - 1), where n = number of inner deciles = 9
    private static final double STANDARD_DEVIATION_SAMPLE_TO_POPULATION_CONVERSION_FACTOR =
            Math.sqrt(VARIANCE_SAMPLE_TO_POPULATION_CONVERSION_FACTOR);

    //parameters for optimizer
    private static final double REL_TOLERANCE = 1E-5;
    private static final double ABS_TOLERANCE = 1E-10;
    private static final int NUM_MAX_EVAL = 1000;
    private static final double DEFAULT_SIMPLEX_STEP = 0.2;

    public static final Logger logger = LogManager.getLogger(TumorHeterogeneityPosteriorData.class);
    private static final MultivariateOptimizer optimizer = new SimplexOptimizer(REL_TOLERANCE, ABS_TOLERANCE);

    private final List<ACNVSegmentPosterior> segmentPosteriors;

    public TumorHeterogeneityPosteriorData(final List<ACNVModeledSegment> segments) {
        Utils.nonNull(segments);
        Utils.validateArg(segments.size() > 0, "Number of segments must be positive.");
        segmentPosteriors = segments.stream().map(ACNVSegmentPosterior::new).collect(Collectors.toList());
    }

    public double logDensity(final int segmentIndex, final double copyRatio, final double minorAlleleFraction) {
        Utils.validateArg(0 <= segmentIndex && segmentIndex < segmentPosteriors.size(), "Segment index is not in valid range.");
        Utils.validateArg(copyRatio >= 0, "Copy ratio must be non-negative.");
        Utils.validateArg(0 <= minorAlleleFraction && minorAlleleFraction <= 0.5, "Minor-allele fraction must be in [0, 0.5].");
        return segmentPosteriors.get(segmentIndex).logDensity(copyRatio, minorAlleleFraction);
    }
    
    private final class ACNVSegmentPosterior {
        private final Function<Double, Double> log2CopyRatioPosteriorLogPDF;
        private final Function<Double, Double> minorAlleleFractionPosteriorLogPDF;

        ACNVSegmentPosterior(final ACNVModeledSegment segment) {
            final List<Double> log2CopyRatioInnerDecilesList = segment.getSegmentMeanPosteriorSummary().getDeciles().getInner();
            final double[] log2CopyRatioInnerDeciles = Doubles.toArray(log2CopyRatioInnerDecilesList);
            logger.info("Fitting normal distribution to inner deciles:\n" + log2CopyRatioInnerDecilesList);
            log2CopyRatioPosteriorLogPDF = fitNormalLogPDFToInnerDeciles(log2CopyRatioInnerDeciles);

            final List<Double> minorAlleleFractionInnerDecilesList = segment.getMinorAlleleFractionPosteriorSummary().getDeciles().getInner();
            final double[] minorAlleleFractionInnerDeciles = Doubles.toArray(minorAlleleFractionInnerDecilesList);
            logger.info("Fitting scaled beta distribution to inner deciles:\n" + minorAlleleFractionInnerDecilesList);
            final boolean isMinorAlleleFractionNaN = Double.isNaN(segment.getMinorAlleleFractionPosteriorSummary().getCenter());
            minorAlleleFractionPosteriorLogPDF = isMinorAlleleFractionNaN ?
                    f -> LN2 :       //flat over minor-allele fraction if NaN (i.e., no hets in segment) = log(1. / 0.5)
                    fitScaledBetaLogPDFToInnerDeciles(minorAlleleFractionInnerDeciles);
        }

        double logDensity(final double copyRatio, final double minorAlleleFraction) {
            final double log2CopyRatio = Math.log(copyRatio) * INV_LN2;
            final double copyRatioPosteriorLogDensity =
                    log2CopyRatioPosteriorLogPDF.apply(log2CopyRatio) - LN_LN2 - Math.log(copyRatio);    //includes Jacobian: p(c) = p(log_2(c)) / (c * ln 2)
            final double minorAlleleFractionPosteriorLogDensity = minorAlleleFractionPosteriorLogPDF.apply(minorAlleleFraction);
            return copyRatioPosteriorLogDensity + minorAlleleFractionPosteriorLogDensity;
        }

        //fit a normal distribution to inner deciles (10th, 20th, ..., 90th percentiles) using least squares
        //and return the log PDF (for use as a posterior for log_2 copy ratio)
        private Function<Double, Double> fitNormalLogPDFToInnerDeciles(final double[] innerDeciles) {
            final ObjectiveFunction innerDecilesL2LossFunction = new ObjectiveFunction(point -> {
                final double mean = point[0];
                final double standardDeviation = Math.abs(point[1]);
                final NormalDistribution normalDistribution = new NormalDistribution(mean, standardDeviation);
                final List<Double> normalInnerDeciles = IntStream.range(1, DecileCollection.NUM_DECILES - 1).boxed()
                        .map(i -> normalDistribution.inverseCumulativeProbability(i / 10.)).collect(Collectors.toList());
                return IntStream.range(0, DecileCollection.NUM_DECILES - 2)
                        .mapToDouble(i -> FastMath.pow(innerDeciles[i] - normalInnerDeciles.get(i), 2)).sum();
            });
            final double meanInitial = new Mean().evaluate(innerDeciles);
            final double standardDeviationInitial = new StandardDeviation().evaluate(innerDeciles) * STANDARD_DEVIATION_SAMPLE_TO_POPULATION_CONVERSION_FACTOR;
            logger.info(String.format("Initial (mean, standard deviation) for normal distribution: (%f, %f)", meanInitial, standardDeviationInitial));
            final PointValuePair optimum = optimizer.optimize(
                            new MaxEval(NUM_MAX_EVAL),
                            innerDecilesL2LossFunction,
                            GoalType.MINIMIZE,
                            new InitialGuess(new double[]{meanInitial, standardDeviationInitial}),
                            new NelderMeadSimplex(new double[]{DEFAULT_SIMPLEX_STEP, DEFAULT_SIMPLEX_STEP}));
            final double mean = optimum.getPoint()[0];
            final double standardDeviation = Math.abs(optimum.getPoint()[1]);
            logger.info(String.format("Final (mean, standard deviation) for normal distribution: (%f, %f)", mean, standardDeviation));
            return log2cr -> new NormalDistribution(mean, standardDeviation).logDensity(log2cr);
        }

        //fit a beta distribution to inner deciles (10th, 20th, ..., 90th percentiles) using least squares
        //and return the log PDF (scaled appropriately for use as a posterior for minor-allele fraction)
        private Function<Double, Double> fitScaledBetaLogPDFToInnerDeciles(final double[] innerDeciles) {
            final double[] scaledInnerDeciles = Arrays.stream(innerDeciles).map(d -> 2. * d).toArray();
            //scale minor-allele fraction deciles to [0, 1] and fit a beta distribution
            final ObjectiveFunction innerDecilesL2LossFunction = new ObjectiveFunction(point -> {
                final double alpha = Math.abs(point[0]);
                final double beta = Math.abs(point[1]);
                final BetaDistribution betaDistribution = new BetaDistribution(alpha, beta);
                final List<Double> betaInnerDeciles = IntStream.range(1, DecileCollection.NUM_DECILES - 1).boxed()
                        .map(i -> betaDistribution.inverseCumulativeProbability(i / 10.)).collect(Collectors.toList());
                return IntStream.range(0, DecileCollection.NUM_DECILES - 2)
                        .mapToDouble(i -> FastMath.pow(scaledInnerDeciles[i] - betaInnerDeciles.get(i), 2)).sum();
            });
            //use moment matching of deciles to initialize
            final double meanInitial = new Mean().evaluate(scaledInnerDeciles);
            final double varianceInitial = new Variance().evaluate(scaledInnerDeciles) * VARIANCE_SAMPLE_TO_POPULATION_CONVERSION_FACTOR;
            final double commonFactor = Math.abs((meanInitial - meanInitial * meanInitial) / varianceInitial - 1.);
            final double alphaInitial = meanInitial * commonFactor;
            final double betaInitial = (1. - meanInitial) * commonFactor;
            logger.info(String.format("Initial (alpha, beta) for scaled beta distribution: (%f, %f)", alphaInitial, betaInitial));
            final PointValuePair optimum = optimizer.optimize(
                    new MaxEval(NUM_MAX_EVAL),
                    innerDecilesL2LossFunction,
                    GoalType.MINIMIZE,
                    new InitialGuess(new double[]{alphaInitial, betaInitial}),
                    new NelderMeadSimplex(new double[]{DEFAULT_SIMPLEX_STEP, DEFAULT_SIMPLEX_STEP}));
            final double alpha = Math.abs(optimum.getPoint()[0]);
            final double beta = Math.abs(optimum.getPoint()[1]);
            logger.info(String.format("Final (alpha, beta) for scaled beta distribution: (%f, %f)", alpha, beta));
            return maf -> new BetaDistribution(alpha, beta).logDensity(2. * maf) + LN2; //scale minor-allele fraction to [0, 1], including Jacobian factor
        }
    }
}