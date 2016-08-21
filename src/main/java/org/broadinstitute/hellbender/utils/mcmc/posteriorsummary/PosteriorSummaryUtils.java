package org.broadinstitute.hellbender.utils.mcmc.posteriorsummary;

import com.google.common.primitives.Doubles;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.SearchInterval;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.mllib.stat.KernelDensity;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Contains methods for computing summary statistics for the posterior of a univariate model parameter.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public class PosteriorSummaryUtils {
    public static final double SILVERMANS_RULE_CONSTANT = 1.06;
    public static final double SILVERMANS_RULE_EXPONENT = -0.2;
    //constants for Brent optimization
    private static final MaxEval BRENT_MAX_EVAL = new MaxEval(100);
    private static final double RELATIVE_TOLERANCE = 0.01;

    private PosteriorSummaryUtils() {
    }

    /**
     * Given a list of posterior samples, returns a PosteriorSummary with central tendency given by the posterior mode
     * (which is estimated using mllib kernel density estimation in {@link KernelDensity},
     * see {@link PosteriorSummaryUtils#calculatePosteriorMode(List, JavaSparkContext)}) and lower/upper credible-interval
     * bounds given by the (1-{@code alpha}) highest-posterior-density interval (i.e., the narrowest interval
     * that contains (1-{@code alpha}) of the samples).  Unimodality is assumed.  If the samples contain
     * {@link Double#NaN}, a {@link PosteriorSummary} with
     * {@link PosteriorSummary#center}, {@link PosteriorSummary#lower}, and {@link PosteriorSummary#upper}
     * all set to {@link Double#NaN} will be returned.
     * @param samples   posterior samples, cannot be {@code null} and number of samples must be greater than 0
     * @param alpha     credible-interval alpha, must be in (0, 1)
     * @param ctx       {@link JavaSparkContext} used by {@link KernelDensity} for mllib kernel density estimation
     */
    public static PosteriorSummary calculateHighestPosteriorDensitySummary(final List<Double> samples,
                                                                           final double alpha,
                                                                           final JavaSparkContext ctx) {
        Utils.nonNull(samples);
        Utils.validateArg(samples.size() > 0, "Number of samples must be greater than zero.");
        Utils.validateArg(0 < alpha && alpha < 1, "Alpha must be in (0, 1).");

        final double central = calculatePosteriorMode(samples, ctx);

        //if samples contain NaN, return NaN for all posterior-summary values
        if (Double.isNaN(central)) {
            return new PosteriorSummary(Double.NaN, Double.NaN, Double.NaN);
        }

        //find highest-posterior-density interval using simple Chen & Shao 1998 procedure
        final List<Double> sortedSamples = new ArrayList<>(samples);
        Collections.sort(sortedSamples);
        final int n = sortedSamples.size();
        double lower = sortedSamples.get(0);
        double upper = sortedSamples.get(n - 1);
        double minIntervalWidth = sortedSamples.get(n - 1) - sortedSamples.get(0);
        final int numSamplesInInterval = (int) Math.floor((1 - alpha) * n);
        for (int i = 0; i < n - numSamplesInInterval; i++) {
            final double intervalWidth = sortedSamples.get(i + numSamplesInInterval) - sortedSamples.get(i);
            if (intervalWidth < minIntervalWidth) {
                minIntervalWidth = intervalWidth;
                lower = sortedSamples.get(i);
                upper = sortedSamples.get(i + numSamplesInInterval);
            }
        }

        return new PosteriorSummary(central, lower, upper);
    }

    /**
     * Given a list of posterior samples, returns a PosteriorSummary with central tendency given by the posterior mode
     * (which is estimated using mllib kernel density estimation in {@link KernelDensity},
     * see {@link PosteriorSummaryUtils#calculatePosteriorMode(List, JavaSparkContext)}), lower/upper credible-interval
     * bounds given by the (1-{@code alpha}) highest-posterior-density interval (i.e., the narrowest interval
     * that contains (1-{@code alpha}) of the samples), and deciles.  Unimodality is assumed.  If the samples contain
     * {@link Double#NaN}, a {@link PosteriorSummary} with
     * {@link PosteriorSummary#center}, {@link PosteriorSummary#lower}, and {@link PosteriorSummary#upper}
     * all set to {@link Double#NaN} will be returned.
     * @param samples   posterior samples, cannot be {@code null} and number of samples must be greater than 0
     * @param alpha     credible-interval alpha, must be in (0, 1)
     * @param ctx       {@link JavaSparkContext} used by {@link KernelDensity} for mllib kernel density estimation
     */
    public static PosteriorSummary calculateHighestPosteriorDensityAndDecilesSummary(final List<Double> samples,
                                                                                     final double alpha,
                                                                                     final JavaSparkContext ctx) {
        Utils.nonNull(samples);
        Utils.validateArg(samples.size() > 0, "Number of samples must be greater than zero.");
        Utils.validateArg(0 < alpha && alpha < 1, "Alpha must be in (0, 1).");

        final PosteriorSummary posteriorSummary = calculateHighestPosteriorDensitySummary(samples, alpha, ctx);
        final DecileCollection deciles = new DecileCollection(samples, DecileCollection.ConstructionMode.SAMPLES);
        posteriorSummary.setDeciles(deciles);
        return posteriorSummary;
    }

    /**
     * Given a list of posterior samples, returns an estimate of the posterior mode (using
     * mllib kernel density estimation in {@link KernelDensity} and {@link BrentOptimizer}).
     * Note that estimate may be poor if number of samples is small (resulting in poor kernel density estimation),
     * or if posterior is not unimodal (or is sufficiently pathological otherwise). If the samples contain
     * {@link Double#NaN}, {@link Double#NaN} will be returned.
     * @param samples   posterior samples, cannot be {@code null} and number of samples must be greater than 0
     * @param ctx       {@link JavaSparkContext} used by {@link KernelDensity} for mllib kernel density estimation
     */
    public static double calculatePosteriorMode(final List<Double> samples, final JavaSparkContext ctx) {
        Utils.nonNull(samples);
        Utils.validateArg(samples.size() > 0, "Number of samples must be greater than zero.");

        //calculate sample min, max, mean, and standard deviation
        final double sampleMin = Collections.min(samples);
        final double sampleMax = Collections.max(samples);
        final double sampleMean = new Mean().evaluate(Doubles.toArray(samples));
        final double sampleStandardDeviation = new StandardDeviation().evaluate(Doubles.toArray(samples));

        //if samples are all the same or contain NaN, can simply return mean
        if (sampleStandardDeviation == 0. || Double.isNaN(sampleMean)) {
            return sampleMean;
        }

        //use Silverman's rule to set bandwidth for kernel density estimation from sample standard deviation
        //see https://en.wikipedia.org/wiki/Kernel_density_estimation#Practical_estimation_of_the_bandwidth
        final double bandwidth =
                SILVERMANS_RULE_CONSTANT * sampleStandardDeviation * Math.pow(samples.size(), SILVERMANS_RULE_EXPONENT);

        //use kernel density estimation to approximate posterior from samples
        final KernelDensity pdf = new KernelDensity().setSample(ctx.parallelize(samples, 1)).setBandwidth(bandwidth);

        //use Brent optimization to find mode (i.e., maximum) of kernel-density-estimated posterior
        final BrentOptimizer optimizer =
                new BrentOptimizer(RELATIVE_TOLERANCE, RELATIVE_TOLERANCE * (sampleMax - sampleMin));
        final UnivariateObjectiveFunction objective =
                new UnivariateObjectiveFunction(f -> pdf.estimate(new double[] {f})[0]);
        //search for mode within sample range, start near sample mean
        final SearchInterval searchInterval = new SearchInterval(sampleMin, sampleMax, sampleMean);
        return optimizer.optimize(objective, GoalType.MAXIMIZE, searchInterval, BRENT_MAX_EVAL).getPoint();
    }
}
