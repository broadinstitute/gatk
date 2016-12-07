package org.broadinstitute.hellbender.tools.exome.pulldown;

import com.google.common.annotations.VisibleForTesting;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.analysis.integration.gauss.GaussIntegrator;
import org.apache.commons.math3.analysis.integration.gauss.GaussIntegratorFactory;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Heterogeneous model prior for heterozygous pileups.
 *
 * This prior is suitable for detecting heterozygous sites from reads that come from tumor or contaminated
 * normal samples.
 *
 * <ul>
 *      <li> The quadrature order can be adaptively chosen based on the pileup size and the accuracy
 *      required for likelihood estimation. In theory, a Gaussian quadrature of order N yields
 *      the exact result for a pileup of size 2N. </li>
 * </ul>
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class HeterogeneousHeterozygousPileupPriorModel extends HeterozygousPileupPriorModel {

    /* these are for building a prior for allele fraction */
    private final double minAbnormalFraction;
    private final double maxAbnormalFraction;
    private final double maxCopyNumber;
    private final double minHetAlleleFraction;
    private final double breakpointHetAlleleFraction;

    /* integration quadrature */
    @VisibleForTesting
    final List<Double> gaussIntegrationWeights = new ArrayList<>();
    private final List<Double> gaussIntegrationLogWeights = new ArrayList<>();
    @VisibleForTesting
    final List<Double> gaussIntegrationAbscissas = new ArrayList<>();

    /* allele fraction prior for Het sites */
    @VisibleForTesting
    final List<Double> alleleFractionPriors = new ArrayList<>();
    private final List<Double> alleleFractionLogPriors = new ArrayList<>();

    /* minimum order of the integration quadrature */
    private static final int MIN_QUADRATURE_ORDER = 50;

    /**
     * Initialize the prior
     * @param minAbnormalFraction estimated minimum fraction of non-germline cells in the sample
     * @param maxAbnormalFraction estimated maximum fraction of non-germline cells in the sample
     * @param maxCopyNumber estimated maximum copy number in non-germline events (note: we use a flat probability
     *                      function for the copy number. using a large value of maxCopyNumber will result in a
     *                      significant spread of the minor allele fraction prior around 1/2. it is recommended not
     *                      to use values > 4).
     * @param quadratureOrder the order of quadrature used in numerical integrations
     */
    public HeterogeneousHeterozygousPileupPriorModel(final double minAbnormalFraction, final double maxAbnormalFraction,
                                                     final int maxCopyNumber, final int quadratureOrder) {
        this.minAbnormalFraction = ParamUtils.inRange(minAbnormalFraction, 0.0, 1.0, "Minimum fraction of abnormal" +
                " cells must be between 0 and 1.");
        this.maxAbnormalFraction = ParamUtils.inRange(maxAbnormalFraction, this.minAbnormalFraction, 1.0, "Maximum fraction of abnormal" +
                " cells must be greater than the provided minimum and less than 1.");
        this.maxCopyNumber = ParamUtils.isPositive(maxCopyNumber, "Maximum copy number must be positive");

        /* auxiliary derived members */
        minHetAlleleFraction = (1 - this.maxAbnormalFraction) / (this.maxCopyNumber * this.maxAbnormalFraction +
                2 * (1 - this.maxAbnormalFraction));
        breakpointHetAlleleFraction = (1 - this.minAbnormalFraction) / (this.maxCopyNumber * this.minAbnormalFraction +
                2 * (1 - this.minAbnormalFraction));

        /* initialize the integration quadrature and calculate the allele fraction prior on the abscissas */
        initializeIntegrationQuadrature(ParamUtils.isPositive(quadratureOrder - MIN_QUADRATURE_ORDER,
                "Quadrature order must be greater than " + MIN_QUADRATURE_ORDER) + MIN_QUADRATURE_ORDER);
        initializeHetAlleleFractionPrior();
    }

    /**
     * Initilizes the quadrature for calculating allele ratio integrals in
     * {@link HeterogeneousHeterozygousPileupPriorModel#getHetLogLikelihood(List)}
     *
     * @param numIntegPoints  number of points in the quadrature
     */
    private void initializeIntegrationQuadrature(final int numIntegPoints) {
        /* get Gauss-Legendre quadrature factory of order @numIntegPoints */
        final GaussIntegratorFactory integratorFactory = new GaussIntegratorFactory();
        final GaussIntegrator gaussIntegrator = integratorFactory.legendre(numIntegPoints,
                minHetAlleleFraction, 1.0 - minHetAlleleFraction);

        /* abscissas */
        gaussIntegrationAbscissas.clear();
        gaussIntegrationAbscissas.addAll(IntStream.range(0, numIntegPoints).
                mapToDouble(gaussIntegrator::getPoint).boxed().collect(Collectors.toList()));

        /* weights */
        gaussIntegrationWeights.clear();
        gaussIntegrationWeights.addAll(IntStream.range(0, numIntegPoints).
                mapToDouble(gaussIntegrator::getWeight).boxed().collect(Collectors.toList()));

        /* log of weights */
        gaussIntegrationLogWeights.clear();
        gaussIntegrationLogWeights.addAll(gaussIntegrationWeights.stream().
                mapToDouble(FastMath::log).boxed().collect(Collectors.toList()));
    }

    /**
     * (advanced) Calculate a simple prior probability distribution for the allele fraction based on
     * (1) purity of the sample, and (2) maximum copy number for abnormal cells.
     * (See CNV-methods.pdf for details.)
     *
     * @param alleleFraction allele fraction to calculate the prior probability distribution on
     * @return allele fraction prior probability distribution
     */
    private double calculateAlleleFractionPriorDistribution(final double alleleFraction) {
        final double minorAlleleFraction = (alleleFraction < 0.5) ? alleleFraction : 1 - alleleFraction;
        if (minorAlleleFraction < minHetAlleleFraction) {
            return 0;
        } else if (minorAlleleFraction < breakpointHetAlleleFraction) {
            final double denom = 2 * FastMath.pow((1 - minorAlleleFraction) * minorAlleleFraction * maxCopyNumber, 2) *
                    maxAbnormalFraction * (maxAbnormalFraction - minAbnormalFraction);
            final double num = (-1 + (-1 + minorAlleleFraction * maxCopyNumber) * maxAbnormalFraction) *
                    (-1 + maxAbnormalFraction + minorAlleleFraction * (2 + (-2 + maxCopyNumber) * maxAbnormalFraction)) +
                    2 * (1 + minorAlleleFraction * (-2 + minorAlleleFraction * maxCopyNumber)) * maxAbnormalFraction *
                            FastMath.log(FastMath.abs(((1 + minorAlleleFraction * (-2 + maxCopyNumber)) * maxAbnormalFraction)) /
                                    (1 - 2 * minorAlleleFraction));
            return num / denom;
        } else { /* breakpointHetAlleleFraction < minorAlleleFraction < 1/2 */
            final double denom = 2 * FastMath.pow((1 - minorAlleleFraction) * minorAlleleFraction * maxCopyNumber, 2) *
                    maxAbnormalFraction * minAbnormalFraction * (maxAbnormalFraction - minAbnormalFraction);
            final double num = (maxAbnormalFraction - minAbnormalFraction) * (-1 + 2 * minorAlleleFraction +
                    (1 + minorAlleleFraction * (-2 + maxCopyNumber)) * (-1 + minorAlleleFraction * maxCopyNumber) *
                            maxAbnormalFraction * minAbnormalFraction) + 2 * (1 + minorAlleleFraction * (-2 + minorAlleleFraction *
                    maxCopyNumber)) * maxAbnormalFraction * minAbnormalFraction * FastMath.log(maxAbnormalFraction / minAbnormalFraction);
            return num / denom;
        }
    }
    /**
     * Calculate the allele fraction distribution function and its log on the abscissas on the integration
     * quadrature
     */
    private void initializeHetAlleleFractionPrior() {
        alleleFractionPriors.clear();
        alleleFractionPriors.addAll(gaussIntegrationAbscissas.stream()
                .mapToDouble(this::calculateAlleleFractionPriorDistribution)
                .boxed().collect(Collectors.toList()));

        /* calculate the log prior */
        alleleFractionLogPriors.clear();
        alleleFractionLogPriors.addAll(alleleFractionPriors.stream()
                .map(FastMath::log).collect(Collectors.toList()));
    }

    /**
     * Marginalize allele fraction by integrating the precomputed allele fraction prior weighted with the
     * likelihoods of the alt/ret portion of reads in the pileup
     * (see CNV-method.pdf for details)
     *
     * @param coeffs list of (alpha_k, beta_k) tuples
     * @return log likelihood
     */
    @Override
    public double getHetLogLikelihood(final Collection<? extends Pair<Double, Double>> coeffs) {
        final ArrayList<Double> refAltLogLikelihoodList = new ArrayList<>(gaussIntegrationAbscissas.size());
        refAltLogLikelihoodList.addAll(
                gaussIntegrationAbscissas.stream()
                        .map(f -> getHetLogLikelihoodFixedAlleleFraction(f, coeffs))
                        .collect(Collectors.toList()));

        final ArrayList<Double> logLikelihoodIntegrandWithPriorAndWeights = new ArrayList<>(gaussIntegrationAbscissas.size());
        logLikelihoodIntegrandWithPriorAndWeights.addAll(
                IntStream.range(0, gaussIntegrationAbscissas.size())
                        .mapToDouble(i -> refAltLogLikelihoodList.get(i) + gaussIntegrationLogWeights.get(i) +
                                alleleFractionLogPriors.get(i))
                        .boxed().collect(Collectors.toList()));

        return GATKProtectedMathUtils.logSumExp(logLikelihoodIntegrandWithPriorAndWeights);
    }
}
