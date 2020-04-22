package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.samtools.util.OverlapDetector;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.mutable.MutableDouble;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.math3.distribution.BinomialDistribution;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.contamination.MinorAlleleFractionRecord;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.NaturalLogUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

public class GermlineFilter extends Mutect2VariantFilter {
    private static final double MIN_ALLELE_FRACTION_FOR_GERMLINE_HOM_ALT = 0.9;
    // numerical precision safeguard in case of bad JVMs inverting the negative log population allele frequency
    private static final double EPSILON = 1.0e-10;

    private final Map<String, OverlapDetector<MinorAlleleFractionRecord>> tumorSegments;

    public GermlineFilter(final List<File> tumorSegmentationTables) {
        tumorSegments = tumorSegmentationTables.stream()
                .map(MinorAlleleFractionRecord::readFromFile)
                .collect(Collectors.toMap(ImmutablePair::getLeft, p -> OverlapDetector.create(p.getRight())));
    }

    @Override
    public ErrorType errorType() { return ErrorType.NON_SOMATIC; }

    @Override
    public double calculateErrorProbability(final VariantContext vc, final Mutect2FilteringEngine filteringEngine, ReferenceContext referenceContext) {
        final double[] somaticLogOdds = Mutect2FilteringEngine.getTumorLogOdds(vc);
        final int maxLodIndex = MathUtils.maxElementIndex(somaticLogOdds);

        final Optional<double[]> normalLogOdds = vc.hasAttribute(GATKVCFConstants.NORMAL_LOG_10_ODDS_KEY) ?
                Optional.of(MathUtils.applyToArrayInPlace(VariantContextGetters.getAttributeAsDoubleArray(vc, GATKVCFConstants.NORMAL_LOG_10_ODDS_KEY), MathUtils::log10ToLog)) : Optional.empty();
        final double[] negativeLog10AlleleFrequencies = VariantContextGetters.getAttributeAsDoubleArray(vc, GATKVCFConstants.POPULATION_AF_KEY);
        final double populationAF = Math.pow(10, -negativeLog10AlleleFrequencies[maxLodIndex]);

        if (populationAF < EPSILON) {
            return 0;
        } else if (populationAF > 1 - EPSILON) {
            return 1;
        }

        // note that this includes the ref
        final int[] alleleCounts = filteringEngine.sumADsOverSamples(vc, true, false);
        final int totalCount = (int) MathUtils.sum(alleleCounts);
        if (totalCount == 0) {  // this could come up in GGA mode
            return 0;
        }
        final int refCount = alleleCounts[0];
        final int altCount = alleleCounts[maxLodIndex + 1];
        final double altAlleleFraction = filteringEngine.weightedAverageOfTumorAFs(vc)[maxLodIndex];

        final double maf = computeMinorAlleleFraction(vc, filteringEngine, alleleCounts);

        // sum of alt minor and alt major possibilities
        final double logGermlineLikelihood = NaturalLogUtils.LOG_ONE_HALF + NaturalLogUtils.logSumExp(
                new BinomialDistribution(null, totalCount, maf).logProbability(altCount),
                new BinomialDistribution(null, totalCount, 1 - maf).logProbability(altCount));

        final double logSomaticLikelihood = filteringEngine.getSomaticClusteringModel().logLikelihoodGivenSomatic(totalCount, altCount);
        // this is \chi in the docs, the correction factor for tumor likelihoods if forced to have maf or 1 - maf
        // as the allele fraction
        final double logOddsOfGermlineHetVsSomatic = logGermlineLikelihood - logSomaticLikelihood;

        // see docs -- basically the tumor likelihood for a germline hom alt is approximately equal to the somatic likelihood
        // as long as the allele fraction is high
        final double logOddsOfGermlineHomAltVsSomatic = altAlleleFraction < MIN_ALLELE_FRACTION_FOR_GERMLINE_HOM_ALT ? Double.NEGATIVE_INFINITY : 0;

        final double normalLod = normalLogOdds.isPresent() ? normalLogOdds.get()[maxLodIndex] : 0;
        // note the minus sign required because Mutect has the convention that this is log odds of allele *NOT* being in the normal
        return germlineProbability(-normalLod, logOddsOfGermlineHetVsSomatic, logOddsOfGermlineHomAltVsSomatic,
                populationAF, filteringEngine.getLogSomaticPrior(vc, maxLodIndex));
    }

    /**
     *
     * @param normalLogOdds log odds of allele in normal as a diploid het or hom var versus not being present in normal.
     * @param logOddsOfGermlineHetVsSomatic  log odds of allele being present in tumor as a germline het versus
     *                                         being present as a somatic variant
     * @param populationAF      frequency of this allele in the population
     * @param logPriorSomatic log prior probability for this allele to be a somatic mutation
     * @return                  probability that this allele exists in the normal sample
     */
    public static double germlineProbability(final double normalLogOdds,
                                             final double logOddsOfGermlineHetVsSomatic,
                                             final double logOddsOfGermlineHomAltVsSomatic,
                                             final double populationAF,
                                             final double logPriorSomatic) {

        final double logPriorNotSomatic = NaturalLogUtils.log1mexp(logPriorSomatic);
        final double logPriorGermlineHet = Math.log(2*populationAF*(1-populationAF));
        final double logPriorGermlineHomAlt = Math.log( MathUtils.square(populationAF));
        final double logPriorNotGermline = Math.log(MathUtils.square(1 - populationAF));

        // the following are unnormalized probabilities
        final double logProbGermlineHet = logPriorGermlineHet + logOddsOfGermlineHetVsSomatic + normalLogOdds + logPriorNotSomatic;
        final double logProbGermlineHomAlt = logPriorGermlineHomAlt + logOddsOfGermlineHomAltVsSomatic + normalLogOdds + logPriorNotSomatic;
        final double logProbGermline = NaturalLogUtils.logSumExp(logProbGermlineHet, logProbGermlineHomAlt);
        final double logProbSomatic = logPriorNotGermline + logPriorSomatic;

        return NaturalLogUtils.normalizeLog(new double[] {logProbGermline, logProbSomatic}, false, true)[0];
    }

    private double computeMinorAlleleFraction(final VariantContext vc, final Mutect2FilteringEngine filteringEngine, final int[] alleleCounts) {
        final MutableDouble weightedSumOfMafs = new MutableDouble(0);
        vc.getGenotypes().stream().filter(filteringEngine::isTumor).forEach(tumorGenotype -> {
            final String sample = tumorGenotype.getSampleName();
            final List<MinorAlleleFractionRecord> segments = tumorSegments.containsKey(sample) ? tumorSegments.get(sample).getOverlaps(vc).stream().collect(Collectors.toList())
                    : Collections.emptyList();

            // minor allele fraction -- we abbreviate the name to make the formulas below less cumbersome
            final double maf = segments.isEmpty() ? 0.5 : segments.get(0).getMinorAlleleFraction();

            weightedSumOfMafs.add(maf * MathUtils.sum(tumorGenotype.getAD()));
        });


        // weighted average of sample minor allele fractions.  This is the expected allele fraction of a germline het in the aggregated read counts
        return weightedSumOfMafs.getValue() / MathUtils.sum(alleleCounts);
    }

    @Override
    public String filterName() {
        return GATKVCFConstants.GERMLINE_RISK_FILTER_NAME;
    }

    @Override
    public Optional<String> phredScaledPosteriorAnnotationName() {
        return Optional.of(GATKVCFConstants.GERMLINE_QUAL_KEY);
    }

    @Override
    protected List<String> requiredInfoAnnotations() {
        return Arrays.asList(GATKVCFConstants.TUMOR_LOG_10_ODDS_KEY, GATKVCFConstants.POPULATION_AF_KEY);
    }

}
