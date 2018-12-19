package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.samtools.util.OverlapDetector;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.mutable.MutableDouble;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.broadinstitute.hellbender.tools.walkers.contamination.MinorAlleleFractionRecord;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

public class GermlineFilter extends Mutect2VariantFilter {
    private static final double MIN_ALLELE_FRACTION_FOR_GERMLINE_HOM_ALT = 0.9;
    // numerical precision safeguard in case of bad JVMs inverting the negative log-10 population allele frequency
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
    public double calculateErrorProbability(final VariantContext vc, final Mutect2FilteringEngine filteringEngine) {
        final double[] somaticLog10Odds = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(vc, GATKVCFConstants.TUMOR_LOD_KEY);
        final int maxLodIndex = MathUtils.maxElementIndex(somaticLog10Odds);

        final Optional<double[]> normalLog10Odds = vc.hasAttribute(GATKVCFConstants.NORMAL_LOD_KEY) ?
                Optional.of(GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(vc, GATKVCFConstants.NORMAL_LOD_KEY)) : Optional.empty();
        final double[] negativeLog10AlleleFrequencies = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(vc, GATKVCFConstants.POPULATION_AF_VCF_ATTRIBUTE);
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
        final double log10GermlineLikelihood = MathUtils.LOG10_ONE_HALF + MathUtils.log10SumLog10(
                MathUtils.log10BinomialProbability(totalCount, altCount, Math.log10(maf)),
                MathUtils.log10BinomialProbability(totalCount, altCount, Math.log10(1 - maf)));

        final double log10SomaticLikelihood = filteringEngine.getSomaticClusteringModel().log10LikelihoodGivenSomatic(totalCount, altCount);
        // this is \chi in the docs, the correction factor for tumor likelihoods if forced to have maf or 1 - maf
        // as the allele fraction
        final double log10OddsOfGermlineHetVsSomatic = log10GermlineLikelihood - log10SomaticLikelihood;

        // see docs -- basically the tumor likelihood for a germline hom alt is approximately equal to the somatic likelihood
        // as long as the allele fraction is high
        final double log10OddsOfGermlineHomAltVsSomatic = altAlleleFraction < MIN_ALLELE_FRACTION_FOR_GERMLINE_HOM_ALT ? Double.NEGATIVE_INFINITY : 0;

        final double normalLod = normalLog10Odds.isPresent() ? normalLog10Odds.get()[maxLodIndex] : 0;
        // note the minus sign required because Mutect has the convention that this is log odds of allele *NOT* being in the normal
        return germlineProbability(-normalLod, log10OddsOfGermlineHetVsSomatic, log10OddsOfGermlineHomAltVsSomatic,
                populationAF, filteringEngine.getLog10PriorOfSomaticVariant(vc, maxLodIndex));
    }

    /**
     *
     * @param normalLog10Odds log10 odds of allele in normal as a diploid het or hom var versus not being present in normal.
     * @param log10OddsOfGermlineHetVsSomatic  log10 odds of allele being present in tumor as a germline het versus
     *                                         being present as a somatic variant
     * @param populationAF      frequency of this allele in the population
     * @param log10PriorSomatic log10 prior probability for this allele to be a somatic mutation
     * @return                  probability that this allele exists in the normal sample
     */
    public static double germlineProbability(final double normalLog10Odds,
                                             final double log10OddsOfGermlineHetVsSomatic,
                                             final double log10OddsOfGermlineHomAltVsSomatic,
                                             final double populationAF,
                                             final double log10PriorSomatic) {

        final double log10PriorNotSomatic = MathUtils.log10OneMinusPow10(log10PriorSomatic);

        final double log10PriorGermlineHet = Math.log10(2*populationAF*(1-populationAF));
        final double log10PriorGermlineHomAlt = Math.log10( MathUtils.square(populationAF));
        final double log10PriorNotGermline = Math.log10(MathUtils.square(1 - populationAF));

        // the following are unnormalized probabilities
        final double log10ProbGermlineHet = log10PriorGermlineHet + log10OddsOfGermlineHetVsSomatic + normalLog10Odds + log10PriorNotSomatic;
        final double log10ProbGermlineHomAlt = log10PriorGermlineHomAlt + log10OddsOfGermlineHomAltVsSomatic + normalLog10Odds + log10PriorNotSomatic;
        final double log10ProbGermline = MathUtils.log10SumLog10(log10ProbGermlineHet, log10ProbGermlineHomAlt);
        final double log10ProbSomatic = log10PriorNotGermline + log10PriorSomatic;

        return MathUtils.normalizeLog10(new double[] {log10ProbGermline, log10ProbSomatic}, false, true)[0];
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
        return Optional.of(GATKVCFConstants.GERMLINE_QUAL_VCF_ATTRIBUTE);
    }

    @Override
    protected List<String> requiredAnnotations() {
        return Arrays.asList(GATKVCFConstants.TUMOR_LOD_KEY, GATKVCFConstants.POPULATION_AF_VCF_ATTRIBUTE);
    }

}
