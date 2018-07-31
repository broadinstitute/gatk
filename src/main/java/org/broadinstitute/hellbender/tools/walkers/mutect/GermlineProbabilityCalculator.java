package org.broadinstitute.hellbender.tools.walkers.mutect;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.ImmutableMap;
import com.google.common.primitives.Doubles;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;
import java.util.stream.IntStream;

/**
 * Created by David Benjamin on 5/4/17.
 */
public class GermlineProbabilityCalculator {

    public static Map<String, Object> getPopulationAFAnnotation(List<VariantContext> germlineResourceVariants,
                                                                final List<Allele> altAlleles,
                                                                final double afOfAllelesNotInGermlineResource) {
        final Optional<VariantContext> germlineVC = germlineResourceVariants.isEmpty() ? Optional.empty()
                : Optional.of(germlineResourceVariants.get(0));  // assume only one VC per site
        final double[] populationAlleleFrequencies = getGermlineAltAlleleFrequencies(altAlleles, germlineVC, afOfAllelesNotInGermlineResource);

        return ImmutableMap.of(GATKVCFConstants.POPULATION_AF_VCF_ATTRIBUTE, populationAlleleFrequencies);
    }

    public static double[] calculateGermlineProbabilities(final double[] populationAlleleFrequencies,
                                                          final double[] log10OddsOfGermlineHetVsSomatic,
                                                          final double[] log10OddsOfGermlineHomAltVsSomatic,
                                                          final Optional<double[]> normalLog10Odds,
                                                          final double log10PriorProbOfSomaticEvent) {
        final int nAltAlleles = populationAlleleFrequencies.length;
        final double[] normalLog10OddsOrFlat = normalLog10Odds.orElseGet(() -> Doubles.toArray(Collections.nCopies(nAltAlleles,0)));
        // note the minus sign required because Mutect has the convention that this is log odds of allele *NOT* being in the normal
        final double[] germlineLog10Posteriors = new IndexRange(0, nAltAlleles).mapToDouble(n ->
                log10PosteriorProbabilityOfGermlineVariant(-normalLog10OddsOrFlat[n], log10OddsOfGermlineHetVsSomatic[n], log10OddsOfGermlineHomAltVsSomatic[n], populationAlleleFrequencies[n], log10PriorProbOfSomaticEvent));

        return germlineLog10Posteriors;
    }

    @VisibleForTesting
    static double[] getGermlineAltAlleleFrequencies(final List<Allele> altAlleles, final Optional<VariantContext> germlineVC, final double afOfAllelesNotInGermlineResource) {
        if (germlineVC.isPresent() && germlineVC.get().hasAttribute(VCFConstants.ALLELE_FREQUENCY_KEY))  {
            final List<Double> germlineAltAFs = germlineVC.get().getAttributeAsDoubleList(VCFConstants.ALLELE_FREQUENCY_KEY, afOfAllelesNotInGermlineResource);
            return altAlleles.stream()
                    .mapToDouble(allele -> {
                        final OptionalInt germlineAltIndex = findAltAlleleIndex(germlineVC.get(), allele);
                        return germlineAltIndex.isPresent() ? germlineAltAFs.get(germlineAltIndex.getAsInt())
                                : afOfAllelesNotInGermlineResource;
                    }).toArray();
        } else {
            return Doubles.toArray(Collections.nCopies(altAlleles.size(), afOfAllelesNotInGermlineResource));
        }
    }

    /**
     * Get the alt allele index (that is 0 is the index of the first alt allele) of an {@link Allele} within a {@link VariantContext}
     */
    private static OptionalInt findAltAlleleIndex(final VariantContext vc, final Allele allele) {
        return IntStream.range(0, vc.getNAlleles() - 1)
                .filter(n -> vc.getAlternateAllele(n).basesMatch(allele))
                .findAny();
    }

    /**
     *
     * @param normalLog10Odds the log10 likelihood ratio between 1) allele being present in normal (as a diploid het or hom var
     *                        and not as an artifact) and 2) not being present.  Since likelihoods are meaningful only up to an
     *                        arbitrary constant factor, we may interpret this as the log10 likelihood that the allele exists
     *                        in the normal provided that the likelihood that it does not exist is log10(1) = 0.
     * @param log10OddsOfGermlineHetVsSomatic  the log10 likelihood ratio between 1) allele being present in tumor as a germline het
     *                                        (using information about the minor allele fraction) and 2)  being present,
     *                                 as a somatic variant as calculated by the {@code SomaticGenotypingEngine}
     * @param populationAlleleFrequency frequency of this allele in the population -- serves as a prior for germline allele
     * @param log10PriorProbOfSomaticEvent the log10 prior probability for this allele to arise de novo in the tumor
     * @return  log10 of the posterior probability that this allele exists in the normal sample
     */
    public static double log10PosteriorProbabilityOfGermlineVariant(final double normalLog10Odds,
                                                                    final double log10OddsOfGermlineHetVsSomatic,
                                                                    final double log10OddsOfGermlineHomAltVsSomatic,
                                                                    final double populationAlleleFrequency,
                                                                    final double log10PriorProbOfSomaticEvent) {
        final double log10OneMinusPriorProbSomatic = MathUtils.log10OneMinusPow10(log10PriorProbOfSomaticEvent);

        final double log10PriorGermlineHet = Math.log10(2*populationAlleleFrequency*(1-populationAlleleFrequency));
        final double log10PriorGermlineHomAlt = Math.log10( MathUtils.square(populationAlleleFrequency));
        final double log10PriorNotGermline = Math.log10(MathUtils.square(1 - populationAlleleFrequency));

        // the following are unnormalized probabilities
        final double log10ProbGermlineHet = log10PriorGermlineHet + log10OddsOfGermlineHetVsSomatic + normalLog10Odds + log10OneMinusPriorProbSomatic;
        final double log10ProbGermlineHomAlt = log10PriorGermlineHomAlt + log10OddsOfGermlineHomAltVsSomatic + normalLog10Odds + log10OneMinusPriorProbSomatic;
        final double log10ProbGermline = MathUtils.log10SumLog10(log10ProbGermlineHet, log10ProbGermlineHomAlt);
        final double log10ProbSomatic = log10PriorNotGermline + log10PriorProbOfSomaticEvent;

        return log10ProbGermline - MathUtils.log10SumLog10(log10ProbGermline, log10ProbSomatic);

    }

    private static double[] getArrayAttribute(final VariantContext vc, final String attribute) {
        return GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(vc, attribute, () -> null, -1);
    }
}
