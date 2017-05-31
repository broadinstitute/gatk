package org.broadinstitute.hellbender.tools.walkers.mutect;

import com.google.cloud.dataflow.sdk.repackaged.com.google.common.primitives.Doubles;
import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.ImmutableMap;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.MathUtils;

import java.util.*;
import java.util.stream.IntStream;

/**
 * Created by David Benjamin on 5/4/17.
 */
public class GermlineProbabilityCalculator {

    public static final String POPULATION_AF_VCF_ATTRIBUTE = "POP_AF";
    public static final String GERMLINE_POSTERIORS_VCF_ATTRIBUTE = "P_GERMLINE";

    public static Map<String, Object> calculateAnnotations(List<VariantContext> germlineResourceVariants,
                                                           final List<Allele> altAlleles,
                                                           final double[] tumorLog10Odds,
                                                           final Optional<double[]> normalLog10Odds,
                                                           final double afOfAllelesNotInGermlineResource,
                                                           final double log10PriorProbOfSomaticEvent) {
        final double[] normalLog10OddsOrFlat = normalLog10Odds.orElseGet(() -> MathUtils.applyToArray(tumorLog10Odds, x -> 0));
        final Optional<VariantContext> germlineVC = germlineResourceVariants.isEmpty() ? Optional.empty()
                : Optional.of(germlineResourceVariants.get(0));  // assume only one VC per site
        final double[] populationAlleleFrequencies = getGermlineAltAlleleFrequencies(altAlleles, germlineVC, afOfAllelesNotInGermlineResource);

        // note the minus sign required because Mutect has the convention that this is log odds of allele *NOT* being in the normal
        final double[] germlineLog10Posteriors = new IndexRange(0, altAlleles.size()).mapToDouble(n ->
                log10PosteriorProbabilityOfGermlineVariant(-normalLog10OddsOrFlat[n], tumorLog10Odds[n], populationAlleleFrequencies[n], log10PriorProbOfSomaticEvent));

        return ImmutableMap.of(POPULATION_AF_VCF_ATTRIBUTE, populationAlleleFrequencies,
                GERMLINE_POSTERIORS_VCF_ATTRIBUTE, germlineLog10Posteriors);
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
     * @param tumorLog10Odds  the log10 likelihood ratio between 1) allele being present in tumor and 2) not being present.
     *                        Since likelihoods are meaningful only up to an arbitrary constant factor, we may interpret this
     *                        as the log10 likelihood that the allele exists in the tumor provided that the likelihood that
     *                        it does not exist is log10(1) = 0.
     * @param populationAlleleFrequency frequency of this allele in the population -- serves as a prior for germline allele
     * @param log10PriorProbOfSomaticEvent the log10 prior probability for this allele to arise de novo in the tumor
     * @return  log10 of the posterior probability that this allele exists in the normal sample
     */
    public static double log10PosteriorProbabilityOfGermlineVariant(final double normalLog10Odds, final double tumorLog10Odds,
                                                                    final double populationAlleleFrequency,
                                                                    final double log10PriorProbOfSomaticEvent) {
        final double log10OneMinusPriorProbSomatic = MathUtils.log10OneMinusPow10(log10PriorProbOfSomaticEvent);

        // the following is log10(p_het + p_homvar)
        final double log10PriorInNormal = Math.log10(2*populationAlleleFrequency*(1-populationAlleleFrequency + MathUtils.square(populationAlleleFrequency)));

        // log10(p_homref)
        final double log10PriorNotInNormal = Math.log10(MathUtils.square(1 - populationAlleleFrequency));

        final double log10UnnormalizedProbInBoth = log10PriorInNormal + normalLog10Odds + tumorLog10Odds + log10OneMinusPriorProbSomatic;
        final double log10UnnormalizedProbTumorOnly = log10PriorNotInNormal + tumorLog10Odds + log10PriorProbOfSomaticEvent;
        final double log10UnnormalizedProbInNeither = log10PriorNotInNormal + log10OneMinusPriorProbSomatic;
        // we neglect the probability that a normal allele is not present in the normal

        // we want log10(10^log10p_both /(10^log10p_both + 10^log10p_tumor + 10^log10p_neither)), which comes out ot
        // log10p_both - log10SumLog10(log10p_both, log10p_tumor, log10p_neither)
        return log10UnnormalizedProbInBoth
                - MathUtils.approximateLog10SumLog10(log10UnnormalizedProbInBoth, log10UnnormalizedProbTumorOnly, log10UnnormalizedProbInNeither);
    }

    private static double[] getArrayAttribute(final VariantContext vc, final String attribute) {
        return GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(vc, attribute, () -> null, -1);
    }
}
