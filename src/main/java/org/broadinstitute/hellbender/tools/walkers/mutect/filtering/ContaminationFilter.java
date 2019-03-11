package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.broadinstitute.hellbender.tools.walkers.contamination.ContaminationRecord;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

public class ContaminationFilter extends Mutect2VariantFilter {
    private final Map<String, Double> contaminationBySample;
    private final double defaultContamination;

    public ContaminationFilter(final List<File> contaminationTables, final double contaminationEstimate) {
        contaminationBySample = contaminationTables.stream()
                .map(file -> ContaminationRecord.readFromFile(file).get(0))
                .collect(Collectors.toMap(rec -> rec.getSample(), rec -> rec.getContamination()));

        defaultContamination = contaminationEstimate;
    }

    @Override
    public ErrorType errorType() { return ErrorType.NON_SOMATIC; }

    @Override
    public double calculateErrorProbability(final VariantContext vc, final Mutect2FilteringEngine filteringEngine) {
        final List<ImmutablePair<Integer, Double>> depthsAndPosteriors = new ArrayList<>();

        for (final Genotype tumorGenotype : vc.getGenotypes()) {
            if (filteringEngine.isNormal(tumorGenotype)) {
                continue;
            }

            final double contamination = contaminationBySample.getOrDefault(tumorGenotype.getSampleName(), defaultContamination);

            final double[] alleleFractions = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(tumorGenotype, VCFConstants.ALLELE_FREQUENCY_KEY,
                    () -> new double[] {1.0}, 1.0);
            final int maxFractionIndex = MathUtils.maxElementIndex(alleleFractions);
            final int[] ADs = tumorGenotype.getAD();
            final int altCount = ADs[maxFractionIndex + 1];   // AD is all alleles, while AF is alts only, hence the +1 offset
            final int depth = (int) MathUtils.sum(ADs);
            final double[] negativeLog10AlleleFrequencies = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(vc,
                    GATKVCFConstants.POPULATION_AF_VCF_ATTRIBUTE, () -> new double[]{Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY}, Double.POSITIVE_INFINITY);
            final double alleleFrequency = MathUtils.applyToArray(negativeLog10AlleleFrequencies, x -> Math.pow(10,-x))[maxFractionIndex];

            final double log10SomaticLikelihood = filteringEngine.getSomaticClusteringModel().log10LikelihoodGivenSomatic(depth, altCount);

            final double singleContaminantLikelihood = 2 * alleleFrequency * (1 - alleleFrequency) * MathUtils.binomialProbability(depth, altCount, contamination /2)
                    + MathUtils.square(alleleFrequency) * MathUtils.binomialProbability(depth, altCount, contamination);
            final double manyContaminantLikelihood = MathUtils.binomialProbability(depth, altCount, contamination * alleleFrequency);
            final double log10ContaminantLikelihood = Math.log10(Math.max(singleContaminantLikelihood, manyContaminantLikelihood));
            final double log10OddsOfRealVsContamination = log10SomaticLikelihood - log10ContaminantLikelihood;
            final double posteriorProbOfContamination = filteringEngine.posteriorProbabilityOfError(vc, log10OddsOfRealVsContamination, maxFractionIndex);

            depthsAndPosteriors.add(ImmutablePair.of(altCount, posteriorProbOfContamination));
        }

        return weightedMedianPosteriorProbability(depthsAndPosteriors);
    }

    @Override
    public String filterName() {
        return GATKVCFConstants.CONTAMINATION_FILTER_NAME;
    }

    @Override
    public Optional<String> phredScaledPosteriorAnnotationName() {
        return Optional.of(GATKVCFConstants.CONTAMINATION_QUAL_ATTRIBUTE);
    }

    @Override
    protected List<String> requiredAnnotations() { return Collections.singletonList(GATKVCFConstants.POPULATION_AF_VCF_ATTRIBUTE); }
}
