package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.contamination.ContaminationRecord;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

public class ContaminationFilter extends Mutect2VariantFilter {
    private final Map<String, Double> contaminationBySample;
    private final double defaultContamination;
    private final double EPSILON = 1.0e-10;

    public ContaminationFilter(final List<File> contaminationTables, final double contaminationEstimate) {
        contaminationBySample = contaminationTables.stream()
                .map(file -> ContaminationRecord.readFromFile(file).get(0))
                .collect(Collectors.toMap(rec -> rec.getSample(), rec -> rec.getContamination()));

        defaultContamination = contaminationEstimate;
    }

    @Override
    public ErrorType errorType() { return ErrorType.NON_SOMATIC; }

    @Override
    public double calculateErrorProbability(final VariantContext vc, final Mutect2FilteringEngine filteringEngine, ReferenceContext referenceContext) {
        final List<ImmutablePair<Integer, Double>> depthsAndPosteriors = new ArrayList<>();

        for (final Genotype tumorGenotype : vc.getGenotypes()) {
            if (filteringEngine.isNormal(tumorGenotype)) {
                continue;
            }

            final double contaminationFromFile = contaminationBySample.getOrDefault(tumorGenotype.getSampleName(), defaultContamination);
            final double contamination = Math.max(0, Math.min(contaminationFromFile, 1 - EPSILON)); // handle file with contamination == 1
            final double[] alleleFractions = VariantContextGetters.getAttributeAsDoubleArray(tumorGenotype, VCFConstants.ALLELE_FREQUENCY_KEY,
                    () -> new double[] {1.0}, 1.0);
            final int maxFractionIndex = MathUtils.maxElementIndex(alleleFractions);
            final int[] ADs = tumorGenotype.getAD();
            final int altCount = ADs[maxFractionIndex + 1];   // AD is all alleles, while AF is alts only, hence the +1 offset
            final int depth = (int) MathUtils.sum(ADs);
            final double[] negativeLog10AlleleFrequencies = VariantContextGetters.getAttributeAsDoubleArray(vc,
                    GATKVCFConstants.POPULATION_AF_KEY, () -> new double[]{Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY}, Double.POSITIVE_INFINITY);
            final double alleleFrequency = MathUtils.applyToArray(negativeLog10AlleleFrequencies, x -> Math.pow(10,-x))[maxFractionIndex];

            final double logSomaticLikelihood = filteringEngine.getSomaticClusteringModel().logLikelihoodGivenSomatic(depth, altCount);

            final double singleContaminantLikelihood = 2 * alleleFrequency * (1 - alleleFrequency) * MathUtils.binomialProbability(depth, altCount, contamination /2)
                    + MathUtils.square(alleleFrequency) * MathUtils.binomialProbability(depth, altCount, contamination);
            final double manyContaminantLikelihood = MathUtils.binomialProbability(depth, altCount, contamination * alleleFrequency);
            final double logContaminantLikelihood = Math.log(Math.max(singleContaminantLikelihood, manyContaminantLikelihood));
            final double logOddsOfRealVsContamination = logSomaticLikelihood - logContaminantLikelihood;
            final double posteriorProbOfContamination = filteringEngine.posteriorProbabilityOfError(vc, logOddsOfRealVsContamination, maxFractionIndex);

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
        return Optional.of(GATKVCFConstants.CONTAMINATION_QUAL_KEY);
    }

    @Override
    protected List<String> requiredAnnotations() { return Collections.singletonList(GATKVCFConstants.POPULATION_AF_KEY); }
}
