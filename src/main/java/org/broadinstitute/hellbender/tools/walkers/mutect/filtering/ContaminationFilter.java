package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.contamination.ContaminationRecord;
import org.broadinstitute.hellbender.tools.walkers.mutect.clustering.SomaticClusteringModel;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

public class ContaminationFilter extends Mutect2AlleleFilter {
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
    public List<Double> calculateErrorProbabilityForAlleles(final VariantContext vc, final Mutect2FilteringEngine filteringEngine, ReferenceContext referenceContext) {
        // for every alt allele, a list of the depth and posterior pair
        final List<List<ImmutablePair<Integer, Double>>> depthsAndPosteriorsPerAllele = new ArrayList<>();
        new IndexRange(0, vc.getNAlleles()-1).forEach(i -> depthsAndPosteriorsPerAllele.add(new ArrayList<>()));

        for (final Genotype tumorGenotype : vc.getGenotypes()) {
            if (filteringEngine.isNormal(tumorGenotype)) {
                continue;
            }

            final double contaminationFromFile = contaminationBySample.getOrDefault(tumorGenotype.getSampleName(), defaultContamination);
            final double contamination = Math.max(0, Math.min(contaminationFromFile, 1 - EPSILON)); // handle file with contamination == 1
            final int[] ADs = tumorGenotype.getAD(); // AD is all alleles, while AF is alts only, hence the +1 offset
            final int totalAD = (int) MathUtils.sum(ADs);
            final int[] altADs = Arrays.copyOfRange(ADs, 1, ADs.length);
            // POPAF has only alt allele data
            final double[] negativeLog10AlleleFrequencies = VariantContextGetters.getAttributeAsDoubleArray(vc,
                    GATKVCFConstants.POPULATION_AF_KEY, () -> new double[]{Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY}, Double.POSITIVE_INFINITY);
            final double[] alleleFrequencies = MathUtils.applyToArray(negativeLog10AlleleFrequencies, x -> Math.pow(10,-x));

            SomaticClusteringModel model = filteringEngine.getSomaticClusteringModel();
            final double[] logSomaticLikelihoodPerAllele = Arrays.stream(altADs).mapToDouble(altCount -> model.logLikelihoodGivenSomatic(totalAD, altCount)).toArray();

            double[] singleContaminantLikelihoodPerAllele = new double[alleleFrequencies.length];
            double[] manyContaminantLikelihoodPerAllele = new double[alleleFrequencies.length];
            double[] logContaminantLikelihoodPerAllele = new double[alleleFrequencies.length];
            double[] logOddsOfRealVsContaminationPerAllele = new double[alleleFrequencies.length];
            double[] posteriorProbOfContaminationPerAllele = new double[alleleFrequencies.length];
            new IndexRange(0,alleleFrequencies.length).forEach(i -> {
                singleContaminantLikelihoodPerAllele[i] = 2 * alleleFrequencies[i] * (1 - alleleFrequencies[i]) * MathUtils.binomialProbability(totalAD,  altADs[i], contamination /2)
                        + MathUtils.square(alleleFrequencies[i]) * MathUtils.binomialProbability(totalAD,  altADs[i], contamination);
                manyContaminantLikelihoodPerAllele[i] = MathUtils.binomialProbability(totalAD, altADs[i], contamination * alleleFrequencies[i]);
                logContaminantLikelihoodPerAllele[i] = Math.log(Math.max(singleContaminantLikelihoodPerAllele[i], manyContaminantLikelihoodPerAllele[i]));
                logOddsOfRealVsContaminationPerAllele[i] = logSomaticLikelihoodPerAllele[i] - logContaminantLikelihoodPerAllele[i];
            });

            new IndexRange(0,alleleFrequencies.length).forEach(i -> {
                posteriorProbOfContaminationPerAllele[i] = filteringEngine.posteriorProbabilityOfError(vc, logOddsOfRealVsContaminationPerAllele[i], i);
                depthsAndPosteriorsPerAllele.get(i).add(ImmutablePair.of(altADs[i], posteriorProbOfContaminationPerAllele[i]));
            });

        }

        return depthsAndPosteriorsPerAllele.stream().map(alleleData -> alleleData.isEmpty() ? Double.NaN : weightedMedianPosteriorProbability(alleleData)).collect(Collectors.toList());
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
