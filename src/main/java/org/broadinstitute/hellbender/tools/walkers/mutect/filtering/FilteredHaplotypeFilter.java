package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;

public class FilteredHaplotypeFilter extends Mutect2VariantFilter {
    private final double maxIntraHaplotypeDistance;

    // for each pgt + pid phasing string, a list of loci-error probability pairs
    private Map<String, List<Pair<Integer, Double>>> accumulatingPhasedProbabilities = new HashMap<>();

    private Map<String, List<Pair<Integer, Double>>> phasedProbabilities = new HashMap<>();

    public FilteredHaplotypeFilter(final double maxIntraHaplotypeDistance) {
        this.maxIntraHaplotypeDistance = maxIntraHaplotypeDistance;
    }

    @Override
    public ErrorType errorType() { return ErrorType.ARTIFACT; }

    @Override
    public double calculateErrorProbability(final VariantContext vc, final Mutect2FilteringEngine filteringEngine) {
        // use phasing of tumor genotype with greatest allele fraction
        final Genotype tumorGenotype = vc.getGenotypes().stream().filter(filteringEngine::isTumor)
                .max(Comparator.comparingDouble(g -> MathUtils.arrayMax(GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(g, VCFConstants.ALLELE_FREQUENCY_KEY,
                        () -> new double[] {0.0}, 0.0)))).get();

        final Optional<String> phasingString = makePhasingString(tumorGenotype);
        if (!phasingString.isPresent()) {
            return 0.0;
        }

        // note that we use the learned probabilities from the previous pass
        final List<Pair<Integer, Double>> phasedProbs = phasedProbabilities.get(phasingString.get());

        if (phasedProbs == null) {
            return 0.0;
        }

        return phasedProbs.stream()
                .filter(pair -> Math.abs(pair.getLeft() - vc.getStart()) <= maxIntraHaplotypeDistance)
                .mapToDouble(Pair::getRight)
                .max().orElse(0.0);
    }

    @Override
    protected void accumulateDataForLearning(final VariantContext vc, final ErrorProbabilities errorProbabilities, final Mutect2FilteringEngine filteringEngine) {
        // we record the maximum non-sequencing artifact that is not this filter itself
        final double artifactProbability = errorProbabilities.getProbabilitiesByFilter().entrySet().stream()
                .filter(e -> e.getKey().errorType() != ErrorType.SEQUENCING)
                .filter(e -> !e.getKey().filterName().equals(filterName()))
                .mapToDouble(e -> e.getValue())
                .max().orElse(0.0);

        for (final Genotype tumorGenotype : vc.getGenotypes()) {
            if (!filteringEngine.isTumor(tumorGenotype)) {
                continue;
            }

            final Optional<String> phasingString = makePhasingString(tumorGenotype);

            if (!phasingString.isPresent()) {
                continue;
            }

            if (!accumulatingPhasedProbabilities.containsKey(phasingString.get())) {
                accumulatingPhasedProbabilities.put(phasingString.get(), new ArrayList<>());
            }

            accumulatingPhasedProbabilities.get(phasingString.get()).add(ImmutablePair.of(vc.getStart(), artifactProbability));
        }
    }


    @Override
    protected void clearAccumulatedData() {
        accumulatingPhasedProbabilities = new HashMap<>();
    }

    @Override
    protected void learnParameters() {
        // move the accumulating probabilities to the learned probabiities
        phasedProbabilities = accumulatingPhasedProbabilities;
    }

    @Override
    public String filterName() {
        return GATKVCFConstants.BAD_HAPLOTYPE_FILTER_NAME;
    }

    @Override
    protected List<String> requiredAnnotations() { return Collections.emptyList(); }

    @Override
    public Optional<String> phredScaledPosteriorAnnotationName() { return Optional.empty(); }

    // concatenate the PGT and PID strings, if present
    private static Optional<String> makePhasingString(final Genotype genotype) {
        final String pgt = (String) genotype.getExtendedAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY, null);
        final String pid = (String) genotype.getExtendedAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY, null);
        return (pgt == null || pid == null) ? Optional.empty() : Optional.of(pgt + pid);
    }
}
