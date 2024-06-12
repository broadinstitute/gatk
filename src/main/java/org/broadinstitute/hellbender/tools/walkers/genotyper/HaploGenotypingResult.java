package org.broadinstitute.hellbender.tools.walkers.genotyper;

import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;

public class HaploGenotypingResult {

    // all fields are log10, the posteriors are normalized
    private final GenotypingLikelihoods<Haplotype> haploGenotypeLikelihoods;
    private final GenotypingLikelihoods<Haplotype> haploGenotypePriors;
    private final GenotypingLikelihoods<Haplotype> haploGenotypePosteriors;

    public HaploGenotypingResult(final GenotypingLikelihoods<Haplotype> haploGenotypeLikelihoods,
                                 final GenotypingLikelihoods<Haplotype> haploGenotypePriors,
                                 final GenotypingLikelihoods<Haplotype> haploGenotypePosteriors) {
        this.haploGenotypeLikelihoods = haploGenotypeLikelihoods;
        this.haploGenotypePriors = haploGenotypePriors;
        this.haploGenotypePosteriors = haploGenotypePosteriors;

        Utils.validateArg(haploGenotypeLikelihoods.asListOfSamples().equals(haploGenotypePriors.asListOfSamples()), () -> "inconsistent samples");
        Utils.validateArg(haploGenotypeLikelihoods.asListOfSamples().equals(haploGenotypePosteriors.asListOfSamples()), () -> "inconsistent samples");

        Utils.validateArg(haploGenotypeLikelihoods.asListOfAlleles().equals(haploGenotypePriors.asListOfAlleles()), () -> "inconsistent haplotypes");
        Utils.validateArg(haploGenotypeLikelihoods.asListOfAlleles().equals(haploGenotypePosteriors.asListOfAlleles()), () -> "inconsistent haplotypes");

        for (int sampleIndex = 0; sampleIndex < haploGenotypeLikelihoods.numberOfSamples(); sampleIndex++) {
            Utils.validateArg(haploGenotypeLikelihoods.samplePloidy(sampleIndex) == haploGenotypePriors.samplePloidy(sampleIndex), () -> "inconsistent ploidies");
            Utils.validateArg(haploGenotypeLikelihoods.samplePloidy(sampleIndex) == haploGenotypePosteriors.samplePloidy(sampleIndex), () -> "inconsistent ploidies");
        }
    }

    public GenotypingLikelihoods<Haplotype> getLikelihoods() {
        return haploGenotypeLikelihoods;
    }

    public GenotypingLikelihoods<Haplotype> getPriors() {
        return haploGenotypePriors;
    }

    public GenotypingLikelihoods<Haplotype> getPosteriors() {
        return haploGenotypePosteriors;
    }


}
