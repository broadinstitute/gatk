package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyState;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.mcmc.Parameter;
import org.broadinstitute.hellbender.utils.mcmc.ParameterizedState;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Represents parameters for the {@link TumorHeterogeneity} model of a mixture of subclones with copy-number variation.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class TumorHeterogeneityState extends ParameterizedState<TumorHeterogeneityParameter> {
    public TumorHeterogeneityState(final double concentration,
                                   final double copyRatioNormalization,
                                   final double copyRatioNoiseConstant,
                                   final double initialPloidy,
                                   final double ploidy,
                                   final PopulationMixture populationMixture) {
        super(Arrays.asList(
                new Parameter<>(TumorHeterogeneityParameter.CONCENTRATION, concentration),
                new Parameter<>(TumorHeterogeneityParameter.COPY_RATIO_NORMALIZATION, copyRatioNormalization),
                new Parameter<>(TumorHeterogeneityParameter.COPY_RATIO_NOISE_CONSTANT, copyRatioNoiseConstant),
                new Parameter<>(TumorHeterogeneityParameter.INITIAL_PLOIDY, initialPloidy),
                new Parameter<>(TumorHeterogeneityParameter.PLOIDY, ploidy),
                new Parameter<>(TumorHeterogeneityParameter.POPULATION_MIXTURE,
                        new PopulationMixture(populationMixture.populationFractions(), populationMixture.variantProfileCollection(), populationMixture.normalPloidyState()))));
        Utils.validateArg(concentration > 0., "Concentration must be positive.");
        Utils.validateArg(copyRatioNormalization > 0., "Copy-ratio normalization must be positive.");
        Utils.validateArg(copyRatioNoiseConstant >= 0., "Copy-ratio noise constant must be non-negative.");
        Utils.validateArg(initialPloidy >= 0, "Initial ploidy must be non-negative.");
        Utils.validateArg(ploidy >= 0, "Ploidy must be non-negative.");
        Utils.nonNull(populationMixture);
    }

    public double concentration() {
        return get(TumorHeterogeneityParameter.CONCENTRATION, Double.class);
    }

    public double copyRatioNormalization() {
        return get(TumorHeterogeneityParameter.COPY_RATIO_NORMALIZATION, Double.class);
    }

    public double copyRatioNoiseConstant() {
        return get(TumorHeterogeneityParameter.COPY_RATIO_NOISE_CONSTANT, Double.class);
    }

    public double initialPloidy() {
        return get(TumorHeterogeneityParameter.INITIAL_PLOIDY, Double.class);
    }

    public double ploidy() {
        return get(TumorHeterogeneityParameter.PLOIDY, Double.class);
    }

    public PopulationMixture populationMixture() {
        return get(TumorHeterogeneityParameter.POPULATION_MIXTURE, PopulationMixture.class);
    }

    /**
     * Initialize state with evenly distributed population fractions and normal variant profiles.
     */
    static TumorHeterogeneityState initializeNormalState(final TumorHeterogeneityPriorCollection priors,
                                                         final int numSegments,
                                                         final int numPopulations) {
        final double concentration = priors.concentrationPriorHyperparameterValues().getAlpha() / priors.concentrationPriorHyperparameterValues().getBeta();
        final double copyRatioNormalization = priors.copyRatioNormalizationPriorHyperparameterValues().getAlpha() / priors.copyRatioNormalizationPriorHyperparameterValues().getBeta();
        final double copyRatioNoiseConstant = priors.copyRatioNoiseConstantPriorHyperparameterValues().getAlpha() / priors.copyRatioNoiseConstantPriorHyperparameterValues().getBeta();
        //initialize population fractions to be evenly distributed
        final PopulationMixture.PopulationFractions populationFractions =
                new PopulationMixture.PopulationFractions(Collections.nCopies(numPopulations, 1. / numPopulations));
        //initialize variant profiles to normal
        final PloidyState normalPloidyState = priors.normalPloidyState();
        final double normalPloidy = normalPloidyState.total();
        final int numVariantPopulations = numPopulations - 1;
        final PopulationMixture.VariantProfileCollection variantProfileCollection =
                initializeNormalProfiles(numVariantPopulations, numSegments, normalPloidyState);
        final PopulationMixture populationMixture = new PopulationMixture(populationFractions, variantProfileCollection, normalPloidyState);
        return new TumorHeterogeneityState(
                concentration, copyRatioNormalization, copyRatioNoiseConstant, normalPloidy, normalPloidy, populationMixture);
    }

    /**
     * Initialize state from a modeller run that assumes one clonal population and one normal population.
     */
    static TumorHeterogeneityState initializeFromClonalState(final TumorHeterogeneityPriorCollection priors,
                                                             final TumorHeterogeneityState clonalState,
                                                             final int maxNumPopulations) {
        //check that the state is clonal
        Utils.validateArg(clonalState.populationMixture().numPopulations() == 2,
                "Clonal state must have two populations (clone + normal).");

        //initialize global parameters
        final double concentration = priors.concentrationPriorHyperparameterValues().getAlpha() / priors.concentrationPriorHyperparameterValues().getBeta();
        final double copyRatioNormalization = clonalState.copyRatioNormalization();
        final double copyRatioNoiseConstant = clonalState.copyRatioNoiseConstant();

        //initialize normal fraction to clonal result
        final double normalFraction = clonalState.populationMixture().populationFractions().normalFraction();
        final double tumorFraction = clonalState.populationMixture().populationFractions().tumorFraction();

        //build new population fractions and variant-profile collection with additional variant populations
        //split clonal fraction evenly among new populations and initialize with clonal profile
        final List<Double> initialFractions = new ArrayList<>();
        final List<PopulationMixture.VariantProfile> initialVariantProfiles = new ArrayList<>();
        final int numVariantPopulations = maxNumPopulations - 1;
        final PopulationMixture.VariantProfile clonalProfile = clonalState.populationMixture().variantProfileCollection().get(0);
        //initialize additional variant profiles
        for (int i = 0; i < numVariantPopulations; i++) {
            initialFractions.add(tumorFraction / numVariantPopulations);
            initialVariantProfiles.add(new PopulationMixture.VariantProfile(clonalProfile));
        }
        //add normal population fraction
        initialFractions.add(normalFraction);
        final PopulationMixture populationMixture = new PopulationMixture(initialFractions, initialVariantProfiles, priors.normalPloidyState());

        final double ploidy = clonalState.ploidy();

        return new TumorHeterogeneityState(
                concentration, copyRatioNormalization, copyRatioNoiseConstant, ploidy, ploidy, populationMixture);
    }

    /**
     * Initialize variant profiles to normal.
     */
    static PopulationMixture.VariantProfileCollection initializeNormalProfiles(final int numVariantPopulations,
                                                                               final int numSegments,
                                                                               final PloidyState normalPloidyState) {
        return new PopulationMixture.VariantProfileCollection(
                IntStream.range(0, numVariantPopulations).boxed().map(i -> initializeNormalProfile(numSegments, normalPloidyState)).collect(Collectors.toList()));
    }

    /**
     * Initialize a variant profile to normal.
     */
    private static PopulationMixture.VariantProfile initializeNormalProfile(final int numSegments,
                                                                            final PloidyState normalPloidyState) {
        final PopulationMixture.VariantProfile ploidyStateIndicators =
                new PopulationMixture.VariantProfile(Collections.nCopies(numSegments, normalPloidyState));
        return new PopulationMixture.VariantProfile(ploidyStateIndicators);
    }
}
