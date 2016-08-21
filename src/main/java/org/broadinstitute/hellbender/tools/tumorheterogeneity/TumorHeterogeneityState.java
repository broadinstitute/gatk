package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyState;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.mcmc.Parameter;
import org.broadinstitute.hellbender.utils.mcmc.ParameterizedState;

import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class TumorHeterogeneityState extends ParameterizedState<TumorHeterogeneityParameter> {
    private static final int NUM_POPULATIONS_CLONAL = 2;

    private final TumorHeterogeneityPriorCollection priors;

    public TumorHeterogeneityState(final double concentration,
                                   final double copyRatioNoiseFloor,
                                   final double copyRatioNoiseFactor,
                                   final double minorAlleleFractionNoiseFactor,
                                   final PopulationMixture populationMixture,
                                   final TumorHeterogeneityPriorCollection priors) {
        super(Arrays.asList(
                new Parameter<>(TumorHeterogeneityParameter.CONCENTRATION, concentration),
                new Parameter<>(TumorHeterogeneityParameter.COPY_RATIO_NOISE_FLOOR, copyRatioNoiseFloor),
                new Parameter<>(TumorHeterogeneityParameter.COPY_RATIO_NOISE_FACTOR, copyRatioNoiseFactor),
                new Parameter<>(TumorHeterogeneityParameter.MINOR_ALLELE_FRACTION_NOISE_FACTOR, minorAlleleFractionNoiseFactor),
                new Parameter<>(TumorHeterogeneityParameter.POPULATION_MIXTURE,
                        new PopulationMixture(populationMixture.populationFractions(), populationMixture.variantProfileCollection(), priors.normalPloidyState()))));
        Utils.validateArg(concentration > 0, "Concentration must be positive.");
        Utils.validateArg(copyRatioNoiseFactor >= 0, "Copy-ratio noise floor must be non-negative.");
        Utils.validateArg(copyRatioNoiseFactor >= 0, "Copy-ratio noise factor must be non-negative.");
        Utils.validateArg(minorAlleleFractionNoiseFactor >= 1, "Minor-allele-fraction noise factor must be >= 1.");
        Utils.nonNull(populationMixture);
        Utils.nonNull(priors);
        this.priors = priors;
    }

    public double concentration() {
        return get(TumorHeterogeneityParameter.CONCENTRATION, Double.class);
    }

    public double copyRatioNoiseFloor() {
        return get(TumorHeterogeneityParameter.COPY_RATIO_NOISE_FLOOR, Double.class);
    }

    public double copyRatioNoiseFactor() {
        return get(TumorHeterogeneityParameter.COPY_RATIO_NOISE_FACTOR, Double.class);
    }

    public double minorAlleleFractionNoiseFactor() {
        return get(TumorHeterogeneityParameter.MINOR_ALLELE_FRACTION_NOISE_FACTOR, Double.class);
    }

    public PopulationMixture populationMixture() {
        return get(TumorHeterogeneityParameter.POPULATION_MIXTURE, PopulationMixture.class);
    }

    public TumorHeterogeneityPriorCollection priors() {
        return priors;
    }

    /**
     * Initialize state with evenly distributed population fractions and normal variant profiles.
     */
    static TumorHeterogeneityState initializeState(final TumorHeterogeneityPriorCollection priors,
                                                   final int numSegments,
                                                   final int numPopulations) {
        final double concentration = priors.concentrationPriorHyperparameterValues().getAlpha() / priors.concentrationPriorHyperparameterValues().getBeta();
        final double copyRatioNoiseFloor = priors.copyRatioNoiseFloorPriorHyperparameterValues().getAlpha() / priors.copyRatioNoiseFloorPriorHyperparameterValues().getBeta();
        final double copyRatioNoiseFactor = 1. + priors.copyRatioNoiseFactorPriorHyperparameterValues().getAlpha() / priors.copyRatioNoiseFactorPriorHyperparameterValues().getBeta();
        final double minorAlleleFractionNoiseFactor = 1. + priors.minorAlleleFractionNoiseFactorPriorHyperparameterValues().getAlpha() / priors.minorAlleleFractionNoiseFactorPriorHyperparameterValues().getBeta();
        //initialize population fractions to be evenly distributed
        final PopulationMixture.PopulationFractions populationFractions =
                new PopulationMixture.PopulationFractions(Collections.nCopies(numPopulations, 1. / numPopulations));
        //initialize variant profiles to normal
        final int numVariantPopulations = numPopulations - 1;
        final PopulationMixture.VariantProfileCollection variantProfileCollection =
                initializeNormalProfiles(numVariantPopulations, numSegments, priors.normalPloidyState());
        final PopulationMixture populationMixture = new PopulationMixture(populationFractions, variantProfileCollection, priors.normalPloidyState());
        return new TumorHeterogeneityState(
                concentration, copyRatioNoiseFloor, copyRatioNoiseFactor, minorAlleleFractionNoiseFactor, populationMixture, priors);
    }

    /**
     * Initialize state from a modeller run that assumes one clonal population and one normal population.
     */
    static TumorHeterogeneityState initializeStateFromClonalResult(final TumorHeterogeneityPriorCollection priors,
                                                                   final TumorHeterogeneityModeller clonalModeller,
                                                                   final int maxNumPopulations) {
        //double check some of the parameters
        Utils.validateArg(clonalModeller.getConcentrationSamples().size() > 0,
                "Clonal modeller must have a non-zero number of samples.");
        Utils.validateArg(clonalModeller.getPopulationFractionsSamples().get(0).size() == 2,
                "Clonal modeller must have two populations (clone + normal).");

        //initialize global parameters to prior mean
        final double concentration = priors.concentrationPriorHyperparameterValues().getAlpha() / priors.concentrationPriorHyperparameterValues().getBeta();
        final double copyRatioNoiseFloor = priors.copyRatioNoiseFloorPriorHyperparameterValues().getAlpha() / priors.copyRatioNoiseFloorPriorHyperparameterValues().getBeta();
        final double copyRatioNoiseFactor = 1. + priors.copyRatioNoiseFactorPriorHyperparameterValues().getAlpha() / priors.copyRatioNoiseFactorPriorHyperparameterValues().getBeta();
        final double minorAlleleFractionNoiseFactor = 1. + priors.minorAlleleFractionNoiseFactorPriorHyperparameterValues().getAlpha() / priors.minorAlleleFractionNoiseFactorPriorHyperparameterValues().getBeta();

        //initialize normal fraction to posterior mean of clonal result
        final double[] normalFractionSamples = clonalModeller.getPopulationFractionsSamples().stream()
                .mapToDouble(pfs -> pfs.get(NUM_POPULATIONS_CLONAL - 1)).toArray();
        final double normalFraction = new Mean().evaluate(normalFractionSamples);

        //initialize one variant profile to posterior mode of clonal result
        final List<PopulationMixture.VariantProfile> clonalProfileSamples = clonalModeller.getVariantProfileCollectionSamples().stream()
                .map(vpc -> vpc.get(0)).collect(Collectors.toList());
        final PopulationMixture.VariantProfile initialClonalProfile = calculateVariantProfilePosteriorMode(clonalProfileSamples);

        //build new population fractions and variant-profile collection with additional variant populations
        //split clonal fraction evenly among new populations and initialize with clonal profile
        final List<Double> initialFractions = new ArrayList<>();
        final List<PopulationMixture.VariantProfile> initialVariantProfiles = new ArrayList<>();
        //add clonal population
        initialFractions.add((1. - normalFraction) / (maxNumPopulations - 1));
        initialVariantProfiles.add(initialClonalProfile);
        //initialize additional variant profiles
        for (int i = 0; i < maxNumPopulations - NUM_POPULATIONS_CLONAL; i++) {
            initialFractions.add((1. - normalFraction) / (maxNumPopulations - 1));
            initialVariantProfiles.add(1, new PopulationMixture.VariantProfile(initialClonalProfile));
        }
        //add normal population fraction
        initialFractions.add(normalFraction);
        final PopulationMixture populationMixture = new PopulationMixture(initialFractions, initialVariantProfiles, priors.normalPloidyState());

        return new TumorHeterogeneityState(
                concentration, copyRatioNoiseFloor, copyRatioNoiseFactor, minorAlleleFractionNoiseFactor, populationMixture, priors);
    }

    /**
     * Calculate the per-segment ploidy-state modes
     */
    private static PopulationMixture.VariantProfile calculateVariantProfilePosteriorMode(final List<PopulationMixture.VariantProfile> variantProfiles) {
        final int numSegments = variantProfiles.get(0).numSegments();
        final List<PloidyState> ploidyStates = new ArrayList<>(numSegments);
        for (int segmentIndex = 0; segmentIndex < numSegments; segmentIndex++) {
            final int si = segmentIndex;
            final List<PloidyState> ploidyStateSamples = variantProfiles.stream()
                    .map(vp -> vp.ploidyState(si)).collect(Collectors.toList());
            //get the mode of the samples; if there is more than one, return the first
            final PloidyState ploidyStateMode = ploidyStateSamples.stream()
                    .collect(Collectors.groupingBy(Function.identity(), Collectors.counting()))
                    .entrySet().stream()
                    .collect(Collectors.groupingBy(Map.Entry::getValue, Collectors.mapping(Map.Entry::getKey, Collectors.toList())))
                    .entrySet().stream().max((o1, o2) -> o1.getKey().compareTo(o2.getKey())).map(Map.Entry::getValue)
                    .get().get(0);
            ploidyStates.add(ploidyStateMode);
        }
        return new PopulationMixture.VariantProfile(ploidyStates);
    }

    /**
     * Initialize variant profiles to normal.
     */
    private static PopulationMixture.VariantProfileCollection initializeNormalProfiles(final int numVariantPopulations,
                                                                                       final int numSegments,
                                                                                       final PloidyState normalPloidyState) {
        return new PopulationMixture.VariantProfileCollection(
                Collections.nCopies(numVariantPopulations, initializeNormalProfile(numSegments, normalPloidyState)));
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
