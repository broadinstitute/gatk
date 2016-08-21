package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import com.google.common.primitives.Doubles;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyState;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.mcmc.GibbsSampler;
import org.broadinstitute.hellbender.utils.mcmc.ParameterizedModel;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class TumorHeterogeneityModeller {
    private static final double EPSILON = TumorHeterogeneityUtils.EPSILON;

    static final double CONCENTRATION_MIN = EPSILON;
    static final double CONCENTRATION_MAX = 10.;

    private static final int NUM_SAMPLES_PER_LOG_ENTRY = 10;
    private static final int NUM_WALKERS = 50;

    private final ParameterizedModel<TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> model;
    private final TumorHeterogeneityData data;
    private final TumorHeterogeneityPriorCollection priors;

    private final List<Double> concentrationSamples = new ArrayList<>();
    private final List<Double> copyRatioNoiseFloorSamples = new ArrayList<>();
    private final List<Double> copyRatioNoiseFactorSamples = new ArrayList<>();
    private final List<Double> minorAlleleFractionNoiseFactorSamples = new ArrayList<>();
    private final List<PopulationMixture> populationMixtureSamples = new ArrayList<>();

    public TumorHeterogeneityModeller(final TumorHeterogeneityData data,
                                      final TumorHeterogeneityPriorCollection priors,
                                      final int numPopulations,
                                      final RandomGenerator rng) {
        this(data, TumorHeterogeneityState.initializeState(priors, data.numSegments(), numPopulations), rng);
    }

    public TumorHeterogeneityModeller(final TumorHeterogeneityData data,
                                      final TumorHeterogeneityState initialState,
                                      final RandomGenerator rng) {
        Utils.nonNull(data);
        Utils.nonNull(initialState);
        Utils.nonNull(rng);

        this.data = data;
        this.priors = initialState.priors();

        final double concentrationSliceSamplingWidth = initialState.concentration();
        final double copyRatioNoiseFloorSliceSamplingWidth = initialState.copyRatioNoiseFloor();
        final double copyRatioNoiseFactorSliceSamplingWidth = initialState.copyRatioNoiseFactor();
        final double minorAlleleFractionNoiseFactorSliceSamplingWidth = initialState.minorAlleleFractionNoiseFactor();
        final int numPopulations = initialState.populationMixture().numPopulations();
        final List<PloidyState> ploidyStates = priors.ploidyStatePrior().ploidyStates();
        final PloidyState normalPloidyState = priors.normalPloidyState();
        final int numWalkers = NUM_WALKERS;

        //define samplers
        final TumorHeterogeneitySamplers.ConcentrationSampler concentrationSampler =
                new TumorHeterogeneitySamplers.ConcentrationSampler(CONCENTRATION_MIN, CONCENTRATION_MAX, concentrationSliceSamplingWidth, numWalkers);
        final TumorHeterogeneitySamplers.CopyRatioNoiseFloorSampler copyRatioNoiseFloorSampler =
                new TumorHeterogeneitySamplers.CopyRatioNoiseFloorSampler(copyRatioNoiseFloorSliceSamplingWidth, numWalkers);
        final TumorHeterogeneitySamplers.CopyRatioNoiseFactorSampler copyRatioNoiseFactorSampler =
                new TumorHeterogeneitySamplers.CopyRatioNoiseFactorSampler(copyRatioNoiseFactorSliceSamplingWidth, numWalkers);
        final TumorHeterogeneitySamplers.MinorAlleleFractionNoiseFactorSampler minorAlleleFractionNoiseFactorSampler =
                new TumorHeterogeneitySamplers.MinorAlleleFractionNoiseFactorSampler(minorAlleleFractionNoiseFactorSliceSamplingWidth, numWalkers);
        final TumorHeterogeneitySamplers.PopulationMixtureSampler populationMixtureSampler =
                new TumorHeterogeneitySamplers.PopulationMixtureSampler(rng, numPopulations, ploidyStates, normalPloidyState, numWalkers);

        model = new ParameterizedModel.GibbsBuilder<>(initialState, data)
                .addParameterSampler(TumorHeterogeneityParameter.CONCENTRATION, concentrationSampler, Double.class)
                .addParameterSampler(TumorHeterogeneityParameter.COPY_RATIO_NOISE_FLOOR, copyRatioNoiseFloorSampler, Double.class)
                .addParameterSampler(TumorHeterogeneityParameter.COPY_RATIO_NOISE_FACTOR, copyRatioNoiseFactorSampler, Double.class)
                .addParameterSampler(TumorHeterogeneityParameter.MINOR_ALLELE_FRACTION_NOISE_FACTOR, minorAlleleFractionNoiseFactorSampler, Double.class)
                .addParameterSampler(TumorHeterogeneityParameter.POPULATION_MIXTURE, populationMixtureSampler, PopulationMixture.class)
                .build();
    }

    /**
     * Adds {@code numSamples - numBurnIn} Markov-Chain Monte-Carlo samples of the parameter posteriors (generated using
     * Gibbs sampling) to the collections held internally.  The current {@link TumorHeterogeneityState} held internally is used
     * to initialize the Markov Chain.
     * @param numSamples    total number of samples per posterior
     * @param numBurnIn     number of burn-in samples to discard
     */
    public void fitMCMC(final int numSamples, final int numBurnIn) {
        Utils.validateArg(numSamples > 0, "Total number of samples must be positive.");
        Utils.validateArg(0 <= numBurnIn && numBurnIn < numSamples,
                "Number of burn-in samples to discard must be non-negative and strictly less than total number of samples.");
        //run MCMC
        final GibbsSampler<TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> gibbsSampler
                = new GibbsSampler<>(NUM_WALKERS * numSamples, model);
        gibbsSampler.setNumSamplesPerLogEntry(NUM_SAMPLES_PER_LOG_ENTRY);
        gibbsSampler.runMCMC();
        //update posterior samples
        concentrationSamples.addAll(gibbsSampler.getSamples(TumorHeterogeneityParameter.CONCENTRATION,
                Double.class, NUM_WALKERS * numBurnIn));
        copyRatioNoiseFloorSamples.addAll(gibbsSampler.getSamples(TumorHeterogeneityParameter.COPY_RATIO_NOISE_FLOOR,
                Double.class, NUM_WALKERS * numBurnIn));
        copyRatioNoiseFactorSamples.addAll(gibbsSampler.getSamples(TumorHeterogeneityParameter.COPY_RATIO_NOISE_FACTOR,
                Double.class, NUM_WALKERS * numBurnIn));
        minorAlleleFractionNoiseFactorSamples.addAll(gibbsSampler.getSamples(TumorHeterogeneityParameter.MINOR_ALLELE_FRACTION_NOISE_FACTOR,
                Double.class, NUM_WALKERS * numBurnIn));
        populationMixtureSamples.addAll(gibbsSampler.getSamples(TumorHeterogeneityParameter.POPULATION_MIXTURE,
                PopulationMixture.class, NUM_WALKERS * numBurnIn));
    }

    /**
     * Returns an unmodifiable view of the list of samples of the concentration posterior.
     * @return  unmodifiable view of the list of samples of the concentration posterior
     */
    public List<Double> getConcentrationSamples() {
        return Collections.unmodifiableList(concentrationSamples);
    }

    /**
     * Returns an unmodifiable view of the list of samples of the copy-ratio noise-floor posterior.
     * @return  unmodifiable view of the list of samples of the copy-ratio noise-floor posterior
     */
    public List<Double> getCopyRatioNoiseFloorSamples() {
        return Collections.unmodifiableList(copyRatioNoiseFloorSamples);
    }

    /**
     * Returns an unmodifiable view of the list of samples of the copy-ratio noise-factor posterior.
     * @return  unmodifiable view of the list of samples of the copy-ratio noise-factor posterior
     */
    public List<Double> getCopyRatioNoiseFactorSamples() {
        return Collections.unmodifiableList(copyRatioNoiseFactorSamples);
    }

    /**
     * Returns an unmodifiable view of the list of samples of the minor-allele-fraction noise-factor posterior.
     * @return  unmodifiable view of the list of samples of the minor-allele-fraction noise-factor posterior
     */
    public List<Double> getMinorAlleleFractionNoiseFactorSamples() {
        return Collections.unmodifiableList(minorAlleleFractionNoiseFactorSamples);
    }

    /**
     * Returns an unmodifiable view of the list of samples of the population-fractions posterior, represented as a list of
     * {@link PopulationMixture.PopulationFractions} objects.
     * @return  unmodifiable view of the list of samples of the population-fractions posterior
     */
    public List<PopulationMixture.PopulationFractions> getPopulationFractionsSamples() {
        return Collections.unmodifiableList(populationMixtureSamples.stream().map(PopulationMixture::populationFractions).collect(Collectors.toList()));
    }

    /**
     * Returns an unmodifiable view of the list of samples of the variant-profile-collection posterior, represented as a list of
     * {@link PopulationMixture.PopulationFractions} objects.
     * @return  unmodifiable view of the list of samples of the variant-profile-collection posterior
     */
    public List<PopulationMixture.VariantProfileCollection> getVariantProfileCollectionSamples() {
        return Collections.unmodifiableList(populationMixtureSamples.stream().map(PopulationMixture::variantProfileCollection).collect(Collectors.toList()));
    }

    public List<Double> getPloidySamples() {
        return Collections.unmodifiableList(populationMixtureSamples.stream().map(pm -> pm.ploidy(data)).collect(Collectors.toList()));
    }

    public void outputSamples(final File outputFile) {
        if (concentrationSamples.size() == 0) {
            throw new IllegalStateException("Cannot output modeller result before samples have been generated.");
        }
        try (final FileWriter writer = new FileWriter(outputFile)) {
            outputFile.createNewFile();

            final List<PopulationMixture.PopulationFractions> populationFractionsSamples = getPopulationFractionsSamples();
            final List<Double> ploidySamples = getPloidySamples();

            final int numPopulations = populationFractionsSamples.get(0).size();

            //column headers
            for (int populationIndex = 0; populationIndex < numPopulations - 1; populationIndex++) {
                writer.write("POPULATION_" + populationIndex + "\t");
            }
            writer.write("POPULATION_NORMAL\tPLOIDY\n");

            //rows
            for (int sampleIndex = 0; sampleIndex < populationFractionsSamples.size(); sampleIndex++) {
                for (int populationIndex = 0; populationIndex < numPopulations; populationIndex++) {
                    final double populationFraction = populationFractionsSamples.get(sampleIndex).get(populationIndex);
                    writer.write(String.format("%.4f\t", populationFraction));
                }
                final double ploidy = ploidySamples.get(sampleIndex);
                writer.write(String.format("%.4f\n", ploidy));
            }
        } catch (final IOException e) {
            throw new GATKException("Error writing modeller samples.");
        }
    }

    public void outputSummary(final File outputFile) {
        if (concentrationSamples.size() == 0) {
            throw new IllegalStateException("Cannot output modeller result before samples have been generated.");
        }
        try (final FileWriter writer = new FileWriter(outputFile)) {
            outputFile.createNewFile();

            //comments
            writePosteriorSummary(writer, "ploidy", getPloidySamples());
            writePosteriorSummary(writer, "concentration", getConcentrationSamples());
            writePosteriorSummary(writer, "CR noise floor", getCopyRatioNoiseFloorSamples());
            writePosteriorSummary(writer, "CR noise factor", getCopyRatioNoiseFactorSamples());
            writePosteriorSummary(writer, "MAF noise factor", getMinorAlleleFractionNoiseFactorSamples());

            final List<PopulationMixture.PopulationFractions> populationFractionsSamples = getPopulationFractionsSamples();
            final List<PopulationMixture.VariantProfileCollection> variantProfileCollectionSamples = getVariantProfileCollectionSamples();

            final int numPopulations = populationFractionsSamples.get(0).size();
            for (int populationIndex = 0; populationIndex < numPopulations; populationIndex++) {
                final int pi = populationIndex;
                final List<Double> populationFractionSamples = populationFractionsSamples.stream()
                        .map(s -> s.get(pi)).collect(Collectors.toList());
                writePosteriorSummary(writer, "population fraction " + populationIndex, populationFractionSamples);
            }

            //column headers
            writer.write("POPULATION_INDEX\tSEGMENT_INDEX\tSEGMENT_INTERVAL\t");
            final List<PloidyState> ploidyStates = priors.ploidyStatePrior().ploidyStates();
            final int numPloidyStates = ploidyStates.size();
            for (int ploidyStateIndex = 0; ploidyStateIndex < numPloidyStates - 1; ploidyStateIndex++) {
                final PloidyState ploidyState = ploidyStates.get(ploidyStateIndex);
                writer.write(ploidyState.m() + "-" + ploidyState.n() + "\t");
            }
            final PloidyState ploidyState = ploidyStates.get(numPloidyStates - 1);
            writer.write(ploidyState.m() + "-" + ploidyState.n());
            writer.write(System.getProperty("line.separator"));

            //rows
            for (int populationIndex = 0; populationIndex < numPopulations; populationIndex++) {
                final int pi = populationIndex;
                final double[] populationFractionSamples = populationFractionsSamples.stream().mapToDouble(s -> s.get(pi)).toArray();
                final double populationFractionPosteriorMean = new Mean().evaluate(populationFractionSamples);

                if (populationIndex != numPopulations - 1 && populationFractionPosteriorMean >= 0.01) {
                    for (int segmentIndex = 0; segmentIndex < data.numSegments(); segmentIndex++) {
                        writer.write(populationIndex + "\t" + segmentIndex + "\t" + data.segments().get(segmentIndex).getInterval() + "\t");

                        final int si = segmentIndex;
                        for (int ploidyStateIndex = 0; ploidyStateIndex < numPloidyStates; ploidyStateIndex++) {
                            final int vpsi = ploidyStateIndex;
                            final double[] isPloidyStateSamples = variantProfileCollectionSamples.stream()
                                    .mapToDouble(vpc -> vpc.get(pi).ploidyState(si).equals(ploidyStates.get(vpsi)) ? 1. : 0)
                                    .toArray();
                            final double ploidyStatePosteriorMean = new Mean().evaluate(isPloidyStateSamples);
                            writer.write(String.format("%.3f", ploidyStatePosteriorMean));
                            if (ploidyStateIndex != numPloidyStates - 1) {
                                writer.write("\t");
                            }
                        }
                        if (!(segmentIndex == data.numSegments() - 1 && populationIndex == numPopulations - 2)) {
                            writer.write(System.getProperty("line.separator"));
                        }
                    }
                }
            }
        } catch (final IOException e) {
            throw new GATKException("Error writing modeller result.");
        }
    }

    private static void writePosteriorSummary(final FileWriter writer,
                                              final String label,
                                              final List<Double> samples) throws IOException {
        final double[] samplesArray = Doubles.toArray(samples);
        final double center = new Mean().evaluate(samplesArray);
        final double width = new StandardDeviation().evaluate(samplesArray);
        writer.write("#" + label + ": " + center + " " + width);
        writer.write(System.getProperty("line.separator"));
    }
}
