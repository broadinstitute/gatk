package org.broadinstitute.hellbender.tools.exome.germlinehmm;

import org.broadinstitute.hellbender.tools.coveragemodel.germline.GermlineCNVCaller;
import org.broadinstitute.hellbender.tools.coveragemodel.germline.IntegerCopyNumberExpectationsCalculator;
import org.broadinstitute.hellbender.tools.coveragemodel.interfaces.TargetLikelihoodCalculator;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.hmm.HMM;

import javax.annotation.Nonnull;
import java.io.Serializable;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * This class represents a hidden Markov model with integer copy number states as its hidden states.
 *
 * It is used for inferring germline copy number variation from coverage data.
 * See {@link IntegerCopyNumberExpectationsCalculator} and {@link GermlineCNVCaller}.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class IntegerCopyNumberHMM<DATA>
        implements HMM<DATA, Target, IntegerCopyNumberState>, Serializable {

    private static final long serialVersionUID = 4013800307709357540L;

    private final TargetLikelihoodCalculator<DATA> emissionProbabilityCalculator;

    private final List<IntegerCopyNumberState> hiddenStates;

    private final IntegerCopyNumberTransitionProbabilityCacheCollection transitionProbabilityCacheCollection;

    private final String sampleSexGenotype;

    /**
     * Public constructor.
     *
     * Note: the copy number transition probability cache collection must be (1) padded, and (2)
     *       include transition matrix data for the sex genotype
     *
     * @param emissionProbabilityCalculator an instance of {@link TargetLikelihoodCalculator}
     * @param transitionProbabilityCacheCollection an instance of {@link IntegerCopyNumberTransitionProbabilityCacheCollection}
     * @param sampleSexGenotype sample sex genotype
     */
    public IntegerCopyNumberHMM(@Nonnull final TargetLikelihoodCalculator<DATA> emissionProbabilityCalculator,
                                @Nonnull final IntegerCopyNumberTransitionProbabilityCacheCollection transitionProbabilityCacheCollection,
                                @Nonnull final String sampleSexGenotype) {
        this.emissionProbabilityCalculator = Utils.nonNull(emissionProbabilityCalculator,
                "The emission probability calculator must be non-null");
        this.transitionProbabilityCacheCollection = Utils.nonNull(transitionProbabilityCacheCollection,
                "The transition probability cache collection must be non-null");
        if (!this.transitionProbabilityCacheCollection.isPadded()) {
            throw new IllegalStateException("The transition probability cache collection must be padded");
        }
        if (!transitionProbabilityCacheCollection.getSexGenotypes().contains(sampleSexGenotype)) {
            throw new IllegalArgumentException(String.format("Sex genotype \"%s\" is not contained in the transition" +
                    " probability cache collection", sampleSexGenotype));
        }
        this.sampleSexGenotype = Utils.nonNull(sampleSexGenotype, "The sample sex genotype must be non-null");
        final int maxCopyNumber = transitionProbabilityCacheCollection.getMaxCopyNumber();
        hiddenStates = Collections.unmodifiableList(IntStream.range(0, maxCopyNumber + 1)
                .mapToObj(IntegerCopyNumberState::new)
                .collect(Collectors.toList()));
    }

    @Override
    public List<IntegerCopyNumberState> hiddenStates() {
        return hiddenStates;
    }

    @Override
    public double logPriorProbability(@Nonnull final IntegerCopyNumberState state,
                                      @Nonnull final Target position) {
        return transitionProbabilityCacheCollection.logStationaryProbability(sampleSexGenotype,
                Utils.nonNull(position, "The target must be non-null").getContig(),
                Utils.nonNull(state, "The copy number state must be non-null"));
    }

    @Override
    public double logTransitionProbability(@Nonnull final IntegerCopyNumberState currentState,
                                           @Nonnull final Target currentPosition,
                                           @Nonnull final IntegerCopyNumberState nextState,
                                           @Nonnull final Target nextPosition) {
        final double distance = Target.calculateDistance(currentPosition, nextPosition);
        if (distance == Double.POSITIVE_INFINITY) {
            return logPriorProbability(nextState, nextPosition);
        } else {
            return transitionProbabilityCacheCollection.logTransitionProbability((int) distance, sampleSexGenotype,
                    currentPosition.getContig(), nextState, currentState);
        }
    }

    @Override
    public double logEmissionProbability(@Nonnull final DATA data,
                                         @Nonnull final IntegerCopyNumberState state,
                                         @Nonnull final Target position) {
        return emissionProbabilityCalculator.logLikelihood(
                Utils.nonNull(data, "Emission data must be non-null"),
                Utils.nonNull(state, "Integer copy number state must be non-null").getCopyNumber(),
                Utils.nonNull(position, "Target must be non-null"));
    }

    public void clearCaches() {
        transitionProbabilityCacheCollection.clearCaches();
    }
}
