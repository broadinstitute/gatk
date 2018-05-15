package org.broadinstitute.hellbender.tools.copynumber.arguments;

import java.util.EnumMap;

/**
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class GermlineCNVHybridADVIArgumentCollection extends HybridADVIArgumentCollection {
    private static final long serialVersionUID = 1L;

    private static final EnumMap<HybridADVIArgument, Object> defaultValue = new EnumMap<>(HybridADVIArgument.class);
    static {
        defaultValue.put(HybridADVIArgument.LEARNING_RATE, 0.01);
        defaultValue.put(HybridADVIArgument.ADAMAX_BETA_1, 0.9);
        defaultValue.put(HybridADVIArgument.ADAMAX_BETA_2, 0.99);
        defaultValue.put(HybridADVIArgument.LOG_EMISSION_SAMPLES_PER_ROUND, 50);
        defaultValue.put(HybridADVIArgument.LOG_EMISSION_SAMPLING_MEDIAN_REL_ERROR, 0.005);
        defaultValue.put(HybridADVIArgument.LOG_EMISSION_SAMPLING_ROUNDS, 10);
        defaultValue.put(HybridADVIArgument.MAX_ADVI_ITER_FIRST_EPOCH, 5000);
        defaultValue.put(HybridADVIArgument.MAX_ADVI_ITER_SUBSEQUENT_EPOCHS, 200);
        defaultValue.put(HybridADVIArgument.MIN_TRAINING_EPOCHS, 10);
        defaultValue.put(HybridADVIArgument.MAX_TRAINING_EPOCHS, 50);
        defaultValue.put(HybridADVIArgument.INITIAL_TEMPERATURE, 1.5);
        defaultValue.put(HybridADVIArgument.NUM_THERMAL_ADVI_ITERS, 2500);
        defaultValue.put(HybridADVIArgument.CONVERGENCE_SNR_AVERAGING_WINDOW, 500);
        defaultValue.put(HybridADVIArgument.CONVERGENCE_SNR_TRIGGER_THRESHOLD, 0.1);
        defaultValue.put(HybridADVIArgument.CONVERGENCE_SNR_COUNTDOWN_WINDOW, 10);
        defaultValue.put(HybridADVIArgument.MAX_CALLING_ITERS, 10);
        defaultValue.put(HybridADVIArgument.CALLER_UPDATE_CONVERGENCE_THRESHOLD, 0.001);
        defaultValue.put(HybridADVIArgument.CALLER_INTERNAL_ADMIXING_RATE, 0.75);
        defaultValue.put(HybridADVIArgument.CALLER_EXTERNAL_ADMIXING_RATE, 1.00);
        defaultValue.put(HybridADVIArgument.DISABLE_ANNEALING, false);
        defaultValue.put(HybridADVIArgument.DISABLE_CALLER, false);
        defaultValue.put(HybridADVIArgument.DISABLE_SAMPLER, false);
    }

    @Override
    public Object getDefaultValue(HybridADVIArgument arg) {
        return defaultValue.get(arg);
    }
}
