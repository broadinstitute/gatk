package org.broadinstitute.hellbender.tools.copynumber.arguments;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public abstract class HybridADVIArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    public static final String LEARNING_RATE_LONG_NAME = "learning-rate";
    public static final String ADAMAX_BETA_1_LONG_NAME = "adamax-beta-1";
    public static final String ADAMAX_BETA_2_LONG_NAME = "adamax-beta-2";
    public static final String LOG_EMISSION_SAMPLES_PER_ROUND_LONG_NAME = "log-emission-samples-per-round";
    public static final String LOG_EMISSION_SAMPLING_MEDIAN_REL_ERROR_LONG_NAME = "log-emission-sampling-median-rel-error";
    public static final String LOG_EMISSION_SAMPLING_ROUNDS_LONG_NAME = "log-emission-sampling-rounds";
    public static final String MAX_ADVI_ITER_FIRST_EPOCH_LONG_NAME = "max-advi-iter-first-epoch";
    public static final String MAX_ADVI_ITER_SUBSEQUENT_EPOCHS_LONG_NAME = "max-advi-iter-subsequent-epochs";
    public static final String MIN_TRAINING_EPOCHS_LONG_NAME = "min-training-epochs";
    public static final String MAX_TRAINING_EPOCHS_LONG_NAME = "max-training-epochs";
    public static final String INITIAL_TEMPERATURE_LONG_NAME = "initial-temperature";
    public static final String NUM_THERMAL_ADVI_ITERS_LONG_NAME = "num-thermal-advi-iters";
    public static final String CONVERGENCE_SNR_AVERAGING_WINDOW_LONG_NAME = "convergence-snr-averaging-window";
    public static final String CONVERGENCE_SNR_TRIGGER_THRESHOLD_LONG_NAME = "convergence-snr-trigger-threshold";
    public static final String CONVERGENCE_SNR_COUNTDOWN_WINDOW_LONG_NAME = "convergence-snr-countdown-window";
    public static final String MAX_CALLING_ITERS_LONG_NAME = "max-calling-iters";
    public static final String CALLER_UPDATE_CONVERGENCE_THRESHOLD_LONG_NAME = "caller-update-convergence-threshold";
    public static final String CALLER_INTERNAL_ADMIXING_RATE_LONG_NAME = "caller-internal-admixing-rate";
    public static final String CALLER_EXTERNAL_ADMIXING_RATE_LONG_NAME = "caller-external-admixing-rate";
    public static final String DISABLE_SAMPLER_LONG_NAME = "disable-sampler";
    public static final String DISABLE_CALLER_LONG_NAME = "disable-caller";
    public static final String DISABLE_ANNEALING_RATE_LONG_NAME = "disable-annealing";

    public enum HybridADVIArgument {
        LEARNING_RATE("learning_rate"),
        ADAMAX_BETA_1("adamax_beta1"),
        ADAMAX_BETA_2("adamax_beta2"),
        LOG_EMISSION_SAMPLES_PER_ROUND("log_emission_samples_per_round"),
        LOG_EMISSION_SAMPLING_MEDIAN_REL_ERROR("log_emission_sampling_median_rel_error"),
        LOG_EMISSION_SAMPLING_ROUNDS("log_emission_sampling_rounds"),
        MAX_ADVI_ITER_FIRST_EPOCH("max_advi_iter_first_epoch"),
        MAX_ADVI_ITER_SUBSEQUENT_EPOCHS("max_advi_iter_subsequent_epochs"),
        MIN_TRAINING_EPOCHS("min_training_epochs"),
        MAX_TRAINING_EPOCHS("max_training_epochs"),
        INITIAL_TEMPERATURE("initial_temperature"),
        NUM_THERMAL_ADVI_ITERS("num_thermal_advi_iters"),
        CONVERGENCE_SNR_AVERAGING_WINDOW("convergence_snr_averaging_window"),
        CONVERGENCE_SNR_TRIGGER_THRESHOLD("convergence_snr_trigger_threshold"),
        CONVERGENCE_SNR_COUNTDOWN_WINDOW("convergence_snr_countdown_window"),
        MAX_CALLING_ITERS("max_calling_iters"),
        CALLER_UPDATE_CONVERGENCE_THRESHOLD("caller_update_convergence_threshold"),
        CALLER_INTERNAL_ADMIXING_RATE("caller_internal_admixing_rate"),
        CALLER_EXTERNAL_ADMIXING_RATE("caller_external_admixing_rate"),
        DISABLE_SAMPLER("disable_sampler"),
        DISABLE_CALLER("disable_caller"),
        DISABLE_ANNEALING("disable_annealing");

        public final String pythonArg;

        HybridADVIArgument(final String pythonArg) {
            this.pythonArg = pythonArg;
        }
    }

    public abstract Object getDefaultValue(final HybridADVIArgument arg);

    @Argument(
            doc="Adamax optimizer learning rate.",
            fullName = LEARNING_RATE_LONG_NAME,
            optional = true,
            minValue = 0.
    )
    private double learningRate =
            (Double)getDefaultValue(HybridADVIArgument.LEARNING_RATE);

    @Argument(
            doc="Adamax optimizer first moment estimation forgetting factor.",
            fullName = ADAMAX_BETA_1_LONG_NAME,
            optional = true,
            minValue = 0.,
            maxValue = 1.0
    )
    private double adamaxBeta1 =
            (Double)getDefaultValue(HybridADVIArgument.ADAMAX_BETA_1);

    @Argument(
            doc="Adamax optimizer second moment estimation forgetting factor.",
            fullName = ADAMAX_BETA_2_LONG_NAME,
            optional = true,
            minValue = 0.,
            maxValue = 1.0
    )
    private double adamaxBeta2 =
            (Double)getDefaultValue(HybridADVIArgument.ADAMAX_BETA_2);

    @Argument(
            doc="Log emission samples drawn per round of sampling.",
            fullName = LOG_EMISSION_SAMPLES_PER_ROUND_LONG_NAME,
            optional = true,
            minValue = 0
    )
    private int logEmissionSamplesPerRound =
            (Integer)getDefaultValue(HybridADVIArgument.LOG_EMISSION_SAMPLES_PER_ROUND);

    @Argument(
            doc="Maximum tolerated median relative error in log emission sampling.",
            fullName = LOG_EMISSION_SAMPLING_MEDIAN_REL_ERROR_LONG_NAME,
            optional = true,
            minValue = 0
    )
    private double logEmissionMedianRelError =
            (Double)getDefaultValue(HybridADVIArgument.LOG_EMISSION_SAMPLING_MEDIAN_REL_ERROR);

    @Argument(
            doc="Log emission maximum sampling rounds.",
            fullName = LOG_EMISSION_SAMPLING_ROUNDS_LONG_NAME,
            optional = true,
            minValue = 0
    )
    private int logEmissionSamplingRounds =
            (Integer)getDefaultValue(HybridADVIArgument.LOG_EMISSION_SAMPLING_ROUNDS);

    @Argument(
            doc="Maximum ADVI iterations in the first epoch.",
            fullName = MAX_ADVI_ITER_FIRST_EPOCH_LONG_NAME,
            optional = true,
            minValue = 0
    )
    private int maxADVIItersFirstEpoch =
            (Integer)getDefaultValue(HybridADVIArgument.MAX_ADVI_ITER_FIRST_EPOCH);

    @Argument(
            doc="Maximum ADVI iterations in subsequent epochs.",
            fullName = MAX_ADVI_ITER_SUBSEQUENT_EPOCHS_LONG_NAME,
            optional = true,
            minValue = 0
    )
    private int maxADVIItersSubsequentEpochs =
            (Integer)getDefaultValue(HybridADVIArgument.MAX_ADVI_ITER_SUBSEQUENT_EPOCHS);

    @Argument(
            doc="Minimum number of training epochs.",
            fullName = MIN_TRAINING_EPOCHS_LONG_NAME,
            optional = true,
            minValue = 0
    )
    private int minTrainingEpochs =
            (Integer)getDefaultValue(HybridADVIArgument.MIN_TRAINING_EPOCHS);

    @Argument(
            doc="Maximum number of training epochs.",
            fullName = MAX_TRAINING_EPOCHS_LONG_NAME,
            optional = true,
            minValue = 0
    )
    private int maxTrainingEpochs =
            (Integer)getDefaultValue(HybridADVIArgument.MAX_TRAINING_EPOCHS);

    @Argument(
            doc="Initial temperature (for DA-ADVI).",
            fullName = INITIAL_TEMPERATURE_LONG_NAME,
            optional = true,
            minValue = 0
    )
    private double initialTemperature =
            (Double)getDefaultValue(HybridADVIArgument.INITIAL_TEMPERATURE);

    @Argument(
            doc="Number of thermal ADVI iterations (for DA-ADVI).",
            fullName = NUM_THERMAL_ADVI_ITERS_LONG_NAME,
            optional = true,
            minValue = 0
    )
    private int numThermalADVIIters =
            (Integer)getDefaultValue(HybridADVIArgument.NUM_THERMAL_ADVI_ITERS);

    @Argument(
            doc="Averaging window for calculating training signal-to-noise ratio (SNR) for convergence checking.",
            fullName = CONVERGENCE_SNR_AVERAGING_WINDOW_LONG_NAME,
            optional = true,
            minValue = 0
    )
    private int convergenceSNRAveragingWindow =
            (Integer)getDefaultValue(HybridADVIArgument.CONVERGENCE_SNR_AVERAGING_WINDOW);

    @Argument(
            doc="The SNR threshold to be reached before triggering the convergence countdown.",
            fullName = CONVERGENCE_SNR_TRIGGER_THRESHOLD_LONG_NAME,
            optional = true,
            minValue = 0
    )
    private double convergenceSNRTriggerThreshold =
            (Double)getDefaultValue(HybridADVIArgument.CONVERGENCE_SNR_TRIGGER_THRESHOLD);

    @Argument(
            doc="The number of ADVI iterations during which the SNR is required to stay below the set threshold " +
                    "for convergence.",
            fullName = CONVERGENCE_SNR_COUNTDOWN_WINDOW_LONG_NAME,
            optional = true,
            minValue = 0
    )
    private int convergenceSNRCountdownWindow =
            (Integer)getDefaultValue(HybridADVIArgument.CONVERGENCE_SNR_COUNTDOWN_WINDOW);

    @Argument(
            doc="Maximum number of internal self-consistency iterations within each calling step.",
            fullName = MAX_CALLING_ITERS_LONG_NAME,
            optional = true,
            minValue = 0
    )
    private int maxCallingIters =
            (Integer)getDefaultValue(HybridADVIArgument.MAX_CALLING_ITERS);

    @Argument(
            doc="Maximum tolerated calling update size for convergence.",
            fullName = CALLER_UPDATE_CONVERGENCE_THRESHOLD_LONG_NAME,
            optional = true,
            minValue = 0
    )
    private double callerUpdateConvergenceThreshold =
            (Double)getDefaultValue(HybridADVIArgument.CALLER_UPDATE_CONVERGENCE_THRESHOLD);

    @Argument(
            doc="Admixing ratio of new and old called posteriors (between 0 and 1; larger values implies using " +
                    "more of the new posterior and less of the old posterior) for internal convergence loops.",
            fullName = CALLER_INTERNAL_ADMIXING_RATE_LONG_NAME,
            optional = true,
            minValue = 0
    )
    private double callerInternalAdmixingRate =
            (Double)getDefaultValue(HybridADVIArgument.CALLER_INTERNAL_ADMIXING_RATE);

    @Argument(
            doc="Admixing ratio of new and old called posteriors (between 0 and 1; larger values implies using " +
                    "more of the new posterior and less of the old posterior) after convergence.",
            fullName = CALLER_EXTERNAL_ADMIXING_RATE_LONG_NAME,
            optional = true,
            minValue = 0
    )
    private double callerExternalAdmixingRate =
            (Double)getDefaultValue(HybridADVIArgument.CALLER_EXTERNAL_ADMIXING_RATE);

    @Argument(
            doc="(advanced) Disable sampler.",
            fullName = DISABLE_SAMPLER_LONG_NAME,
            optional = true
    )
    private boolean disableSampler =
            (Boolean)getDefaultValue(HybridADVIArgument.DISABLE_SAMPLER);

    @Argument(
            doc="(advanced) Disable caller.",
            fullName = DISABLE_CALLER_LONG_NAME,
            optional = true
    )
    private boolean disableCaller =
            (Boolean)getDefaultValue(HybridADVIArgument.DISABLE_CALLER);

    @Argument(
            doc="(advanced) Disable annealing.",
            fullName = DISABLE_ANNEALING_RATE_LONG_NAME,
            optional = true
    )
    private boolean disableAnnealing =
            (Boolean)getDefaultValue(HybridADVIArgument.DISABLE_ANNEALING);

    public void validate() {
        Utils.validateArg(maxTrainingEpochs >= minTrainingEpochs,
                "Maximum number of training epochs must be greater or equal to minimum number of training " +
                        "epochs.");
    }

    public List<String> generatePythonArguments() {
        final List<String> arguments = new ArrayList<>(Arrays.asList(
                String.format("--" + HybridADVIArgument.LEARNING_RATE.pythonArg + "=%e", learningRate),
                String.format("--" + HybridADVIArgument.ADAMAX_BETA_1.pythonArg + "=%e", adamaxBeta1),
                String.format("--" + HybridADVIArgument.ADAMAX_BETA_2.pythonArg + "=%e", adamaxBeta2),
                String.format("--" + HybridADVIArgument.LOG_EMISSION_SAMPLES_PER_ROUND.pythonArg + "=%d",
                        logEmissionSamplesPerRound),
                String.format("--" + HybridADVIArgument.LOG_EMISSION_SAMPLING_ROUNDS.pythonArg + "=%d",
                        logEmissionSamplingRounds),
                String.format("--" + HybridADVIArgument.LOG_EMISSION_SAMPLING_MEDIAN_REL_ERROR.pythonArg + "=%e",
                        logEmissionMedianRelError),
                String.format("--" + HybridADVIArgument.MAX_ADVI_ITER_FIRST_EPOCH.pythonArg + "=%d",
                        maxADVIItersFirstEpoch),
                String.format("--" + HybridADVIArgument.MAX_ADVI_ITER_SUBSEQUENT_EPOCHS.pythonArg + "=%d",
                        maxADVIItersSubsequentEpochs),
                String.format("--" + HybridADVIArgument.MIN_TRAINING_EPOCHS.pythonArg + "=%d",
                        minTrainingEpochs),
                String.format("--" + HybridADVIArgument.MAX_TRAINING_EPOCHS.pythonArg + "=%d",
                        maxTrainingEpochs),
                String.format("--" + HybridADVIArgument.INITIAL_TEMPERATURE.pythonArg + "=%e",
                        initialTemperature),
                String.format("--" + HybridADVIArgument.NUM_THERMAL_ADVI_ITERS.pythonArg + "=%d",
                        numThermalADVIIters),
                String.format("--" + HybridADVIArgument.CONVERGENCE_SNR_AVERAGING_WINDOW.pythonArg + "=%d",
                        convergenceSNRAveragingWindow),
                String.format("--" + HybridADVIArgument.CONVERGENCE_SNR_TRIGGER_THRESHOLD.pythonArg + "=%e",
                        convergenceSNRTriggerThreshold),
                String.format("--" + HybridADVIArgument.CONVERGENCE_SNR_COUNTDOWN_WINDOW.pythonArg + "=%d",
                        convergenceSNRCountdownWindow),
                String.format("--" + HybridADVIArgument.MAX_CALLING_ITERS.pythonArg + "=%d",
                        maxCallingIters),
                String.format("--" + HybridADVIArgument.CALLER_UPDATE_CONVERGENCE_THRESHOLD.pythonArg + "=%e",
                        callerUpdateConvergenceThreshold),
                String.format("--" + HybridADVIArgument.CALLER_INTERNAL_ADMIXING_RATE.pythonArg + "=%e",
                        callerInternalAdmixingRate),
                String.format("--" + HybridADVIArgument.CALLER_EXTERNAL_ADMIXING_RATE.pythonArg + "=%e",
                        callerExternalAdmixingRate)));
        if (disableCaller) {
            arguments.add("--" + HybridADVIArgument.DISABLE_CALLER.pythonArg + "=true");
        } else {
            arguments.add("--" + HybridADVIArgument.DISABLE_CALLER.pythonArg + "=false");
        }
        if (disableSampler) {
            arguments.add("--" + HybridADVIArgument.DISABLE_SAMPLER.pythonArg + "=true");
        } else {
            arguments.add("--" + HybridADVIArgument.DISABLE_SAMPLER.pythonArg + "=false");
        }
        if (disableAnnealing) {
            arguments.add("--" + HybridADVIArgument.DISABLE_ANNEALING.pythonArg + "=true");
        } else {
            arguments.add("--" + HybridADVIArgument.DISABLE_ANNEALING.pythonArg + "=false");
        }
        return arguments;
    }
}
