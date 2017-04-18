package org.broadinstitute.hellbender.tools.coveragemodel;

import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.util.FastMath;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.hmm.interfaces.AlleleMetadataProducer;
import org.broadinstitute.hellbender.utils.hmm.interfaces.CallStringProducer;
import org.broadinstitute.hellbender.utils.hmm.interfaces.ScalarProducer;

import javax.annotation.Nonnull;
import java.io.File;
import java.util.Arrays;
import java.util.Collections;
import java.util.function.Function;
import java.util.function.Supplier;

/**
 * This class implements the EM algorithm for the GATK coverage model. It is used for:
 *
 * (1) performing max likelihood estimation of the model parameters AND calculating posterior expectations
 *     via calling {@link CoverageModelEMAlgorithm {@link #runExpectationMaximization()}}
 *
 * (2) Calculating posterior expectations for the model parameters already existing in the workspace
 *     via calling {@link CoverageModelEMAlgorithm {@link #runExpectation()}}
 *
 * (see CNV-methods.pdf for technical details).
 *
 * This class does not store or perform any of the calculations. Rather, it controls the flow of the EM
 * algorithm subroutines (E- and M- steps) and makes calls to various methods in the provided
 * {@link CoverageModelEMWorkspace} in order to perform the actual computations.
 *
 * @param <STATE> copy ratio (or copy number) hidden state type
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */

public final class CoverageModelEMAlgorithm<STATE extends AlleleMetadataProducer & CallStringProducer & ScalarProducer> {
    private static final Logger logger = LogManager.getLogger(CoverageModelEMAlgorithm.class);
    private final CoverageModelArgumentCollection config;
    private final CoverageModelEMWorkspace<STATE> workspace;
    private EMAlgorithmStatus status;

    /**
     * Extract iteration count from a {@link SubroutineSignal}
     */
    private static Function<SubroutineSignal, String> ITERATION_EXTRACTOR = s ->
            "iters: " + s.<Integer>get(StandardSubroutineSignals.ITERATIONS);

    /**
     * This supplies a "N/A" string no matter what; used when nothing needs to be reported from
     * the {@link SubroutineSignal}
     */
    private static Function<SubroutineSignal, String> NOT_APPLICABLE_EXTRACTOR = s -> "N/A";

    /**
     * Public constructor.
     *
     * @param config configuration of the EM algorithm
     * @param workspace workspace for actual calculations
     */
    public CoverageModelEMAlgorithm(@Nonnull final CoverageModelArgumentCollection config,
                                    @Nonnull CoverageModelEMWorkspace<STATE> workspace) {
        this.config = Utils.nonNull(config, "Coverage model EM algorithm configuration must be non-null");
        config.validate();
        this.workspace = Utils.nonNull(workspace, "The coverage modeller EM algorithm workspace must be non-null");
        this.status = EMAlgorithmStatus.TBD;
        logger.info("EM algorithm initialized.");
    }

    /**
     * Log the header for iteration info
     */
    private void showIterationHeader() {
        final String header = String.format("%-15s%-20s%-20s%-20s%-20s", "Iterations", "Type", "Log Likelihood", "Update Size", "Misc.");
        logger.info(header);
        logger.info(StringUtils.repeat("=", header.length()));
    }

    /**
     * Log the iteration info
     *
     * @param iter iteration number
     * @param title iteration title
     * @param logLikelihood current log likelihood
     * @param updateSize current update size
     * @param misc other miscellaneous info
     */
    private void showIterationInfo(final int iter, final String title, final double logLikelihood,
                                   final double updateSize, final String misc) {
        final String header = String.format("%-15d%-20s%-20.6e%-20.6e%-20s", iter, title, logLikelihood, updateSize, misc);
        logger.info(header);
    }

    /**
     * A generic runner of a step in the EM algorithm
     *
     * @param routine a routine to run
     * @param miscFactory a function from the routine signal to a miscellaneous iteration info string
     * @param title title of the routine
     * @param iterInfo the iteration info class (will be updated)
     */
    private void runRoutine(@Nonnull final Supplier<SubroutineSignal> routine,
                            @Nonnull final Function<SubroutineSignal, String> miscFactory,
                            @Nonnull final String title,
                            @Nonnull EMAlgorithmIterationInfo iterInfo) {
        final SubroutineSignal sig = routine.get();
        final String misc = miscFactory.apply(sig);
        iterInfo.errorNorm = sig.<Double>get(StandardSubroutineSignals.RESIDUAL_ERROR_NORM);
        iterInfo.logLikelihood = getLogLikelihood();
        finalizeIteration();
        showIterationInfo(iterInfo.getIterationCount(), title, iterInfo.getLogLikelihood(), iterInfo.getErrorNorm(), misc);
    }

    /**
     * Run the full EM algorithm stack for variational inference and model parameter optimization
     *
     * @return the exit status of the EM algorithm
     */
    public EMAlgorithmStatus runExpectationMaximization() {
        config.validate();
        /* If copy ratio posterior calling is enabled, the first few iterations need to be robust.
         * Therefore, we disable the target-resolved unexplained variance estimation (if enabled) and
         * enforce isotropic unexplained variance */
        CoverageModelArgumentCollection.TargetSpecificVarianceUpdateMode currentTargetSpecificVarianceUpdateMode =
                config.getTargetSpecificVarianceUpdateMode();
        if (config.adaptiveTargetSpecificVarianceSolverModeSwitchingEnabled() &&
                !currentTargetSpecificVarianceUpdateMode.equals(CoverageModelArgumentCollection.TargetSpecificVarianceUpdateMode.ISOTROPIC)) {
            currentTargetSpecificVarianceUpdateMode = CoverageModelArgumentCollection.TargetSpecificVarianceUpdateMode.ISOTROPIC;
            logger.info("Overriding the requested unexplained variance solver to " +
                    CoverageModelArgumentCollection.TargetSpecificVarianceUpdateMode.ISOTROPIC.name());
        }

        showIterationHeader();

        double prevMStepLikelihood = Double.NEGATIVE_INFINITY;
        double latestMStepLikelihood = Double.NEGATIVE_INFINITY;
        final EMAlgorithmIterationInfo iterInfo = new EMAlgorithmIterationInfo(Double.NEGATIVE_INFINITY, 0, 1);
        final boolean updateBiasCovariates = config.getNumLatents() > 0;

        /* initial states -- these will change over the course of the algorithm adaptively */
        boolean updateCopyRatioPosteriors = false;
        boolean updateARDCoefficients = false;
        boolean paramEstimationConverged = false;
        boolean performMStep = true;

        while (iterInfo.iter <= config.getMaxEMIterations()) {

            /* cycle through E-step mean-field equations until they are satisfied to the desired degree */
            int iterEStep = 0;
            double maxPosteriorErrorNorm = 2 * config.getPosteriorAbsTol();
            while (iterEStep < config.getMaxEStepCycles() && maxPosteriorErrorNorm > config.getPosteriorAbsTol()) {
                double posteriorErrorNormReadDepth;
                double posteriorErrorNormSampleUnexplainedVariance;
                double posteriorErrorNormBias;
                double posteriorErrorNormCopyRatio;
                double posteriorErrorNormBiasCovariates;

                runRoutine(this::updateReadDepthLatentPosteriorExpectations, NOT_APPLICABLE_EXTRACTOR, "E_STEP_D", iterInfo);
                posteriorErrorNormReadDepth = iterInfo.errorNorm;

                /* initially, we just want read depth estimate */
                if (iterInfo.iter == 1) {
                    break;
                }

                if (updateBiasCovariates) {
                    runRoutine(this::updateBiasLatentPosteriorExpectations, NOT_APPLICABLE_EXTRACTOR, "E_STEP_Z", iterInfo);
                    posteriorErrorNormBias = iterInfo.errorNorm;
                    runRoutine(this::updateBiasCovariates, NOT_APPLICABLE_EXTRACTOR, "E_STEP_W", iterInfo);
                    posteriorErrorNormBiasCovariates = iterInfo.errorNorm;
                } else {
                    posteriorErrorNormBias = 0;
                    posteriorErrorNormBiasCovariates = 0;
                }

                if (config.sampleSpecificVarianceUpdateEnabled()) {
                    runRoutine(this::updateSampleUnexplainedVariance, ITERATION_EXTRACTOR, "E_STEP_GAMMA", iterInfo);
                    posteriorErrorNormSampleUnexplainedVariance = iterInfo.errorNorm;
                } else {
                    posteriorErrorNormSampleUnexplainedVariance = 0;
                }

                if (updateCopyRatioPosteriors) {
                    runRoutine(this::updateCopyRatioLatentPosteriorExpectations, NOT_APPLICABLE_EXTRACTOR, "E_STEP_C", iterInfo);
                    posteriorErrorNormCopyRatio = iterInfo.errorNorm;
                } else {
                    posteriorErrorNormCopyRatio = 0;
                }

                /* calculate the maximum change of posteriors in this cycle */
                maxPosteriorErrorNorm = Collections.max(Arrays.asList(
                        posteriorErrorNormReadDepth,
                        posteriorErrorNormBiasCovariates,
                        posteriorErrorNormSampleUnexplainedVariance,
                        posteriorErrorNormBias,
                        posteriorErrorNormCopyRatio));

                iterEStep++;
            }

            if (maxPosteriorErrorNorm > config.getPosteriorAbsTol()) {
                logger.warn("E-step cycles did not fully converge. Increase the maximum number of E-step cycles." +
                        " Continuing...");
            }

            /* parameter estimation */
            if (performMStep && !paramEstimationConverged) {
                int iterMStep = 0;
                double maxParamErrorNorm;

                /* sequential M-steps */
                while (iterMStep < config.getMaxMStepCycles()) {
                    double errorNormMeanTargetBias;
                    double errorNormUnexplainedVariance;
                    double errorNormARDCoefficients;

                    /* neglect the contribution from principal components in the first iteration */
                    final boolean ignoreBiasCovariates = iterInfo.iter == 1;
                    runRoutine(() -> updateTargetMeanBias(ignoreBiasCovariates), NOT_APPLICABLE_EXTRACTOR, "M_STEP_M", iterInfo);
                    errorNormMeanTargetBias = iterInfo.errorNorm;

                    if (config.targetSpecificVarianceUpdateEnabled()) {
                        final CoverageModelArgumentCollection.TargetSpecificVarianceUpdateMode finalizedPsiUpdateMode = currentTargetSpecificVarianceUpdateMode;
                        runRoutine(() -> updateTargetUnexplainedVariance(finalizedPsiUpdateMode), ITERATION_EXTRACTOR, "M_STEP_PSI", iterInfo);
                        errorNormUnexplainedVariance = iterInfo.errorNorm;
                    } else {
                        errorNormUnexplainedVariance = 0;
                    }

                    if (updateARDCoefficients) {
                        runRoutine(this::updateBiasCovariatesARDCoefficients, NOT_APPLICABLE_EXTRACTOR, "M_STEP_ALPHA", iterInfo);
                        errorNormARDCoefficients = iterInfo.errorNorm;
                    } else {
                        errorNormARDCoefficients = 0;
                    }

                    maxParamErrorNorm = Collections.max(Arrays.asList(
                            errorNormMeanTargetBias,
                            errorNormUnexplainedVariance,
                            errorNormARDCoefficients));

                    /* check convergence of parameter estimation */
                    if (updateCopyRatioPosteriors && maxParamErrorNorm < config.getParameterEstimationAbsoluteTolerance()) {
                        status = EMAlgorithmStatus.SUCCESS_PARAMS_TOL;
                        paramEstimationConverged = true;
                    }
                    iterMStep++;
                }

                prevMStepLikelihood = latestMStepLikelihood;
                latestMStepLikelihood = iterInfo.logLikelihood;

                /* if partially converged, start updating copy number posteriors */
                final boolean partiallyConvergedForCRCalling = FastMath.abs(latestMStepLikelihood - prevMStepLikelihood) <
                        config.getLogLikelihoodTolThresholdCRCalling();
                if (partiallyConvergedForCRCalling && config.copyRatioUpdateEnabled() && !updateCopyRatioPosteriors) {
                    updateCopyRatioPosteriors = true;
                    logger.info("Partial convergence achieved; starting to update copy ratio posteriors in" +
                            " the next iteration");
                }

                /* if partially converged, switch to target-resolved psi update mode (if enabled) */
                final boolean partiallyConvergedForPsiSwitching = FastMath.abs(latestMStepLikelihood - prevMStepLikelihood) <
                        config.getLogLikelihoodTolThresholdTargetSpecificVarianceSwitching();
                if (partiallyConvergedForPsiSwitching && config.targetSpecificVarianceUpdateEnabled() &&
                        config.adaptiveTargetSpecificVarianceSolverModeSwitchingEnabled() &&
                        !currentTargetSpecificVarianceUpdateMode.equals(CoverageModelArgumentCollection.TargetSpecificVarianceUpdateMode.TARGET_RESOLVED)) {
                    currentTargetSpecificVarianceUpdateMode = CoverageModelArgumentCollection.TargetSpecificVarianceUpdateMode.TARGET_RESOLVED;
                    logger.info("Partial convergence achieved; switching to target-specific unexplained variance in" +
                            " the next iteration");
                }

                /* if partially converged, switch to target-resolved psi update mode (if enabled) */
                final boolean partiallyConvergedForARDUpdate = FastMath.abs(latestMStepLikelihood - prevMStepLikelihood) <
                        config.getLogLikelihoodTolThresholdARDUpdate();
                if (partiallyConvergedForARDUpdate && config.isARDEnabled() && !updateARDCoefficients) {
                    updateARDCoefficients = true;
                    logger.info("Partial convergence achieved; enabling update of ARD coefficients");
                    runRoutine(this::updateBiasCovariatesARDCoefficients, NOT_APPLICABLE_EXTRACTOR, "M_STEP_ALPHA", iterInfo);
                }
            }

            /* check convergence in log likelihood change */
            if (FastMath.abs(latestMStepLikelihood - prevMStepLikelihood) < config.getLogLikelihoodTolerance()) {
                /* make sure that we have either already called copy ratio posteriors, or we are not required to */
                if (!config.copyRatioUpdateEnabled() || updateCopyRatioPosteriors) {
                    status = EMAlgorithmStatus.SUCCESS_LIKELIHOOD_TOL;
                    break;
                }
            } else if (iterInfo.iter == config.getMaxEMIterations() - 2) {
                performMStep = false; /* so that we end with an E-step */
            }

            iterInfo.incrementIterationCount();

            /* checkpointing */
            if (config.isRunCheckpointingEnabled() && iterInfo.iter % config.getRunCheckpointingInterval() == 0) {
                final String modelOutputAbsolutePath = new File(config.getRunCheckpointingPath(),
                        String.format("%s_iter_%d", CoverageModelGlobalConstants.MODEL_CHECKPOINT_PATH_PREFIX, iterInfo.iter)).getAbsolutePath();
                final String posteriorOutputAbsolutePath = new File(config.getRunCheckpointingPath(),
                        String.format("%s_iter_%d", CoverageModelGlobalConstants.POSTERIOR_CHECKPOINT_PATH_PREFIX, iterInfo.iter)).getAbsolutePath();
                saveModel(modelOutputAbsolutePath);
                savePosteriors(posteriorOutputAbsolutePath, CoverageModelEMWorkspace.PosteriorVerbosityLevel.BASIC);
            }
        }

        if (iterInfo.iter == config.getMaxEMIterations()) {
            status = EMAlgorithmStatus.FAILURE_MAX_ITERS_REACHED;
        }

        performPostEMOperations();

        return status;
    }

    /**
     * Perform variational inference for given model parameters
     *
     * @return the exit status
     */
    public EMAlgorithmStatus runExpectation() {
        config.validate();
        showIterationHeader();

        double prevEStepLikelihood;
        double latestEStepLikelihood = Double.NEGATIVE_INFINITY;
        final EMAlgorithmIterationInfo iterInfo = new EMAlgorithmIterationInfo(Double.NEGATIVE_INFINITY, 0, 1);
        /* disable copy ratio posterior calculation until bias estimation is stabilized */
        boolean updateCopyRatioPosteriors = false;

        while (iterInfo.iter <= config.getMaxEMIterations()) {

            /* cycle through E-step mean-field equations until they are satisfied to the desired degree */
            double maxPosteriorErrorNorm;
            double posteriorErrorNormReadDepth;
            double posteriorErrorNormSampleUnexplainedVariance;
            double posteriorErrorNormBias;
            double posteriorErrorNormCopyRatio;

            runRoutine(this::updateReadDepthLatentPosteriorExpectations, NOT_APPLICABLE_EXTRACTOR, "E_STEP_D", iterInfo);
            posteriorErrorNormReadDepth = iterInfo.errorNorm;

            if (config.getNumLatents() > 0) {
                runRoutine(this::updateBiasLatentPosteriorExpectations, NOT_APPLICABLE_EXTRACTOR, "E_STEP_Z", iterInfo);
                posteriorErrorNormBias = iterInfo.errorNorm;
            } else {
                posteriorErrorNormBias = 0;
            }

            if (config.sampleSpecificVarianceUpdateEnabled()) {
                runRoutine(this::updateSampleUnexplainedVariance, ITERATION_EXTRACTOR, "E_STEP_GAMMA", iterInfo);
                posteriorErrorNormSampleUnexplainedVariance = iterInfo.errorNorm;
            } else {
                posteriorErrorNormSampleUnexplainedVariance = 0;
            }

            if (updateCopyRatioPosteriors) {
                runRoutine(this::updateCopyRatioLatentPosteriorExpectations, NOT_APPLICABLE_EXTRACTOR, "E_STEP_C", iterInfo);
                posteriorErrorNormCopyRatio = iterInfo.errorNorm;
            } else {
                posteriorErrorNormCopyRatio = 0;
            }

            /* calculate the maximum change of posteriors in this cycle */
            maxPosteriorErrorNorm = Collections.max(Arrays.asList(
                    posteriorErrorNormReadDepth,
                    posteriorErrorNormSampleUnexplainedVariance,
                    posteriorErrorNormBias,
                    posteriorErrorNormCopyRatio));

            prevEStepLikelihood = latestEStepLikelihood;
            latestEStepLikelihood = iterInfo.getLogLikelihood();

            /* check convergence of the E-step */
            if (maxPosteriorErrorNorm < config.getPosteriorAbsTol() &&
                    FastMath.abs(latestEStepLikelihood - prevEStepLikelihood) < config.getLogLikelihoodTolerance()) {
                if (!config.copyRatioUpdateEnabled() || updateCopyRatioPosteriors) {
                    status = EMAlgorithmStatus.SUCCESS_POSTERIOR_CONVERGENCE;
                    break;
                }
            }

            /* if the change in log likelihood is small enough, start updating copy number posteriors */
            if (config.copyRatioUpdateEnabled() && !updateCopyRatioPosteriors &&
                    FastMath.abs(latestEStepLikelihood - prevEStepLikelihood) < config.getLogLikelihoodTolThresholdCRCalling()) {
                updateCopyRatioPosteriors = true;
                logger.info("Partial convergence achieved; starting to update copy ratio posteriors in the next iteration");
            }

            iterInfo.incrementIterationCount();

            if (config.isRunCheckpointingEnabled() && iterInfo.iter % config.getRunCheckpointingInterval() == 0) {
                final String posteriorOutputAbsolutePath = new File(config.getRunCheckpointingPath(),
                        String.format("%s_iter_%d", CoverageModelGlobalConstants.POSTERIOR_CHECKPOINT_PATH_PREFIX, iterInfo.iter)).getAbsolutePath();
                /* the following will automatically create the directory if it doesn't exist */
                savePosteriors(posteriorOutputAbsolutePath, CoverageModelEMWorkspace.PosteriorVerbosityLevel.BASIC);
            }
        }

        if (iterInfo.iter == config.getMaxEMIterations()) {
            status = EMAlgorithmStatus.FAILURE_MAX_ITERS_REACHED;
        }

        performPostEMOperations();

        return status;
    }

    private void performPostEMOperations() {
        logger.info("EM algorithm status: " + status.getMessage());
    }

    /**
     * E-step -- Update E[z] and E[z z^T] for all samples using the current estimate of model parameters
     */
    public SubroutineSignal updateBiasLatentPosteriorExpectations() {
        return workspace.updateBiasLatentPosteriorExpectations();
    }

    /**
     * E-step -- Update E[log(d_s)] and E[log(d_s)^2]
     */
    public SubroutineSignal updateReadDepthLatentPosteriorExpectations() {
        return workspace.updateReadDepthPosteriorExpectations();
    }

    /**
     * E-step -- Update gamma_s
     */
    public SubroutineSignal updateSampleUnexplainedVariance() {
        return workspace.updateSampleUnexplainedVariance();
    }

    /**
     * E-step -- Update E[log(c_{st})] and E[log(c_{st})^2]
     */
    public SubroutineSignal updateCopyRatioLatentPosteriorExpectations() {
        return workspace.updateCopyRatioPosteriorExpectations();
    }

    /**
     * E-step -- Update W
     */
    public SubroutineSignal updateBiasCovariates() {
        return workspace.updateBiasCovariates();
    }

    /**
     * M-step -- Update m
     */
    public SubroutineSignal updateTargetMeanBias(final boolean neglectBiasCovariates) {
        return workspace.updateMeanLogBias(neglectBiasCovariates);
    }

    /**
     * M-step -- Update Psi
     */
    public SubroutineSignal updateTargetUnexplainedVariance(final CoverageModelArgumentCollection.TargetSpecificVarianceUpdateMode targetSpecificVarianceUpdateMode) {
        return workspace.updateTargetUnexplainedVariance(targetSpecificVarianceUpdateMode);
    }

    /**
     * M-step -- Update ARD coefficients
     */
    public SubroutineSignal updateBiasCovariatesARDCoefficients() {
        return workspace.updateBiasCovariatesARDCoefficients();
    }

    /**
     * This routine is called after each EM iteration
     */
    private void finalizeIteration() {
        workspace.performGarbageCollection();
    }

    /**
     * Calculate the log likelihood
     *
     * @return log likelihood
     */
    public double getLogLikelihood() {
        return workspace.getLogLikelihood();
    }

    /**
     * Query the latest status of the EM algorithm (reflecting the data in the workspace)
     *
     * @return an instance of {@link EMAlgorithmStatus}
     */
    public EMAlgorithmStatus getStatus() { return status; }

    /**
     * Save model
     *
     * @param modelOutputPath where to save
     */
    public void saveModel(final String modelOutputPath) {
        workspace.writeModel(modelOutputPath);
    }

    /**
     * Save posteriors
     *
     * @param posteriorOutputPath where to save
     * @param verbosity the verbosity level
     */
    public void savePosteriors(final String posteriorOutputPath, final CoverageModelEMWorkspace.PosteriorVerbosityLevel verbosity) {
        workspace.writePosteriors(posteriorOutputPath, verbosity);
    }

    /**
     * This enum represents the status of the EM algorithm
     */
    public enum EMAlgorithmStatus {
        TBD(false, "Status is not determined yet."),
        SUCCESS_LIKELIHOOD_TOL(true, "Success -- converged in likelihood change tolerance."),
        SUCCESS_PARAMS_TOL(true, "Success -- converged in parameters change tolerance."),
        FAILURE_MAX_ITERS_REACHED(false, "Failure -- maximum iterations reached."),
        SUCCESS_POSTERIOR_CONVERGENCE(true, "Success -- converged in posterior and likelihood change tolerance.");

        private final boolean success;
        private final String message;

        EMAlgorithmStatus(final boolean success, final String message) {
            this.success = success;
            this.message = message;
        }

        public String getMessage() {
            return message;
        }

        public boolean isSuccessful() {
            return success;
        }
    }

    /**
     * This class stores basic info about each iteration of the EM algorithm
     */
    private final class EMAlgorithmIterationInfo {
        private double logLikelihood;
        private double errorNorm;
        private int iter;

        EMAlgorithmIterationInfo(final double logLikelihood, final double errorNorm, final int iter) {
            this.logLikelihood = logLikelihood;
            this.errorNorm = errorNorm;
            this.iter = iter;
        }

        public double getLogLikelihood() {
            return logLikelihood;
        }

        public double getErrorNorm() {
            return errorNorm;
        }

        public int getIterationCount() {
            return iter;
        }

        void incrementIterationCount() {
            iter++;
        }
    }
}
