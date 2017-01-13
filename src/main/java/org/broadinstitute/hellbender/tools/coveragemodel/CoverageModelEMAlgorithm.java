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
 * This class implements the EM algorithm for the GATK Bayesian coverage model. It is used for:
 *
 * (1) performing max likelihood estimation of the model parameters AND calculating posterior expectations
 *     via calling {@link CoverageModelEMAlgorithm {@link #runExpectationMaximization()}}
 *
 * (2) Calculating posterior expectations for the model parameters already existing in the workspace
 *     via calling {@link CoverageModelEMAlgorithm {@link #runExpectation()}}

 * (see CNV-methods.pdf for technical details).
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */

public final class CoverageModelEMAlgorithm<S extends AlleleMetadataProducer & CallStringProducer & ScalarProducer> {

    protected final Logger logger = LogManager.getLogger(CoverageModelEMAlgorithm.class);
    protected final CoverageModelEMParams params;
    protected final CoverageModelEMWorkspace<S> workspace;
    protected EMAlgorithmStatus status;

    /**
     * Public constructor.
     *
     * @param params specification of the hyper-parameters and other details of the algorithm
     * @param workspace workspace for actual calculations
     */
    public CoverageModelEMAlgorithm(@Nonnull final CoverageModelEMParams params,
                                    @Nonnull CoverageModelEMWorkspace<S> workspace) {
        this.params = Utils.nonNull(params, "Target coverage EM algorithm parameters must be non-null");
        params.validate();
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
        iterInfo.errorNorm = sig.getDouble("error_norm");
        iterInfo.logLikelihood = getLogLikelihood();
        showIterationInfo(iterInfo.iter, title, iterInfo.logLikelihood, iterInfo.errorNorm, misc);
    }

    /**
     * Run the full EM algorithm stack
     *
     * @return the exit status of the EM algorithm
     */
    public EMAlgorithmStatus runExpectationMaximization() {
        params.validate();
        if (params.adaptivePsiSolverModeSwitchingEnabled()) {
            /* if copy ratio posterior calling is enabled, the first few iterations need to be robust.
             * Therefore, we disable the target-resolved unexplained variance estimation (if enabled) and
             * enforce isotropic unexplained variance */
            if (params.adaptivePsiSolverModeSwitchingEnabled() &&
                    !this.params.getPsiUpdateMode().equals(CoverageModelEMParams.PsiUpdateMode.PSI_ISOTROPIC)) {
                this.params.setPsiPsiolverType(CoverageModelEMParams.PsiUpdateMode.PSI_ISOTROPIC);
                logger.info("Overriding the requested unexplained variance solver to " +
                        CoverageModelEMParams.PsiUpdateMode.PSI_ISOTROPIC.name());
            }
        }

        showIterationHeader();

        double prevEStepLikelihood = Double.NEGATIVE_INFINITY;
        double latestEStepLikelihood = Double.NEGATIVE_INFINITY;
        double prevMStepLikelihood = Double.NEGATIVE_INFINITY;
        double latestMStepLikelihood = Double.NEGATIVE_INFINITY;
        final EMAlgorithmIterationInfo iterInfo = new EMAlgorithmIterationInfo(Double.NEGATIVE_INFINITY, 0, 1);
        boolean updateCopyRatioPosteriors = false;
        boolean paramEstimationConverged = false;
        boolean performMStep = true;

        while (iterInfo.iter <= params.getMaxEMIterations()) {

            /* cycle through E-step mean-field equations until they are satisfied to the desired degree */
            double maxPosteriorErrorNorm = 0;
            int iterEStep = 0;

            while (iterEStep < params.getMaxEStepCycles()) {
                double posteriorErrorNormReadDepth, posteriorErrorNormSampleUnexplainedVariance,
                        posteriorErrorNormBias, posteriorErrorNormCopyRatio;

                runRoutine(this::updateReadDepthLatentPosteriorExpectations, s -> "N/A", "E_STEP_D", iterInfo);
                posteriorErrorNormReadDepth = iterInfo.errorNorm;

                runRoutine(this::updateBiasLatentPosteriorExpectations, s -> "N/A", "E_STEP_Z", iterInfo);
                posteriorErrorNormBias = iterInfo.errorNorm;

                if (params.gammaUpdateEnabled()) {
                    runRoutine(this::updateSampleUnexplainedVariance,
                            s -> "iters: " + s.getInteger("iterations"), "E_STEP_GAMMA", iterInfo);
                    posteriorErrorNormSampleUnexplainedVariance = iterInfo.errorNorm;
                } else {
                    posteriorErrorNormSampleUnexplainedVariance = 0;
                }

                if (updateCopyRatioPosteriors) {
                    runRoutine(this::updateCopyRatioLatentPosteriorExpectations, s -> "N/A", "E_STEP_C", iterInfo);
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

                /* check convergence of the E-step */
                if (maxPosteriorErrorNorm < params.getPosteriorAbsTol()) {
                    break;
                }

                iterEStep++;
            }

            if (maxPosteriorErrorNorm > params.getPosteriorAbsTol()) {
                logger.warn("E-step cycles did not fully converge. Increase the maximum number of E-step cycles." +
                        " Continuing...");
            }

            prevEStepLikelihood = latestEStepLikelihood;
            latestEStepLikelihood = iterInfo.logLikelihood;

            /* parameter estimation */
            if (performMStep && !paramEstimationConverged) {
                int iterMStep = 0;
                double maxParamErrorNorm;

                /* sequential M-steps */
                while (iterMStep < params.getMaxMStepCycles()) {
                    double errorNormMeanTargetBias = 0, errorNormUnexplainedVariance = 0, errorNormPrincipalMap = 0;

                    if (iterInfo.iter == 0) {
                        /* neglect the contribution from principal components in the first iteration */
                        runRoutine(() -> updateTargetMeanBias(true), s -> "N/A", "M_STEP_M", iterInfo);
                    } else {
                        runRoutine(() -> updateTargetMeanBias(false), s -> "N/A", "M_STEP_M", iterInfo);
                    }
                    errorNormMeanTargetBias = iterInfo.errorNorm;

                    if (params.psiUpdateEnabled()) {
                        runRoutine(this::updateTargetUnexplainedVariance, s -> "iters: " + s.getInteger("iterations"),
                                "M_STEP_PSI", iterInfo);
                        errorNormUnexplainedVariance = iterInfo.errorNorm;
                    }

                    if (params.fourierRegularizationEnabled()) {
                        runRoutine(this::updateBiasCovariates, s -> "iters: " + s.getInteger("iterations"),
                                "M_STEP_W", iterInfo);
                    } else {
                        runRoutine(this::updateBiasCovariates, s -> "N/A",
                                "M_STEP_W", iterInfo);
                    }
                    errorNormPrincipalMap = iterInfo.errorNorm;

                    maxParamErrorNorm = Collections.max(Arrays.asList(errorNormMeanTargetBias,
                            errorNormUnexplainedVariance, errorNormPrincipalMap));

                    /* check convergence of parameter estimation */
                    if (updateCopyRatioPosteriors && maxParamErrorNorm < params.getParameterEstimationAbsoluteTolerance()) {
                        status = EMAlgorithmStatus.SUCCESS_PARAMS_TOL;
                        paramEstimationConverged = true;
                    }
                    iterMStep++;
                }

                prevMStepLikelihood = latestMStepLikelihood;
                latestMStepLikelihood = iterInfo.logLikelihood;

                /* if the likelihood has increased and the increment is small, start updating copy number posteriors */
                if (params.copyRatioUpdateEnabled() &&
                        !updateCopyRatioPosteriors &&
                        (latestMStepLikelihood - prevMStepLikelihood) > 0 &&
                        (latestMStepLikelihood - prevMStepLikelihood) < params.getLogLikelihoodTolThresholdCRCalling()) {
                    updateCopyRatioPosteriors = true;
                    logger.info("Partial convergence achieved; starting to update copy ratio posteriors (if enabled)" +
                            " and sample-specific unexplained variance (if enabled) in the next iteration");
                    if (params.adaptivePsiSolverModeSwitchingEnabled()) {
                        params.setPsiPsiolverType(CoverageModelEMParams.PsiUpdateMode.PSI_TARGET_RESOLVED);
                    }
                }
            }

            /* check convergence in log likelihood change */
            if (FastMath.abs(latestMStepLikelihood - prevMStepLikelihood) < params.getLogLikelihoodTolerance()) {
                /* make sure that we have either already called copy ratio posteriors, or we are not required to */
                if (!params.copyRatioUpdateEnabled() || updateCopyRatioPosteriors) {
                    status = EMAlgorithmStatus.SUCCESS_LIKELIHOOD_TOL;
                    break;
                }
            } else if (iterInfo.iter == params.getMaxEMIterations() - 2) {
                performMStep = false; /* so that we end with an E-step */
            }

            iterInfo.increaseIterationCount();

            finalizeIteration();

            /* checkpointing */
            if (params.isRunCheckpointingEnabled() && iterInfo.iter % params.getRunCheckpointingInterval() == 0) {
                final String modelOutputAbsolutePath = new File(params.getRunCheckpointingPath(),
                        String.format("%s_iter_%d", CoverageModelGlobalConstants.MODEL_CHECKPOINT_PATH_PREFIX, iterInfo.iter)).getAbsolutePath();
                final String posteriorOutputAbsolutePath = new File(params.getRunCheckpointingPath(),
                        String.format("%s_iter_%d", CoverageModelGlobalConstants.POSTERIOR_CHECKPOINT_PATH_PREFIX, iterInfo.iter)).getAbsolutePath();
                saveModel(modelOutputAbsolutePath);
                savePosteriors(posteriorOutputAbsolutePath, PosteriorVerbosityLevel.BASIC);
            }
        }

        if (iterInfo.iter == params.getMaxEMIterations()) {
            status = EMAlgorithmStatus.FAILURE_MAX_ITERS_REACHED;
        }

        performPostEMOperations();

        return status;
    }

    /**
     *
     */
    public EMAlgorithmStatus runExpectation() {
        params.validate();
        showIterationHeader();

        double prevEStepLikelihood = Double.NEGATIVE_INFINITY;
        double latestEStepLikelihood = Double.NEGATIVE_INFINITY;
        final EMAlgorithmIterationInfo iterInfo = new EMAlgorithmIterationInfo(Double.NEGATIVE_INFINITY, 0, 1);
        /* disable copy ratio posterior calculation until bias estimation is stabilized */
        boolean updateCopyRatioPosteriors = false;

        while (iterInfo.iter <= params.getMaxEMIterations()) {

            /* cycle through E-step mean-field equations until they are satisfied to the desired degree */
            double maxPosteriorErrorNorm;
            double posteriorErrorNormReadDepth, posteriorErrorNormSampleUnexplainedVariance,
                    posteriorErrorNormBias, posteriorErrorNormCopyRatio;

            runRoutine(this::updateReadDepthLatentPosteriorExpectations, s -> "N/A", "E_STEP_D", iterInfo);
            posteriorErrorNormReadDepth = iterInfo.errorNorm;

            runRoutine(this::updateBiasLatentPosteriorExpectations, s -> "N/A", "E_STEP_Z", iterInfo);
            posteriorErrorNormBias = iterInfo.errorNorm;

            if (params.gammaUpdateEnabled()) {
                runRoutine(this::updateSampleUnexplainedVariance,
                        s -> "iters: " + s.getInteger("iterations"), "E_STEP_GAMMA", iterInfo);
                posteriorErrorNormSampleUnexplainedVariance = iterInfo.errorNorm;
            } else {
                posteriorErrorNormSampleUnexplainedVariance = 0;
            }

            if (updateCopyRatioPosteriors) {
                runRoutine(this::updateCopyRatioLatentPosteriorExpectations, s -> "N/A", "E_STEP_C", iterInfo);
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
            latestEStepLikelihood = iterInfo.logLikelihood;

            /* check convergence of the E-step */
            if (maxPosteriorErrorNorm < params.getPosteriorAbsTol() &&
                    FastMath.abs(latestEStepLikelihood - prevEStepLikelihood) < params.getLogLikelihoodTolerance()) {
                if (!params.copyRatioUpdateEnabled() || updateCopyRatioPosteriors) {
                    status = EMAlgorithmStatus.SUCCESS_POSTERIOR_CONVERGENCE;
                    break;
                }
            }

            /* if the likelihood has increased and the increment is small, start updating copy number posteriors */
            if (params.copyRatioUpdateEnabled() &&
                    !updateCopyRatioPosteriors &&
                    (latestEStepLikelihood - prevEStepLikelihood) > 0 &&
                    (latestEStepLikelihood - prevEStepLikelihood) < params.getLogLikelihoodTolThresholdCRCalling()) {
                updateCopyRatioPosteriors = true;
                logger.info("Partial convergence achieved; will start calling copy ratio posteriors");
            }

            iterInfo.increaseIterationCount();

            if (params.isRunCheckpointingEnabled() && iterInfo.iter % params.getRunCheckpointingInterval() == 0) {
                final String posteriorOutputAbsolutePath = new File(params.getRunCheckpointingPath(),
                        String.format("%s_iter_%d", CoverageModelGlobalConstants.POSTERIOR_CHECKPOINT_PATH_PREFIX, iterInfo.iter)).getAbsolutePath();
                /* the following will automatically create the directory if it doesn't exist */
                savePosteriors(posteriorOutputAbsolutePath, PosteriorVerbosityLevel.BASIC);
            }
        }

        if (iterInfo.iter == params.getMaxEMIterations()) {
            status = EMAlgorithmStatus.FAILURE_MAX_ITERS_REACHED;
        }

        performPostEMOperations();

        return status;
    }

    private void performPostEMOperations() {
        logger.info("EM algorithm status: " + status.message);
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
     * M-step -- Update m
     */
    public SubroutineSignal updateTargetMeanBias(final boolean neglectBiasCovariates) {
        return workspace.updateMeanLogBias(neglectBiasCovariates);
    }

    /**
     * M-step -- Update Psi
     */
    public SubroutineSignal updateTargetUnexplainedVariance() {
        return workspace.updateTargetUnexplainedVariance();
    }

    /**
     * M-step -- Update W
     */
    public SubroutineSignal updateBiasCovariates() {
        return workspace.updateBiasCovariates();
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
        workspace.saveModel(modelOutputPath);
    }

    /**
     * Save posteriors
     *
     * @param posteriorOutputPath where to save
     * @param verbosity the verbosity level
     */
    public void savePosteriors(final String posteriorOutputPath, final PosteriorVerbosityLevel verbosity) {
        workspace.savePosteriors(posteriorOutputPath, verbosity);
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

        final boolean success;
        final String message;

        EMAlgorithmStatus(final boolean success, final String message) {
            this.success = success;
            this.message = message;
        }
    }

    /**
     * This class stores basic info about each iteration of the EM algorithm
     */
    private final class EMAlgorithmIterationInfo {
        double logLikelihood, errorNorm;
        int iter;

        EMAlgorithmIterationInfo(final double logLikelihood, final double errorNorm, final int iter) {
            this.logLikelihood = logLikelihood;
            this.errorNorm = errorNorm;
            this.iter = iter;
        }

        void increaseIterationCount() {
            iter++;
        }
    }

}
