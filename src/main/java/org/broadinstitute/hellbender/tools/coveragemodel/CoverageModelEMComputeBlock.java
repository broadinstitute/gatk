package org.broadinstitute.hellbender.tools.coveragemodel;

import com.google.common.annotations.VisibleForTesting;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.ImmutableTriple;
import org.apache.commons.math3.exception.NoBracketingException;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.tools.coveragemodel.cachemanager.Duplicable;
import org.broadinstitute.hellbender.tools.coveragemodel.cachemanager.DuplicableNDArray;
import org.broadinstitute.hellbender.tools.coveragemodel.cachemanager.DuplicableNumber;
import org.broadinstitute.hellbender.tools.coveragemodel.cachemanager.ImmutableComputableGraph;
import org.broadinstitute.hellbender.tools.coveragemodel.math.RobustBrentSolver;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;
import org.nd4j.linalg.indexing.NDArrayIndex;
import org.nd4j.linalg.ops.transforms.Transforms;

import javax.annotation.Nonnull;
import javax.annotation.Nullable;
import java.io.Serializable;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ForkJoinPool;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * This class represents an immutable block of data containers, query methods and cloners corresponding to
 * a partition of the target space.
 *
 * TODO github/gatk-protected issue #853 -- logging in spark mode (log4j is not serializable)
 * TODO github/gatk-protected issue #853 -- use instrumentation to measure memory consumption

 * @implNote Methods that manipulate INDArrays must make sure to leave queried values from
 * {@link CoverageModelEMComputeBlock#icg} unchanged. For example, to calculate the A.B.C + D,
 * the correct template is first fetch the data from ICG as follows:
 *
 * // query cached values
 * final INDArray A = getINDArrayFromCache(CoverageModelICGCacheNode.A);
 * final INDArray B = getINDArrayFromCache(CoverageModelICGCacheNode.B);
 * final INDArray C = getINDArrayFromCache(CoverageModelICGCacheNode.C);
 * final INDArray D = getINDArrayFromCache(CoverageModelICGCacheNode.D);
 *
 * and then call:
 *
 * // calculations
 * final INDArray result_1 = A.mul(B).muli(C).addi(D); // OK
 *                           ---1--- (a copy is instantiated; call it tmp)
 *                                ----2---  (tmp is modified in-place)
 *                                    -----3---- (tmp is modified in-place)
 *
 * (note the difference between mul(...) and muli(...)
 *
 * Here, A and multiplied by B to create a new INDArray; subsequent operations can be performed
 * in-place to avoid unnecessary memory allocation/deallocation and garbage generation.
 *
 * It is, however, illegal to perform the following operation:
 *
 * final INDArray result_3 = A.muli(B).muli(C).addi(D); // NOT OK
 *                           ---1--- (A is modified in-place!)
 *
 * Here, A is modified in-place which corrupts the data in {@link CoverageModelEMComputeBlock#icg}
 *
 * The code reviewer must CAREFULLY check all of the matrix computations to ensure that the such
 * case does not occur, and suggest optimizations if noted.
 *
 * Also, since we do not have an immutable INDArray (or a wrapper for that),
 * it is the user's responsibility to ensure that the queried INDArray's are not mutated.
 * The safest approach is to call dup() on cached queries, however, it leads to
 * poor memory performance.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class CoverageModelEMComputeBlock {

    /**
     * An instance of {@link ImmutableComputableGraph} for automatic evaluation, caching, and bookkeeping
     * of most evaluations
     */
    private final ImmutableComputableGraph icg;

    /**
     * The index range of targets included in this compute block
     */
    private final LinearSpaceBlock targetBlock;

    private final int numSamples, numLatents, numTargets;

    private final boolean ardEnabled, biasCovariatesEnabled;

    /**
     * The latest signal from calling an M-step subroutine
     */
    private final SubroutineSignal latestMStepSignal;

    /**
     * Immutable computable graph cache nodes
     */
    public enum CoverageModelICGCacheNode {
        n_st("Read count"),
        M_st("Learning mask (only 0 and 1 values)"),
        err_st("Mapping error probability"),
        m_t("Mean log bias"),
        Psi_t("Target-specific unexplained variance of log bias"),
        W_tl("Mean of bias covariates (principal components)"),
        alpha_l("Precision of bias covariates"),
        log_c_st("Posterior mean log copy ratio (or copy number)"),
        var_log_c_st("Posterior variance of log copy ratio (or copy number)"),
        log_d_s("Posterior mean of log read depth"),
        var_log_d_s("Posterior variance of log read depth"),
        gamma_s("Posterior mean of sample-specific unexplained variance of the log bias"),
        z_sl("Posterior mean of log bias continuous latent variables"),
        zz_sll("Posterior mean of pair-product of log bias continuous latent variables"),
        F_W_tl("Fourier transform of bias covariates in the target space"),
        log_n_st("Log read count"),
        Sigma_st("Poisson noise"),
        sum_M_t("Sample-summed learning mask"),
        sum_M_s("Target-summed learning mask"),
        Wz_st("W_{t\\mu} E[z_{s\\mu}]"),
        WzzWT_st("W_{t\\mu} W_{t\\nu} E[z_{s\\mu} z_{s\\nu}]"),
        tot_Psi_st("The sum of target-specific and sample-specific unexplained variance of the log bias"),
        Delta_st("\\Delta_{st} (refer to the technical white paper)"),
        Delta_PCA_st("\\Delta_PCA_{st} (refer to the technical white paper)"),
        loglike_normalization_s("Additive normalization factor of the variational log likelihood"),
        M_Psi_inv_st("Masked precision of log bias"),
        v_tl("v_{t\\mu}} (refer to the technical white paper)"),
        Q_tll("Q_{t\\mu\\nu} (refer to the technical white paper)"),
        sum_Q_ll("\\sum_t Q_{t\\mu\\nu} (refer to the technical white paper)"),
        B_st("B_{st} (refer to the technical white paper)"),
        loglike_unreg("Contribution of this block of targets to the model log likelihood (w/o Fourier regularizer)"),
        loglike_reg("Contribution of this block of targets to the model log likelihood (w/ Fourier regularizer)");

        public final String description;

        CoverageModelICGCacheNode(final String description) {
            this.description = description;
        }
    }

    /**
     * Immutable computable graph cache tags
     */
    public enum CoverageModelICGCacheTag {
        E_STEP_D("Cache nodes to be updated for the E-step for read depth"),
        E_STEP_C("Cache nodes to be updated for the E-step for copy ratio (or copy number)"),
        E_STEP_Z("Cache nodes to be updated for the E-step for bias continuous latent variables"),
        E_STEP_GAMMA("Cache nodes to be updated for the E-step for sample-specific unexplained variance of the log bias"),
        E_STEP_W_REG("Cache nodes to be updated for the E-step for bias covariates (w/ regularization)"),
        E_STEP_W_UNREG("Cache nodes to be updated for the E-step for bias covariates (w/o regularization)"),
        M_STEP_M("Cache nodes to be updated for the M-step for mean log bias"),
        M_STEP_PSI("Cache nodes to be updated for the M-step for target-specific unexplained variance"),
        LOGLIKE_UNREG("Cache nodes to be updated for log likelihood calculation (w/o regularization)"),
        LOGLIKE_REG("Cache nodes to be updated for log likelihood calculation (w/ regularization)");
        public final String description;

        CoverageModelICGCacheTag(final String description) {
            this.description = description;
        }
    }

    /**
     * This method creates an empty instance of {@link ImmutableComputableGraph} for worker-side
     * calculations.
     *
     * @param biasCovariatesEnabled whether or not create nodes relating to bias covariates
     * @param ardEnabled whether or not create nodes for ARD calculation of bias bias covariates
     * @return an instance of {@link ImmutableComputableGraph}
     */
    private static ImmutableComputableGraph createEmptyCacheGraph(final boolean biasCovariatesEnabled,
                                                                  final boolean ardEnabled) {
        Utils.validateArg(!ardEnabled || biasCovariatesEnabled, "If ARD is enabled, bias covariates must be" +
                " enabled as well");
        final ImmutableComputableGraph.ImmutableComputableGraphBuilder cgbuilder =
                ImmutableComputableGraph.builder();
        /*
         * Data nodes
         */
        cgbuilder
                /* raw read counts */
                .addPrimitiveNode(CoverageModelICGCacheNode.n_st.name(),
                        new String[]{},
                        new DuplicableNDArray())
                /* mask */
                .addPrimitiveNode(CoverageModelICGCacheNode.M_st.name(),
                        new String[]{},
                        new DuplicableNDArray())
                /* mapping error probability */
                .addPrimitiveNode(CoverageModelICGCacheNode.err_st.name(),
                        new String[]{},
                        new DuplicableNDArray());

        /*
         * Model parameters
         */
        cgbuilder
                /* mean log bias */
                .addPrimitiveNode(CoverageModelICGCacheNode.m_t.name(),
                        new String[]{},
                        new DuplicableNDArray())
                /* unexplained variance */
                .addPrimitiveNode(CoverageModelICGCacheNode.Psi_t.name(),
                        new String[]{},
                        new DuplicableNDArray());

        /* if ARD is enabled, add a node for ARD coefficients */
        if (ardEnabled) {
            cgbuilder
                    /* precision of bias covariates */
                    .addPrimitiveNode(CoverageModelICGCacheNode.alpha_l.name(),
                            new String[]{},
                            new DuplicableNDArray());
        }

        /*
         * Externally determined computable nodes (all of latent posterior expectations + etc)
         */
        cgbuilder
                /* E[log(c_{st})] */
                .addComputableNode(CoverageModelICGCacheNode.log_c_st.name(),
                        new String[]{},
                        new String[]{},
                        null, true)
                /* var[log(c_{st})] */
                .addComputableNode(CoverageModelICGCacheNode.var_log_c_st.name(),
                        new String[]{},
                        new String[]{},
                        null, true)
                /* E[log(d_s)] */
                .addComputableNode(CoverageModelICGCacheNode.log_d_s.name(),
                        new String[]{},
                        new String[]{},
                        null, true)
                /* var[log(d_s)] */
                .addComputableNode(CoverageModelICGCacheNode.var_log_d_s.name(),
                        new String[]{},
                        new String[]{},
                        null, true)
                /* E[\gamma_s] */
                .addComputableNode(CoverageModelICGCacheNode.gamma_s.name(),
                        new String[]{},
                        new String[]{},
                        null, true);

        if (biasCovariatesEnabled) {
            cgbuilder
                    /* mean of bias covariates */
                    .addComputableNode(CoverageModelICGCacheNode.W_tl.name(),
                            new String[]{},
                            new String[]{},
                            null, true)
                    /* E[z_{sm}] */
                    .addComputableNode(CoverageModelICGCacheNode.z_sl.name(),
                            new String[]{},
                            new String[]{},
                            null, true)
                    /* E[z_{sm} z_{sn}] */
                    .addComputableNode(CoverageModelICGCacheNode.zz_sll.name(),
                            new String[]{},
                            new String[]{},
                            null, true);
        }

        /*
         * Automatically computable nodes
         */
        cgbuilder
                /* log read counts */
                .addComputableNode(CoverageModelICGCacheNode.log_n_st.name(),
                        new String[]{
                                CoverageModelICGCacheTag.M_STEP_M.name(),
                                CoverageModelICGCacheTag.E_STEP_D.name()},
                        new String[]{
                                CoverageModelICGCacheNode.n_st.name(),
                                CoverageModelICGCacheNode.M_st.name()},
                        calculate_log_n_st, true)
                /* Poisson noise */
                .addComputableNode(CoverageModelICGCacheNode.Sigma_st.name(),
                        new String[]{},
                        new String[]{
                                CoverageModelICGCacheNode.n_st.name(),
                                CoverageModelICGCacheNode.M_st.name()},
                        calculate_Sigma_st, true)
                /* \sum_s M_{st} */
                .addComputableNode(CoverageModelICGCacheNode.sum_M_t.name(),
                        new String[]{},
                        new String[]{CoverageModelICGCacheNode.M_st.name()},
                        calculate_sum_M_t, true)
                /* \sum_t M_{st} */
                .addComputableNode(CoverageModelICGCacheNode.sum_M_s.name(),
                        new String[]{
                                CoverageModelICGCacheTag.LOGLIKE_UNREG.name(),
                                CoverageModelICGCacheTag.LOGLIKE_REG.name()},
                        new String[]{CoverageModelICGCacheNode.M_st.name()},
                        calculate_sum_M_s, true)
                /* \Psi_{st} = \Psi_t + \Sigma_{st} + E[\gamma_s] */
                .addComputableNode(CoverageModelICGCacheNode.tot_Psi_st.name(),
                        new String[]{},
                        new String[]{
                                CoverageModelICGCacheNode.Sigma_st.name(),
                                CoverageModelICGCacheNode.Psi_t.name(),
                                CoverageModelICGCacheNode.gamma_s.name()},
                        calculate_tot_Psi_st, true)
                /* log(n_{st}) - E[log(c_{st})] - E[log(d_s)] - m_t */
                .addComputableNode(CoverageModelICGCacheNode.Delta_st.name(),
                        new String[]{CoverageModelICGCacheTag.E_STEP_Z.name()},
                        new String[]{
                                CoverageModelICGCacheNode.log_n_st.name(),
                                CoverageModelICGCacheNode.log_c_st.name(),
                                CoverageModelICGCacheNode.log_d_s.name(),
                                CoverageModelICGCacheNode.m_t.name()},
                        calculate_Delta_st, true)
                /* \sum_{t} M_{st} \log(\Psi_{st}) */
                .addComputableNode(CoverageModelICGCacheNode.loglike_normalization_s.name(),
                        new String[]{
                                CoverageModelICGCacheTag.LOGLIKE_REG.name(),
                                CoverageModelICGCacheTag.LOGLIKE_UNREG.name()},
                        new String[]{
                                CoverageModelICGCacheNode.M_st.name(),
                                CoverageModelICGCacheNode.tot_Psi_st.name(),
                                CoverageModelICGCacheNode.log_n_st.name()},
                        calculate_loglike_normalization_s, true)
                /* M_{st} \Psi_{st}^{-1} */
                .addComputableNode(CoverageModelICGCacheNode.M_Psi_inv_st.name(),
                        new String[]{
                                CoverageModelICGCacheTag.E_STEP_W_UNREG.name(),
                                CoverageModelICGCacheTag.E_STEP_W_REG.name(),
                                CoverageModelICGCacheTag.M_STEP_M.name(),
                                CoverageModelICGCacheTag.E_STEP_Z.name(),
                                CoverageModelICGCacheTag.E_STEP_D.name(),
                                CoverageModelICGCacheTag.E_STEP_C.name()},
                        new String[]{
                                CoverageModelICGCacheNode.M_st.name(),
                                CoverageModelICGCacheNode.tot_Psi_st.name()},
                        calculate_M_Psi_inv_st, true)
                /* log likelihood (w/o regularization) */
                .addComputableNode(CoverageModelICGCacheNode.loglike_unreg.name(),
                        new String[]{CoverageModelICGCacheTag.LOGLIKE_UNREG.name()},
                        new String[]{
                                CoverageModelICGCacheNode.B_st.name(),
                                CoverageModelICGCacheNode.M_Psi_inv_st.name(),
                                CoverageModelICGCacheNode.loglike_normalization_s.name()},
                        calculate_loglike_unreg, true);

        /* nodes specific to bias covariates */
        if (biasCovariatesEnabled) {
            cgbuilder
                    /* E[W] E[z_s] */
                    .addComputableNode(CoverageModelICGCacheNode.Wz_st.name(),
                            new String[]{CoverageModelICGCacheTag.M_STEP_M.name(),
                                    CoverageModelICGCacheTag.E_STEP_D.name(),
                                    CoverageModelICGCacheTag.E_STEP_C.name()},
                            new String[]{
                                    CoverageModelICGCacheNode.W_tl.name(),
                                    CoverageModelICGCacheNode.z_sl.name()},
                            calculate_Wz_st, true)
                    /* (E[z_s z_s^T] E[W W^T])_{tt} */
                    .addComputableNode(CoverageModelICGCacheNode.WzzWT_st.name(),
                            new String[]{},
                            new String[]{
                                    CoverageModelICGCacheNode.W_tl.name(),
                                    CoverageModelICGCacheNode.zz_sll.name()},
                            calculate_WzzWT_st, true)
                    /* v_{t\mu} */
                    .addComputableNode(CoverageModelICGCacheNode.v_tl.name(),
                            new String[]{
                                    CoverageModelICGCacheTag.E_STEP_W_REG.name(),
                                    CoverageModelICGCacheTag.E_STEP_W_UNREG.name()},
                            new String[]{
                                    CoverageModelICGCacheNode.M_Psi_inv_st.name(),
                                    CoverageModelICGCacheNode.Delta_st.name(),
                                    CoverageModelICGCacheNode.z_sl.name()},
                            calculate_v_tl, false)
                    /* Q_{t\mu\nu} */
                    .addComputableNode(CoverageModelICGCacheNode.Q_tll.name(),
                            new String[]{
                                    CoverageModelICGCacheTag.E_STEP_W_REG.name(),
                                    CoverageModelICGCacheTag.E_STEP_W_UNREG.name()},
                            new String[]{
                                    CoverageModelICGCacheNode.M_Psi_inv_st.name(),
                                    CoverageModelICGCacheNode.zz_sll.name()},
                            calculate_Q_tll, false)
                    /* B_{st} */
                    .addComputableNode(CoverageModelICGCacheNode.B_st.name(),
                            new String[]{
                                    CoverageModelICGCacheTag.M_STEP_PSI.name(),
                                    CoverageModelICGCacheTag.E_STEP_GAMMA.name()},
                            new String[]{
                                    CoverageModelICGCacheNode.Delta_st.name(),
                                    CoverageModelICGCacheNode.var_log_c_st.name(),
                                    CoverageModelICGCacheNode.var_log_d_s.name(),
                                    CoverageModelICGCacheNode.WzzWT_st.name(),
                                    CoverageModelICGCacheNode.Wz_st.name()},
                            calculate_B_st_with_bias_covariates, true)
                    /* M_{st} . (log(n_{st}) - E[log(c_{st})] - E[log(d_s)]) - mean --- externally computed */
                    .addComputableNode(CoverageModelICGCacheNode.Delta_PCA_st.name(),
                            new String[]{},
                            new String[]{
                                    CoverageModelICGCacheNode.log_n_st.name(),
                                    CoverageModelICGCacheNode.log_c_st.name(),
                                    CoverageModelICGCacheNode.log_d_s.name()},
                            null, true);

        } else { /* no bias covariates */
            cgbuilder
                    /* B_{st} */
                    .addComputableNode(CoverageModelICGCacheNode.B_st.name(),
                            new String[]{
                                    CoverageModelICGCacheTag.M_STEP_PSI.name(),
                                    CoverageModelICGCacheTag.E_STEP_GAMMA.name()},
                            new String[]{
                                    CoverageModelICGCacheNode.Delta_st.name(),
                                    CoverageModelICGCacheNode.var_log_c_st.name(),
                                    CoverageModelICGCacheNode.var_log_d_s.name()},
                            calculate_B_st_without_bias_covariates, true);
        }

        // TODO github/gatk-protected issue #701 -- this class is part of the upcoming CNV-avoiding regularizer
        //
        //            /* FFT[W] */
        //            .addComputableNode(CoverageModelICGCacheNode.F_W_tl.name(),
        //                   new String[]{},
        //                   new String[]{CoverageModelICGCacheNode.W_tl.name()},
        //                   null, true);
        //
        //            /* log likelihood (w/ regularization) */
        //            .addComputableNode(CoverageModelICGCacheNode.loglike_reg.name(),
        //                    new String[]{CoverageModelICGCacheTag.LOGLIKE_REG.name()},
        //                    new String[]{
        //                            CoverageModelICGCacheNode.B_st.name(),
        //                            CoverageModelICGCacheNode.M_Psi_inv_st.name(),
        //                            CoverageModelICGCacheNode.loglike_normalization_s.name(),
        //                            CoverageModelICGCacheNode.W_tl.name(),
        //                            CoverageModelICGCacheNode.F_W_tl.name(),
        //                            CoverageModelICGCacheNode.zz_sll.name()},
        //                    calculate_loglike_reg, true);
        //            /* \sum_t Q_{t\mu\nu} */
        //            .addComputableNode(CoverageModelICGCacheNode.sum_Q_ll.name(),
        //        new String[]{CoverageModelICGCacheTag.E_STEP_W_REG.name()},
        //        new String[]{CoverageModelICGCacheNode.Q_tll.name()},
        //        calculate_sum_Q_ll, true);

        return cgbuilder.build();
    }

    /**
     * Private constructor (for cloners)
     *
     * @param targetBlock target space block
     * @param numSamples number of samples
     * @param numLatents dimension of bias latent space
     * @param ardEnabled enable/disable ARD for bias covariates
     * @param icg an instance of {@link ImmutableComputableGraph}
     */
    private CoverageModelEMComputeBlock(@Nonnull final LinearSpaceBlock targetBlock,
                                        final int numSamples, final int numLatents,
                                        final boolean ardEnabled,
                                        @Nonnull final ImmutableComputableGraph icg,
                                        @Nullable SubroutineSignal latestMStepSignal) {
        this.numSamples = ParamUtils.isPositive(numSamples, "Number of samples must be positive.");
        this.numLatents = ParamUtils.isPositiveOrZero(numLatents, "Dimension of the latent space must be non-negative.");
        this.targetBlock = Utils.nonNull(targetBlock, "Target space block identifier can not be null.");
        this.icg = Utils.nonNull(icg, "The immutable computable graph can not be null.");
        this.ardEnabled = ardEnabled;
        this.biasCovariatesEnabled = numLatents > 0;
        this.latestMStepSignal = latestMStepSignal;
        this.numTargets = targetBlock.getNumTargets();
    }

    /**
     * Public constructor.
     *
     * @param targetBlock target space block
     * @param numSamples number of samples
     * @param numLatents dimension of the bias latent space
     * @param ardEnabled enable/disable ARD for bias covariates
     */
    public CoverageModelEMComputeBlock(@Nonnull final LinearSpaceBlock targetBlock,
                                       final int numSamples, final int numLatents,
                                       final boolean ardEnabled) {
        this(targetBlock, numSamples, numLatents, ardEnabled,
                createEmptyCacheGraph(numLatents > 0, ardEnabled),
                SubroutineSignal.builder().build());
    }

    public LinearSpaceBlock getTargetSpaceBlock() {
        return targetBlock;
    }

    public SubroutineSignal getLatestMStepSignal() { return latestMStepSignal; }

    /**
     * Fetch the value of a cache node and cast it to an INDArray
     *
     * Note: up-casting is not checked
     *
     * @param key key of the cache node
     * @return an INDArray
     */
    public INDArray getINDArrayFromCache(final CoverageModelICGCacheNode key) {
        return ((DuplicableNDArray)icg.getValueWithRequiredEvaluations(key.name())).value();
    }

    /**
     * Fetch the value of a cache node and cast it to a primitive double
     *
     * Note: up-casting is not checked

     * @param key key of the cache
     * @return a double value
     */
    public double getDoubleFromCache(final CoverageModelICGCacheNode key) {
        return ((DuplicableNumber) icg.getValueWithRequiredEvaluations(key.name())).value().doubleValue();
    }

    private void assertBiasCovariatesEnabled() {
        Utils.validateArg(biasCovariatesEnabled, "Bias covariates are disabled but a method was invoked that" +
                " is related to bias covariates");
    }

    private void assertARDEnabled() {
        Utils.validateArg(ardEnabled, "ARD is disabled but a method was invoked that is related to ARD");
    }

    /**
     * Calculates the contribution of the target-space block represented by this compute node
     * to the E-step for $z_{s\mu}$ w/o regularization. The final calculation is performed on the
     * driver node in {@link CoverageModelEMWorkspace}.
     *
     * The result is the following pair of INDArrays:
     *
     *      contribGMatrix_{s\mu,\nu} = \sum_{t \in targetBlock} E[W_{\mu t} W_{t\nu}] [diag(M_st \Psi_st)]
     *
     *      contribZ_{\mu,s} = \sum_{t \in targetBlock} E[W_{\mu t}] [diag(M_st E[\Psi_st])] \Delta_st
     *      (note the order of indices)
     *
     * @return an {@link ImmutablePair} of contribGMatrix (left) and contribZ (right)
     */
    public ImmutablePair<INDArray, INDArray> getBiasLatentPosteriorDataUnregularized() {
        assertBiasCovariatesEnabled();
        /* fetch data from cache */
        final INDArray M_Psi_inv_st = getINDArrayFromCache(CoverageModelICGCacheNode.M_Psi_inv_st);
        final INDArray Delta_st = getINDArrayFromCache(CoverageModelICGCacheNode.Delta_st);
        final INDArray W_tl = getINDArrayFromCache(CoverageModelICGCacheNode.W_tl);

        /* calculate the required quantities */
        final INDArray contribGMatrix = Nd4j.create(new int[] {numSamples, numLatents, numLatents});
        final INDArray contribZ = Nd4j.create(new int[] {numLatents, numSamples});
        for (int si = 0; si < numSamples; si++) {
            final INDArray M_psi_int_row_t = M_Psi_inv_st.getRow(si);
            final INDArray M_psi_int_row_t_trans = M_psi_int_row_t.transpose();
            final INDArray Delta_row_t = Delta_st.getRow(si);
            final INDArray contribZ_col_l = contribZ.getColumn(si);
            contribZ_col_l.assign(W_tl.transpose().mmul(M_psi_int_row_t.mul(Delta_row_t).transpose()));

            /* mean_W_part_{mn} = \sum_t E[W_{tm}] E[W_{tn}] M_{st} Psi^{-1}_{st} */
            final INDArray mean_W_part_ll = W_tl.transpose().mmul(W_tl.mulColumnVector(M_psi_int_row_t_trans));
            contribGMatrix.get(NDArrayIndex.point(si), NDArrayIndex.all(), NDArrayIndex.all())
                    .assign(mean_W_part_ll);
        }

        return new ImmutablePair<>(contribGMatrix, contribZ);
    }

    /**
     * Calculates the contribution of the target-space block represented by this compute node
     * to the E-step for $z_{s\mu}$ w/ regularization. The result is a triple
     * of INDArrays (contribGMatrix, contribZ, contribFilter). The first two were defined
     * in {@link CoverageModelEMComputeBlock#getBiasLatentPosteriorDataUnregularized()}. The third
     * INDArray is:
     *
     *      contribFilter = [W]^T [F.W]
     *
     * @return an {@link ImmutableTriple} of contribGMatrix (left), contribZ (middle), contribFilter (right)
     */
    public ImmutableTriple<INDArray, INDArray, INDArray> getBiasLatentPosteriorDataRegularized() {
        assertBiasCovariatesEnabled();
        final ImmutablePair<INDArray, INDArray> unreg = getBiasLatentPosteriorDataUnregularized();
        final INDArray contribFilter = getINDArrayFromCache(CoverageModelICGCacheNode.W_tl)
                .transpose().mmul(getINDArrayFromCache(CoverageModelICGCacheNode.F_W_tl));
        return new ImmutableTriple<>(unreg.left, unreg.right, contribFilter);
    }

    /**
     * Calculates the contribution of the target-space block represented by this compute node
     * to the E-step for sample read depth estimation. The final calculation is
     * performed on the driver node by {@link CoverageModelEMWorkspace}. The output of this
     * method is the pair of INDArrays (A_s, B_s) where:
     *
     *     A_s = \sum_{t \in block} M_{st} [\log(n_{st}/c_{st}) - b_{st}]
     *     B_s = \sum_{t \in block} M_{st}
     *
     * @return an {@link ImmutablePair} of (A_s, B_s)
     */
    public ImmutablePair<INDArray, INDArray> getReadDepthLatentPosteriorData(final boolean neglectBiasCovariates) {
        /* fetch data from cache */
        final INDArray log_n_st = getINDArrayFromCache(CoverageModelICGCacheNode.log_n_st);
        final INDArray log_c_st = getINDArrayFromCache(CoverageModelICGCacheNode.log_c_st);
        final INDArray m_t = getINDArrayFromCache(CoverageModelICGCacheNode.m_t);
        final INDArray M_Psi_inv_st = getINDArrayFromCache(CoverageModelICGCacheNode.M_Psi_inv_st);

        /* calculate the required quantities */
        final INDArray denominator = M_Psi_inv_st.sum(1);
        final INDArray numerator;
        if (biasCovariatesEnabled && !neglectBiasCovariates) {
            final INDArray Wz_st = getINDArrayFromCache(CoverageModelICGCacheNode.Wz_st);
            numerator = log_n_st.sub(log_c_st).subi(Wz_st).subiRowVector(m_t).muli(M_Psi_inv_st).sum(1);
        } else {
            numerator = log_n_st.sub(log_c_st).subiRowVector(m_t).muli(M_Psi_inv_st).sum(1);
        }
        return ImmutablePair.of(numerator, denominator);
    }

    /**
     * Calculates the contribution of the target-space block represented by this compute node
     * to the function subject to root finding in order to estimate the sample-specific log bias
     * unexplained variance (done on the driver done via {@link CoverageModelEMWorkspace}).
     *
     * This methods queries the objective function for multiple samples simultaneously.
     *
     * @param sampleIndices indices of samples
     * @param gamma_s the trial unexplained variance for the samples in question
     * @return the objective function as an {@link INDArray} with shape [sampleIndices.length, 1]
     */
    public INDArray calculateGammaObjectiveFunctionMultiSample(final int[] sampleIndices,
                                                               final INDArray gamma_s) {
        final int numQueries = sampleIndices.length;
        Utils.validateArg(numQueries > 0, "The objective function for at least one sample must be queried");
        Utils.validateArg(gamma_s.length() == numQueries, "The INDArray of trial values for gamma must have" +
                " the same length as the corresponding list of samples");

        final INDArray M_st = getINDArrayFromCache(CoverageModelICGCacheNode.M_st);
        final INDArray Psi_t = getINDArrayFromCache(CoverageModelICGCacheNode.Psi_t);
        final INDArray Sigma_st = getINDArrayFromCache(CoverageModelICGCacheNode.Sigma_st);
        final INDArray B_st = getINDArrayFromCache(CoverageModelICGCacheNode.B_st);

        final INDArray assembled_M_st = Nd4j.create(numQueries, numTargets);
        final INDArray assembled_Sigma_st = Nd4j.create(numQueries, numTargets);
        final INDArray assembled_B_st = Nd4j.create(numQueries, numTargets);
        for (int i = 0; i < numQueries; i++) {
            assembled_M_st.getRow(i).assign(M_st.getRow(sampleIndices[i]));
            assembled_Sigma_st.getRow(i).assign(Sigma_st.getRow(sampleIndices[i]));
            assembled_B_st.getRow(i).assign(B_st.getRow(sampleIndices[i]));
        }
        final INDArray totalMaskedPsiInverse = assembled_M_st.div(
                assembled_Sigma_st.addRowVector(Psi_t).addiColumnVector(gamma_s));
        final INDArray totalMaskedPsiInverseMulB = assembled_B_st.mul(totalMaskedPsiInverse);
        return totalMaskedPsiInverse.mul(assembled_M_st.sub(totalMaskedPsiInverseMulB)).sum(1);
    }

    /**
     * Calculates the contribution of the target-space block represented by this compute node
     * to the E-step for copy ratio (or copy number) (done on the driver done via
     * {@link CoverageModelEMWorkspace}).
     *
     * @return a double-list of emission data (first for samples, second for targets)
     */
    public List<List<CoverageModelCopyRatioEmissionData>> getSampleCopyRatioLatentPosteriorData() {
        /* fetch data from cache */
        final INDArray n_st = getINDArrayFromCache(CoverageModelICGCacheNode.n_st);
        final INDArray Psi_t = getINDArrayFromCache(CoverageModelICGCacheNode.Psi_t);
        final INDArray gamma_s = getINDArrayFromCache(CoverageModelICGCacheNode.gamma_s);
        final INDArray m_t = getINDArrayFromCache(CoverageModelICGCacheNode.m_t);
        final INDArray err_st = getINDArrayFromCache(CoverageModelICGCacheNode.err_st);

        /* calculate the required quantities */
        final INDArray mu_st;
        if (biasCovariatesEnabled) {
            final INDArray Wz_st = getINDArrayFromCache(CoverageModelICGCacheNode.Wz_st);
             mu_st = Wz_st.addRowVector(m_t);
        } else {
             mu_st = Nd4j.vstack(Collections.nCopies(numSamples, m_t));
        }

        final double[] psiArray = Psi_t.dup().data().asDouble();
        final double[] gammaArray = gamma_s.dup().data().asDouble();

        return IntStream.range(0, numSamples)
                .mapToObj(si -> {
                    final double[] currentSampleReadCountArray = n_st.getRow(si).dup().data().asDouble();
                    final double[] currentSampleMuArray = mu_st.getRow(si).dup().data().asDouble();
                    final double[] currentSampleMappingErrorArray = err_st.getRow(si).dup().data().asDouble();
                    return IntStream.range(0, targetBlock.getNumTargets())
                            .mapToObj(ti -> new CoverageModelCopyRatioEmissionData(
                                    currentSampleMuArray[ti], /* log bias */
                                    psiArray[ti] + gammaArray[si], /* total unexplained variance */
                                    (int) currentSampleReadCountArray[ti], /* raw read count */
                                    currentSampleMappingErrorArray[ti]) /* mapping error probability */
                            ).collect(Collectors.toList());
                }).collect(Collectors.toList());
    }

    /**
     * Returns the partially target-summed posterior second moment of W:
     *
     *     \sum_{t in target-space block} E[W_{tl} W_{tl}]
     *
     * We currently treat W in the max likelihood sense and neglect its covariance. Therefore,
     * E[W_{tl} W_{tl}] = E[W_{tl}] E[W_{tl}].
     *
     * @return an {@link INDArray}
     */
    public INDArray getBiasCovariatesSecondMomentPosteriorsPartialTargetSum() {
        assertBiasCovariatesEnabled();
        assertARDEnabled();
        final INDArray W_tl = getINDArrayFromCache(CoverageModelICGCacheNode.W_tl);
        return Transforms.pow(W_tl, 2, true).sum(0);
    }

    private double calculatePsiSolverObjectiveFunction(final int targetIndex,
                                              final double psi,
                                              @Nonnull final INDArray M_st,
                                              @Nonnull final INDArray Sigma_st,
                                              @Nonnull final INDArray gamma_s,
                                              @Nonnull final INDArray B_st) {
        final INDArray M_s = M_st.get(NDArrayIndex.all(), NDArrayIndex.point(targetIndex));
        final INDArray Sigma_s = Sigma_st.get(NDArrayIndex.all(), NDArrayIndex.point(targetIndex));
        final INDArray B_s = B_st.get(NDArrayIndex.all(), NDArrayIndex.point(targetIndex));

        final INDArray totalMaskedPsiInverse = M_s.div(Sigma_s.add(gamma_s).addi(psi));
        final INDArray totalMaskedPsiInverseSqrMulB_s = B_s.mul(totalMaskedPsiInverse).muli(totalMaskedPsiInverse);

        return totalMaskedPsiInverse.subi(totalMaskedPsiInverseSqrMulB_s).sumNumber().doubleValue();
    }

    private double calculatePsiSolverMeritFunction(final int targetIndex,
                                                   final double psi,
                                                   @Nonnull final INDArray M_st,
                                                   @Nonnull final INDArray Sigma_st,
                                                   @Nonnull final INDArray gamma_s,
                                                   @Nonnull final INDArray B_st) {
        final INDArray M_s = M_st.get(NDArrayIndex.all(), NDArrayIndex.point(targetIndex));
        final INDArray Sigma_s = Sigma_st.get(NDArrayIndex.all(), NDArrayIndex.point(targetIndex));
        final INDArray B_s = B_st.get(NDArrayIndex.all(), NDArrayIndex.point(targetIndex));
        final INDArray totalPsi = Sigma_s.add(gamma_s).addi(psi);
        final INDArray totalMaskedPsiInverse = M_s.div(totalPsi);
        final INDArray maskedLogTotalPsi = Transforms.log(totalPsi, false).muli(M_s); /* this mutates totalPsi */

        return maskedLogTotalPsi.addi(B_s.mul(totalMaskedPsiInverse)).muli(-0.5).sumNumber().doubleValue();
    }

    /**
     * Calculates the contribution of this target-space block to the object function subject
     * to root finding in the M-step update of isotropic unexplained variance.
     *
     * @param psi trial value of isotropic unexplained variance
     * @return a double
     */
    public double calculateSampleTargetSummedPsiObjectiveFunction(final double psi) {
        final INDArray M_st = getINDArrayFromCache(CoverageModelICGCacheNode.M_st);
        final INDArray Sigma_st = getINDArrayFromCache(CoverageModelICGCacheNode.Sigma_st);
        final INDArray B_st = getINDArrayFromCache(CoverageModelICGCacheNode.B_st);
        final INDArray gamma_s = getINDArrayFromCache(CoverageModelICGCacheNode.gamma_s);
        final INDArray totalMaskedPsiInverse = M_st.div(Sigma_st.addColumnVector(gamma_s).addi(psi));
        final INDArray totalMaskedPsiInverseMulB = B_st.mul(totalMaskedPsiInverse);
        return totalMaskedPsiInverse.mul(M_st.sub(totalMaskedPsiInverseMulB)).sumNumber().doubleValue();
    }

    /**
     * Calculates the contribution of this target-space block to the merit function for choosing
     * the best solution in the M-step update of isotropic unexplained variance.
     *
     * @param psi trial value of isotropic unexplained variance
     * @return a double
     */
    public double calculateSampleTargetSummedPsiMeritFunction(final double psi) {
        final INDArray M_st = getINDArrayFromCache(CoverageModelICGCacheNode.M_st);
        final INDArray Sigma_st = getINDArrayFromCache(CoverageModelICGCacheNode.Sigma_st);
        final INDArray B_st = getINDArrayFromCache(CoverageModelICGCacheNode.B_st);
        final INDArray gamma_s = getINDArrayFromCache(CoverageModelICGCacheNode.gamma_s);

        final INDArray totalPsi = Sigma_st.addColumnVector(gamma_s).addi(psi);
        final INDArray totalMaskedPsiInverse = M_st.div(totalPsi);
        final INDArray maskedLogTotalPsi = Transforms.log(totalPsi, false).muli(M_st); /* this mutates totalPsi */

        return maskedLogTotalPsi.addi(B_st.mul(totalMaskedPsiInverse)).muli(-0.5).sumNumber().doubleValue();
    }

    /**
     * Calculates the contribution of this target-space block to log model evidence for given ARD coefficients.
     *
     * @param alpha_l ARD coefficients (as a 1 x D row vector)
     * @return
     */
    public double calculateBiasCovariatesLogEvidencePartialTargetSum(final INDArray alpha_l) {
        assertARDEnabled();
        /* fetch the required caches */
        final INDArray Q_tll = getINDArrayFromCache(CoverageModelICGCacheNode.Q_tll);
        final INDArray v_tl = getINDArrayFromCache(CoverageModelICGCacheNode.v_tl);

        final int numTargets = v_tl.shape()[0];
        final INDArray alpha_ll = Nd4j.diag(alpha_l);

        /* calculate log|A| */
        double logDetATerm = numTargets * CoverageModelEMWorkspaceMathUtils.logdet(alpha_ll);

        /* calculate log|A + Q| */
        double logDetXTerm = 0;
        for (int ti = 0; ti < numTargets; ti++) {
            final INDArray X_ll = Q_tll.get(NDArrayIndex.point(ti), NDArrayIndex.all(), NDArrayIndex.all())
                    .add(alpha_ll);
            logDetXTerm += CoverageModelEMWorkspaceMathUtils.logdet(X_ll);
        }

        /* calculate v^T (Q + A)^{-1} v */
        double v_X_v = 0;
        final double LINSOLVE_SINGULARITY_THRESHOLD = Double.MIN_VALUE;
        for (int ti = 0; ti < numTargets; ti++) {
            final INDArray v_l = v_tl.get(NDArrayIndex.point(ti), NDArrayIndex.all());
            final INDArray X_inv_v_l = CoverageModelEMWorkspaceMathUtils.linsolve(
                    Q_tll.get(NDArrayIndex.point(ti), NDArrayIndex.all(), NDArrayIndex.all())
                            .add(alpha_ll),
                    v_tl.get(NDArrayIndex.point(ti), NDArrayIndex.all()),
                    LINSOLVE_SINGULARITY_THRESHOLD);
            v_X_v += v_l.mul(X_inv_v_l.transpose()).sumNumber().doubleValue();
        }

        return 0.5 * (logDetATerm - logDetXTerm + v_X_v);
    }

    /**
     * Calculates the contribution of this target-space block to the target covariance matrix.
     * This method is used in {@link CoverageModelEMWorkspace#initializeWorkersWithPCA}.
     *
     * @return an {@link INDArray}
     */
    public INDArray calculateTargetCovarianceMatrixForPCAInitialization() {
        assertBiasCovariatesEnabled();
        final INDArray Delta_PCA_st = getINDArrayFromCache(CoverageModelICGCacheNode.Delta_PCA_st);
        return Delta_PCA_st.mmul(Delta_PCA_st.transpose());
    }

    /**
     * Takes an arbitrary INDArray {@code arr} and a mask array {@code mask} with the same shape
     * and replace all entries in {@code arr} on which {@code mask} is 0 to {@code value}
     *
     * TODO github/gatk-protected issue #853 -- this could be optimized perhaps using the existing
     *      nd4j native Ops. I haven't found a way yet.
     *
     * @param arr an arbitrary INDArray
     * @param mask a mask array
     * @param value replacement value on masked entries
     * @return a copy of {@code arr} with replaced values
     */
    @VisibleForTesting
    public static INDArray replaceMaskedEntries(final INDArray arr, final INDArray mask, final double value) {
        Utils.validateArg(Arrays.equals(arr.shape(), mask.shape()), "The value and mask arrays have different shapes");
        Utils.validateArg(mask.rank() <= 2, "The mask array must be at most rank 2");
        final int rowDim = mask.shape()[0];
        final int colDim = mask.shape()[1];
        final INDArray res = arr.dup();
        for (int i = 0; i < rowDim; i++) {
            for (int j = 0; j < colDim; j++) {
                if (mask.getInt(i, j) == 0) {
                    res.putScalar(i, j, value);
                }
            }
        }
        return res;
    }

    /**
     * Clone the compute block with initialized data
     *
     * @param initialDataBlock initial data block
     * @return a new instance of {@link CoverageModelEMComputeBlock}
     */
    public CoverageModelEMComputeBlock cloneWithInitializedData(@Nonnull final InitialDataBlock initialDataBlock) {
        Utils.nonNull(initialDataBlock, "The initial data block must be non-null");
        initialDataBlock.assertDataBlockSize(targetBlock.getNumTargets(), numSamples);

        /* create NDArrays from the data blocks */
        final INDArray readCounts = Nd4j.create(
                Arrays.stream(initialDataBlock.readCountBlock).mapToDouble(n -> (double)n).toArray(),
                new int[] {numSamples, numTargets}, 'f');
        final INDArray mask = Nd4j.create(
                Arrays.stream(initialDataBlock.maskBlock).mapToDouble(n -> (double)n).toArray(),
                new int[] {numSamples, numTargets}, 'f');
        final INDArray mappingErrorRate = Nd4j.create(initialDataBlock.mappingErrorRateBlock,
                new int[] {numSamples, numTargets}, 'f');

        return cloneWithUpdatedPrimitive(CoverageModelICGCacheNode.n_st, readCounts)
                .cloneWithUpdatedPrimitive(CoverageModelICGCacheNode.M_st, mask)
                .cloneWithUpdatedPrimitive(CoverageModelICGCacheNode.err_st, mappingErrorRate);
    }

    /**
     * Clone the compute block with updated copy ratio posteriors
     *
     * @param log_c_st posterior mean of log copy ratios
     * @param var_log_c_st posterior variance of log copy ratios
     * @param admixingRatio the admixing ratio of old and new posterior expectations
     * @return the new instance of compute block
     */
    public CoverageModelEMComputeBlock cloneWithUpdatedCopyRatioPosteriors(@Nonnull final INDArray log_c_st,
                                                                           @Nonnull final INDArray var_log_c_st,
                                                                           final double admixingRatio) {
        Utils.nonNull(log_c_st, "Log copy ratio posterior means must be non-null");
        Utils.nonNull(var_log_c_st, "Log copy ratio posterior variances must be non-null");
        Utils.validateArg(admixingRatio >= 0 && admixingRatio <= 1, "The admixing ratio must be between 0 and 1");

        /* fetch required quantities from cache */
        final INDArray M_st = getINDArrayFromCache(CoverageModelICGCacheNode.M_st);
        final INDArray old_log_c_st = getINDArrayFromCache(CoverageModelICGCacheNode.log_c_st);
        final INDArray old_var_log_c_st = getINDArrayFromCache(CoverageModelICGCacheNode.var_log_c_st);

        /* replace values on masked targets with default values */
        final INDArray rectified_log_c_st = replaceMaskedEntries(log_c_st, M_st,
                CoverageModelGlobalConstants.MEAN_LOG_COPY_RATIO_ON_MASKED_TARGETS);
        final INDArray rectified_var_log_c_st = replaceMaskedEntries(var_log_c_st, M_st,
                CoverageModelGlobalConstants.VAR_LOG_COPY_RATIO_ON_MASKED_TARGETS);

        /* admix */
        final INDArray admixed_log_c_st = rectified_log_c_st
                .muli(admixingRatio)
                .addi(old_log_c_st.mul(1.0 - admixingRatio));
        final INDArray admixed_var_log_c_st = rectified_var_log_c_st
                .muli(admixingRatio)
                .addi(old_var_log_c_st.mul(1.0 - admixingRatio));

        final double errNormInfinity = CoverageModelEMWorkspaceMathUtils.getINDArrayNormInfinity(
                old_log_c_st.sub(admixed_log_c_st));
        return cloneWithUpdatedPrimitive(CoverageModelICGCacheNode.log_c_st, admixed_log_c_st)
                .cloneWithUpdatedPrimitive(CoverageModelICGCacheNode.var_log_c_st, admixed_var_log_c_st)
                .cloneWithUpdatedSignal(SubroutineSignal.builder().put("error_norm", errNormInfinity).build());
    }

    /**
     * Clones this compute block with updated prior mean and variance of log copy ratios.
     *
     * @param log_c_st prior mean of log copy ratio
     * @param var_log_c_st prior variance of log copy ratio
     * @return an instance of {@link CoverageModelEMComputeBlock}
     */
    public CoverageModelEMComputeBlock cloneWithUpdateCopyRatioPriors(@Nonnull final INDArray log_c_st,
                                                                      @Nonnull final INDArray var_log_c_st) {
        Utils.nonNull(log_c_st, "Log copy ratio posterior means must be non-null");
        Utils.nonNull(var_log_c_st, "Log copy ratio posterior variances must be non-null");

        /* fetch required quantities from cache */
        final INDArray M_st = getINDArrayFromCache(CoverageModelICGCacheNode.M_st);

        /* replace values on masked targets with default values */
        final INDArray rectified_log_c_st = replaceMaskedEntries(log_c_st, M_st,
                CoverageModelGlobalConstants.MEAN_LOG_COPY_RATIO_ON_MASKED_TARGETS);
        final INDArray rectified_var_log_c_st = replaceMaskedEntries(var_log_c_st, M_st,
                CoverageModelGlobalConstants.VAR_LOG_COPY_RATIO_ON_MASKED_TARGETS);

        return cloneWithUpdatedPrimitive(CoverageModelICGCacheNode.log_c_st, rectified_log_c_st)
                .cloneWithUpdatedPrimitive(CoverageModelICGCacheNode.var_log_c_st, rectified_var_log_c_st);
    }

    /**
     * Performs the M-step for target-specific unexplained variance and clones the compute block
     * with the updated value.
     *
     * @param maxIters maximum number of iterations
     * @param psiUpperLimit upper limit for the unexplained variance
     * @param absTol absolute error tolerance (used in root finding)
     * @param relTol relative error tolerance (used in root finding)
     * @param numBisections number of bisections (used in root finding)
     * @param refinementDepth depth of search (used in root finding)
     *
     * @return a new instance of {@link CoverageModelEMComputeBlock}
     */
    public CoverageModelEMComputeBlock cloneWithUpdatedTargetUnexplainedVarianceTargetResolved(final int maxIters,
                                                                                               final double psiUpperLimit,
                                                                                               final double absTol,
                                                                                               final double relTol,
                                                                                               final int numBisections,
                                                                                               final int refinementDepth,
                                                                                               final int numThreads) {
        Utils.validateArg(maxIters > 0, "At least one iteration is required");
        Utils.validateArg(psiUpperLimit >= 0, "The upper limit must be non-negative");
        Utils.validateArg(absTol >= 0, "The absolute error tolerance must be non-negative");
        Utils.validateArg(relTol >= 0, "The relative error tolerance must be non-negative");
        Utils.validateArg(numBisections >= 0, "The number of bisections must be non-negative");
        Utils.validateArg(refinementDepth >= 0, "The refinement depth must be non-negative");
        Utils.validateArg(numThreads > 0, "Number of execution threads must be positive");

        /* fetch the required caches */
        final INDArray Psi_t = getINDArrayFromCache(CoverageModelICGCacheNode.Psi_t);
        final INDArray M_st = getINDArrayFromCache(CoverageModelICGCacheNode.M_st);
        final INDArray Sigma_st = getINDArrayFromCache(CoverageModelICGCacheNode.Sigma_st);
        final INDArray gamma_s = getINDArrayFromCache(CoverageModelICGCacheNode.gamma_s);
        final INDArray B_st = getINDArrayFromCache(CoverageModelICGCacheNode.B_st);

        final ForkJoinPool forkJoinPool = new ForkJoinPool(numThreads);
        final List<ImmutablePair<Double, Integer>> res;
        try {
            res = forkJoinPool.submit(() -> {
                return IntStream.range(0, numTargets)
                        .parallel()
                        .mapToObj(ti -> {
                            final RobustBrentSolver solver = new RobustBrentSolver(relTol, absTol,
                                    CoverageModelGlobalConstants.DEFAULT_FUNCTION_EVALUATION_ACCURACY);
                            double newPsi;
                            try {
                                newPsi = solver.solve(maxIters,
                                        psi -> calculatePsiSolverObjectiveFunction(ti, psi, M_st, Sigma_st, gamma_s, B_st),
                                        psi -> calculatePsiSolverMeritFunction(ti, psi, M_st, Sigma_st, gamma_s, B_st),
                                        null, 0, psiUpperLimit, numBisections, refinementDepth);
                            } catch (NoBracketingException | TooManyEvaluationsException e) {
                            /* if a solution can not be found, set Psi to its old value */
                                newPsi = Psi_t.getDouble(ti);
                            }
                            return new ImmutablePair<>(newPsi, solver.getEvaluations());
                        })
                        .collect(Collectors.toList());
            }).get();
        } catch (InterruptedException | ExecutionException ex) {
            throw new RuntimeException("Failure in concurrent update of target-specific unexplained variance");
        }

        final INDArray newPsi_t = Nd4j.create(res.stream().mapToDouble(p -> p.left).toArray(), Psi_t.shape());
        final int maxIterations = Collections.max(res.stream().mapToInt(p -> p.right).boxed().collect(Collectors.toList()));
        final double errNormInfinity = CoverageModelEMWorkspaceMathUtils.getINDArrayNormInfinity(newPsi_t.sub(Psi_t));
        return cloneWithUpdatedPrimitiveAndSignal(CoverageModelICGCacheNode.Psi_t, newPsi_t, SubroutineSignal.builder()
                .put("error_norm", errNormInfinity).put("iterations", maxIterations).build());
    }

    /**
     * Performs the M-step for updating the bias covariates and clones the compute with updated values
     *
     * @return a new instance of {@link CoverageModelEMComputeBlock}
     */
    public CoverageModelEMComputeBlock cloneWithUpdatedBiasCovariatesUnregularized(final double admixingRatio) {
        assertBiasCovariatesEnabled();
        /* fetch the required caches */
        final INDArray Q_tll = getINDArrayFromCache(CoverageModelICGCacheNode.Q_tll);
        final INDArray v_tl = getINDArrayFromCache(CoverageModelICGCacheNode.v_tl);
        final INDArray W_tl = getINDArrayFromCache(CoverageModelICGCacheNode.W_tl);

        final int numTargets = v_tl.shape()[0];
        final int numLatents = v_tl.shape()[1];
        final INDArray new_W_tl = Nd4j.create(numTargets, numLatents);

        final INDArray alpha_ll;
        if (ardEnabled) {
            final INDArray alpha_l = getINDArrayFromCache(CoverageModelICGCacheNode.alpha_l);
            alpha_ll = Nd4j.diag(alpha_l);
        } else {
            alpha_ll = Nd4j.zeros(new int[] {numLatents, numLatents});
        }

        /* calculate E[W_{tm}] */
        final double LINSOLVE_SINGULARITY_THRESHOLD = Double.MIN_VALUE;
        for (int ti = 0; ti < numTargets; ti++) {
            new_W_tl.get(NDArrayIndex.point(ti), NDArrayIndex.all()).assign(
                    CoverageModelEMWorkspaceMathUtils.linsolve(
                            Q_tll.get(NDArrayIndex.point(ti), NDArrayIndex.all(), NDArrayIndex.all())
                                    .add(alpha_ll),
                            v_tl.get(NDArrayIndex.point(ti), NDArrayIndex.all()), LINSOLVE_SINGULARITY_THRESHOLD));
        }

        /* admix with old posteriors */
        final INDArray new_W_tl_admixed = new_W_tl.mul(admixingRatio).addi(W_tl.mul(1.0 - admixingRatio));

        /* calculate error norm only based on the change in E[W_{tl}] */
        final double errNormInfinity = CoverageModelEMWorkspaceMathUtils
                .getINDArrayNormInfinity(new_W_tl.sub(getINDArrayFromCache(CoverageModelICGCacheNode.W_tl)));
        final SubroutineSignal sig = SubroutineSignal.builder().put("error_norm", errNormInfinity).build();

        return cloneWithUpdatedPrimitiveAndSignal(CoverageModelICGCacheNode.W_tl, new_W_tl_admixed, sig);
    }

    /**
     * Performs the M-step for updating the mean log bias and clones the compute block with the updated value
     *
     * @param neglectBiasCovariates if true, the contribution of bias covariates will be neglected
     * @return a new instance of {@link CoverageModelEMComputeBlock}
     */
    public CoverageModelEMComputeBlock cloneWithUpdatedMeanLogBias(final boolean neglectBiasCovariates) {
        /* fetch the required caches */
        final INDArray log_n_st = getINDArrayFromCache(CoverageModelICGCacheNode.log_n_st);
        final INDArray log_c_st = getINDArrayFromCache(CoverageModelICGCacheNode.log_c_st);
        final INDArray log_d_s = getINDArrayFromCache(CoverageModelICGCacheNode.log_d_s);
        final INDArray M_Psi_inv_st = getINDArrayFromCache(CoverageModelICGCacheNode.M_Psi_inv_st);

        final INDArray numerator;
        if (neglectBiasCovariates || !biasCovariatesEnabled) {
            numerator = M_Psi_inv_st.mul(log_n_st.sub(log_c_st).subiColumnVector(log_d_s)).sum(0);
        } else {
            final INDArray Wz_st = getINDArrayFromCache(CoverageModelICGCacheNode.Wz_st);
            numerator = M_Psi_inv_st.mul(log_n_st.sub(log_c_st).subiColumnVector(log_d_s).subi(Wz_st)).sum(0);
        }
        final INDArray denominator = M_Psi_inv_st.sum(0);
        final INDArray newTargetMeanBias = numerator.divi(denominator);

        double errNormInfinity = CoverageModelEMWorkspaceMathUtils.getINDArrayNormInfinity(
                getINDArrayFromCache(CoverageModelICGCacheNode.m_t).sub(newTargetMeanBias));

        return cloneWithUpdatedPrimitiveAndSignal(CoverageModelICGCacheNode.m_t, newTargetMeanBias,
                SubroutineSignal.builder().put("error_norm", errNormInfinity).build());
    }

    /**
     * Clones this compute block with updated cache node {@link CoverageModelICGCacheNode#Delta_PCA_st}.
     * Used for PCA initialization of the model parameters in {@link CoverageModelEMWorkspace#initializeWorkersWithPCA()}
     *
     * @implNote The design choice for defining {@link CoverageModelICGCacheNode#Delta_PCA_st} as externally mutable
     *           computable node is to be able to nullify it later on to save memory.
     *           See {{@link #cloneWithRemovedPCAInitializationData()}}.
     * @return an updated instance of {@link CoverageModelEMComputeBlock}
     */
    public CoverageModelEMComputeBlock cloneWithPCAInitializationData(final int minIncludedReadCount,
                                                                      final int maxIncludedReadCount) {
        /* fetch the required INDArrays */
        final INDArray M_st = getINDArrayFromCache(CoverageModelICGCacheNode.M_st);
        final INDArray n_st = getINDArrayFromCache(CoverageModelICGCacheNode.n_st);
        final INDArray log_n_st = getINDArrayFromCache(CoverageModelICGCacheNode.log_n_st);
        final INDArray log_c_st = getINDArrayFromCache(CoverageModelICGCacheNode.log_c_st);
        final INDArray log_d_s = getINDArrayFromCache(CoverageModelICGCacheNode.log_d_s);

        /* create the second mask based on the provided min and max read counts */
        final INDArray M_extra_st = M_st.dup();
        for (int si = 0; si < numSamples; si++) {
            for (int ti = 0; ti < numTargets; ti++) {
                if (M_st.getInt(si, ti) == 0) { /* if previously masked */
                    continue;
                }
                final int nextReadCount = n_st.getInt(si, ti);
                if (nextReadCount >= minIncludedReadCount && nextReadCount <= maxIncludedReadCount) {
                    M_extra_st.put(si, ti, 1);
                } else {
                    M_extra_st.put(si, ti, 0);
                }
            }
        }

        final INDArray M_Delta_st = log_n_st.sub(log_c_st).subiColumnVector(log_d_s).muli(M_extra_st);
        final int numSamples = M_Delta_st.shape()[0];
        final INDArray M_Delta_mean_t = M_Delta_st.sum(0).divi(numSamples);
        final INDArray log_bias_var_st = M_Delta_st.subiRowVector(M_Delta_mean_t);

        return cloneWithUpdatedPrimitive(CoverageModelICGCacheNode.Delta_PCA_st, log_bias_var_st);
    }

    /**
     * Clone this compute block with nullified {@link CoverageModelICGCacheNode#Delta_PCA_st}. This method is
     * invoked after the PCA initialization is complete.
     *
     * @return an updated instance of {@link CoverageModelEMComputeBlock}
     */
    public CoverageModelEMComputeBlock cloneWithRemovedPCAInitializationData() {
        return cloneWithUpdatedPrimitive(CoverageModelICGCacheNode.Delta_PCA_st, null);
    }


    /**
     * Creates a new instance of this compute block with an updated primitive node
     *
     * @param key the key for the node to update
     * @param value the new value
     * @return a new instance of {@link CoverageModelEMComputeBlock}
     */
    public CoverageModelEMComputeBlock cloneWithUpdatedPrimitive(@Nonnull final CoverageModelICGCacheNode key,
                                                                 @Nullable final INDArray value) {
        if (value == null) {
            return new CoverageModelEMComputeBlock(targetBlock, numSamples, numLatents, ardEnabled,
                    icg.nullifyNode(key.name()), latestMStepSignal);
        } else {
            return new CoverageModelEMComputeBlock(targetBlock, numSamples, numLatents, ardEnabled,
                    icg.setValue(key.name(), new DuplicableNDArray(value)), latestMStepSignal);
        }
    }

    /**
     * Creates a new instance of this compute block with an updated primitive node and a subroutine
     * signal calculated externally (i.e. on the driver node)
     *
     * @param key the key for the node to update
     * @param value the new value
     * @param latestMStepSignal an externally calculated M-step signal
     * @return a new instance of {@link CoverageModelEMComputeBlock}
     */
    public CoverageModelEMComputeBlock cloneWithUpdatedPrimitiveAndSignal(@Nonnull final CoverageModelICGCacheNode key,
                                                                          @Nullable final INDArray value,
                                                                          @Nonnull final SubroutineSignal latestMStepSignal) {
        if (value == null) {
            return new CoverageModelEMComputeBlock(targetBlock, numSamples, numLatents, ardEnabled,
                    icg.setValue(key.name(), new DuplicableNDArray()), latestMStepSignal);
        } else {
            return new CoverageModelEMComputeBlock(targetBlock, numSamples, numLatents, ardEnabled,
                    icg.setValue(key.name(), new DuplicableNDArray(value)), latestMStepSignal);
        }
    }

    /**
     * Creates a new instance of this compute block with an updated subroutine signal calculated externally
     * (i.e. on the driver node)
     *
     * @param latestMStepSignal an externally calculated M-step signal
     * @return a new instance of {@link CoverageModelEMComputeBlock}
     */
    public CoverageModelEMComputeBlock cloneWithUpdatedSignal(@Nonnull final SubroutineSignal latestMStepSignal) {
        return new CoverageModelEMComputeBlock(targetBlock, numSamples, numLatents, ardEnabled, icg, latestMStepSignal);
    }

    /**
     * Creates a new instance of this compute block with updated caches associated to a cache tag
     *
     * @param tag the cache tag
     * @return a new instance of {@link CoverageModelEMComputeBlock}
     */
    public CoverageModelEMComputeBlock cloneWithUpdatedCachesByTag(final CoverageModelICGCacheTag tag) {
        return new CoverageModelEMComputeBlock(targetBlock, numSamples, numLatents, ardEnabled,
                icg.updateCachesForTag(tag.name()), latestMStepSignal);
    }

    /**
     * Invoke garbage collection.
     *
     * @implNote The driver node must call this method periodically to perform garbage collection
     * on the worker nodes. Since INDArray data buffers are stored off-heap (with only pointers
     * living on-heap), JVM fails to perform garbage collection in a timely fashion. In theory,
     * JavaCPP, the off-head memory allocation native backend of Nd4j, in meant to properly
     * keep track of memory allocation, call gc() when necessary, and free up the off-heap memory
     * via deallocation finalizers. However, its gc() calling strategy is very unreliable, i.e.
     * it calls gc() when the off-heap memory is nearly full (and not concurrently or scheduled).
     *
     * -- Mehrtash Babadi (Jan 2, 2017)
     */
    public void performGarbageCollection() {
        System.gc();
    }

    /**
     * Retrieves the value of a cache node from a map and casts it to an {@link INDArray}
     *
     * Note: up-casting and non-null conditions are not checked
     *
     * @param key the key for the cache node
     * @param map a map from keys to values
     * @return an {@link INDArray} instance
     */
    @SuppressWarnings("unchecked")
    private static INDArray getINDArrayFromMap(@Nonnull final CoverageModelICGCacheNode key,
                                               @Nonnull Map<String, ? extends Duplicable> map) {
        return ((DuplicableNDArray)map.get(key.name())).value();
    }

    /**
     * Retrieves the value of a cache node from a map and casts it to a primitive double
     *
     * Note: up-casting and non-null conditions are not checked
     *
     * @param key the key for the cache node
     * @param map a map from keys to values
     * @return a double value
     */
    @SuppressWarnings("unchecked")
    private static double getDoubleFromMap(@Nonnull final CoverageModelICGCacheNode key,
                                           @Nonnull Map<String, ? extends Duplicable> map) {
        return ((DuplicableNumber<Double>)map.get(key.name())).value();
    }

    /* graphical computation functions */

    /* dependents: [M_st] */
    private static final Function<Map<String, ? extends Duplicable>, ? extends Duplicable> calculate_sum_M_t =
            new Function<Map<String, ? extends Duplicable>, Duplicable>() {
                @Override
                public Duplicable apply(Map<String, ? extends Duplicable> dat) {
                    return new DuplicableNDArray(getINDArrayFromMap(CoverageModelICGCacheNode.M_st, dat).sum(0));
                }
            };

    /* dependents: [M_st] */
    private static final Function<Map<String, ? extends Duplicable>, ? extends Duplicable> calculate_sum_M_s =
            new Function<Map<String, ? extends Duplicable>, Duplicable>() {
                @Override
                public Duplicable apply(Map<String, ? extends Duplicable> dat) {
                    return new DuplicableNDArray(getINDArrayFromMap(CoverageModelICGCacheNode.M_st, dat).sum(1));
                }
            };

    /* dependents: [n_st, M_st] */
    private static final Function<Map<String, ? extends Duplicable>, ? extends Duplicable> calculate_Sigma_st =
            new Function<Map<String, ? extends Duplicable>, Duplicable>() {
                @Override
                public Duplicable apply(Map<String, ? extends Duplicable> dat) {
                    final INDArray n_st = getINDArrayFromMap(CoverageModelICGCacheNode.n_st, dat);
                    final INDArray M_st = getINDArrayFromMap(CoverageModelICGCacheNode.M_st, dat);
                    final INDArray Sigma_st = replaceMaskedEntries(
                            Nd4j.ones(n_st.shape()).divi(n_st),
                            M_st,
                            CoverageModelGlobalConstants.POISSON_STATISTICAL_VARIANCE_ON_MASKED_TARGETS);
                    return new DuplicableNDArray(Sigma_st);
                }
            };

    /* dependents: [n_st, M_st] */
    private static final Function<Map<String, ? extends Duplicable>, ? extends Duplicable> calculate_log_n_st =
            new Function<Map<String, ? extends Duplicable>, Duplicable>() {
                @Override
                public Duplicable apply(Map<String, ? extends Duplicable> dat) {
                    final INDArray n_st = getINDArrayFromMap(CoverageModelICGCacheNode.n_st, dat);
                    final INDArray M_st = getINDArrayFromMap(CoverageModelICGCacheNode.M_st, dat);
                    final INDArray log_n_st = replaceMaskedEntries(
                            Transforms.log(n_st, true),
                            M_st,
                            CoverageModelGlobalConstants.LOG_READ_COUNT_ON_MASKED_TARGETS);
                    return new DuplicableNDArray(log_n_st);
                }
            };

    private static final Function<Map<String, ? extends Duplicable>, ? extends Duplicable> calculate_Delta_st =
            new Function<Map<String, ? extends Duplicable>, Duplicable>() {
                @Override
                public Duplicable apply(Map<String, ? extends Duplicable> dat) {
                    final INDArray log_n_st = getINDArrayFromMap(CoverageModelICGCacheNode.log_n_st, dat);
                    final INDArray log_c_st = getINDArrayFromMap(CoverageModelICGCacheNode.log_c_st, dat);
                    final INDArray log_d_s = getINDArrayFromMap(CoverageModelICGCacheNode.log_d_s, dat);
                    final INDArray m_t = getINDArrayFromMap(CoverageModelICGCacheNode.m_t, dat);
                    return new DuplicableNDArray(log_n_st.sub(log_c_st).subiColumnVector(log_d_s).subiRowVector(m_t));
                }
            };

    /* dependents: [W_tl, z_sl] */
    private static final Function<Map<String, ? extends Duplicable>, ? extends Duplicable> calculate_Wz_st =
            new Function<Map<String, ? extends Duplicable>, Duplicable>() {
                @Override
                public Duplicable apply(Map<String, ? extends Duplicable> dat) {
                    final INDArray W_tl = getINDArrayFromMap(CoverageModelICGCacheNode.W_tl, dat);
                    final INDArray z_sl = getINDArrayFromMap(CoverageModelICGCacheNode.z_sl, dat);
                    return new DuplicableNDArray(W_tl.mmul(z_sl.transpose()).transpose());
                }
            };

    /* dependents: [W_tl, zz_sll] */
    private static final Function<Map<String, ? extends Duplicable>, ? extends Duplicable> calculate_WzzWT_st =
            new Function<Map<String, ? extends Duplicable>, Duplicable>() {
                @Override
                public Duplicable apply(Map<String, ? extends Duplicable> dat) {
                    final INDArray W_tl = getINDArrayFromMap(CoverageModelICGCacheNode.W_tl, dat);
                    final INDArray zz_sll = getINDArrayFromMap(CoverageModelICGCacheNode.zz_sll, dat);
                    final int numSamples = zz_sll.shape()[0];
                    final int numTargets = W_tl.shape()[0];
                    final INDArray WzzWT_st = Nd4j.create(numSamples, numTargets);
                    for (int si = 0; si < numSamples; si++) {
                        final INDArray zz_ll = zz_sll.get(NDArrayIndex.point(si), NDArrayIndex.all(), NDArrayIndex.all());
                        /* mean_W_contrib_t = \sum_{m,n} E[W_{tm}] E[W_{tn}] E[z_{sm} z_{sn}] */
                        final INDArray mean_W_contrib_t = W_tl.mmul(zz_ll).muli(W_tl).sum(1).transpose();
                        WzzWT_st.get(NDArrayIndex.point(si), NDArrayIndex.all()).assign(mean_W_contrib_t);
                    }
                    return new DuplicableNDArray(WzzWT_st);
                }
            };

    /* dependents: ["Sigma_st", "Psi_t"] */
    private static final Function<Map<String, ? extends Duplicable>, ? extends Duplicable> calculate_tot_Psi_st =
            new Function<Map<String, ? extends Duplicable>, Duplicable>() {
                @Override
                public Duplicable apply(Map<String, ? extends Duplicable> dat) {
                    final INDArray Sigma_st = getINDArrayFromMap(CoverageModelICGCacheNode.Sigma_st, dat);
                    final INDArray Psi_t = getINDArrayFromMap(CoverageModelICGCacheNode.Psi_t, dat);
                    final INDArray gamma_s = getINDArrayFromMap(CoverageModelICGCacheNode.gamma_s, dat);
                    return new DuplicableNDArray(Sigma_st.addRowVector(Psi_t).addiColumnVector(gamma_s));
                }
            };

    /* dependents: ["M_st", "tot_Psi_st", "log_n_st"] */
    private static final Function<Map<String, ? extends Duplicable>, ? extends Duplicable> calculate_loglike_normalization_s =
            new Function<Map<String, ? extends Duplicable>, Duplicable>() {
                @Override
                public Duplicable apply(Map<String, ? extends Duplicable> dat) {
                    final INDArray M_st = getINDArrayFromMap(CoverageModelICGCacheNode.M_st, dat);
                    final INDArray tot_Psi_st = getINDArrayFromMap(CoverageModelICGCacheNode.tot_Psi_st, dat);
                    final INDArray log_n_st = getINDArrayFromMap(CoverageModelICGCacheNode.log_n_st, dat);
                    final INDArray loglike_normalization_s = Transforms.log(tot_Psi_st, true)
                            .addi(FastMath.log(2 * FastMath.PI))
                            .muli(-0.5)
                            .subi(log_n_st)
                            .muli(M_st)
                            .sum(1);
                    return new DuplicableNDArray(loglike_normalization_s);
                }
            };

    /* dependents: ["M_st", "tot_Psi_st"] */
    private static final Function<Map<String, ? extends Duplicable>, ? extends Duplicable> calculate_M_Psi_inv_st =
            new Function<Map<String, ? extends Duplicable>, Duplicable>() {
                @Override
                public Duplicable apply(Map<String, ? extends Duplicable> dat) {
                    return new DuplicableNDArray(getINDArrayFromMap(CoverageModelICGCacheNode.M_st, dat).div(
                            getINDArrayFromMap(CoverageModelICGCacheNode.tot_Psi_st, dat)));
                }
            };

    /* dependents: ["M_Psi_inv_st", "Delta_st", "z_sl"] */
    private static final Function<Map<String, ? extends Duplicable>, ? extends Duplicable> calculate_v_tl =
            new Function<Map<String, ? extends Duplicable>, Duplicable>() {
                @Override
                public Duplicable apply(Map<String, ? extends Duplicable> dat) {
                    final INDArray M_Psi_inv_st = getINDArrayFromMap(CoverageModelICGCacheNode.M_Psi_inv_st, dat);
                    final INDArray Delta_st = getINDArrayFromMap(CoverageModelICGCacheNode.Delta_st, dat);
                    final INDArray z_sl = getINDArrayFromMap(CoverageModelICGCacheNode.z_sl, dat);
                    return new DuplicableNDArray(M_Psi_inv_st.mul(Delta_st).transpose().mmul(z_sl));
                }
            };

    /* dependents: ["M_Psi_inv_st", "zz_sll"] */
    private static final Function<Map<String, ? extends Duplicable>, ? extends Duplicable> calculate_Q_tll =
            new Function<Map<String, ? extends Duplicable>, Duplicable>() {
                @Override
                public Duplicable apply(Map<String, ? extends Duplicable> dat) {
                    final INDArray M_Psi_inv_st_trans = getINDArrayFromMap(CoverageModelICGCacheNode.M_Psi_inv_st, dat).transpose();
                    final INDArray zz_sll = getINDArrayFromMap(CoverageModelICGCacheNode.zz_sll, dat);
                    final int numTargets = M_Psi_inv_st_trans.shape()[0];
                    final int numLatents = zz_sll.shape()[1];
                    final INDArray res = Nd4j.create(numTargets, numLatents, numLatents);
                    for (int li = 0; li < numLatents; li++) {
                        res.get(NDArrayIndex.all(), NDArrayIndex.all(), NDArrayIndex.point(li)).assign(
                                    M_Psi_inv_st_trans.mmul(zz_sll.get(NDArrayIndex.all(), NDArrayIndex.all(), NDArrayIndex.point(li))));
                    }
                    return new DuplicableNDArray(res);
                }
            };

    /* dependents: ["Q_tll"] */
    private static final Function<Map<String, ? extends Duplicable>, ? extends Duplicable> calculate_sum_Q_ll =
            new Function<Map<String, ? extends Duplicable>, Duplicable>() {
                @Override
                public Duplicable apply(Map<String, ? extends Duplicable> dat) {
                    return new DuplicableNDArray(getINDArrayFromMap(CoverageModelICGCacheNode.Q_tll, dat).sum(0));
                }
            };

    /* dependents: ["Delta_st", "var_log_c_st", "var_log_d_s", "WzzWT_st", "Wz_st"] */
    private static final Function<Map<String, ? extends Duplicable>, ? extends Duplicable> calculate_B_st_with_bias_covariates =
            new Function<Map<String, ? extends Duplicable>, Duplicable>() {
                @Override
                public Duplicable apply(Map<String, ? extends Duplicable> dat) {
                    final INDArray Delta_st = getINDArrayFromMap(CoverageModelICGCacheNode.Delta_st, dat);
                    final INDArray var_log_c_st = getINDArrayFromMap(CoverageModelICGCacheNode.var_log_c_st, dat);
                    final INDArray var_log_d_s = getINDArrayFromMap(CoverageModelICGCacheNode.var_log_d_s, dat);
                    final INDArray WzzWT_st = getINDArrayFromMap(CoverageModelICGCacheNode.WzzWT_st, dat);
                    final INDArray Wz_st = getINDArrayFromMap(CoverageModelICGCacheNode.Wz_st, dat);
                    return new DuplicableNDArray(
                            Delta_st.mul(Delta_st)
                                    .addi(var_log_c_st)
                                    .addiColumnVector(var_log_d_s)
                                    .addi(WzzWT_st)
                                    .subi(Delta_st.mul(Wz_st.mul(2))));
                }
            };

    /* dependents: ["Delta_st", "var_log_c_st", "var_log_d_s"] */
    private static final Function<Map<String, ? extends Duplicable>, ? extends Duplicable> calculate_B_st_without_bias_covariates =
            new Function<Map<String, ? extends Duplicable>, Duplicable>() {
                @Override
                public Duplicable apply(Map<String, ? extends Duplicable> dat) {
                    final INDArray Delta_st = getINDArrayFromMap(CoverageModelICGCacheNode.Delta_st, dat);
                    final INDArray var_log_c_st = getINDArrayFromMap(CoverageModelICGCacheNode.var_log_c_st, dat);
                    final INDArray var_log_d_s = getINDArrayFromMap(CoverageModelICGCacheNode.var_log_d_s, dat);
                    return new DuplicableNDArray(Delta_st.mul(Delta_st).addi(var_log_c_st).addiColumnVector(var_log_d_s));
                }
            };

    /************************
     * log likelihood nodes *
     ************************/

    /* dependents: ["B_st", "M_Psi_inv_st", "loglike_normalization_s"] */
    private static final Function<Map<String, ? extends Duplicable>, ? extends Duplicable> calculate_loglike_unreg =
            new Function<Map<String, ? extends Duplicable>, Duplicable>() {
                @Override
                public Duplicable apply(Map<String, ? extends Duplicable> dat) {
                    /* fetch */
                    final INDArray B_st = getINDArrayFromMap(CoverageModelICGCacheNode.B_st, dat);
                    final INDArray M_Psi_inv_st = getINDArrayFromMap(CoverageModelICGCacheNode.M_Psi_inv_st, dat);
                    final INDArray loglike_normalization_s = getINDArrayFromMap(CoverageModelICGCacheNode.loglike_normalization_s, dat);

                    return new DuplicableNDArray(B_st.mul(M_Psi_inv_st).sum(1)
                            .muli(-0.5)
                            .addi(loglike_normalization_s));
                }
            };

    /* TODO github/gatk-protected issue #853 -- the filter contribution part could be cached */
    /* dependents: ["B_st", "M_Psi_inv_st", "loglike_normalization_s", "W_tl", "F_W_tl", "zz_sll"] */
    private static final Function<Map<String, ? extends Duplicable>, ? extends Duplicable> calculate_loglike_reg =
            new Function<Map<String, ? extends Duplicable>, Duplicable>() {
                @Override
                public Duplicable apply(Map<String, ? extends Duplicable> dat) {
                    final INDArray B_st = getINDArrayFromMap(CoverageModelICGCacheNode.B_st, dat);
                    final INDArray M_Psi_inv_st = getINDArrayFromMap(CoverageModelICGCacheNode.M_Psi_inv_st, dat);
                    final INDArray loglike_normalization_s = getINDArrayFromMap(CoverageModelICGCacheNode.loglike_normalization_s, dat);
                    final INDArray regularPart = B_st.mul(M_Psi_inv_st).sum(1)
                            .muli(-0.5)
                            .addi(loglike_normalization_s);

                    final INDArray W_tl = getINDArrayFromMap(CoverageModelICGCacheNode.W_tl, dat);
                    final INDArray F_W_tl = getINDArrayFromMap(CoverageModelICGCacheNode.F_W_tl, dat);
                    final INDArray zz_sll = getINDArrayFromMap(CoverageModelICGCacheNode.zz_sll, dat);
                    final INDArray WFW_ll = W_tl.transpose().mmul(F_W_tl).muli(-0.5);
                    final int numSamples = B_st.shape()[0];
                    final INDArray filterPart = Nd4j.create(new int[]{numSamples, 1},
                            IntStream.range(0, numSamples).mapToDouble(si ->
                                    zz_sll.get(NDArrayIndex.point(si)).mul(WFW_ll).sumNumber().doubleValue()).toArray());

                    return new DuplicableNDArray(regularPart.addi(filterPart));
                }
            };

    /**
     * This class represents the initial data passed to a compute block
     */
    public static final class InitialDataBlock implements Serializable {

        private static final long serialVersionUID = -1843744574954246765L;

        public final int[] readCountBlock;
        public final int[] maskBlock;
        public final double[] mappingErrorRateBlock;

        /**
         * Public constructor.
         *
         * @param readCountBlock read counts raveled in Fortran (column major) order
         * @param maskBlock learning mask raveled in Fortran (column major) order
         * @param mappingErrorRateBlock mapping error rate raveled in Fortran (column major) order
         */
        InitialDataBlock(@Nonnull final int[] readCountBlock,
                         @Nonnull final int[] maskBlock,
                         @Nonnull final double[] mappingErrorRateBlock) {
            this.readCountBlock = Utils.nonNull(readCountBlock, "Read count data block must be non-null");
            this.maskBlock = Utils.nonNull(maskBlock, "Mask data block must be non-null");
            this.mappingErrorRateBlock = Utils.nonNull(mappingErrorRateBlock, "Mapping error rate data block must" +
                    " be non-null");
        }

        void assertDataBlockSize(final int numTargets, final int numSamples) {
            final int size = numTargets * numSamples;
            Utils.validateArg(readCountBlock.length == size, "The read count data block has the wrong length");
            Utils.validateArg(maskBlock.length == size, "The learning mask data block has the wrong length");
            Utils.validateArg(mappingErrorRateBlock.length == size, "The mapping error rate data block has the" +
                    " wrong length");
        }
    }
}
