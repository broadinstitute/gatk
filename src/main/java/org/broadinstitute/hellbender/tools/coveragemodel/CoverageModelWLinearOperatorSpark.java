package org.broadinstitute.hellbender.tools.coveragemodel;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.tools.coveragemodel.linalg.FourierLinearOperatorNDArray;
import org.broadinstitute.hellbender.tools.coveragemodel.linalg.GeneralLinearOperator;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;
import org.nd4j.linalg.indexing.NDArrayIndex;

import javax.annotation.Nonnull;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Implementes the (real and symmetric) linear operator,
 *
 *      Q_{t}^{\mu\nu}\delta_{t,t'} + F_{t,t'}Z^{\mu\nu}.
 *
 * It acts on W_{t'}^{\nu} and returns an INDArray of similar shape.
 *
 * It uses Spark for calculating the first term.
 *
 * TODO github/gatk-protected issue #701 -- this class is part of the CNV-avoiding regularizer and will undergo
 * significant changes soon. To the reviewer: let this go through as is for the time being.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class CoverageModelWLinearOperatorSpark extends GeneralLinearOperator<INDArray> {

    private final Logger logger = LogManager.getLogger(CoverageModelWLinearOperatorSpark.class);

    private final int numLatents, numTargets;
    private final INDArray Z_ll;
    private final FourierLinearOperatorNDArray F_tt;

    /* sparky stuff */
    private final JavaSparkContext ctx;
    private final List<LinearSpaceBlock> targetSpaceBlocks;
    private final JavaPairRDD<LinearSpaceBlock, CoverageModelEMComputeBlock> computeRDD;

    public CoverageModelWLinearOperatorSpark(@Nonnull final INDArray Z_ll,
                                             @Nonnull final FourierLinearOperatorNDArray F_tt,
                                             final int numTargets,
                                             @Nonnull final JavaSparkContext ctx,
                                             @Nonnull final JavaPairRDD<LinearSpaceBlock, CoverageModelEMComputeBlock> computeRDD,
                                             @Nonnull final List<LinearSpaceBlock> targetSpaceBlocks) {
        this.numTargets = numTargets;
        numLatents = Z_ll.shape()[0];
        if (Z_ll.shape()[0] != numLatents || Z_ll.shape()[1] != numLatents)
            throw new IllegalArgumentException("Malformed Z_ll.");
        if (F_tt.getRowDimension() != numTargets || F_tt.getColumnDimension() != numTargets)
            throw new IllegalArgumentException("Malformed F_tt.");
        this.Z_ll = Z_ll;
        this.F_tt = F_tt;

        /* sparky stuff */
        this.ctx = ctx;
        this.targetSpaceBlocks = targetSpaceBlocks;
        this.computeRDD = computeRDD;
    }

    @Override
    public int getRowDimension() {
        return numLatents * numTargets;
    }

    @Override
    public int getColumnDimension() {
        return numLatents * numTargets;
    }

    @Override
    public INDArray operate(@Nonnull final INDArray W_tl) throws DimensionMismatchException {
        if (W_tl.rank() != 2 || W_tl.shape()[0] != numTargets || W_tl.shape()[1] != numLatents)
            throw new DimensionMismatchException(W_tl.length(), numTargets * numLatents);

        /* Z F W */
        final long startTimeZFW = System.nanoTime();
        final INDArray Z_F_W_tl = Nd4j.create(numTargets, numLatents);
        IntStream.range(0, numLatents).parallel().forEach(li ->
                Z_F_W_tl.get(NDArrayIndex.all(), NDArrayIndex.point(li)).assign(
                        F_tt.operate(W_tl.get(NDArrayIndex.all(), NDArrayIndex.point(li)))));
        Z_F_W_tl.assign(Nd4j.gemm(Z_F_W_tl, Z_ll, false, false));
        final long endTimeZFW = System.nanoTime();

        /* perform a broadcast hash join */
        final long startTimeQW = System.nanoTime();
        final Map<LinearSpaceBlock, INDArray> W_tl_map = CoverageModelSparkUtils.partitionINDArrayToAMap(targetSpaceBlocks, W_tl);
        final Broadcast<Map<LinearSpaceBlock, INDArray>> W_tl_bc = ctx.broadcast(W_tl_map);
        final INDArray Q_W_tl = CoverageModelSparkUtils.assembleINDArrayBlocksFromRDD(
                computeRDD.mapValues(cb -> {
                    final INDArray W_tl_chunk = W_tl_bc.value().get(cb.getTargetSpaceBlock());
                    final INDArray Q_tll_chunk = cb.getINDArrayFromCache(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.Q_tll);
                    final Collection<INDArray> W_Q_chunk = IntStream.range(0, cb.getTargetSpaceBlock().getNumTargets()).parallel()
                            .mapToObj(ti -> Q_tll_chunk.get(NDArrayIndex.point(ti))
                                    .mmul(W_tl_chunk.get(NDArrayIndex.point(ti)).transpose()))
                            .collect(Collectors.toList());
                    return Nd4j.vstack(W_Q_chunk);
                }), 0);
        W_tl_bc.destroy();
//        final JavaPairRDD<LinearSpaceBlock, INDArray> W_tl_RDD = CoverageModelSparkUtils.rddFromINDArray(W_tl,
//                targetSpaceBlocks, ctx, true);
//        final INDArray Q_W_tl = CoverageModelSparkUtils.assembleINDArrayBlocks(
//                computeRDD.join(W_tl_RDD).mapValues(p -> {
//                    final CoverageModelEMComputeBlock cb = p._1;
//                    final INDArray W_tl_chunk = p._2;
//                    final INDArray Q_tll_chunk = cb.getINDArrayFromCache("Q_tll");
//                    return Nd4j.vstack(IntStream.range(0, cb.getTargetSpaceBlock().getNumTargets()).parallel()
//                            .mapToObj(ti -> Q_tll_chunk.get(NDArrayIndex.point(ti)).mmul(W_tl_chunk.get(NDArrayIndex.point(ti)).transpose()))
//                            .collect(Collectors.toList()));
//                }), false);
//        W_tl_RDD.unpersist();
        final long endTimeQW = System.nanoTime();

        logger.debug("Local [Z] [F] [W] timing: " + (endTimeZFW - startTimeZFW)/1000000 + " ms");
        logger.debug("Spark [Q] [W] timing: " + (endTimeQW - startTimeQW)/1000000 + " ms");

        return Q_W_tl.addi(Z_F_W_tl);
    }

    @Override
    public INDArray operateTranspose(@Nonnull final INDArray x)
            throws DimensionMismatchException, UnsupportedOperationException {
        return operate(x);
    }

    @Override
    public boolean isTransposable() {
        return true;
    }
}
