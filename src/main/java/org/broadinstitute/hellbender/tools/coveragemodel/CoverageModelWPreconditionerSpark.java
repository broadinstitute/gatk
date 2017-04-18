package org.broadinstitute.hellbender.tools.coveragemodel;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.tools.coveragemodel.linalg.FourierLinearOperatorNDArray;
import org.broadinstitute.hellbender.tools.coveragemodel.linalg.GeneralLinearOperator;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;
import org.nd4j.linalg.indexing.NDArrayIndex;
import scala.Tuple2;

import javax.annotation.Nonnull;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * TODO github/gatk-protected issue #701 -- this class is part of the CNV-avoiding regularizer and will undergo
 * significant changes soon. To the reviewer: let this go through as is for the time being.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class CoverageModelWPreconditionerSpark implements GeneralLinearOperator<INDArray> {

    private static final Logger logger = LogManager.getLogger(CoverageModelWPreconditionerSpark.class);

    private final int numLatents, numTargets, fftSize;
    private final FourierLinearOperatorNDArray F_tt;
    private final INDArray orderedFourierFactors;

    /* sparky stuff */
    private final JavaSparkContext ctx;
    private final List<LinearlySpacedIndexBlock> fourierSpaceBlocks;
    private final JavaPairRDD<LinearlySpacedIndexBlock, INDArray> linOpPairRDD;

    public CoverageModelWPreconditionerSpark(@Nonnull final INDArray Q_ll,
                                             @Nonnull final INDArray Z_ll,
                                             @Nonnull final FourierLinearOperatorNDArray F_tt,
                                             final int numTargets,
                                             @Nonnull JavaSparkContext ctx,
                                             final int numPartitions) {
        this.numTargets = ParamUtils.isPositive(numTargets, "Number of target must be positive");
        this.numLatents = Q_ll.shape()[0];
        this.fftSize = F_tt.getFFTSize();
        if (Q_ll.shape()[1] != numLatents)
            throw new IllegalArgumentException("Malformed Q_ll.");
        if (Z_ll.shape()[0] != numLatents || Z_ll.shape()[1] != numLatents)
            throw new IllegalArgumentException("Malformed Z_ll.");
        if (F_tt.getRowDimension() != numTargets || F_tt.getColumnDimension() != numTargets)
            throw new IllegalArgumentException("Malformed F_tt.");
        this.F_tt = F_tt;
        orderedFourierFactors = Nd4j.create(F_tt.getOrderedFourierFactors(), new int[]{fftSize, 1});

        /* sparky stuff */
        this.ctx = ctx;

        fourierSpaceBlocks = CoverageModelSparkUtils.createLinearlySpacedIndexBlocks(fftSize, numPartitions, 1);

        final INDArray linOp = Nd4j.create(fftSize, numLatents, numLatents);
        IntStream.range(0, fftSize).parallel().forEach(k -> linOp.get(NDArrayIndex.point(k)).assign(
                Z_ll.mul(orderedFourierFactors.getDouble(k)).addi(Q_ll)));
//        linOpPairRDD = CoverageModelSparkUtils.rddFromINDArray(linOp, fourierSpaceBlocks, ctx, true);
        /* for a broadcast hash join, repartitioning is unncessary */
        linOpPairRDD = ctx.parallelizePairs(
                CoverageModelSparkUtils.partitionINDArrayToList(fourierSpaceBlocks, linOp), fourierSpaceBlocks.size());
        linOpPairRDD.cache();
    }

    @Override
    public INDArray operate(@Nonnull final INDArray W_tl)
            throws DimensionMismatchException {
        if (W_tl.rank() != 2 || W_tl.shape()[0] != numTargets || W_tl.shape()[1] != numLatents) {
            throw new DimensionMismatchException(W_tl.length(), numTargets * numLatents);
        }
        long startTimeRFFT = System.nanoTime();
        /* forward rfft */
        final INDArray W_kl = Nd4j.create(fftSize, numLatents);
        IntStream.range(0, numLatents).parallel().forEach(li ->
                W_kl.get(NDArrayIndex.all(), NDArrayIndex.point(li)).assign(
                        Nd4j.create(F_tt.getForwardFFT(W_tl.get(NDArrayIndex.all(), NDArrayIndex.point(li))),
                                new int[]{fftSize, 1})));
        long endTimeRFFT = System.nanoTime();

        /* apply the preconditioner in the Fourier space */
        long startTimePrecond = System.nanoTime();
        final Map<LinearlySpacedIndexBlock, INDArray> W_kl_map = CoverageModelSparkUtils.partitionINDArrayToMap(fourierSpaceBlocks, W_kl);
        final Broadcast<Map<LinearlySpacedIndexBlock, INDArray>> W_kl_bc = ctx.broadcast(W_kl_map);
        final JavaPairRDD<LinearlySpacedIndexBlock, INDArray> preconditionedWRDD = linOpPairRDD
                .mapToPair(p -> {
                    final INDArray W_kl_chuck = W_kl_bc.value().get(p._1);
                    final INDArray linOp_chunk = p._2;
                    final int blockSize = linOp_chunk.shape()[0];
                    final List<INDArray> linOpWList = IntStream.range(0, blockSize).parallel()
                            .mapToObj(k -> CoverageModelEMWorkspaceMathUtils.linsolve(linOp_chunk.get(NDArrayIndex.point(k)),
                                    W_kl_chuck.get(NDArrayIndex.point(k))))
                            .collect(Collectors.toList());
                    return new Tuple2<>(p._1, Nd4j.vstack(linOpWList));
                });
        W_kl.assign(CoverageModelSparkUtils.assembleINDArrayBlocksFromRDD(preconditionedWRDD, 0));
        W_kl_bc.destroy();
//        final JavaPairRDD<LinearlySpacedIndexBlock, INDArray> W_kl_RDD = CoverageModelSparkUtils.rddFromINDArray(W_kl,
//                fourierSpaceBlocks, ctx, true);
//        W_kl.assign(CoverageModelSparkUtils.assembleINDArrayBlocks(linOpPairRDD.join((W_kl_RDD))
//                .mapValues(p -> {
//                    final INDArray linOp = p._1;
//                    final INDArray W = p._2;
//                    final int blockSize = linOp.shape()[0];
//                    final List<INDArray> linOpWList = IntStream.range(0, blockSize).parallel().mapToObj(k ->
//                            CoverageModelEMWorkspaceMathUtils.linsolve(linOp.get(NDArrayIndex.point(k)),
//                                    W.get(NDArrayIndex.point(k))))
//                            .collect(Collectors.toList());
//                    return Nd4j.vstack(linOpWList);
//                }), false));
//        W_kl_RDD.unpersist();
        long endTimePrecond = System.nanoTime();

        /* irfft */
        long startTimeIRFFT = System.nanoTime();
        final INDArray res = Nd4j.create(numTargets, numLatents);
        IntStream.range(0, numLatents).parallel().forEach(li ->
                res.get(NDArrayIndex.all(), NDArrayIndex.point(li)).assign(
                        F_tt.getInverseFFT(W_kl.get(NDArrayIndex.all(), NDArrayIndex.point(li)))));
        long endTimeIRFFT = System.nanoTime();

        logger.debug("Local FFT timing: " + (endTimeRFFT - startTimeRFFT + endTimeIRFFT - startTimeIRFFT)/1000000 + " ms");
        logger.debug("Spark preconditioner application timing: " + (endTimePrecond - startTimePrecond)/1000000 + " ms");

        return res;
    }

    @Override
    public int getRowDimension() {
        return numTargets * numLatents;
    }

    @Override
    public int getColumnDimension() {
        return numTargets * numLatents;
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

    /**
     * Unpersists the caches
     */
    @Override
    public void cleanupAfter() {
        linOpPairRDD.unpersist();
    }
}
