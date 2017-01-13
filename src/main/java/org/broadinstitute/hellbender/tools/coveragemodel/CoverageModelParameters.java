package org.broadinstitute.hellbender.tools.coveragemodel;

import com.google.cloud.dataflow.sdk.repackaged.com.google.common.collect.Sets;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.math3.util.FastMath;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.coveragemodel.nd4jutils.Nd4jIOUtils;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.tools.exome.TargetTableReader;
import org.broadinstitute.hellbender.tools.exome.TargetWriter;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;
import org.nd4j.linalg.indexing.NDArrayIndex;

import javax.annotation.Nonnull;
import java.io.*;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * This class reads, writes, and stores the coverage model parameters.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class CoverageModelParameters implements Serializable {

    private static final long serialVersionUID = -4350342293001054849L;

    private final List<Target> targetList;

    /* 1 x T */
    private final INDArray targetMeanLogBias;

    /* 1 x T */
    private final INDArray targetUnexplainedVariance;

    /* T x L */
    private final INDArray biasCovariates;

    private final int numTargets, numLatents;

    /**
     * Public constructor.
     *
     * @param targetMeanLogBias target-specific mean log bias (m_t)
     * @param targetUnexplainedVariance target-specific unexplained bias variance (\Psi_t)
     * @param biasCovariates bias covariates matrix (W_{t\mu})
     */
    public CoverageModelParameters(@Nonnull final List<Target> targetList,
                                   @Nonnull final INDArray targetMeanLogBias,
                                   @Nonnull final INDArray targetUnexplainedVariance,
                                   @Nonnull final INDArray biasCovariates) {
        this.targetList = Utils.nonNull(targetList, "Target list must be non-null");
        this.targetMeanLogBias = Utils.nonNull(targetMeanLogBias, "Target-specific mean log bias must be non-null");
        this.targetUnexplainedVariance = Utils.nonNull(targetUnexplainedVariance, "Target-specific unexplained variance" +
                " must be non-null");
        this.biasCovariates = Utils.nonNull(biasCovariates, "Bias covariates matrix must be non-null");
        this.numTargets = targetList.size();
        this.numLatents = biasCovariates.size(1);

        /* check the dimensions */
        Utils.validateArg(numTargets == targetMeanLogBias.size(1) && numTargets == targetUnexplainedVariance.size(1) &&
                numTargets == biasCovariates.size(0), "The dimension of target space (expected: " + numTargets + ")" +
                " does not agree with the dimensions of the provided model data: " +
                    "target mean bias = " + targetMeanLogBias.shapeInfoToString() + ", " +
                    "target unexplained variance = " + targetUnexplainedVariance.shapeInfoToString() + ", " +
                    "principal latent to target map = " + biasCovariates.shapeInfoToString());
    }

    public INDArray getTargetMeanLogBias() {
        return targetMeanLogBias;
    }

    public INDArray getTargetUnexplainedVariance() {
        return targetUnexplainedVariance;
    }

    public INDArray getBiasCovariates() {
        return biasCovariates;
    }

    public INDArray getTargetMeanBiasOnTargetBlock(@Nonnull final LinearSpaceBlock tb) {
        checkTargetBlock(tb);
        return targetMeanLogBias.get(NDArrayIndex.all(), NDArrayIndex.interval(tb.getBegIndex(), tb.getEndIndex()));
    }

    private void checkTargetBlock(@Nonnull LinearSpaceBlock tb) {
        ParamUtils.inRange(tb.getBegIndex(), 0, numTargets, "The begin index of target block is out of range");
        ParamUtils.inRange(tb.getEndIndex(), 0, numTargets, "The begin index of target block is out of range");
    }

    public INDArray getTargetUnexplainedVarianceOnTargetBlock(@Nonnull final LinearSpaceBlock tb) {
        checkTargetBlock(tb);
        return targetUnexplainedVariance.get(NDArrayIndex.all(), NDArrayIndex.interval(tb.getBegIndex(), tb.getEndIndex()));
    }

    public INDArray getBiasCovariatesOnTargetBlock(@Nonnull final LinearSpaceBlock tb) {
        checkTargetBlock(tb);
        return biasCovariates.get(NDArrayIndex.interval(tb.getBegIndex(), tb.getEndIndex()), NDArrayIndex.all());
    }

    private static void createOutputPath(final String outputPath) {
        final File outputPathFile = new File(outputPath);
        if (!outputPathFile.exists()) {
            if (!outputPathFile.mkdirs()) {
                throw new UserException.CouldNotCreateOutputFile(outputPathFile, "Could not create the output directory");
            }
        }
    }

    public int getNumTargets() {
        return numTargets;
    }

    public int getNumLatents() {
        return numLatents;
    }

    public List<Target> getTargetList() {
        return targetList;
    }

    /**
     * Check whether the bias covariates are orthogonal to each other or not.
     *
     * @param tol orthogonality tolerance
     * @param logResults print results or not
     */
    public boolean checkOrthogonalityOfBiasCovariates(final double tol, final boolean logResults,
                                                      @Nonnull final Logger logger) {
        ParamUtils.isPositive(tol, "Orthogonality tolerance must be positive");
        boolean orthogonal = true;
        for (int mu = 0; mu < numLatents; mu++) {
            for (int nu = mu; nu < numLatents; nu++) {
                final double innerProd = biasCovariates.get(NDArrayIndex.all(),
                        NDArrayIndex.point(mu)).mul(biasCovariates.get(NDArrayIndex.all(),
                        NDArrayIndex.point(nu))).sumNumber().doubleValue();
                if (mu != nu && FastMath.abs(innerProd) > tol) {
                    orthogonal = false;
                    logger.info("Inner product test failed on (" + mu + ", " + nu + ")");
                }
                if (logResults) {
                    logger.info("Inner product of (" + mu + ", " + nu + "): " + innerProd);
                } /* no need to continue */
                if (!logResults && !orthogonal) {
                    break;
                }
            }
        }
        return orthogonal;
    }

    /**
     * Reads the model from disk
     *
     * @param modelPath input model path
     * @return an instance of {@link CoverageModelParameters}
     */
    public static CoverageModelParameters read(@Nonnull final String modelPath) {
        final File modelPathFile = new File(Utils.nonNull(modelPath, "The input model path must be non-null"));
        Utils.validateArg(modelPathFile.exists(), "The model path does not exist: " + modelPathFile.getAbsolutePath());

        final File targetListFile = new File(modelPath, CoverageModelGlobalConstants.TARGET_LIST_OUTPUT_FILE);
        final List<Target> targetList;
        try (final Reader reader = new FileReader(targetListFile)) {
            targetList = TargetTableReader.readTargetFromReader(targetListFile.getAbsolutePath(), reader);
        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile(targetListFile, "Could not read targets interval list");
        }

        final File targetMeanLogBiasFile = new File(modelPath, CoverageModelGlobalConstants.TARGET_MEAN_LOG_BIAS_OUTPUT_FILE);
        final INDArray targetMeanLogBias = Nd4jIOUtils.readNDArrayFromTextFile(targetMeanLogBiasFile);

        final File targetUnexplainedVarianceFile = new File(modelPath, CoverageModelGlobalConstants.TARGET_UNEXPLAINED_VARIANCE_OUTPUT_FILE);
        final INDArray targetUnexplainedVariance = Nd4jIOUtils.readNDArrayFromTextFile(targetUnexplainedVarianceFile);

        final File biasCovariatesFile = new File(modelPath, CoverageModelGlobalConstants.BIAS_COVARIATES_OUTPUT_FILE);
        final INDArray biasCovariates = Nd4jIOUtils.readNDArrayFromTextFile(biasCovariatesFile);

        return new CoverageModelParameters(targetList, targetMeanLogBias, targetUnexplainedVariance, biasCovariates);
    }

    /**
     * Writes the model to disk
     *
     * @param outputPath model output path
     */
    public static void write(@Nonnull CoverageModelParameters model, @Nonnull final String outputPath) {
        /* create output directory if it doesn't exist */
        createOutputPath(Utils.nonNull(outputPath, "The output path string must be non-null"));

        /* write targets list */
        final File targetListFile = new File(outputPath, CoverageModelGlobalConstants.TARGET_LIST_OUTPUT_FILE);
        TargetWriter.writeTargetsToFile(targetListFile, model.getTargetList());

        final List<String> targetNames = model.getTargetList().stream()
                .map(Target::getName).collect(Collectors.toList());

        /* write target mean bias to file */
        final File targetMeanBiasFile = new File(outputPath, CoverageModelGlobalConstants.TARGET_MEAN_LOG_BIAS_OUTPUT_FILE);
        Nd4jIOUtils.writeNDArrayToTextFile(model.getTargetMeanLogBias(), targetMeanBiasFile, null, targetNames);

        /* write target unexplained variance to file */
        final File targetUnexplainedVarianceFile = new File(outputPath, CoverageModelGlobalConstants.TARGET_UNEXPLAINED_VARIANCE_OUTPUT_FILE);
        Nd4jIOUtils.writeNDArrayToTextFile(model.getTargetUnexplainedVariance(), targetUnexplainedVarianceFile,
                null, targetNames);

        /* write bias covariates to file */
        final List<String> biasCovariatesNames = IntStream.range(0, model.getNumLatents())
                .mapToObj(li -> String.format("BC_%d", li)).collect(Collectors.toList());
        final File biasCovariatesFile = new File(outputPath, CoverageModelGlobalConstants.BIAS_COVARIATES_OUTPUT_FILE);
        Nd4jIOUtils.writeNDArrayToTextFile(model.getBiasCovariates(), biasCovariatesFile, targetNames, biasCovariatesNames);

        /* write norm_2 of bias covariates to file */
        final double[] biasCovariatesNorm2 = new double[model.numLatents];
        final INDArray WTW = model.getBiasCovariates().transpose().mmul(model.getBiasCovariates());
        for (int li = 0; li < model.getNumLatents(); li++) {
            biasCovariatesNorm2[li] = WTW.getDouble(li, li);
        }
        final File biasCovariatesNorm2File = new File(outputPath, CoverageModelGlobalConstants.BIAS_COVARIATES_NORM2_OUTPUT_FILE);
        Nd4jIOUtils.writeNDArrayToTextFile(Nd4j.create(biasCovariatesNorm2, new int[] {1, model.getNumLatents()}),
                biasCovariatesNorm2File, Collections.singletonList("NORM_2"), biasCovariatesNames);
    }

    /**
     * This method "adapts" a model to a read count collection in the following sense:
     *
     *     - removes targets that are not included in the model from the read counts collection
     *     - removes targets that are in the read count collection from the model
     *     - rearranges model targets in the same order as read count collection targets
     *
     * The modifications are not done in-plane and the original input parameters remain intact.
     *
     * @param model a model
     * @param readCounts a read count collection
     * @return a pair of model and read count collection
     */
    public static ImmutablePair<CoverageModelParameters, ReadCountCollection> adaptModelToReadCountCollection(
            @Nonnull final CoverageModelParameters model, @Nonnull final ReadCountCollection readCounts,
            @Nonnull final Logger logger) {
        logger.info("Adapting model to read counts...");
        Utils.nonNull(model, "The model parameters must be non-null");
        Utils.nonNull(readCounts, "The read count collection must be non-null");
        Utils.nonNull(logger, "The logger must be non-null");

        final List<Target> modelTargetList = model.getTargetList();
        final List<Target> readCountsTargetList = readCounts.targets();
        final Set<Target> mutualTargetList = Sets.intersection(new HashSet<>(modelTargetList),
                new HashSet<>(readCountsTargetList));
        final List<Target> finalTargetList = readCountsTargetList.stream()
                .filter(mutualTargetList::contains)
                .collect(Collectors.toList());
        final Set<Target> finalTargetsSet = new LinkedHashSet<>(finalTargetList);

        logger.info("Number of mutual targets: " + finalTargetList.size());
        Utils.validateArg(finalTargetList.size() > 0, "The intersection between model targets and targets from read count" +
                    " collection is empty. Please check there the model is compatible with the given read count" +
                    " collection.");

        if (modelTargetList.size() > finalTargetList.size()) {
            logger.info("The following targets dropped from the model: " + Sets.difference(new HashSet<>(modelTargetList),
                    finalTargetsSet).stream().map(Target::getName).collect(Collectors.joining(", ", "[", "]")));
        }

        if (readCountsTargetList.size() > finalTargetList.size()) {
            logger.info("The following targets dropped from read counts: " + Sets.difference(new HashSet<>(readCountsTargetList),
                    finalTargetsSet).stream().map(Target::getName).collect(Collectors.joining(", ", "[", "]")));
        }

        /* the targets in {@code subsetReadCounts} follow the original order of targets in {@code readCounts} */
        final ReadCountCollection subsetReadCounts = readCounts.subsetTargets(finalTargetsSet);

        /* fetch original model parameters */
        final INDArray originalModelTargetMeanBias = model.getTargetMeanLogBias();
        final INDArray originalModelTargetUnexplainedVariance = model.getTargetUnexplainedVariance();
        final INDArray originalModelPrincipalLatentToTargetMap = model.getBiasCovariates();

        /* re-arrange targets */
        final int[] newTargetIndicesInOriginalModel = finalTargetList.stream()
                .mapToInt(modelTargetList::indexOf)
                .toArray();
        final INDArray newModelTargetMeanBias = Nd4j.create(new int[] {1, finalTargetList.size()});
        final INDArray newModelTargetUnexplainedVariance = Nd4j.create(new int[] {1, finalTargetList.size()});
        final INDArray newModelPrincipalLatentToTargetMap = Nd4j.create(new int[] {finalTargetList.size(),
                model.getNumLatents()});
        IntStream.range(0, finalTargetList.size())
                .forEach(ti -> {
                    newModelTargetMeanBias.put(0, ti,
                            originalModelTargetMeanBias.getDouble(0, newTargetIndicesInOriginalModel[ti]));
                    newModelTargetUnexplainedVariance.put(0, ti,
                            originalModelTargetUnexplainedVariance.getDouble(0, newTargetIndicesInOriginalModel[ti]));
                    newModelPrincipalLatentToTargetMap.get(NDArrayIndex.point(ti), NDArrayIndex.all())
                            .assign(originalModelPrincipalLatentToTargetMap.get(NDArrayIndex.point(newTargetIndicesInOriginalModel[ti]),
                                    NDArrayIndex.all()));
                });

        return ImmutablePair.of(new CoverageModelParameters(finalTargetList, newModelTargetMeanBias,
                newModelTargetUnexplainedVariance, newModelPrincipalLatentToTargetMap), subsetReadCounts);
    }

}
