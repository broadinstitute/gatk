package org.broadinstitute.hellbender.tools.coveragemodel.germline;

import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.coveragemodel.*;
import org.broadinstitute.hellbender.tools.coveragemodel.interfaces.CopyRatioExpectationsCalculator;
import org.broadinstitute.hellbender.tools.coveragemodel.nd4jutils.Nd4jIOUtils;
import org.broadinstitute.hellbender.tools.exome.HashedListTargetCollection;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.tools.exome.TargetCollection;
import org.broadinstitute.hellbender.tools.exome.germlinehmm.IntegerCopyNumberState;
import org.broadinstitute.hellbender.tools.exome.sexgenotyper.GermlinePloidyAnnotatedTargetCollection;
import org.broadinstitute.hellbender.tools.exome.sexgenotyper.SexGenotypeData;
import org.broadinstitute.hellbender.tools.exome.sexgenotyper.SexGenotypeDataCollection;
import org.broadinstitute.hellbender.utils.hmm.segmentation.HiddenMarkovModelSegmentProcessor;
import org.broadinstitute.hellbender.utils.hmm.segmentation.HiddenStateSegmentRecordWriter;

import javax.annotation.Nonnull;
import javax.annotation.Nullable;
import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.function.BiFunction;
import java.util.stream.Collectors;

/**
 * This class represents the driver-node workspace for EM algorithm calculations of the coverage model
 * specific for germline copy number states.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class CoverageModelEMWorkspaceGermline extends CoverageModelEMWorkspace<IntegerCopyNumberState> {

    /**
     * Public constructor
     *
     * @param rawReadCounts an instance of {@link ReadCountCollection} containing raw read counts
     * @param germlinePloidyAnnotatedTargetCollection an instance of {@link GermlinePloidyAnnotatedTargetCollection}
     *                                                for obtaining target ploidies for different sex genotypes
     * @param sexGenotypeDataCollection an instance of {@link SexGenotypeDataCollection} for obtaining sample sex genotypes
     * @param params EM algorithm parameters
     * @param ctx the Spark context
     * @param copyRatioExpectationsCalculator an implementation of {@link CopyRatioExpectationsCalculator}
     */
    public CoverageModelEMWorkspaceGermline(@Nonnull final ReadCountCollection rawReadCounts,
                                            @Nonnull final GermlinePloidyAnnotatedTargetCollection germlinePloidyAnnotatedTargetCollection,
                                            @Nonnull final SexGenotypeDataCollection sexGenotypeDataCollection,
                                            @Nonnull final IntegerCopyNumberExpectationsCalculator copyRatioExpectationsCalculator,
                                            @Nonnull final CoverageModelEMParams params,
                                            @Nullable final CoverageModelParameters model,
                                            @Nullable final JavaSparkContext ctx) {
        super(rawReadCounts, germlinePloidyAnnotatedTargetCollection, sexGenotypeDataCollection,
                copyRatioExpectationsCalculator, params, model, ctx);
    }

    /**
     * Save copy ratio (here, integer copy number) posteriors to disk
     *
     * @param outputPath the output path
     */
    protected void saveCopyRatioPosteriors(final String outputPath) {
        /* fetch forward-backward and Viterbi results */
        final List<CopyRatioHiddenMarkovModelResults<CoverageModelCopyRatioEmissionData,
                        IntegerCopyNumberState>> copyRatioHMMResult = getCopyRatioHiddenMarkovModelResults();

        /* perform segmentation */
        final BiFunction<SexGenotypeData, Target, IntegerCopyNumberState> referenceStateFactory = (sexGenotypeData, target) ->
                new IntegerCopyNumberState(germlinePloidyAnnotatedTargetCollection.getTargetGermlinePloidyByGenotype(target,
                        sexGenotypeData.getSexGenotype()));
        final List<TargetCollection<Target>> sampleTargetCollection = Collections.nCopies(numSamples,
                new HashedListTargetCollection<>(processedTargetList));
        final HiddenMarkovModelSegmentProcessor<CoverageModelCopyRatioEmissionData, IntegerCopyNumberState, Target>
                copyRatioSegmenter = new HiddenMarkovModelSegmentProcessor<>(processedSampleNameList,
                processedSampleSexGenotypeData, referenceStateFactory, sampleTargetCollection,
                copyRatioHMMResult.stream()
                                .map(CopyRatioHiddenMarkovModelResults::getForwardBackwardResult)
                                .collect(Collectors.toList()),
                        copyRatioHMMResult.stream()
                                .map(CopyRatioHiddenMarkovModelResults::getViterbiResult)
                                .collect(Collectors.toList()));

        /* save copy number segments to disk */
        final File segmentsFile = new File(outputPath, CoverageModelGlobalConstants.COPY_RATIO_SEGMENTS_FILENAME);
        try (final HiddenStateSegmentRecordWriter<IntegerCopyNumberState, Target> segWriter =
                     new HiddenStateSegmentRecordWriter<>(segmentsFile)) {
            copyRatioSegmenter.writeSegmentsToTableWriter(segWriter);
        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile(segmentsFile, "Could not create copy ratio segments file");
        }

        /* save the Viterbi copy number chains as a sample-target matrix */
        if (params.extendedPosteriorOutputEnabled()) {
            final File copyRatioViterbiFile = new File(outputPath, CoverageModelGlobalConstants.COPY_RATIO_VITERBI_FILENAME);
            final List<String> targetNames = processedReadCounts.targets().stream()
                    .map(Target::getName).collect(Collectors.toList());
            Nd4jIOUtils.writeNDArrayToTextFile(getViterbiAsNDArray(copyRatioHMMResult, sampleTargetCollection),
                    copyRatioViterbiFile, processedSampleNameList, targetNames);
        }
    }
}
