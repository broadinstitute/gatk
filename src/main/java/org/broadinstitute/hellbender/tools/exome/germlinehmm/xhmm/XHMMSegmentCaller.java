package org.broadinstitute.hellbender.tools.exome.germlinehmm.xhmm;

import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.NormalizeSomaticReadCounts;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.tools.exome.TargetCollection;
import org.broadinstitute.hellbender.tools.exome.germlinehmm.CopyNumberTriState;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.hmm.segmentation.HiddenMarkovModelPostProcessor;
import org.broadinstitute.hellbender.utils.hmm.segmentation.HiddenStateSegmentRecordWriter;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * Discovers potential copy-number segments in input read-counts as the most likely segment sequences
 * based on {@link XHMMModel} HMM model.
 *
 * <p>
 *     You normally want to run this tool on values normalized using {@link NormalizeSomaticReadCounts}
 *     which should get rid of systematic biases due to sequencing and capture technology
 * </p>
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        programGroup = CopyNumberProgramGroup.class,
        summary = "Find possible locations for rare copy number variation events in germline samples using a HMM",
        oneLineSummary = "Discover possible locations of copy number variation"
)
public final class XHMMSegmentCaller extends XHMMSegmentCallerBase {

    private HiddenStateSegmentRecordWriter<CopyNumberTriState, Target> outputWriter;

    @Override
    protected void openOutput(final File outputFile, final XHMMModel model,
                              final TargetCollection<Target> targets, final ReadCountCollection inputCounts) {
        Utils.nonNull(outputFile);
        Utils.nonNull(model);
        Utils.nonNull(targets);
        Utils.nonNull(inputCounts);
        try {
            outputWriter = new HiddenStateSegmentRecordWriter<>(outputFile);
        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, ex);
        }
    }

    @Override
    protected void closeOutput(final File outputFile) {
        try {
            outputWriter.close();
        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, ex);
        }
    }

    @Override
    protected void makeCalls(XHMMModel model, TargetCollection<Target> targets, ReadCountCollection inputCounts) {
        /* perform segmentation and write calls to outputWriter */
        final HiddenMarkovModelPostProcessor<XHMMEmissionData, CopyNumberTriState, Target> processor =
                new HiddenMarkovModelPostProcessor<>(inputCounts.columnNames(),
                        Collections.nCopies(inputCounts.columnNames().size(), targets),
                        sampleForwardBackwardResults, sampleBestPaths,  CopyNumberTriState.NEUTRAL);
        processor.writeSegmentsToTableWriter(outputWriter);
    }
}
