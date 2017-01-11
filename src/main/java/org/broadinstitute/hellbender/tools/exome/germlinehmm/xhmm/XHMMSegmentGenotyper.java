package org.broadinstitute.hellbender.tools.exome.germlinehmm.xhmm;

import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.tools.exome.TargetCollection;
import org.broadinstitute.hellbender.tools.exome.germlinehmm.CopyNumberTriState;
import org.broadinstitute.hellbender.utils.hmm.segmentation.HiddenMarkovModelPostProcessor;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.io.File;
import java.util.*;

/**
 * Genotype predetermined segments passed in the inputs together with the targets and their coverage per sample.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        programGroup = CopyNumberProgramGroup.class,
        summary = "Genotype locations for copy number variation in germline samples using a HMM",
        oneLineSummary = "Genotype location for copy number variation"
)
public final class XHMMSegmentGenotyper extends XHMMSegmentCallerBase {

    /**
     * VCF header keys
     */
    public static final String DISCOVERY_KEY = HiddenMarkovModelPostProcessor.DISCOVERY_KEY;
    public static final String NUMBER_OF_TARGETS_KEY = HiddenMarkovModelPostProcessor.NUMBER_OF_POOLED_TARGETS_KEY;
    public static final String SOME_QUALITY_KEY = HiddenMarkovModelPostProcessor.SOME_QUALITY_KEY;
    public static final String START_QUALITY_KEY = HiddenMarkovModelPostProcessor.START_QUALITY_KEY;
    public static final String END_QUALITY_KEY = HiddenMarkovModelPostProcessor.END_QUALITY_KEY;
    public static final String DISCOVERY_TRUE = HiddenMarkovModelPostProcessor.DISCOVERY_TRUE;
    public static final String DISCOVERY_FALSE = HiddenMarkovModelPostProcessor.DISCOVERY_FALSE;

    public static final int MAX_GQ = HiddenMarkovModelPostProcessor.MAX_GQ;
    public static final double MAX_LOG_PROB = HiddenMarkovModelPostProcessor.MAX_LOG_PROB;
    public static final double MAX_QUAL_SCORE = HiddenMarkovModelPostProcessor.MAX_QUAL_SCORE;
    public static final double PHRED_SCORE_PRECISION = HiddenMarkovModelPostProcessor.PHRED_SCORE_PRECISION;

    public static final String DISCOVERY_FILE_SHORT_NAME = "segments";
    public static final String DISCOVERY_FILE_FULL_NAME = "segmentsFile";

    @Argument(
            doc = "Discovered segments file",
            fullName = DISCOVERY_FILE_FULL_NAME,
            shortName = DISCOVERY_FILE_SHORT_NAME
    )
    protected File segmentsFile;

    private VariantContextWriter outputWriter;

    @Override
    protected void openOutput(final File outputFile, final XHMMModel model, final TargetCollection<Target> targets,
                              ReadCountCollection inputCounts) {
        outputWriter = GATKVariantContextUtils.createVCFWriter(outputFile, null, false);
    }

    @Override
    protected void closeOutput(final File outputFile) {
        outputWriter.close();
    }

    @Override
    protected void makeCalls(final XHMMModel model, final TargetCollection<Target> targets,
                             final ReadCountCollection inputCounts) {
        /* perform segmentation and write calls to outputWriter */
        final HiddenMarkovModelPostProcessor<XHMMEmissionData, CopyNumberTriState, Target> processor =
                new HiddenMarkovModelPostProcessor<>(inputCounts.columnNames(),
                        Collections.nCopies(inputCounts.columnNames().size(), targets),
                        sampleForwardBackwardResults, sampleBestPaths,  CopyNumberTriState.NEUTRAL);
        processor.writeVariantsToVCFWriterOnGivenSegments(segmentsFile, CopyNumberTriState::fromCallString,
                outputWriter, "CNV", this.getCommandLine());
    }
}
