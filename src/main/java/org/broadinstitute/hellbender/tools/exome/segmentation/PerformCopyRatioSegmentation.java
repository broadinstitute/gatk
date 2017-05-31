package org.broadinstitute.hellbender.tools.exome.segmentation;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.ModeledSegment;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollectionUtils;
import org.broadinstitute.hellbender.tools.exome.SegmentUtils;

import java.io.File;
import java.io.IOException;
import java.util.List;

/**
 * Groups contiguous regions of constant copy ratio per sample.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Segment genomic data into regions of constant copy ratio.  Only supports one sample input.",
        oneLineSummary = "(Experimental) Segment genomic data into regions of constant copy ratio",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
public final class PerformCopyRatioSegmentation extends CommandLineProgram {
    protected static final String INITIAL_NUM_STATES_LONG_NAME = "initialNumberOfStates";
    protected static final String INITIAL_NUM_STATES_SHORT_NAME = "initialNumStates";

    @Argument(
            doc = "Tangent-normalized log2 read counts file",
            fullName = ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME,
            optional = false
    )
    protected String coverageFile;

    @Argument(
            doc = "Initial number of hidden copy-ratio states",
            fullName = INITIAL_NUM_STATES_LONG_NAME,
            shortName = INITIAL_NUM_STATES_SHORT_NAME,
            optional = false
    )
    protected int initialNumStates;

    @Argument(
            doc = "Output file for copy-ratio segments.  Output is not logged.",
            fullName = ExomeStandardArgumentDefinitions.SEGMENT_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME,
            optional = false
    )
    protected File outputSegmentsFile;

    @Override
    public Object doWork() {
        final String sampleName = ReadCountCollectionUtils.getSampleNameForCLIsFromReadCountsFile(new File(coverageFile));
        final ReadCountCollection rcc;
        try {
            rcc = ReadCountCollectionUtils.parse(new File(coverageFile));
        } catch (final IOException ex) {
            throw new UserException.BadInput("could not read input file");
        }

        final CopyRatioSegmenter segmenter = new CopyRatioSegmenter(initialNumStates, rcc);
        final List<ModeledSegment> segments = segmenter.getModeledSegments();
        SegmentUtils.writeModeledSegmentFile(outputSegmentsFile, segments, sampleName, false);

        return "SUCCESS";
    }
}
