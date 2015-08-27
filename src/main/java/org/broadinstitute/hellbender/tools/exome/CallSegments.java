package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ExomeAnalysisProgramGroup;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.io.File;
import java.util.List;

/**
 * Calls segments as amplified, deleted, or copy number neutral given files containing tangent-normalized
 * read counts by target and a list of segments
 *
 * @author David Benjamin
 */
@CommandLineProgramProperties(
        summary = "Call segments as amplified, deleted, or copy number neutral given files containing tangent-normalized" +
                " read counts by target and a list of segments",
        oneLineSummary = "Call segments as amplified, deleted, or copy number neutral",
        programGroup = ExomeAnalysisProgramGroup.class
)
public final class CallSegments extends CommandLineProgram{

    protected static final String SEGFILE_SHORT_NAME = "S";
    protected static final String SEGFILE_LONG_NAME = "segments";

    protected static final String TARGET_FILE_SHORT_NAME = "T";
    protected static final String TARGET_FILE_LONG_NAME = "targets";

    protected static final String OUTPUT_SHORT_NAME = StandardArgumentDefinitions.OUTPUT_SHORT_NAME;
    protected static final String OUTPUT_LONG_NAME = StandardArgumentDefinitions.OUTPUT_LONG_NAME;

    protected static final String Z_THRESHOLD_SHORT_NAME = "Z";
    protected static final String Z_THRESHOLD_LONG_NAME = "threshold";

    protected static final String SAMPLE_LONG_NAME = "sample";

    @Argument(
            doc = "normalized read counts input file.",
            shortName = TARGET_FILE_SHORT_NAME,
            fullName = TARGET_FILE_LONG_NAME,
            optional = false
    )
    protected File targetsFile;

    @Argument(
            doc = "segments files",
            shortName = SEGFILE_SHORT_NAME,
            fullName = SEGFILE_LONG_NAME,
            optional = false
    )
    protected File segmentsFile;

    @Argument(
            doc = "Called segments output",
            shortName = OUTPUT_SHORT_NAME,
            fullName = OUTPUT_LONG_NAME,
            optional = false
    )
    protected File outFile;

    @Argument(
            doc = "Sample",
            fullName = SAMPLE_LONG_NAME,
            optional = false
    )
    protected String sample;

    @Argument(
            doc = "(Advanced) Number of standard deviations of targets' coverage a segment mean must deviate from copy neutral"
            + " to be considered an amplification or deletion.  This parameter controls the trade-off between"
            + " sensitivity and specificity, with smaller values favoring sensitivity.",
            shortName = Z_THRESHOLD_SHORT_NAME,
            fullName = Z_THRESHOLD_LONG_NAME,
            optional = true
    )
    protected double zThreshold = ReCapSegCaller.DEFAULT_Z_SCORE_THRESHOLD;

    @Override
    protected Object doWork() {


        final TargetCollection<TargetCoverage> targets = TargetCoverageUtils.readModeledTargetFileIntoTargetCollection(targetsFile);

        final List<ModeledSegment> segments = SegmentUtils.readModeledSegmentsFromSegfile(segmentsFile);

        //add calls to segments in-place
        ReCapSegCaller.makeCalls(targets, segments, zThreshold);

        SegmentUtils.writeModeledSegmentsToSegfile(outFile, segments, sample);

        return "SUCCESS";
    }
}
