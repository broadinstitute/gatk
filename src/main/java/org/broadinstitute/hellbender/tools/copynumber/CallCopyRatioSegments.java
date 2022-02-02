package org.broadinstitute.hellbender.tools.copynumber;

import com.google.common.annotations.VisibleForTesting;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberArgumentValidationUtils;
import org.broadinstitute.hellbender.tools.copynumber.caller.SimpleCopyRatioCaller;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CalledCopyRatioSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CalledLegacySegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CopyRatioSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CalledCopyRatioSegment;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CalledLegacySegment;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyRatioSegment;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Calls copy-ratio segments as amplified, deleted, or copy-number neutral.
 *
 * <p>
 *     This is a relatively naive caller that takes the modeled-segments output of {@link ModelSegments} and
 *     performs a simple statistical test on the segmented log2 copy ratios to call amplifications and deletions,
 *     given a specified range for determining copy-number neutral segments.  This caller is
 *     based on the calling functionality of ReCapSeg.
 * </p>
 *
 * <p>
 *     If provided {@link ModelSegments} results that incorporate allele-fraction data, i.e. data with presumably improved segmentation,
 *     the statistical test performed by CallCopyRatioSegments ignores the modeled minor-allele fractions when making calls.
  * </p>
 *
 * <h3>Input</h3>
 *
 * <ul>
 *     <li>
 *         Copy-ratio-segments .cr.seg file from {@link ModelSegments}.
 *         This is a tab-separated values (TSV) file with a SAM-style header containing a read group sample name, a sequence dictionary,
 *         a row specifying the column headers contained in {@link CopyRatioSegmentCollection.CopyRatioSegmentTableColumn},
 *         and the corresponding entry rows.
 *     </li>
 * </ul>
 *
 * <h3>Outputs</h3>
 *
 * <ul>
 *     <li>
 *         Called copy-ratio-segments file.
 *         This is a tab-separated values (TSV) file with a SAM-style header containing a read group sample name, a sequence dictionary,
 *         a row specifying the column headers contained in {@link CalledCopyRatioSegmentCollection.CalledCopyRatioSegmentTableColumn},
 *         and the corresponding entry rows.
 *     </li>
 *     <li>
 *         CBS-formatted .igv.seg file (compatible with IGV).
 *         This is a tab-separated values (TSV) file with CBS-format column headers
 *         (see <a href="http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#CBS">
 *             http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#CBS</a>)
 *         and the corresponding entry rows that can be plotted using IGV (see
 *         <a href="https://software.broadinstitute.org/software/igv/SEG">
 *             https://software.broadinstitute.org/software/igv/SEG</a>).
 *         The mean log2 copy ratio is given in the SEGMENT_MEAN column.
 *     </li>
 * </ul>
 *
 * <h3>Usage examples</h3>
 *
 * <pre>
 *     gatk CallCopyRatioSegments \
 *          -I tumor.cr.seg \
 *          -O tumor.called.seg
 * </pre>
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Calls copy-ratio segments as amplified, deleted, or copy-number neutral",
        oneLineSummary = "Calls copy-ratio segments as amplified, deleted, or copy-number neutral",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
public final class CallCopyRatioSegments extends CommandLineProgram {
    public static final String NEUTRAL_SEGMENT_COPY_RATIO_LOWER_BOUND_LONG_NAME = "neutral-segment-copy-ratio-lower-bound";
    public static final String NEUTRAL_SEGMENT_COPY_RATIO_UPPER_BOUND_LONG_NAME = "neutral-segment-copy-ratio-upper-bound";
    public static final String OUTLIER_NEUTRAL_SEGMENT_COPY_RATIO_Z_SCORE_THRESHOLD_LONG_NAME = "outlier-neutral-segment-copy-ratio-z-score-threshold";
    public static final String CALLING_COPY_RATIO_Z_SCORE_THRESHOLD_LONG_NAME = "calling-copy-ratio-z-score-threshold";
    public static final String LEGACY_SEGMENTS_FILE_SUFFIX = ".igv.seg";

    @Argument(
            doc = "Input file containing copy-ratio segments (.cr.seg output of ModelSegments).",
            fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME
    )
    private File inputCopyRatioSegmentsFile;

    @Argument(
            doc = "Output file for called copy-ratio segments.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private File outputCalledCopyRatioSegmentsFile;

    @Argument(
            doc = "Lower bound on non-log2 copy ratio used for determining copy-neutral segments.",
            fullName = NEUTRAL_SEGMENT_COPY_RATIO_LOWER_BOUND_LONG_NAME,
            optional = true,
            minValue = 0.
    )
    private double neutralSegmentCopyRatioLowerBound = 0.9;

    @Argument(
            doc = "Upper bound on non-log2 copy ratio used for determining copy-neutral segments.",
            fullName = NEUTRAL_SEGMENT_COPY_RATIO_UPPER_BOUND_LONG_NAME,
            optional = true,
            minValue = 0.
    )
    private double neutralSegmentCopyRatioUpperBound = 1.1;

    @Argument(
            doc = "Threshold on z-score of non-log2 copy ratio used for determining outlier copy-neutral segments.  " +
                    "If non-log2 copy ratio z-score is above this threshold for a copy-neutral segment, " +
                    "it is considered an outlier and not used in the calculation of the length-weighted mean and standard deviation " +
                    "used for calling.",
            fullName = OUTLIER_NEUTRAL_SEGMENT_COPY_RATIO_Z_SCORE_THRESHOLD_LONG_NAME,
            optional = true,
            minValue = 0.
    )
    private double outlierNeutralSegmentCopyRatioZScoreThreshold = 2.;

    @Argument(
            doc = "Threshold on z-score of non-log2 copy ratio used for calling segments.",
            fullName = CALLING_COPY_RATIO_Z_SCORE_THRESHOLD_LONG_NAME,
            optional = true,
            minValue = 0.
    )
    private double callingCopyRatioZScoreThreshold = 2.;

    private File outputCalledLegacySegmentsFile;

    @Override
    protected Object doWork() {
        validateArguments();

        final CopyRatioSegmentCollection copyRatioSegments = new CopyRatioSegmentCollection(inputCopyRatioSegmentsFile);

        // BEGIN sato
        final double[] logCopyRatios = copyRatioSegments.getRecords().stream()
                .mapToDouble(r -> r.getMeanLog2CopyRatio()).toArray();

        // Not the overall medians, by the median of segment means.
        double medianCopyRatio = Math.pow(2, new Median().evaluate(logCopyRatios));
        double delta = 0.1;
        neutralSegmentCopyRatioLowerBound = medianCopyRatio - delta;
        neutralSegmentCopyRatioUpperBound = medianCopyRatio + delta;
        // END sato

        final CalledCopyRatioSegmentCollection calledCopyRatioSegments =
                new SimpleCopyRatioCaller(copyRatioSegments,
                        neutralSegmentCopyRatioLowerBound, neutralSegmentCopyRatioUpperBound,
                        outlierNeutralSegmentCopyRatioZScoreThreshold, callingCopyRatioZScoreThreshold)
                        .makeCalls();

        logger.info(String.format("Writing called segments to %s...", outputCalledCopyRatioSegmentsFile.getAbsolutePath()));
        calledCopyRatioSegments.write(outputCalledCopyRatioSegmentsFile);

        // Write an IGV compatible collection
        final CalledLegacySegmentCollection legacySegmentCollection = createCalledLegacySegmentCollection(calledCopyRatioSegments);
        logger.info(String.format("Writing called segments in IGV-compatible format to %s...", outputCalledLegacySegmentsFile.getAbsolutePath()));
        legacySegmentCollection.write(outputCalledLegacySegmentsFile);

        logger.info(String.format("%s complete.", getClass().getSimpleName()));

        return null;
    }

    private void validateArguments() {
        Utils.validateArg(neutralSegmentCopyRatioLowerBound < neutralSegmentCopyRatioUpperBound,
                String.format("Copy-neutral lower bound (%s) must be less than upper bound (%s).",
                        NEUTRAL_SEGMENT_COPY_RATIO_LOWER_BOUND_LONG_NAME,
                        NEUTRAL_SEGMENT_COPY_RATIO_UPPER_BOUND_LONG_NAME));

        CopyNumberArgumentValidationUtils.validateInputs(inputCopyRatioSegmentsFile);
        outputCalledLegacySegmentsFile = createCalledLegacySegmentsFile(outputCalledCopyRatioSegmentsFile);
        CopyNumberArgumentValidationUtils.validateOutputFiles(
                outputCalledCopyRatioSegmentsFile,
                outputCalledLegacySegmentsFile);
    }

    @VisibleForTesting
    static File createCalledLegacySegmentsFile(final File calledCopyRatioSegmentsFile) {
        return new File(FilenameUtils.removeExtension(calledCopyRatioSegmentsFile.getAbsolutePath()) + LEGACY_SEGMENTS_FILE_SUFFIX);
    }

    private static CalledLegacySegmentCollection createCalledLegacySegmentCollection(final CalledCopyRatioSegmentCollection segments) {
        return new CalledLegacySegmentCollection(segments.getMetadata(), segments.getRecords().stream()
                .map(r -> convert(r, segments.getMetadata().getSampleName())).collect(Collectors.toList()));
    }

    private static CalledLegacySegment convert(final CalledCopyRatioSegment segment, final String sampleName) {
        return new CalledLegacySegment(sampleName, segment.getInterval(), segment.getNumPoints(), segment.getMeanLog2CopyRatio(), segment.getCall());
    }

}
