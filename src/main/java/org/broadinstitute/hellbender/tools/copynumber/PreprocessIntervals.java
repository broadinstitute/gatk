package org.broadinstitute.hellbender.tools.copynumber;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.io.File;
import java.util.List;

/**
 * This tool takes in intervals via the standard arguments of {@link GATKTool}, prepares them for binning, and creates
 * bins that cover the processed intervals. The intervals are first padded, uniquified, then the overlapping
 * intervals are merged. Finally, the bins are created.
 *
 * <p>Standard GATK engine arguments include -L and -XL. For example, for the -L argument, the tool accepts GATK-style
 * intervals (.list or .intervals), BED files and VCF files. If no intervals are given, each contig will be assumed
 * to be a single interval (whole genome sequencing). </p>
 *
 * <p>Using the -P flag, the user can specify the amount of padding (in bp) added to each side of the intervals.
 * This padding is in addition to the padding added by the -ip (or --interval_padding) argument of
 * IntervalArgumentCollection. However, we encourage using only the -P flag.
 *
 * <p>The user can also specify the length of the bins (in bp) using the -BL option. If this is not commensurate with
 * the length of the padded intervals, then the last bin will be of different length than the others. </p>
 *
 * <p> The -O argument specifies a filename for the output bins, stored as a Picard interval list. </p>
 *
 * <h3>Example</h3>
 * <pre>
 * gatk-launch --javaOptions "-Xmx4g" PreprocessIntervals \
 *   -R ref_fasta.fa \
 *   -L intervals.list \
 *   -BL 10000 \
 *   -P 500 \
 *   -O preprocessed-intervals.interval_list
 * </pre>
 *
 * @author Marton Kanasz-Nagy &lt;mkanaszn@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Prepare intervals for binning, then create bins that cover them. "
                + "The intervals are first padded, the overlapping ones are merged "
                + "and bins are created that cover the intervals. "
                + "If the bin length is incommensurate with the length of an interval, "
                + "the last bin of that interval will be of a different length. "
                + "The length of the padding regions at both sides of the intervals (-P) and "
                + "the length of the bins (-BL) can be specified. ",
        oneLineSummary = "Prepare intervals for binning then create bins that cover them.",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
public final class PreprocessIntervals extends GATKTool {
    public static final String BIN_LENGTH_LONG_NAME = "binLength";
    public static final String BIN_LENGTH_SHORT_NAME = "BL";

    public static final String PADDING_LONG_NAME = "padding";
    public static final String PADDING_SHORT_NAME = "P";

    @Argument(
            doc = "Length (in bp) of the bins.",
            fullName = BIN_LENGTH_LONG_NAME,
            shortName = BIN_LENGTH_SHORT_NAME,
            optional = true,
            minValue = 1
    )
    private int binLength = 1000;

    @Argument(
            doc = "Length (in bp) of the padding regions on each side of the intervals.",
            fullName = PADDING_LONG_NAME,
            shortName = PADDING_SHORT_NAME,
            optional = true,
            minValue = 0
    )
    private int padding = 0;

    @Argument(
            doc = "Output Picard interval-list file containing the preprocessed intervals.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private File outputFile;

    @Override
    public void onTraversalStart() {
        final SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();

        // if the user didn't add any intervals, we assume that they wanted to do whole genome sequencing
        final List<SimpleInterval> inputIntervals = hasIntervals() ? intervalArgumentCollection.getIntervals(sequenceDictionary)
                : IntervalUtils.getAllIntervalsForReference(sequenceDictionary);

        logger.info("Padding and merging intervals...");
        final IntervalList preparedIntervalList = padAndMergeIntervals(inputIntervals, padding, sequenceDictionary);

        logger.info("Generating bins...");
        final IntervalList bins = generateBins(preparedIntervalList, binLength, sequenceDictionary);

        logger.info(String.format("Writing bins to %s...", outputFile));
        bins.write(outputFile);
    }

    private static IntervalList padAndMergeIntervals(final List<SimpleInterval> inputIntervals, final int padding, final SAMSequenceDictionary sequenceDictionary) {
        final IntervalList inputIntervalList = new IntervalList(sequenceDictionary);
        inputIntervals.stream().map(si -> new Interval(si.getContig(), si.getStart(), si.getEnd())).forEach(inputIntervalList::add);

        final IntervalList uniquedIntervalList = inputIntervalList.uniqued();
        final IntervalList paddedIntervalList = uniquedIntervalList.padded(padding, padding);
        final IntervalList mergedIntervalList = IntervalList.intersection(paddedIntervalList, paddedIntervalList);

        return mergedIntervalList;
    }

    private static IntervalList generateBins(final IntervalList preparedIntervalList, final int binLength, final SAMSequenceDictionary sequenceDictionary) {
        final IntervalList bins = new IntervalList(sequenceDictionary);
        for (final Interval interval : preparedIntervalList) {
            for (int binStart = interval.getStart(); binStart <= interval.getEnd(); binStart += binLength) {
                final int binEnd = FastMath.min(binStart + binLength - 1, interval.getEnd());
                bins.add(new Interval(interval.getContig(), binStart, binEnd));
            }
        }
        return bins;
    }

    @Override
    public void traverse() {}  // no traversal for this tool
}