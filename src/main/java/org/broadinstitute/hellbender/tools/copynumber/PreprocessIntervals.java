package org.broadinstitute.hellbender.tools.copynumber;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberArgumentValidationUtils;
import org.broadinstitute.hellbender.utils.*;
import picard.cmdline.programgroups.IntervalsManipulationProgramGroup;

import java.io.File;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Prepares bins for coverage collection.
 *
 * <p>
 *     The input intervals are first checked for overlapping intervals, which are merged.
 *     The resulting intervals are then padded.  The padded intervals are then split into bins.
 *     Finally, bins that contain only Ns are filtered out.
 * </p>
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         Reference FASTA file
 *     </li>
 *     <li>
 *         Intervals to be preprocessed.
 *         The argument {@code interval-merging-rule} must be set to {@link IntervalMergingRule#OVERLAPPING_ONLY}
 *         and all other common arguments for interval padding or merging must be set to their defaults.
 *         If no intervals are specified, then each contig will be assumed to be a single interval and binned accordingly;
 *         this produces bins appropriate for whole genome sequencing analyses.
 *     </li>
 *     <li>
 *         Padding length (in bp).
 *         Use {@code padding} to specify the size of each of the regions added to both ends of the intervals that result
 *         after overlapping intervals have been merged.  Do not use the common {@code interval-padding} argument.
 *         Intervals that would overlap after padding by the specified amount are instead only
 *         padded until they are adjacent.
 *     </li>
 *     <li>
 *         Bin length (in bp).
 *         If this length is not commensurate with the length of a padded interval, then the last bin will be of
 *         different length than the others in that interval.  If zero is specified, then no binning will be performed;
 *         this is generally appropriate for targeted analyses.
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         Preprocessed Picard interval-list file.
 *     </li>
 * </ul>
 *
 * <h3>Usage examples</h3>
 *
 * To pad intervals by 250 bases and disable binning (e.g., for targeted analyses):
 *
 * <pre>
 *     gatk PreprocessIntervals \
 *          -R reference.fa \
 *          -L intervals.interval_list \
 *          --bin-length 0 \
 *          --padding 250 \
 *          -O preprocessed_intervals.interval_list
 * </pre>
 *
 * To generate consecutive bins of 1000 bases from the reference (e.g., for whole genome sequencing analyses):
 *
 * <pre>
 *     gatk PreprocessIntervals \
 *          -R reference.fa \
 *          --bin-length 1000 \
 *          --padding 0 \
 *          -O preprocessed_intervals.interval_list
 * </pre>
 *
 * @author Marton Kanasz-Nagy &lt;mkanaszn@broadinstitute.org&gt;
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Prepares bins for coverage collection",
        oneLineSummary = "Prepares bins for coverage collection",
        programGroup = IntervalsManipulationProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public final class PreprocessIntervals extends GATKTool {
    public static final String BIN_LENGTH_LONG_NAME = "bin-length";
    public static final String PADDING_LONG_NAME = "padding";

    @Argument(
            doc = "Length (in bp) of the bins.  If zero, no binning will be performed.",
            fullName = BIN_LENGTH_LONG_NAME,
            optional = true,
            minValue = 0
    )
    private int binLength = 1000;

    @Argument(
            doc = "Length (in bp) of the padding regions on each side of the intervals.",
            fullName = PADDING_LONG_NAME,
            optional = true,
            minValue = 0
    )
    private int padding = 250;

    @Argument(
            doc = "Output Picard interval-list file containing the preprocessed intervals.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private File outputFile;

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public void onTraversalStart() {
        final SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();

        final List<SimpleInterval> inputIntervals;
        if (hasUserSuppliedIntervals()) {
            CopyNumberArgumentValidationUtils.validateIntervalArgumentCollection(intervalArgumentCollection);
            inputIntervals = intervalArgumentCollection.getIntervals(sequenceDictionary);
        } else {
            // if the user didn't add any intervals, we assume that they wanted to do whole genome sequencing
            inputIntervals = IntervalUtils.getAllIntervalsForReference(sequenceDictionary);
        }

        logger.info("Padding intervals...");
        final IntervalList paddedIntervalList = padIntervals(inputIntervals, padding, sequenceDictionary);

        logger.info("Generating bins...");
        final IntervalList unfilteredBins = generateBins(paddedIntervalList, binLength, sequenceDictionary);

        logger.info("Filtering bins containing only Ns...");
        final ReferenceDataSource reference = ReferenceDataSource.of(referenceArguments.getReferencePath());
        final IntervalList bins = filterBinsContainingOnlyNs(unfilteredBins, reference);

        logger.info(String.format("Writing bins to %s...", outputFile));
        bins.write(outputFile);
    }

    private static IntervalList padIntervals(final List<SimpleInterval> inputIntervals, final int padding, final SAMSequenceDictionary sequenceDictionary) {
        final List<SimpleInterval> paddedIntervals = inputIntervals.stream()
                .map(i -> new SimpleInterval(
                        i.getContig(), 
                        Math.max(1, i.getStart() - padding), 
                        Math.min(i.getEnd() + padding, sequenceDictionary.getSequence(i.getContig()).getSequenceLength())))
                .collect(Collectors.toList());

        // alter the padded intervals in place to eliminate overlaps
        for (int i = 0; i < paddedIntervals.size() - 1; i++) {
            final SimpleInterval thisInterval = paddedIntervals.get(i);
            final SimpleInterval nextInterval = paddedIntervals.get(i + 1);
            if (thisInterval.overlaps(nextInterval)) {
                final int originalThisEnd = inputIntervals.get(i).getEnd();
                final int originalNextStart = inputIntervals.get(i + 1).getStart();

                final int newThisEnd = (originalThisEnd + originalNextStart) / 2;
                final int newNextStart = newThisEnd + 1;

                paddedIntervals.set(i, new SimpleInterval(thisInterval.getContig(), thisInterval.getStart(), newThisEnd));
                paddedIntervals.set(i + 1, new SimpleInterval(nextInterval.getContig(), newNextStart, nextInterval.getEnd()));
            }
        }

        final IntervalList paddedIntervalList = new IntervalList(sequenceDictionary);
        paddedIntervals.forEach(i -> paddedIntervalList.add(new Interval(i.getContig(), i.getStart(), i.getEnd())));
        return paddedIntervalList;
    }

    private static IntervalList generateBins(final IntervalList preparedIntervalList, final int binLength, final SAMSequenceDictionary sequenceDictionary) {
        if (binLength == 0) {
            return IntervalList.copyOf(preparedIntervalList);
        }
        final IntervalList bins = new IntervalList(sequenceDictionary);
        for (final Interval interval : preparedIntervalList) {
            for (int binStart = interval.getStart(); binStart <= interval.getEnd(); binStart += binLength) {
                final int binEnd = FastMath.min(binStart + binLength - 1, interval.getEnd());
                bins.add(new Interval(interval.getContig(), binStart, binEnd));
            }
        }
        return bins;
    }

    private static IntervalList filterBinsContainingOnlyNs(final IntervalList unfilteredBins, final ReferenceDataSource reference) {
        final IntervalList bins = new IntervalList(reference.getSequenceDictionary());
        for (final Interval unfilteredBin : unfilteredBins) {
            if (!Utils.stream(reference.query(new SimpleInterval(unfilteredBin))).allMatch(b -> Nucleotide.decode(b) == Nucleotide.N)) {
                bins.add(unfilteredBin);
            }
        }
        return bins;
    }

    @Override
    public void traverse() {}  // no traversal for this tool
}