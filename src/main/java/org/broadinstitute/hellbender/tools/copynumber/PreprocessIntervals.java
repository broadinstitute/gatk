package org.broadinstitute.hellbender.tools.copynumber;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberArgumentValidationUtils;
import org.broadinstitute.hellbender.utils.*;
import picard.cmdline.programgroups.IntervalsManipulationProgramGroup;

import java.io.File;
import java.util.Comparator;
import java.util.List;
import java.util.Spliterator;
import java.util.function.Consumer;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

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
 *     <li>
 *         Minimum bin length (in bp).
 *         With {@code min-bin-length} you can specify the minimum size for any given bin in case that these need to be
 *         truncated to accommodate the size of the enclosing interval or contig. Bin shorter than that will be excluded
 *         from the output.
 *         By default this argument is set to 1, so that any non-empty bin will be emitted.
 *     </li>
 *     <li>
 *         Gridded output bins.
 *         Using the {@code grid} flag you can request that the bins start at specific position with in the interval regardless
 *         of the start of the interval itself. These "grid points" position would be {@code n * bin-length + 1} where {@code n}
 *         is any integer equal or greater than 0. As a consequence the first and last bin would be truncated or
 *         skipped (see {@code min-bin-length}). By default bins are not gridded.
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
public final class PreprocessIntervals extends GATKTool {
    public static final String BIN_LENGTH_LONG_NAME = "bin-length";
    public static final String PADDING_LONG_NAME = "padding";
    public static final String GRID_LONG_NAME = "grid";
    public static final String MINIMUM_BIN_LENGTH_LONG_NAME = "min-bin-length";

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

    @Argument(
            doc = "Intervals are layout in a grid so that bin start at fixed positions (grid points) 1 + n * bin_length. When intervals " +
                    "do not exactly align with a grid starting position, the first bin in that interval will be smaller or simply skipped (see min-bin-length) " +
                    "so that the rest of the  bins will start at a grid point",
            fullName = GRID_LONG_NAME,
            optional = true
    )
    private boolean grid = false;

    @Argument(
            doc = "The minimum bin length. Any truncated bin that are smaller than this length will be excluded, " +
                    "when set to 1 means that no  bin will be excluded",
            fullName = MINIMUM_BIN_LENGTH_LONG_NAME,
            minValue = 1,
            optional = true
    )
    private int minimumBinLength = 1;

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public void onTraversalStart() {
        final SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();

        if (binLength > 0 && minimumBinLength > binLength) {
            throw new UserException.BadInput(
                    String.format("%s (%s) must be less or equal than %s (%s)", MINIMUM_BIN_LENGTH_LONG_NAME,
                            minimumBinLength, BIN_LENGTH_LONG_NAME, binLength));
        } else if (grid && binLength <= 0) {
            throw new UserException.BadInput("when requested gridded bin you must speficy a bin-length larger than 0");
        }

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
        final SAMFileHeader header = paddedIntervalList.getHeader();
        final Stream<Interval> unfilteredBins = generateBinStream(paddedIntervalList, binLength, grid, minimumBinLength);

        logger.info("Filtering bins containing only Ns...");
        final ReferenceDataSource reference = ReferenceDataSource.of(referenceArguments.getReferencePath());
        final Stream<Interval> bins = filterBinsContainingOnlyNs(unfilteredBins, reference);

        logger.info(String.format("Writing bins to %s...", outputFile));
        IntervalUtils.write(header, bins, outputFile.toString());
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

    private static Stream<Interval> generateBinStream(final IntervalList preparedIntervalList, final int binLength, final boolean grid,
                                             final int minimumBinLength) {
        if (binLength == 0) {
            return Utils.stream(preparedIntervalList);
        }
        return Utils.stream(preparedIntervalList)
                .flatMap(interval ->
                        StreamSupport.stream(
                                new IntervalSpliterator(interval, binLength, minimumBinLength, grid), false));
    }

    private static Stream<Interval> filterBinsContainingOnlyNs(final Stream<Interval> unfilteredBins, final ReferenceDataSource reference) {
        return unfilteredBins.filter(interval -> {
            final byte[] bases = reference.queryAndPrefetch(interval.getContig(), interval.getStart(), interval.getEnd()).getBases();
            for (byte b : bases) {
                if (Nucleotide.decode(b) != Nucleotide.N) {
                    return true;
                }
            }
            return false;
        });
    }

    @Override
    public void traverse() {}  // no traversal for this tool

    private static class IntervalSpliterator implements Spliterator<Interval> {

        private int binStart;
        private final Interval interval;
        private final int binLength;
        private final boolean grid;
        private final int minBinLength;

        IntervalSpliterator(final Interval interval,
                            final int binLength,
                            final int minBinLength,
                            final boolean grid) {
            this.binStart = interval.getStart();
            this.interval = interval;
            this.binLength = binLength;
            this.minBinLength = minBinLength;
            this.grid = grid;
        }

        @Override
        public boolean tryAdvance(final Consumer<? super Interval> action) {
             while (true) {
                 int gridOffset;
                 // no enough bases left for a minimum bin?
                 if (binStart > interval.getEnd() - minBinLength + 1) {
                    return false;
                    // no gridded bins
                 } else if (!grid || (gridOffset = binStart % binLength) == 1) {
                    int binEnd = FastMath.min(binStart + binLength - 1, interval.getEnd());
                    action.accept(new Interval(interval.getContig(), binStart, binEnd));
                    binStart = binEnd + 1;
                    return true;
                    // gridded bins:
                 } else { // grid == true && not on grid.
                    // just 1bp before the next grid start!
                    if (gridOffset == 0) {
                        if (minBinLength <= 1) {
                            action.accept(new Interval(interval.getContig(), binStart, binStart++));
                            return true;
                        } else {
                            binStart++;
                            // continue:
                        }
                    } else {
                        int binEnd = FastMath.min(binStart + binLength - gridOffset, interval.getEnd());
                        if (minBinLength <= binEnd - binStart + 1) { // gridOffset > 1
                            action.accept(new Interval(interval.getContig(), binStart, binEnd));
                            binStart = binEnd + 1;
                            return true;
                        } else { // we skip to next grid start and we try again:
                            binStart = binEnd + 1;
                            // continue;
                        }
                    }
                }
            }
        }

        @Override
        public Spliterator<Interval> trySplit() {
            return null;
        }

        @Override
        public long estimateSize() {
            return (interval.getLengthOnReference() / binLength) + 1;
        }

        @Override
        public int characteristics() {
            return Spliterator.NONNULL | Spliterator.DISTINCT | Spliterator.SORTED;
        }

        @Override
        public Comparator<Interval> getComparator() {
            return Comparator.comparingInt(Interval::getStart).thenComparing(Interval::getEnd);
        }

    }
}