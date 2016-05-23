package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMReadGroupRecord;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.TargetWalker;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.function.BiFunction;
import java.util.function.IntBinaryOperator;
import java.util.function.ToIntFunction;
import java.util.stream.*;

/**
 * Calculates the coverage for each target at each input sample.
 * <h4>Syntax example:</h4>
 * <pre>
 *     $ java -jar hellbender.jar CalculateTargetCoverage -targets mytargets.tab \
 *          -I sample1.bam -I sample2.bam -unit AVERAGE_DEPTH \
 *          -max 5000 -minBQ 20 -minMQ 30 -o coverage.tab
 * </pre>
 *
 * </h4>
 * <h4>Input and Output</h4>
 * <p>
 * The user must indicate the target table file using the argument {@link TargetArgumentCollection#targetsFile normalizedCoveragedFile}.
 * That file must follow the format described in {@link TargetTableReader}.
 * </p>
 * <p>
 * Alignment files must also be provided using {@link StandardArgumentDefinitions#INPUT_LONG_NAME} argument.
 * </p>
 *
 * <p>
 *     The output follows the same format as the input target file including an additional column per each sample
 *     with its coverage.
 * </p>
 *
 * <h4>Qualifying reads</h4>
 * <p>Only evidence contained in qualifying reads will be taken into account to
 * calculate coverage.</p>
 * <p>
 *     For a read to qualify it must match the following criteria:
 *     <ul>
 *         <li>must be <i>well-formed</i> as defined in {@link WellformedReadFilter},</li>
 *         <li>must be mapped (unmapped reads with mapped mates are not taken into account),</li>
 *         <li>must align with at least one base in the reference (all insertion reads are discarded) and</li>
 *         <li>must not be marked as duplicated.</li>
 *     </ul>
 * </p>
 *
 * <p>In addition, reads that have a poor mapping quality will be discarded.
 *    The user can indicate the minimum mapping quality with the {@link #minimumMappingQuality} argument.
 *    By default this limit is set to {@value #MINIMUM_MAPPING_QUALITY_DEFAULT}.</p>
 *
 * <h4>Supported coverage units</h4>
 * <p>The user has a choice of two different way to count coverage
 * through the {@link #coverageUnit} argument:</p>
 * <dl>
 *     <dt>{@link CoverageUnit#OVERLAPPING_READ OVERLAPPING_READ}</dt>
 *     <dd>Number of reads that overlap the target at least by one base. <i>This is the default.</i></dd>
 *     <dt>{@link CoverageUnit#AVERAGE_READ_DEPTH AVERAGE_DEPTH}</dt>
 *     <dd>Average number of base-calls per base-pair in the target.</dd>
 * </dl>
 *
 * <h4>Qualifying base-calls</h4>
 * <p>When the user selects {@link CoverageUnit#AVERAGE_READ_DEPTH} as the
 * coverage unit, only qualifying base-calls will be taken into account.</p>
 * A base is qualifying iff:
 * <ul>
 *     <li>is enclosed in a qualifying read,</li>
 *     <li>aligned with a base in the reference (i.e. not part for an insertion nor a clip) and</li>
 *     <li>its not a poor quality call.</li>
 * </ul>
 * <p>The user can indicate the minimum base call quality to qualify
 * using the {@link #minimumBaseQuality} argument. By default this argument is
 * set to {@value #MINIMUM_BASE_QUALITY_DEFAULT}.</p>
 * <p>Notice that overlapping deletion "bases" in a qualifying read
 * do not count towards coverage even if they would contribute to the depth
 * of the pile-up on that site in the reference.</p>
 *
 * <h4>Maximum coverage value</h4>
 *
 * <p>The user can indicate a maximum coverage reportable value using the
 * {@link #maximumCoverage} argument. By default the reported coverage
 * is unbound.
 * </p>
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        oneLineSummary = "Calculate Target Coverage",
        summary = "Calculate Target Coverage",
        programGroup = CopyNumberProgramGroup.class
)
public class CalculateTargetBaseCallCoverage extends TargetWalker {

    /**
     * Short name for the {@link #maximumCoverage} argument.
     */
    public static final String MAXIMUM_COVERAGE_SHORT_NAME = "max";

    /**
     * Long name for the {@link #maximumCoverage} argument.
     */
    public static final String MAXIMUM_COVERAGE_FULL_NAME = "maximumCoverage";

    /**
     * Default value for the {@link #maximumCoverage} argument.
     */
    public static final double MAXIMUM_COVERAGE_DEFAULT = Double.POSITIVE_INFINITY;

    /**
     * Short name for the {@link #minimumMappingQuality} argument.
     */
    public static final String MINIMUM_MAPPING_QUALITY_SHORT_NAME = "minMQ";

    /**
     * Long name for the {@link #minimumMappingQuality} argument.
     */
    public static final String MINIMUM_MAPPING_QUALITY_FULL_NAME = "minimumMappingQuality";

    /**
     * Short name for the {@link #minimumBaseQuality} argument.
     */
    public static final String MINIMUM_BASE_QUALITY_SHORT_NAME = "minBQ";

    /**
     * Long name for the {@link #minimumBaseQuality} argument.
     */
    public static final String MINIMUM_BASE_QUALITY_FULL_NAME = "minimumBaseQuality";

    /**
     * Default value for the {@link #minimumMappingQuality} argument.
     */
    public static final int MINIMUM_MAPPING_QUALITY_DEFAULT = 0;

    /**
     * Default value for the {@link #minimumBaseQuality} argument.
     */
    public static final int MINIMUM_BASE_QUALITY_DEFAULT = 0;

    /**
     * Short name for the {@link #coverageUnit} argument.
     */
    public static final String COVERAGE_UNIT_SHORT_NAME = "unit";

    /**
     * Long name for the {@link #coverageUnit} argument.
     */
    public static final String COVERAGE_UNIT_FULL_NAME = "coverageUnit";

    /**
     * Default value for the {@link #coverageUnit} argument.
     */
    public static final CoverageUnit DEFAULT_COVERAGE_UNIT = CoverageUnit.AVERAGE_READ_DEPTH;

    @Argument(
            doc = "Output file",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName  = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional  = false
    )
    protected File outputFile;

    @Argument(
            doc = "Maximum coverage per target and coverage group",
            shortName = MAXIMUM_COVERAGE_SHORT_NAME,
            fullName = MAXIMUM_COVERAGE_FULL_NAME,
            optional = true
    )
    protected double maximumCoverage = MAXIMUM_COVERAGE_DEFAULT;

    @Argument(
            doc = "Minimum mapping quality",
            shortName = MINIMUM_MAPPING_QUALITY_SHORT_NAME,
            fullName  = MINIMUM_MAPPING_QUALITY_FULL_NAME,
            optional = true
    )
    protected int minimumMappingQuality = MINIMUM_MAPPING_QUALITY_DEFAULT;

    @Argument(
            doc = "Coverage unit",
            shortName = COVERAGE_UNIT_SHORT_NAME,
            fullName = COVERAGE_UNIT_FULL_NAME,
            optional = true
    )
    protected CoverageUnit coverageUnit = DEFAULT_COVERAGE_UNIT;

    @Argument(
            doc = "Minimum base quality; base calls with lower quality won't be taken into account",
            shortName = MINIMUM_BASE_QUALITY_SHORT_NAME,
            fullName = MINIMUM_BASE_QUALITY_FULL_NAME,
            optional = true
    )
    protected int minimumBaseQuality = MINIMUM_BASE_QUALITY_DEFAULT;

    protected TableWriter<ReadCountRecord> outputTableWriter;

    protected ToIntFunction<GATKRead> readToColumn;

    protected CountingReadFilter readFilter;

    protected int countColumnCount = 0;

    /**
     * Reference to the counter code.
     */
    private BiFunction<Target, Stream<GATKRead>, long[]> counter;

    @Override
    public boolean requiresReads() {
        return true;
    }

    @Override
    public void onTraversalStart() {

        checkThresholdArgumentValues();
        final List<String> sampleList = getHeaderForReads().getReadGroups().stream()
                .map(SAMReadGroupRecord::getSample)
                .filter(Objects::nonNull)
                .sorted()
                .distinct()
                .collect(Collectors.toList());
        final Map<String, Integer> readGroupToIndex = getHeaderForReads().getReadGroups().stream()
                .filter(rg -> rg.getSample() != null)
                .collect(Collectors.toMap(SAMReadGroupRecord::getId, rg -> sampleList.indexOf(rg.getSample())));
        counter = composeCoverageCounter();
        readToColumn = (read) -> readGroupToIndex.getOrDefault(read.getReadGroup(), -1);
        countColumnCount = sampleList.size();
        readFilter = makeReadFilter();

        try  {
            final Writer outputWriter = createOutputWriter(outputFile);
            outputTableWriter = ReadCountCollectionUtils.writerWithIntervals(outputWriter, sampleList);
        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, ex);
        }
    }

    /**
     * Checks that the values of threshold arguments are valid.
     *
     * <p>
     *     In addition it warns the user if he sets a value for the {@link #minimumBaseQuality}
     *     argument yet base-call qualities are not relevant given the {@link #coverageUnit}
     *     selected.
     * </p>
     *
     * @throws UserException.BadArgumentValue if any of the argument has a bad value.
     */
    private void checkThresholdArgumentValues() {

        if (minimumBaseQuality < 0) {
            throw new UserException.BadArgumentValue(MINIMUM_BASE_QUALITY_FULL_NAME,
                    String.valueOf(minimumBaseQuality),
                    "Its value must be 0 or greater");
        } else if (minimumMappingQuality < 0) {
            throw new UserException.BadArgumentValue(MINIMUM_MAPPING_QUALITY_FULL_NAME,
                    String.valueOf(minimumMappingQuality),
                    "Its value most be 0 or greater");
        } else if (maximumCoverage < 0 || Double.isNaN(maximumCoverage)) {
            throw new UserException.BadArgumentValue(MAXIMUM_COVERAGE_FULL_NAME,
                    String.valueOf(maximumCoverage),
                    "Its value must be 0 or greater.");
        }

        // Warn the user if he specifies a value for the minBQ but
        // the Coverage reported is not based on base-calls.
        if (coverageUnit.datum == CoverageDatum.BASE_CALL && minimumBaseQuality != MINIMUM_BASE_QUALITY_DEFAULT) {
            logger.warn(
                    String.format("The value provided for the argument '%s' " +
                            "won't have any effect as coverage does " +
                            "not depend on individual base-calls",
                            MINIMUM_BASE_QUALITY_FULL_NAME));
        }
    }

    /**
     * Chooses the appropriate coverage calculator based on user arguments.
     * @return never {@code null}.
     */
    private BiFunction<Target, Stream<GATKRead>, long[]> composeCoverageCounter() {

        final BiFunction<Target, Stream<GATKRead>, long[]> unitCounter;
        switch (coverageUnit.datum) {
            case READ:
                unitCounter = new ReadsCounter();
                break;
            case BASE_CALL:
                unitCounter = new BaseCallsCounter();
                break;
            default:
                // This exception means that a new datum type was added to {@link CoverageDatum}.
                // Yet this switch was not updated to handle it.
                // Coverage Note: please ignore lack of coverage on this line as it is virtually impossible to cover.
                throw new GATKException("Unsupported coverage datum: " + coverageUnit.datum);
        }

        if (coverageUnit.basedOnFragments) {
            return (target, readContext) -> unitCounter.apply(target, condenseReadsIntoFragments(readContext));
        } else {
            return unitCounter;
        }

    }

    private Stream<GATKRead> condenseReadsIntoFragments(final Stream<GATKRead> readContext) {
        final Map<String, List<GATKRead>> readsByFragmentName = readContext
                .filter(read -> !read.isUnmapped())
                .collect(Collectors.groupingBy(GATKRead::getName));
        return readsByFragmentName.values().stream()
                .map(this::mergeReads);
    }

    private GATKRead mergeReads(final List<GATKRead> reads) {
        if (reads.size() == 1) {
            return reads.get(0);
        } else if (reads.size() == 2) {
            final GATKRead first = reads.get(0);
            final GATKRead second = reads.get(1);
            return MergedGATKReadPair.mergeReadPair(first, second);
        } else {
            throw new IllegalArgumentException("reads must either have one or two elements");
        }
    }

    @Override
    public Object onTraversalSuccess() {
        try {
            outputTableWriter.close();
            return null;
        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, "problem closing the output:"  + ex.getMessage());
        } finally {
            outputTableWriter = null;
        }
    }

    /**
     * Creates the output file writer.
     * @param outputFile the destination output file assumed no to be {@code null}.
     * @return never {@code null}.
     * @throws UserException.CouldNotCreateOutputFile if a {@link IOException} is thrown when creating the writer.
     */
    private static Writer createOutputWriter(final File outputFile) {
        try {
            return new FileWriter(outputFile);
        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, "could not open output file", ex);
        }
    }

    /**
     * Creates the read filter.
     *
     * <p>A read will fail the filter iff it does not qualify to count towards coverage</p>.
     * @return never {@code null}.
     */
    private CountingReadFilter makeReadFilter() {
        final CountingReadFilter baseFilter = new CountingReadFilter("Wellformed", new WellformedReadFilter(getHeaderForReads()))
                .and(new CountingReadFilter("Mapped", ReadFilterLibrary.MAPPED))
                .and(new CountingReadFilter("Primary", ReadFilterLibrary.PRIMARY_ALIGNMENT))
                .and(new CountingReadFilter("VendorQuality", ReadFilterLibrary.PASSES_VENDOR_QUALITY_CHECK))
                .and(new CountingReadFilter("Not_Duplicate", ReadFilterLibrary.NOT_DUPLICATE))
                .and(new CountingReadFilter("Non_Zero_Reference_Length", ReadFilterLibrary.NON_ZERO_REFERENCE_LENGTH_ALIGNMENT));

        // Adds the MQ filter only if it applies, i.e minMQ > 0:
        return minimumMappingQuality <= 0 ? baseFilter :
                baseFilter.and(new CountingReadFilter("MinMQ_" + minimumMappingQuality, read -> read.getMappingQuality() >= minimumMappingQuality));
    }

    @Override
    public void apply(final Target target, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {

        final long[] absoluteCounts = counter.apply(target, StreamSupport.stream(readsContext.spliterator(), false).filter(readFilter));
        final ReadCountRecord record = createReadCountRecord(target, absoluteCounts);
        try {
            outputTableWriter.writeRecord(record);
        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, ex);
        }
    }

    /**
     * Creates a read-count-record given the absolute counts of the requested coverage unit's datum.
     * <p>
     *     This method applies the transformation required by the requested coverage unit.
     * </p>
     * @param target the corresponding target.
     * @param absoluteCounts the absolute counts for that target.
     * @return never {@code null}.
     */
    private ReadCountRecord createReadCountRecord(final Target target, final long[] absoluteCounts) {
        final LongStream absoluteCountsStream = LongStream.of(absoluteCounts);

        final DoubleStream transformedCountsStream;
        switch (coverageUnit.transformation) {
            case NONE:
                transformedCountsStream = absoluteCountsStream.mapToDouble(l -> l);
                break;
            case AVERAGE_PER_BP:
                final double targetWidth = target.getInterval().getEnd() - target.getInterval().getStart() + 1;
                transformedCountsStream = absoluteCountsStream.mapToDouble(l -> l / targetWidth);
                break;
            default:
                // This exception mean that a new transformation was added to
                // {@link CoverageTransformation} but this switch was not updated to handle it.
                // Coverage Note: please ignore lack of coverage on this line as it is virtually impossible to cover.
                throw new GATKException("Unsupported coverage-unit transformation: " + coverageUnit.transformation);
        }

        // The final coverage counts stream applies the maximumCoverage cap if it applies:
        final DoubleStream finalCountsStream = (maximumCoverage == Double.POSITIVE_INFINITY) ? transformedCountsStream : transformedCountsStream.map(d -> Math.min(d, maximumCoverage));

        return new ReadCountRecord(target, finalCountsStream.toArray());
    }

    ///////////////////////////////////////
    // Code for the different coverage data counters at our disposal.
    ///////////////////////////////////////

    /**
     * Each read is considered an independent coverage unit on the target.
     */
    private class ReadsCounter implements BiFunction<Target, Stream<GATKRead>, long[]> {

        @Override
        public long[] apply(final Target target, final Stream<GATKRead> readsContext) {
            final long[] counts = new long[countColumnCount];
            readsContext
                    .mapToInt(readToColumn)
                    .filter(i -> i >= 0)
                    .forEach(i -> counts[i]++);
            return counts;
        }
    }

    /**
     * Each base call is counts toward coverage where the final output is the total number of overlapping
     * based divided by the target width, thus the average base call depth per bp.
     */
    private class BaseCallsCounter implements BiFunction<Target, Stream<GATKRead>, long[]> {

        @Override
        public long[] apply(final Target target, final Stream<GATKRead> readsContext) {
            final long[] counts = new long[countColumnCount];
            final SimpleInterval interval = target.getInterval();
            readsContext
                    .forEach(read -> {
                        final int columnIndex = readToColumn.applyAsInt(read);
                        if (columnIndex >= 0) {
                            counts[columnIndex] += countOverlappingQualifyingBaseCalls(read, interval);
                        }
                    });
            return counts;
        }

        /**
         * Given a read and an interval, this method calculates how many base calls in that read fall within the
         * interval and qualify to be counted
         * @param read the query read.
         * @param interval the target interval.
         * @return 0 or greater.
         */
        private long countOverlappingQualifyingBaseCalls(final GATKRead read, final SimpleInterval interval) {
            final int intervalStart = interval.getStart();
            final int intervalEnd = interval.getEnd();
            int referencePosition = read.getStart();
            int readPosition = 0;
            long overlap = 0;
            final IntBinaryOperator poorBaseCallsCounter = composePoorBaseCallCalculator(read);

            // We iterate thru each cigar element updating readPosition and referencePosition to reflect
            // the offset of the alignment block on the read and reference contig respectively.
            // We use these offsets, the interval boundaries (intervalStart, intervalEnd) and poor
            // base call count calculator to calculate the number of qualifying base calls at each
            // alignment blocks.
            // These are added in the overlap variable.
            for (final CigarElement element : read.getCigar().getCigarElements()) {
                final CigarOperator operator = element.getOperator();
                final int length = element.getLength();
                final boolean consumesReadBases = operator.consumesReadBases();
                final boolean consumesReferenceBases = operator.consumesReferenceBases();
                // If the block passes the start of the target interval and consumes reference and
                // read bases then it might contain some overlapping base-calls:
                if (consumesReadBases && consumesReferenceBases && (referencePosition + length) > intervalStart) {
                    overlap += calculateOverlapIncrement(intervalStart, intervalEnd, referencePosition, readPosition,
                            length, poorBaseCallsCounter);
                }
                // Update reference and read positions:
                if ((referencePosition = consumesReferenceBases ? referencePosition + length : referencePosition) > intervalEnd) {
                   break; // early termination when the alignment goes beyond the end of the target interval.
                }
                readPosition = consumesReadBases ? readPosition + length : readPosition;
            }
            return overlap;
        }
    }

    /**
     * Calculates the number of qualify base calls from a read sub-region to a reference sub-region within an interval.
     * @param intervalStart the enclosing interval start position.
     * @param intervalEnd the enclosing interval end position.
     * @param referenceStart the offset of the sub-region on the reference contig.
     * @param readStart the offset of the sub-region on the read
     * @param poorBaseCallsCounter function that calculates the number of non-qualifying base calls in a read sub-region.
     * @param length the length of the sub-region.
     * @return 0 or greater.
     */
    private int calculateOverlapIncrement(final int intervalStart, final int intervalEnd, final int referenceStart,
                                          final int readStart, final int length, final IntBinaryOperator poorBaseCallsCounter) {
        final int overlapStart = Math.max(intervalStart, referenceStart);
        final int overlapEnd = Math.min(intervalEnd, referenceStart + length - 1);
        final int overlapLength = overlapEnd - overlapStart + 1;
        final int readOverlapStart = readStart + (overlapStart - referenceStart);
        final int poorBaseCalls = poorBaseCallsCounter.applyAsInt(readOverlapStart, overlapLength);
        return overlapLength - poorBaseCalls;
    }

    // Creates a poor base-call counting function for sub-region on a particular read.
    // This function takes on two int arguments.
    // The first one is the first position within the read or (offset).
    // and the second is the length of the sub-region in the read to consider.
    private IntBinaryOperator composePoorBaseCallCalculator(final GATKRead read) {
        if (minimumBaseQuality == 0) {
            return (a, b) -> 0;
        } else {
            final byte[] qualities = read.getBaseQualities();
            // If there is no base-qualities, no bases are considered poor (optimistic assumption).
            if (qualities == null || qualities.length == 0) {
                return (a, b) -> 0;
            } else {
                return (offset, length) -> {
                    int result = 0;
                    final int to = offset + length;
                    for (int i = offset; i < to; i++) {
                        if (qualities[i] < minimumBaseQuality) {
                            result++;
                        }
                    }
                    return result;
                };
            }
        }
    }
}
