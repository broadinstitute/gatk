package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import org.apache.commons.collections4.ListUtils;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.tsv.*;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.*;
import java.util.function.BiConsumer;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class SegmentUtils {
    private SegmentUtils() {}

    /**
     * Returns the mean of coverages for targets in a given collection that overlap a given interval.
     */
    public static double meanTargetCoverage(final Locatable interval,
                                            final TargetCollection<TargetCoverage> targets) {
        Utils.nonNull(interval, "Can't get mean coverage of null interval.");
        Utils.nonNull(targets, "The collection of target coverages cannot be null.");

        final List<TargetCoverage> myTargets = targets.targets(interval);

        if (myTargets.size() == 0) {
            throw new IllegalStateException("Empty segment -- no overlapping targets.");
        }
        return myTargets.stream().mapToDouble(TargetCoverage::getCoverage).average().getAsDouble();
    }

    /**
     * Get difference between a segment mean and the overlapping targets.
     * This will select the overlapping targets from the input collection.
     *
     * @param segment   segment to use in calculating all difference of target CRs from the mean
     * @param targets   targets.  Note that this method will use only the targets overlapping the given segment.
     * @return          never {@code null}.  List of doubles of the difference, possibly empty.  Modifiable.
     */
    public static List<Double> segmentMeanTargetDifference(final ModeledSegment segment,
                                                           final TargetCollection<TargetCoverage> targets) {
        Utils.nonNull(segment, "Can't count targets of null segment.");
        Utils.nonNull(targets, "Counting targets requires non-null targets collection.");

        final double segMean = segment.getSegmentMeanInCRSpace();
        final List<TargetCoverage> myTargets = targets.targets(segment);
        return myTargets.stream().map(t -> (Math.pow(2, t.getCoverage()) - segMean)).collect(Collectors.toList());
    }

    /*===============================================================================================================*
     * PUBLIC METHODS FOR SEGMENT UNION (COMBINING TARGET AND SNP SEGMENTS)                                          *
     *===============================================================================================================*/

    /**
     * Target segments from CBS (i.e., {@link PerformSegmentation}) are specified by target-end--target-end intervals.
     * This method converts these to the corresponding target-start--target-end intervals, returning a new,
     * modifiable list of segments.
     */
    public static List<SimpleInterval> fixTargetSegmentStarts(final List<SimpleInterval> targetSegments,
                                                              final TargetCollection<? extends Locatable> targets) {
        Utils.nonNull(targetSegments, "The list of target segments cannot be null.");
        Utils.nonNull(targets, "The collection of targets cannot be null.");

        final List<SimpleInterval> targetSegmentsFixed = new ArrayList<>();
        for (final SimpleInterval segment : targetSegments) {
            try {
                final Locatable firstTarget = targets.targets(segment).get(0);
                if (firstTarget.getEnd() == segment.getStart()) {
                    targetSegmentsFixed.add(
                            new SimpleInterval(segment.getContig(), firstTarget.getStart(), segment.getEnd()));
                } else {
                    throw new IllegalArgumentException("List of targets and segments do not correspond.");
                }
            } catch (final IndexOutOfBoundsException ex) {
                throw new IllegalArgumentException("List of targets and segments do not correspond.");
            }
        }
        return targetSegmentsFixed;
    }

    /**
     * Given an interval and collections of targets and SNPs, returns a trimmed interval produced by removing the empty
     * portions at the start and the end of the original interval that do not overlap the targets and SNPs that overlap
     * with the original interval.  If this procedure would remove the entire interval, the original interval is
     * returned instead.  Note that this method will not expand an interval to the start of the first overlapping target
     * and the end of the last overlapping target; it will only shrink the interval or leave it alone.  This is to
     * avoid overlapping segments (which would occur if a SNP breakpoint fell within a target and the interval
     * were expanded, for example).
     */
    public static SimpleInterval trimInterval(final Locatable interval,
                                              final TargetCollection<? extends Locatable> targets,
                                              final TargetCollection<? extends Locatable> snps) {
        Utils.nonNull(interval, "The interval cannot be null.");
        Utils.nonNull(targets, "The collection of targets cannot be null.");
        Utils.nonNull(snps, "The collection of SNPs cannot be null.");

        final IndexRange targetRange = targets.indexRange(interval);
        final IndexRange snpRange = snps.indexRange(interval);

        final int numTargetsInInterval = targetRange.size();
        final int numSNPsInInterval = snpRange.size();

        int start = interval.getStart();
        int end = interval.getEnd();

        if (numTargetsInInterval == 0 && numSNPsInInterval > 0) {
            //if there are no targets overlapping interval, use SNPs to determine trimmed interval
            start = snps.target(snpRange.from).getStart();
            end = snps.target(snpRange.to - 1).getEnd();
        } else if (numTargetsInInterval > 0) {
            //if interval start does not fall within first target, use start of first target as start of trimmed interval
            start = Math.max(start, targets.target(targetRange.from).getStart());
            //if interval end does not fall within last target, use end of last target as end of trimmed interval
            end = Math.min(end, targets.target(targetRange.to - 1).getEnd());

            if (numSNPsInInterval > 0) {
                //if there are also SNPs within interval, check to see if they give a larger trimmed interval
                start = Math.min(start, snps.target(snpRange.from).getStart());
                end = Math.max(end, snps.target(snpRange.to - 1).getEnd());
            }
        }
        if (start < interval.getStart() || end > interval.getEnd() || end < start) {
            throw new GATKException.ShouldNeverReachHereException("Something went wrong in trimming interval.");
        }
        return new SimpleInterval(interval.getContig(), start, end);
    }

    /**
     * Returns the union of target and SNP segments.  First, all breakpoints from both sets of segments are combined
     * to form new segments. Spurious segments (i.e., segments containing only targets that are created by
     * SNP-segment breakpoints and are not present in the original set of target segments)
     * at the starts and ends of the original target segments are then remerged to the right and left, respectively;
     * spurious segments introduced within the original target segments are merged with adjacent segments
     * according to the logic in {@link org.broadinstitute.hellbender.tools.exome.SegmentMergeUtils.SmallSegments}.
     * Finally, the segments are trimmed by {@link SegmentUtils#trimInterval}.
     */
    public static List<SimpleInterval> unionSegments(final List<SimpleInterval> targetSegments,
                                                     final List<SimpleInterval> snpSegments,
                                                     final Genome genome) {
        Utils.nonNull(targetSegments, "The list of target segments cannot be null.");
        Utils.nonNull(snpSegments, "The list of SNP segments cannot be null.");
        Utils.nonNull(genome, "The genome cannot be null.");

        final SortedMap<String, List<Breakpoint>> breakpointsByContig =
                collectBreakpointsByContig(targetSegments, snpSegments);
        List<SimpleInterval> unionedSegments = constructUntrimmedSegments(genome, breakpointsByContig);
        unionedSegments = SegmentUtils.mergeSpuriousStartsAndEnds(unionedSegments, targetSegments, genome.getSNPs());
        unionedSegments = SegmentMergeUtils.mergeSpuriousMiddles(unionedSegments, targetSegments, genome);
        return unionedSegments.stream().map(s -> trimInterval(s, genome.getTargets(), genome.getSNPs()))
                .collect(Collectors.toList());
    }

    /*===============================================================================================================*
     * PUBLIC METHODS FOR READING/WRITING SEGMENT FILES                                                              *
     *===============================================================================================================*/

    /**
     * Reads a list of intervals from a segment file with header:
     * <p>
     *      Sample  Chromosome  Start   End     ...
     * </p>
     *
     * Other columns are optional and are ignored.
     */
    public static List<SimpleInterval> readIntervalsFromSegmentFile(final File segmentsFile) {
        return readSegmentFile(segmentsFile, SegmentTableColumns.INTERVAL_COLUMN_NAME_ARRAY, SegmentUtils::toInterval);
    }

    /**
     * Reads a list of segments as represented by {@link ModeledSegment}
     * (with calls, if possible) from a file with header:
     * <p>
     *      Sample  Chromosome  Start   End     Num_Probes  Segment_Mean    Segment_Call (optional)     ...
     * </p>
     *
     * Segment means in the input file are assumed to be in CR space (not log2).  Please note that
     * {@link ModeledSegment} stores the segment mean as log2CR.
     *
     * Note that Segment_Call is optional.  Other columns are optional and are ignored.
     */
    public static List<ModeledSegment> readModeledSegmentsFromSegmentFile(final File segmentsFile) {
        return readSegmentFile(segmentsFile, SegmentTableColumns.NO_CALL_COLUMN_NAME_ARRAY,
                dataLine -> toModeledSegment(dataLine, true));
    }

    /**
     * Reads a list of segments as represented by {@link ModeledSegment}
     * (with calls, if possible) from a legacy file with header:
     * <p>
     *      Sample  Chromosome  Start   End     Num_Probes  Segment_Mean    Segment_Call (optional)     ...
     * </p>
     *
     * Segment means in the input file are assumed to be in log2CR space.  Please note that
     * {@link ModeledSegment} stores the segment mean as log2CR.
     *
     * Note that Segment_Call is optional.  Other columns are optional and are ignored.
     */
    public static List<ModeledSegment> readModeledSegmentsFromLegacySegmentFile(final File segmentsFile) {
        return readSegmentFile(segmentsFile, SegmentTableColumns.NO_CALL_COLUMN_NAME_ARRAY,
                dataLine -> toModeledSegment(dataLine, false));
    }

    /**
     * Writes a list of segments with calls represented by {@link ModeledSegment} to a file with header:
     * <p>
     *      Sample  Chromosome  Start   End     Num_Probes  Segment_Mean    Segment_Call
     * </p>
     *
     * Note that Segment_Mean is output in CR space (not log2).
     *
     * If there is no call, the corresponding field is populated with a blank value.
     */
    public static void writeModeledSegmentFile(final File outFile,
                                               final List<ModeledSegment> segments,
                                               final String sampleName) {
        writeSegmentFile(outFile, segments, sampleName, SegmentTableColumns.MEAN_AND_CALL_COLUMN_NAME_ARRAY,
                //lambda for appending segment fields to a DataLine
                (segment, dataLine) ->
                    dataLine.append(segment.getContig()).append(segment.getStart(), segment.getEnd())
                            .append(segment.getOriginalProbeCount())
                            .append(segment.getSegmentMeanInCRSpace())
                            .append(segment.getCall()));
    }

    /**
     * Writes a list of segments with number of targets and SNPs in each segment to a file with header:
     * <p>
     *      Sample  Chromosome  Start   End     Num_Targets     Num_SNPs
     * </p>
     */
    public static void writeSegmentFileWithNumTargetsAndNumSNPs(final File outFile,
                                                                final List<? extends Locatable> segments,
                                                                final Genome genome) {
        Utils.nonNull(genome, "The genome cannot be null.");

        final TargetCollection<TargetCoverage> targets = genome.getTargets();
        final TargetCollection<AllelicCount> snps = genome.getSNPs();
        final String sampleName = genome.getSampleName();

        writeSegmentFile(outFile, segments, sampleName, SegmentTableColumns.NUM_TARGETS_AND_SNPS_COLUMN_NAME_ARRAY,
                //lambda for appending segment fields to a DataLine
                (segment, dataLine) ->
                        dataLine.append(segment.getContig()).append(segment.getStart(), segment.getEnd())
                                .append(targets.targetCount(segment)).append(snps.targetCount(segment)));
    }

    /**
     * Writes a list of segments represented by {@link ACNVModeledSegment} to a file with header:
     * <p>
     *      Sample  Chromosome  Start   End     Num_Targets     Num_SNPs
     *      Segment_Mean_Post_Mean  Segment_Mean_Post_Std   MAF_Post_Mean   MAF_Post_Std
     * </p>
     */
    public static void writeACNVModeledSegmentFile(final File outFile,
                                                   final List<ACNVModeledSegment> segments,
                                                   final Genome genome) {
        Utils.nonNull(genome, "The genome cannot be null.");

        final TargetCollection<TargetCoverage> targets = genome.getTargets();
        final TargetCollection<AllelicCount> snps = genome.getSNPs();
        final String sampleName = genome.getSampleName();

        writeSegmentFile(outFile, segments, sampleName, SegmentTableColumns.ACNV_MODELED_SEGMENT_COLUMN_NAME_ARRAY,
                //lambda for appending segment fields to a DataLine
                (segment, dataLine) ->
                        dataLine.append(segment.getContig()).append(segment.getStart(), segment.getEnd())
                                .append(targets.targetCount(segment)).append(snps.targetCount(segment))
                                .append(String.format("%6.4f", segment.getSegmentMeanPosteriorSummary().mean()))
                                .append(String.format("%6.4f", segment.getSegmentMeanPosteriorSummary().standardDeviation()))
                                .append(String.format("%6.4f", segment.getMinorAlleleFractionPosteriorSummary().mean()))
                                .append(String.format("%6.4f", segment.getMinorAlleleFractionPosteriorSummary().standardDeviation())));
    }

    /*===============================================================================================================*
     * PRIVATE METHODS FOR SEGMENT UNION (COMBINING TARGET AND SNP SEGMENTS)                                         *
     *===============================================================================================================*/

    private enum BreakpointType {
        START, END
    }

    private static final class Breakpoint extends Interval {
        final BreakpointType type;

        public Breakpoint(final BreakpointType type, final String contig, final int site) {
            super(contig, site, site);
            this.type = type;
        }
        public BreakpointType getType() {
            return type;
        }

        public int getSite() {
            return getStart();
        }
    }

    //returns a map of (contig -> combined breakpoints from target and SNP segments on that contig)
    private static SortedMap<String, List<Breakpoint>> collectBreakpointsByContig(final List<SimpleInterval> targetSegments,
                                                                                  final List<SimpleInterval> snpSegments) {
        return ListUtils.union(targetSegments, snpSegments).stream()
                .map(s -> Arrays.asList(
                        new Breakpoint(BreakpointType.START, s.getContig(), s.getStart()),
                        new Breakpoint(BreakpointType.END, s.getContig(), s.getEnd())))
                .flatMap(Collection::stream)
                .sorted()
                .collect(Collectors.groupingBy(Interval::getContig, TreeMap::new, Collectors.toList()));
    }

    //given breakpointsByContig map, constructs list of untrimmed segments
    //(i.e., segments that are directly adjacent to each other; segment boundaries at each breakpoint are determined
    //by whether that breakpoint was a start or an end for the corresponding original target/SNP segment),
    //then returns only those untrimmed segments that are non-empty
    private static List<SimpleInterval> constructUntrimmedSegments(final Genome genome,
                                                                   final SortedMap<String, List<Breakpoint>> breakpointsByContig) {
        final List<SimpleInterval> segments = new ArrayList<>();
        for (final String contig : breakpointsByContig.keySet()) {
            final List<Breakpoint> breakpoints = breakpointsByContig.get(contig);
            int start = breakpoints.get(0).getStart();
            for (int i = 1; i < breakpoints.size(); i++) {
                final Breakpoint breakpoint = breakpoints.get(i);
                //adjust segment boundaries according to breakpoint type; this prevents, e.g., SNPs that originally
                //started SNP segments from being stranded
                final int end = breakpoint.getSite() - (breakpoint.getType() == BreakpointType.START ? 1 : 0);

                if (end < start) {  //this could happen if there are adjacent breakpoints
                    continue;
                }

                final SimpleInterval segment = new SimpleInterval(contig, start, end);
                if (genome.getTargets().targetCount(segment) > 0 || genome.getSNPs().targetCount(segment) > 0) {
                    segments.add(segment);
                }
                start = segment.getEnd() + 1;
            }
        }
        return segments;
    }

    //remerge spurious target-only segments introduced at starts and ends of the
    //original target-coverage segments by union of breakpoints
    private static List<SimpleInterval> mergeSpuriousStartsAndEnds(final List<SimpleInterval> segments,
                                                                   final List<SimpleInterval> targetSegments,
                                                                   final TargetCollection<AllelicCount> snps) {
        //get original target-segment starts and ends
        final Set<SimpleInterval> targetSegmentStarts =
                targetSegments.stream().map(s -> new SimpleInterval(s.getContig(), s.getStart(), s.getStart()))
                        .collect(Collectors.toSet());
        final Set<SimpleInterval> targetSegmentEnds =
                targetSegments.stream().map(s -> new SimpleInterval(s.getContig(), s.getEnd(), s.getEnd()))
                        .collect(Collectors.toSet());

        final List<SimpleInterval> mergedSegments = new ArrayList<>();
        final ListIterator<SimpleInterval> segmentsIter = segments.listIterator();
        while (segmentsIter.hasNext()) {
            final SimpleInterval segment = segmentsIter.next();
            //do not remerge non-target-only segments (i.e., those containing SNPs)
            if (snps.targetCount(segment) > 0) {
                mergedSegments.add(segment);
                continue;
            }
            final SimpleInterval segmentStart =
                    new SimpleInterval(segment.getContig(), segment.getStart(), segment.getStart());
            final SimpleInterval segmentEnd =
                    new SimpleInterval(segment.getContig(), segment.getEnd(), segment.getEnd());
            if (targetSegmentStarts.contains(segmentStart) && !targetSegmentEnds.contains(segmentEnd)) {
                //remerge segments introduced at starts to the right
                final SimpleInterval nextSegment = segmentsIter.next();
                mergedSegments.add(SegmentMergeUtils.mergeSegments(segment, nextSegment));
            } else if (!targetSegmentStarts.contains(segmentStart) && targetSegmentEnds.contains(segmentEnd)) {
                //remerge segments introduced at ends to the left
                final int previousIndex = mergedSegments.size() - 1;
                final SimpleInterval previousSegment = mergedSegments.get(previousIndex);
                mergedSegments.set(previousIndex, SegmentMergeUtils.mergeSegments(previousSegment, segment));
            } else {
                //do not merge otherwise; although some spurious segments remain, they will be merged in a later step
                mergedSegments.add(segment);
            }
        }
        return mergedSegments;
    }

    /*===============================================================================================================*
     * PRIVATE METHODS FOR READING/WRITING SEGMENT FILES                                                             *
     *===============================================================================================================*/

    private static SimpleInterval toInterval(final DataLine dataLine) {
        return new SimpleInterval(dataLine.get(SegmentTableColumns.CONTIG.toString()),
                dataLine.getInt(SegmentTableColumns.START.toString()), dataLine.getInt(SegmentTableColumns.END.toString()));

    }

    private static ModeledSegment toModeledSegment(final DataLine dataLine, final boolean takeLog2) {
        final double mean = dataLine.getDouble(SegmentTableColumns.MEAN.toString());

        return new ModeledSegment(
                toInterval(dataLine),
                dataLine.get(SegmentTableColumns.CALL.toString(), ModeledSegment.NO_CALL),
                dataLine.getLong(SegmentTableColumns.NUM_PROBES.toString()),
                takeLog2 ? ParamUtils.log2(mean) : mean);
    }

    private static <T extends Locatable> List<T> readSegmentFile(final File segmentsFile,
                                                                 final String[] mandatoryColumns,
                                                                 final Function<DataLine, T> dataLineToSegmentFunction) {
        try (final TableReader<T> reader = TableUtils.reader(segmentsFile,
                (columns, formatExceptionFactory) -> {
                    if (!columns.containsAll(mandatoryColumns)) {
                        throw formatExceptionFactory.apply("Bad header in segment file.");
                    }
                    //return the lambda to translate dataLines into called segments
                    return dataLineToSegmentFunction;
                })) {
            return reader.stream().collect(Collectors.toList());
        } catch (final IOException | UncheckedIOException e) {
            throw new UserException.CouldNotReadInputFile(segmentsFile, e);
        }
    }

    private static <T extends Locatable> void writeSegmentFile(final File outFile,
                                                               final List<T> segments,
                                                               final String sampleName,
                                                               final String[] columnArray,
                                                               final BiConsumer<T, DataLine> segmentToDataLineFunction) {
        Utils.nonNull(segments, "The list of segments cannot be null.");
        try (final TableWriter<T> writer =
                     TableUtils.writer(outFile, new TableColumnCollection(columnArray),
                             //lambda for creating DataLine with sampleName and segment fields
                             (segment, dataLine) -> {
                                 segmentToDataLineFunction.accept(segment, dataLine.append(sampleName));
                             }
                     )) {
            for (final T segment : segments) {
                writer.writeRecord(Utils.nonNull(segment, "List of segments contains a null."));
            }
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(outFile, e);
        }
    }
}
