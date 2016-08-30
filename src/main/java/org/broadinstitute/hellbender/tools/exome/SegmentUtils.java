package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import org.apache.commons.collections4.ListUtils;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.tools.exome.cnlohcaller.CNLOHCall;
import org.broadinstitute.hellbender.tools.exome.samplenamefinder.SampleNameFinder;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.mcmc.Decile;
import org.broadinstitute.hellbender.utils.mcmc.DecileCollection;
import org.broadinstitute.hellbender.utils.mcmc.PosteriorSummary;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.tsv.*;

import java.io.*;
import java.util.*;
import java.util.function.BiConsumer;
import java.util.function.Function;
import java.util.stream.Collectors;

import static org.broadinstitute.hellbender.tools.exome.ACNVModeller.ACNV_DOUBLE_FORMAT;

/**
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class SegmentUtils {
    private SegmentUtils() {}

    /**
     * Get difference between a segment mean and the overlapping targets.
     * This will select the overlapping targets from the input collection.
     *
     * @param segment   segment to use in calculating all difference of target CRs from the mean
     * @param targets   targets.  Note that this method will use only the targets overlapping the given segment.
     * @return          never {@code null}.  List of doubles of the difference, possibly empty.  Modifiable.
     */
    public static List<Double> segmentMeanTargetDifference(final ModeledSegment segment,
                                                           final TargetCollection<ReadCountRecord.SingleSampleRecord> targets) {
        Utils.nonNull(segment, "Can't count targets of null segment.");
        Utils.nonNull(targets, "Counting targets requires non-null targets collection.");

        final double segMean = segment.getSegmentMeanInCRSpace();
        return targets.targets(segment).stream().map(t -> (Math.pow(2, t.getCount()) - segMean)).collect(Collectors.toList());
    }

    /*===============================================================================================================*
     * PUBLIC METHODS FOR SEGMENT UNION (COMBINING TARGET AND SNP SEGMENTS)                                          *
     *===============================================================================================================*/

    /**
     * Legacy target segments from CBS (i.e., from legacy versions of {@link PerformSegmentation})
     * are specified by target-end--target-end intervals.  This method converts these to the corresponding
     * target-start--target-end intervals, returning a new, modifiable list of segments.  Non-legacy segments that are
     * already specified by target-start--target-end intervals are also returned as a new, modifiable list of segments
     * with no changes.
     */
    public static List<SimpleInterval> fixLegacyTargetSegmentStarts(final List<SimpleInterval> targetSegments,
                                                                    final TargetCollection<? extends Locatable> targets) {
        Utils.nonNull(targetSegments, "The list of target segments cannot be null.");
        Utils.nonNull(targets, "The collection of targets cannot be null.");

        final List<SimpleInterval> targetSegmentsFixed = new ArrayList<>();
        for (final SimpleInterval segment : targetSegments) {
            try {
                final Locatable firstTarget = targets.targets(segment).get(0);
                if (firstTarget.getEnd() == segment.getStart() || firstTarget.getStart() == segment.getStart()) {
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
        return readSegmentFile(segmentsFile, SegmentTableColumn.INTERVAL_COLUMNS, SegmentUtils::toInterval);
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
        return readSegmentFile(segmentsFile, SegmentTableColumn.MEAN_AND_NO_CALL_COLUMNS,
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
        return readSegmentFile(segmentsFile, SegmentTableColumn.MEAN_AND_NO_CALL_COLUMNS,
                dataLine -> toModeledSegment(dataLine, false));
    }

    /**
     * Read only the sample names from a segment file.
     *
     * @param segmentFile Never {@code null}
     * @return a list of the unique sample names in the seg file.
     */
    public static List<String> readSampleNamesFromSegmentFile(final File segmentFile) {
        Utils.nonNull(segmentFile);
        Utils.regularReadableUserFile(segmentFile);
        try (final Reader segmentFileReader = new FileReader(segmentFile)) {
            return readSampleNamesFromSegmentReader(segmentFileReader, segmentFile.getAbsolutePath());
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile("Could not read the segment file: " + segmentFile.getAbsolutePath());
        }
    }

    /**
     * Read only the sample names from a segment reader.
     *
     * @param segmentReader Never {@code null}; an instance of {@link Reader} for input segments
     * @param segmentSourceName Never {@code null}; a string identifier for the input segments reader (used in error messages)
     * @return a list of the unique sample names in the seg file.
     */
    public static List<String> readSampleNamesFromSegmentReader(final Reader segmentReader, final String segmentSourceName) {
        Utils.nonNull(segmentReader);
        Utils.nonNull(segmentSourceName);
        try (final TableReader<String> reader = TableUtils.reader(segmentReader,
                (columns, formatExceptionFactory) -> {
                    if (!columns.contains(SegmentTableColumn.SAMPLE.toString())) {
                        throw formatExceptionFactory.apply("No sample name in segment file.");
                    }
                    //return the lambda to translate dataLines into called segments
                    return SegmentUtils::toSampleName;
                })) {

            // Return the unique sample name values in a list.
            return new ArrayList<>(reader.stream().collect(Collectors.toSet()));
        } catch (final IOException | UncheckedIOException e) {
            throw new UserException.CouldNotReadInputFile(segmentSourceName, e);
        }
    }

    /**
     * Writes a list of segments with calls represented by {@link ModeledSegment} to a file with header:
     * <p>
     *      Sample  Chromosome  Start   End     Num_Probes  Segment_Mean    Segment_Call
     * </p>
     *
     * Note that Segment_Mean is output in CR space (not log2) if isLog2Output is false
     *
     * If there is no call, the corresponding field is populated with a blank value.
     */
    public static <S extends ModeledSegment> void writeModeledSegmentFile(final File outFile,
                                               final List<S> segments,
                                               final String sampleName, final boolean isLog2Output) {
        writeSegmentFile(outFile, segments, sampleName, SegmentTableColumn.MEAN_AND_CALL_COLUMNS,
                //lambda for appending segment fields to a DataLine
                (segment, dataLine) ->
                        dataLine.append(segment.getContig()).append(segment.getStart(), segment.getEnd())
                                .append(segment.getTargetCount())
                                .append(isLog2Output ? segment.getSegmentMean() : segment.getSegmentMeanInCRSpace())
                                .append(segment.getCall()));
    }

    /**
     * Writes a list of segments with calls represented by {@link ModeledSegment} to a file with header:
     * <p>
     *      Sample  Chromosome  Start   End     Num_Probes  Segment_Mean    Segment_Call
     * </p>
     *
     * Note that Segment_Mean is output in CR space (not log2)
     *
     * If there is no call, the corresponding field is populated with a blank value.
     */
    public static <S extends ModeledSegment> void writeModeledSegmentFile(final File outFile,
                                                                          final List<S> segments,
                                                                          final String sampleName) {
        writeModeledSegmentFile(outFile, segments, sampleName, false);
    }

    /**
     * Search a segment file and extract all sample names.
     *
     * Throws an exception if more than one sample name is found.
     *
     * @param segmentFile segment file generated by one of the GATK tools, such as PerformSegmentation.
     * @return sample name
     */
    public static String getSampleNameForCLIsFromSegmentFile(final File segmentFile) {
        String sampleName;
        final List<String> sampleNames = SampleNameFinder.determineSampleNamesFromSegmentFile(segmentFile);
        if (sampleNames.size() == 1) {
            sampleName = sampleNames.get(0);
        } else if (sampleNames.size() > 1) {
            throw new UserException.BadInput("Input file must contain data for only one sample.  Found samples: " + StringUtils.join(sampleNames, ", "));
        } else {
            throw new UserException.BadInput("Input file must contain data for only one sample.  Could not find any sample information.");
        }
        return sampleName;
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

        final TargetCollection<ReadCountRecord.SingleSampleRecord> targets = genome.getTargets();
        final TargetCollection<AllelicCount> snps = genome.getSNPs();
        final String sampleName = genome.getSampleName();

        writeSegmentFile(outFile, segments, sampleName, SegmentTableColumn.NUM_TARGETS_AND_SNPS_COLUMNS,
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

        final TargetCollection<ReadCountRecord.SingleSampleRecord> targets = genome.getTargets();
        final TargetCollection<AllelicCount> snps = genome.getSNPs();
        final String sampleName = genome.getSampleName();

        writeSegmentFile(outFile, segments, sampleName, SegmentTableColumn.ACNV_MODELED_SEGMENT_COLUMNS,
                createDataLineFromACNVModeledSegmentFunction(targets, snps));
    }

    /**
     * Writes a list of segments represented by {@link ACNVModeledSegment} to a file with header:
     * <p>
     *      Sample  Chromosome  Start   End     Num_Targets     Num_SNPs
     *      Segment_Mean_Post_Mean  Segment_Mean_Post_Std   MAF_Post_Mean   MAF_Post_Std
     * </p>
     */
    public static void writeCnLoHACNVModeledSegmentFile(final File outFile,
                                                   final List<CNLOHCall> calls,
                                                   final Genome genome) {
        Utils.nonNull(genome, "The genome cannot be null.");

        final TargetCollection<ReadCountRecord.SingleSampleRecord> targets = genome.getTargets();
        final TargetCollection<AllelicCount> snps = genome.getSNPs();
        final String sampleName = genome.getSampleName();

        writeSegmentFile(outFile, calls, sampleName, SegmentTableColumn.CNLOH_ACNV_MODELED_SEGMENT_COLUMNS,
                createDataLineFromCNLoHCallsFunction(targets, snps));
    }

    /**
     * Read ACNV file into ACNV modeled segments.
     *
     * @param acnvSegFile File in the ACNV format (sim-*.seg).  Must be readable file and never {@code null}
     * @return Never {@code null}.  Can be empty list if no entries exist in the input file.
     */
    public static List<ACNVModeledSegment> readACNVModeledSegmentFile(final File acnvSegFile) {
        Utils.nonNull(acnvSegFile);
        Utils.regularReadableUserFile(acnvSegFile);
        return readSegmentFile(acnvSegFile, SegmentTableColumn.ACNV_MODELED_SEGMENT_COLUMNS, SegmentUtils::toACNVModeledSegment);
    }

    /** Simple conversion to an interval given standard segment columns
     *
     * @param dataLine line in a file
     * @return interval in the dataLine.  Never {@code null}
     */
    public static SimpleInterval toInterval(final DataLine dataLine) {
        Utils.nonNull(dataLine);
        return toInterval(dataLine, SegmentTableColumn.CONTIG.toString(), SegmentTableColumn.START.toString(),
                SegmentTableColumn.END.toString());

    }

    /** Simple conversion to an interval given standard segment columns
     *
     * @param dataLine line in a file
     * @param contigColumn name of the contig column
     * @param startColumn column name for the start position
     * @param endColumn column name for the end position
     * @return interval in the dataLine.  Never {@code null}
     */
    public static SimpleInterval toInterval(final DataLine dataLine, final String contigColumn, final String startColumn,
                                            final String endColumn) {
        Utils.nonNull(dataLine);
        Utils.nonNull(contigColumn);
        Utils.nonNull(startColumn);
        Utils.nonNull(endColumn);
        return new SimpleInterval(dataLine.get(contigColumn), dataLine.getInt(startColumn), dataLine.getInt(endColumn));

    }

    public static <T extends Locatable> List<T> readSegmentFile(final File segmentsFile,
                                                                final TableColumnCollection mandatoryColumns,
                                                                final Function<DataLine, T> dataLineToSegmentFunction) {
        Utils.nonNull(segmentsFile);
        Utils.regularReadableUserFile(segmentsFile);
        try (final TableReader<T> reader = TableUtils.reader(segmentsFile,
                (columns, formatExceptionFactory) -> {
                    TableUtils.checkMandatoryColumns(columns, mandatoryColumns, formatExceptionFactory);
                    //return the lambda to translate dataLines into called segments
                    return dataLineToSegmentFunction;
                })) {
            return reader.stream().collect(Collectors.toList());
        } catch (final IOException | UncheckedIOException e) {
            throw new UserException.CouldNotReadInputFile(segmentsFile, e);
        }
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


    private static ModeledSegment toModeledSegment(final DataLine dataLine, final boolean takeLog2) {
        final double mean = dataLine.getDouble(SegmentTableColumn.MEAN.toString());

        return new ModeledSegment(
                toInterval(dataLine),
                dataLine.get(SegmentTableColumn.CALL.toString(), ModeledSegment.NO_CALL),
                dataLine.getLong(SegmentTableColumn.NUM_PROBES.toString()),
                takeLog2 ? ParamUtils.log2(mean) : mean);
    }

    private static ACNVModeledSegment toACNVModeledSegment(final DataLine dataLine) {
        final PosteriorSummary segmentMeanPosteriorSummary = new PosteriorSummary(dataLine.getDouble(SegmentTableColumn.SEGMENT_MEAN_POSTERIOR_MODE.toString()),
                dataLine.getDouble(SegmentTableColumn.SEGMENT_MEAN_POSTERIOR_LOWER.toString()),
                dataLine.getDouble(SegmentTableColumn.SEGMENT_MEAN_POSTERIOR_UPPER.toString()));
        final PosteriorSummary minorAlleleFractionPosteriorSummary = new PosteriorSummary(dataLine.getDouble(SegmentTableColumn.MINOR_ALLELE_FRACTION_POSTERIOR_MODE.toString()),
                dataLine.getDouble(SegmentTableColumn.MINOR_ALLELE_FRACTION_POSTERIOR_LOWER.toString()),
                dataLine.getDouble(SegmentTableColumn.MINOR_ALLELE_FRACTION_POSTERIOR_UPPER.toString()));

        final DecileCollection mafDecileCollection = new DecileCollection(
                SegmentTableColumn.ACNV_MODELED_SEGMENT_MAF_DECILES_SUMMARY_COLUMNS.names().stream().map(dataLine::getDouble).collect(Collectors.toList()));

        final DecileCollection segmentMeanDecileCollection = new DecileCollection(
                SegmentTableColumn.ACNV_MODELED_SEGMENT_MEAN_DECILES_SUMMARY_COLUMNS.names().stream().map(dataLine::getDouble).collect(Collectors.toList()));

        minorAlleleFractionPosteriorSummary.setDeciles(mafDecileCollection);
        segmentMeanPosteriorSummary.setDeciles(segmentMeanDecileCollection);

        return new ACNVModeledSegment(
                toInterval(dataLine),
                segmentMeanPosteriorSummary,
                minorAlleleFractionPosteriorSummary);
    }

    private static String toSampleName(final DataLine dataLine) {
        return dataLine.get(SegmentTableColumn.SAMPLE.toString());
    }

    private static <T extends Locatable> void writeSegmentFile(final File outFile,
                                                               final List<T> segments,
                                                               final String sampleName,
                                                               final TableColumnCollection columns,
                                                               final BiConsumer<T, DataLine> segmentToDataLineFunction) {
        Utils.nonNull(segments, "The list of segments cannot be null.");
        try (final TableWriter<T> writer =
                     TableUtils.writer(outFile, columns,
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

    private static BiConsumer<CNLOHCall, DataLine> createDataLineFromCNLoHCallsFunction
            (final TargetCollection<ReadCountRecord.SingleSampleRecord> targets,
             final TargetCollection<AllelicCount> snps) {
        return (cnlohCall, dataLine) -> {
            createDataLineFromACNVModeledSegmentFunction(targets,snps).accept(cnlohCall.getAcnvSegment(), dataLine);
            dataLine.append(cnlohCall.getBalancedCall().toString())
                    .append(cnlohCall.getCnlohCall().toString())
                    .append(cnlohCall.getRho())
                    .append(cnlohCall.getM())
                    .append(cnlohCall.getN())
                    .append(cnlohCall.getfCr())
                    .append(cnlohCall.getfMaf());
        };
    }

    private static BiConsumer<ACNVModeledSegment, DataLine> createDataLineFromACNVModeledSegmentFunction
            (final TargetCollection<ReadCountRecord.SingleSampleRecord> targets,
             final TargetCollection<AllelicCount> snps) {

        //lambda for appending segment fields to a DataLine
        return (segment, dataLine) ->
                dataLine.append(segment.getContig()).append(segment.getStart(), segment.getEnd())
                        .append(targets.targetCount(segment)).append(snps.targetCount(segment))
                        .append(String.format(ACNV_DOUBLE_FORMAT, segment.getSegmentMeanPosteriorSummary().getCenter()))
                        .append(String.format(ACNV_DOUBLE_FORMAT, segment.getSegmentMeanPosteriorSummary().getLower()))
                        .append(String.format(ACNV_DOUBLE_FORMAT, segment.getSegmentMeanPosteriorSummary().getUpper()))
                        .append(String.format(ACNV_DOUBLE_FORMAT, segment.getSegmentMeanPosteriorSummary().getDeciles().get(Decile.DECILE_10)))
                        .append(String.format(ACNV_DOUBLE_FORMAT, segment.getSegmentMeanPosteriorSummary().getDeciles().get(Decile.DECILE_20)))
                        .append(String.format(ACNV_DOUBLE_FORMAT, segment.getSegmentMeanPosteriorSummary().getDeciles().get(Decile.DECILE_30)))
                        .append(String.format(ACNV_DOUBLE_FORMAT, segment.getSegmentMeanPosteriorSummary().getDeciles().get(Decile.DECILE_40)))
                        .append(String.format(ACNV_DOUBLE_FORMAT, segment.getSegmentMeanPosteriorSummary().getDeciles().get(Decile.DECILE_50)))
                        .append(String.format(ACNV_DOUBLE_FORMAT, segment.getSegmentMeanPosteriorSummary().getDeciles().get(Decile.DECILE_60)))
                        .append(String.format(ACNV_DOUBLE_FORMAT, segment.getSegmentMeanPosteriorSummary().getDeciles().get(Decile.DECILE_70)))
                        .append(String.format(ACNV_DOUBLE_FORMAT, segment.getSegmentMeanPosteriorSummary().getDeciles().get(Decile.DECILE_80)))
                        .append(String.format(ACNV_DOUBLE_FORMAT, segment.getSegmentMeanPosteriorSummary().getDeciles().get(Decile.DECILE_90)))
                        .append(String.format(ACNV_DOUBLE_FORMAT, segment.getMinorAlleleFractionPosteriorSummary().getCenter()).replaceAll("\\s+", "")) //removes whitespace for NaN entries
                        .append(String.format(ACNV_DOUBLE_FORMAT, segment.getMinorAlleleFractionPosteriorSummary().getLower()).replaceAll("\\s+", ""))
                        .append(String.format(ACNV_DOUBLE_FORMAT, segment.getMinorAlleleFractionPosteriorSummary().getUpper()).replaceAll("\\s+", ""))
                        .append(String.format(ACNV_DOUBLE_FORMAT, segment.getMinorAlleleFractionPosteriorSummary().getDeciles().get(Decile.DECILE_10)).replaceAll("\\s+", ""))
                        .append(String.format(ACNV_DOUBLE_FORMAT, segment.getMinorAlleleFractionPosteriorSummary().getDeciles().get(Decile.DECILE_20)).replaceAll("\\s+", ""))
                        .append(String.format(ACNV_DOUBLE_FORMAT, segment.getMinorAlleleFractionPosteriorSummary().getDeciles().get(Decile.DECILE_30)).replaceAll("\\s+", ""))
                        .append(String.format(ACNV_DOUBLE_FORMAT, segment.getMinorAlleleFractionPosteriorSummary().getDeciles().get(Decile.DECILE_40)).replaceAll("\\s+", ""))
                        .append(String.format(ACNV_DOUBLE_FORMAT, segment.getMinorAlleleFractionPosteriorSummary().getDeciles().get(Decile.DECILE_50)).replaceAll("\\s+", ""))
                        .append(String.format(ACNV_DOUBLE_FORMAT, segment.getMinorAlleleFractionPosteriorSummary().getDeciles().get(Decile.DECILE_60)).replaceAll("\\s+", ""))
                        .append(String.format(ACNV_DOUBLE_FORMAT, segment.getMinorAlleleFractionPosteriorSummary().getDeciles().get(Decile.DECILE_70)).replaceAll("\\s+", ""))
                        .append(String.format(ACNV_DOUBLE_FORMAT, segment.getMinorAlleleFractionPosteriorSummary().getDeciles().get(Decile.DECILE_80)).replaceAll("\\s+", ""))
                        .append(String.format(ACNV_DOUBLE_FORMAT, segment.getMinorAlleleFractionPosteriorSummary().getDeciles().get(Decile.DECILE_90)).replaceAll("\\s+", ""));
    }
}
