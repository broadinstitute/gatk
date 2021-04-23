package org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import scala.Tuple2;

import javax.annotation.Nonnull;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigAlignmentsArgumentCollection.GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY;

/**
 * Locally assembled contig:
 * its name
 * its sequence as produced by the assembler (no reverse complement like in the SAM record if it maps to '-' strand), and
 * its stripped-down alignment information.
 */
@DefaultSerializer(AlignedContig.Serializer.class)
public final class AlignedContig {

    /**
     * parameters to be passed to {@link AlignedContig#removeNonUniqueMappings(GoodAndBadMappings, int, int)}
     * for dropping alignments that offer either low read uniqueness.
     */
    public static final int ALIGNMENT_LOW_READ_UNIQUENESS_THRESHOLD = 10;

    /**
     * <p>
     * A filter that is used to remove contigs upfront which don't meet the following criteria
     * either:
     * has only 1 mapping, with MQ strictly above this threshold
     * or:
     * has more than 1 mapping, but only 1 mapping has MQ strictly above this threshold and it has a large gap in it.
     * </p>
     *
     * <p>
     * This value is also used for filtering out alignments that has low mapping quality in
     * {@link #removeDueToLowMQ(List, int, List, List)}
     * </p>
     */
    static final int ALIGNMENT_MQ_THRESHOLD = 20;

    /**
     * Default value to be passed to
     * {@link AlignedContig#filterSecondaryConfigurationsByMappingQualityThreshold(List, int)}
     */
    static final int SECONDARY_CONFIGURATION_MQ_FILTER_THRESHOLD = 0;

    // normalization factor by which all alignments will be divided to compute a MQ-based weight.
    static final double COVERAGE_MQ_NORMALIZATION_CONST = 60.0;

    /**
     * A filter to boost configuration scoring implemented here:
     * if the configuration has more than 10 mappings, then
     * any mappings in such configuration with MQ
     * not strictly above this threshold is classified as bad and filtered.
     */
    static final int ALIGNMENT_MQ_THRESHOLD_FOR_SPEED_BOOST = 10;

    private final String contigName;
    private final byte[] contigSequence;
    private final List<AlignmentInterval> alignmentIntervals;

    // throws if alignment interval is null
    public AlignedContig( final String contigName, final byte[] contigSequence, final List<AlignmentInterval> alignmentIntervals ) {
        if ( alignmentIntervals == null ) {
            throw new IllegalArgumentException("AlignedContig being constructed with null alignments: " + contigName);
        }
        this.contigName = contigName;
        this.contigSequence = contigSequence;
        this.alignmentIntervals = alignmentIntervals.stream().sorted(getAlignmentIntervalComparator()).collect(Collectors.toList());
    }

    AlignedContig( final Kryo kryo, final Input input ) {

        contigName = input.readString();

        final int nBases = input.readInt();
        contigSequence = new byte[nBases];
        for ( int b = 0; b < nBases; ++b ) {
            contigSequence[b] = input.readByte();
        }

        final int nAlignments = input.readInt();
        alignmentIntervals = new ArrayList<>(nAlignments);
        for ( int i = 0; i < nAlignments; ++i ) {
            alignmentIntervals.add(new AlignmentInterval(kryo, input));
        }
    }

    /**
     * Iterates through the input {@code noSecondaryAlignments}, which are assumed to contain no secondary alignment (i.e. records with "XA" tag),
     * converts to custom {@link AlignmentInterval} format and
     * split the records when the gap in the alignment reaches the specified {@code sensitivity}.
     * The size of the returned iterable of {@link AlignmentInterval}'s is guaranteed to be no lower than that of the input iterable.
     */
    @VisibleForTesting
    public static AlignedContig parseReadsAndOptionallySplitGappedAlignments(final Iterable<SAMRecord> noSecondaryAlignments,
                                                                             final int gapSplitSensitivity,
                                                                             final boolean splitGapped) {

        Utils.validateArg(noSecondaryAlignments.iterator().hasNext(), "input collection of GATK reads is empty");

        final SAMRecord primaryAlignment
                = Utils.stream(noSecondaryAlignments).filter(sam -> !sam.getSupplementaryAlignmentFlag())
                .findFirst()
                .orElseThrow(() -> new GATKException("no primary alignment for read " + noSecondaryAlignments.iterator().next().getReadName()));

        Utils.validate(!primaryAlignment.getCigar().containsOperator(CigarOperator.H),
                "assumption that primary alignment does not contain hard clipping is invalid for read: " + primaryAlignment.toString());

        final byte[] contigSequence = primaryAlignment.getReadBases().clone();
        final List<AlignmentInterval> parsedAlignments;
        if ( primaryAlignment.getReadUnmappedFlag() ) { // the Cigar
            parsedAlignments = Collections.emptyList();
        } else {
            if (primaryAlignment.getReadNegativeStrandFlag()) {
                SequenceUtil.reverseComplement(contigSequence);
            }

            final Stream<AlignmentInterval> unSplitAIList = Utils.stream(noSecondaryAlignments).map(AlignmentInterval::new);
            if (splitGapped) {
                final int unClippedContigLength = primaryAlignment.getReadLength();
                parsedAlignments = unSplitAIList.map(ar ->
                        ContigAlignmentsModifier.splitGappedAlignment(ar, gapSplitSensitivity, unClippedContigLength))
                        .flatMap(Utils::stream).collect(Collectors.toList());
            } else {
                parsedAlignments = unSplitAIList.collect(Collectors.toList());
            }
        }
        return new AlignedContig(primaryAlignment.getReadName(), contigSequence, parsedAlignments);
    }

    boolean hasOnly2Alignments() {
        return alignmentIntervals.size() == 2;
    }

    /**
     * @return first alignment of the contig
     */
    public AlignmentInterval getHeadAlignment() {
        return alignmentIntervals.get(0);
    }

    /**
     * @return last alignment of the contig
     */
    public AlignmentInterval getTailAlignment() {
        return alignmentIntervals.get(alignmentIntervals.size() - 1);
    }

    public String getContigName() {
        return contigName;
    }

    public byte[] getContigSequence() {
        return contigSequence;
    }

    public boolean isUnmapped() {
        return alignmentIntervals.isEmpty();
    }

    public List<AlignmentInterval> getAlignments() {
        return alignmentIntervals;
    }


    /**
     * Idea is to keep mapped contig that
     * either has at least two alignments over {@link #ALIGNMENT_MQ_THRESHOLD},
     * or in the case of a single alignment, it must be MQ > {@link #ALIGNMENT_MQ_THRESHOLD}.
     * Note that we are not simply filtering out contigs with only 1 alignment because
     * they might contain large (> 50) gaps hence should be kept.
     * <p>
     * a point that could use improvements:
     * the current implementation exhaustively checks the power set of all possible alignments of each assembly contig,
     * which is computationally impossible for contigs having many-but-barely-any-good alignments, yet bringing in no value,
     * hence this primitive filtering step to get rid of these bad assembly contigs.
     */
    public boolean hasGoodMQ() {
        if ( alignmentIntervals.size() < 2 ) {
            return (!alignmentIntervals.isEmpty()) && alignmentIntervals.get(0).mapQual > ALIGNMENT_MQ_THRESHOLD;
        } else {
            int notBadMappingsCount = 0; // either more than 1 non-bad mappings, or at least 1 non-bad mapping containing large gap
            for ( final AlignmentInterval alignment : alignmentIntervals ) {
                if ( alignment.mapQual > ALIGNMENT_MQ_THRESHOLD ) {
                    if ( alignment.containsGapOfEqualOrLargerSize(GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY) ) {
                        return true;// early return when a not-bad one contains a large gap
                    } else {
                        ++notBadMappingsCount;
                    }
                }
            }
            return notBadMappingsCount > 1;
        }
    }

    /**
     * Reconstructs (possibly more than one) {@link AssemblyContigWithFineTunedAlignments} based on
     * the best-scored configuration(s).
     * <p>
     * Essentially, it
     * <ul>
     *     <li>
     *         retrieves a list of best-scoring configurations
     *     </li>
     *     <li>
     *         splits alignments that contain large gaps into multiple child alignments
     *     </li>
     *     <li>
     *         performs a final round of cleanup by reclassifying alignments in
     *         {@link GoodAndBadMappings#getGoodMappings()} that offers low uniqueness as bad mappings
     *         via {@link #removeNonUniqueMappings(GoodAndBadMappings, int, int)}.
     *         See {@link #removeNonUniqueMappings(GoodAndBadMappings, int, int)} for meaning of "low uniqueness"
     *     </li>
     * </ul>
     */
    @VisibleForTesting
    public List<AssemblyContigWithFineTunedAlignments> reconstructContigFromBestConfiguration(
            final Set<String> canonicalChromosomes,
            final double scoreDiffTolerance ) {
            final List<GoodAndBadMappings> bestConfigurations =
                    pickAndFilterConfigurations(canonicalChromosomes, scoreDiffTolerance);
        if ( bestConfigurations.size() > 1 ) { // more than one best configuration
            return bestConfigurations.stream()
                    .map(mappings -> splitGaps(mappings, false))
                    .map(mappings -> removeNonUniqueMappings(mappings, ALIGNMENT_MQ_THRESHOLD, ALIGNMENT_LOW_READ_UNIQUENESS_THRESHOLD))
                    .filter(mappings -> alignmentShouldNotBeStitchedTogether(mappings.getGoodMappings()))
                    .map(mappings -> createContigGivenClassifiedAlignments(contigName, contigSequence, mappings, true))
                    .sorted(getConfigurationComparator()).collect(Collectors.toList());
        } else {
            final GoodAndBadMappings intermediate = splitGaps(bestConfigurations.get(0), false);
            final GoodAndBadMappings result = removeNonUniqueMappings(intermediate, ALIGNMENT_MQ_THRESHOLD, ALIGNMENT_LOW_READ_UNIQUENESS_THRESHOLD);
            if ( alignmentShouldNotBeStitchedTogether(result.getGoodMappings()) )
                return Collections.singletonList(createContigGivenClassifiedAlignments(contigName, contigSequence, result, false));
            else {
                return Collections.emptyList();
            }
        }
    }

    private static boolean alignmentShouldNotBeStitchedTogether( final List<AlignmentInterval> alignments ) {
        return alignments.size() != 2
                ||
                !simpleChimeraWithStichableAlignments(alignments.get(0), alignments.get(1));
    }

    private static AssemblyContigWithFineTunedAlignments createContigGivenClassifiedAlignments( final String contigName, final byte[] contigSeq,
                                                                                                final GoodAndBadMappings goodAndBadMappings,
                                                                                                final boolean setResultContigAsAmbiguous ) {

        return new AssemblyContigWithFineTunedAlignments(
                new AlignedContig(contigName, contigSeq, goodAndBadMappings.getGoodMappings()),
                goodAndBadMappings.getBadMappings().stream().map(AlignmentInterval::toPackedString).collect(Collectors.toList()),
                setResultContigAsAmbiguous,
                goodAndBadMappings.getMayBeNullGoodMappingToNonCanonicalChromosome());
    }

    /**
     * when two configurations are the same,
     * implement ordering that
     * prefers the configuration with less alignments,
     * and prefers the configuration with a lower number of summed mismatches in case of a tie
     */
    @VisibleForTesting
    static Comparator<AssemblyContigWithFineTunedAlignments> getConfigurationComparator() {
        final Comparator<AssemblyContigWithFineTunedAlignments> numFirst
                = Comparator.comparingInt(tig -> tig.getAlignments().size());
        final Comparator<AssemblyContigWithFineTunedAlignments> mismatchSecond
                = ( AssemblyContigWithFineTunedAlignments x, AssemblyContigWithFineTunedAlignments y )
                -> Integer.compare(x.getAlignments().stream().mapToInt(ai -> ai.mismatches).sum(),
                y.getAlignments().stream().mapToInt(ai -> ai.mismatches).sum());
        return numFirst.thenComparing(mismatchSecond);
    }

    /**
     * Split the good alignments stored in {@code configuration}.
     * <p>
     * Note that an edge case is possible where for some gapped alignment,
     * after gap-split, one or more child alignments might be contained by
     * another original alignment.
     * <p>
     * Example:
     * picture: ref span by gapped original alignment   -------------          --------
     * ref span by overlapping other alignment                  ---------------------------
     * <p>
     * Here the method offers two options:
     * 1) keep doing the gap split, and compare if the gapped original alignment offers
     * better read coverage/breadth compared to the other, if so, keep the children alignments from the gap split,
     * and ditch the other alignment, otherwise ditch the alignment gapped original alignment,
     * the aim is to keep the split children alignments together
     * 2) drop the gap-split child alignments whose read spans are contained in other alignment spans
     * TODO: based on evaluation done on 2018-06-30, making either choice has small effect on the final call set; one could further evaluate options when making improvements.
     */
    @VisibleForTesting
    static GoodAndBadMappings splitGaps( final GoodAndBadMappings configuration, final boolean keepSplitChildrenTogether ) {
        return keepSplitChildrenTogether ? splitGapsAndKeepChildrenTogether(configuration) : splitGapsAndDropAlignmentContainedByOtherOnRead(configuration);
    }

    /**
     * See {@link #splitGaps(GoodAndBadMappings, boolean)}.
     * This implementation keeps the split children together
     */
    @VisibleForTesting
    static GoodAndBadMappings splitGapsAndKeepChildrenTogether( final GoodAndBadMappings configuration ) {

        final List<AlignmentInterval> originalBadMappings = configuration.getBadMappings();
        final List<AlignmentInterval> scan = configuration.getGoodMappings();

        // 1st pass, split gapped alignments when available, and all defaults to good
        final List<Tuple2<Boolean, Iterable<AlignmentInterval>>> alignmentSplitChildren =
                scan.stream()
                        .map(alignment -> {
                            final Iterable<AlignmentInterval> split;
                            if ( alignment.containsGapOfEqualOrLargerSize(GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY) ) {
                                split = ContigAlignmentsModifier.splitGappedAlignment(alignment, GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY,
                                        CigarUtils.countUnclippedReadBases(alignment.cigarAlong5to3DirectionOfContig));
                            } else {
                                split = Collections.singletonList(alignment);
                            }
                            return new Tuple2<>(true, split);
                        }).collect(Collectors.toList());

        // 2nd pass make a choice between gapped and overlapping alignment (alignments that are not favored has its "2nd" set to false)
        final int count = scan.size();
        for ( int i = 0; i < count; ++i ) {
            final AlignmentInterval alignment = scan.get(i);
            final Tuple2<Boolean, Iterable<AlignmentInterval>> split = alignmentSplitChildren.get(i);
            if ( split._1 && Iterables.size(split._2) != 1 ) { // the split could be already marked bad (i.e. to be filtered out), or could contain no large gaps (i.e. should not check it)
                for ( int j = 0; j < count; ++j ) {
                    final AlignmentInterval other = scan.get(j);
                    if ( j == i || AlignmentInterval.overlapOnContig(alignment, other) == 0 )
                        continue;

                    if ( Utils.stream(split._2).anyMatch(other::containsOnRead) ) {
                        if ( gappedAlignmentOffersBetterCoverage(alignment, other) ) {
                            final Iterable<AlignmentInterval> copy = alignmentSplitChildren.get(j)._2;
                            alignmentSplitChildren.set(j, new Tuple2<>(false, copy));
                        } else {
                            final Iterable<AlignmentInterval> copy = alignmentSplitChildren.get(i)._2;
                            alignmentSplitChildren.set(i, new Tuple2<>(false, copy));
                        }
                    }
                }
            }
        }

        final List<AlignmentInterval> good = new ArrayList<>();
        final List<AlignmentInterval> bad = new ArrayList<>(originalBadMappings);
        for ( final Tuple2<Boolean, Iterable<AlignmentInterval>> pair : alignmentSplitChildren ) {
            if ( pair._1 ) {
                good.addAll(Lists.newArrayList(pair._2));
            } else {
                bad.addAll(Lists.newArrayList(pair._2));
            }
        }
        good.sort(getAlignmentIntervalComparator());
        return new GoodAndBadMappings(good, bad, configuration.getMayBeNullGoodMappingToNonCanonicalChromosome());
    }

    /**
     * See {@link #splitGaps(GoodAndBadMappings, boolean)}.
     * This implementation drops alignment that are contained by others, in terms of their read span.
     */
    @VisibleForTesting
    static GoodAndBadMappings splitGapsAndDropAlignmentContainedByOtherOnRead( final GoodAndBadMappings configuration ) {
        final List<AlignmentInterval> originalGoodMappings = configuration.getGoodMappings();
        final List<AlignmentInterval> gapSplit = new ArrayList<>(originalGoodMappings.size());
        originalGoodMappings.forEach(alignment -> {
            if ( alignment.containsGapOfEqualOrLargerSize(GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY) ) {
                ContigAlignmentsModifier.splitGappedAlignment(alignment, GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY,
                        CigarUtils.countUnclippedReadBases(alignment.cigarAlong5to3DirectionOfContig))
                        .forEach(gapSplit::add);
            } else {
                gapSplit.add(alignment);
            }
        });
        final int count = gapSplit.size();
        final List<AlignmentInterval> bad = new ArrayList<>(configuration.getBadMappings());
        for ( int i = 0; i < count; ++i ) {
            final AlignmentInterval one = gapSplit.get(i);
            for ( int j = i + 1; j < count; ++j ) {
                final AlignmentInterval two = gapSplit.get(j);
                if ( one.containsOnRead(two) ) {
                    bad.add(two);
                } else if ( two.containsOnRead(one) ) {
                    bad.add(one);
                }
            }
        }
        gapSplit.removeAll(bad);
        gapSplit.sort(getAlignmentIntervalComparator());
        return new GoodAndBadMappings(gapSplit, bad, configuration.getMayBeNullGoodMappingToNonCanonicalChromosome());
    }

    @VisibleForTesting
    static boolean gappedAlignmentOffersBetterCoverage( final AlignmentInterval gapped,
                                                        final AlignmentInterval overlappingNonGapped ) {
        final int diff = gapped.getSizeOnRead() - overlappingNonGapped.getSizeOnRead();
        if ( diff == 0 ) {
            return gapped.alnScore > overlappingNonGapped.alnScore;
        } else {
            return diff > 0;
        }
    }

    /**
     * Process provided {@code originalConfiguration} of an assembly contig and split between good and bad alignments.
     *
     * <p>
     * What is considered good and bad?
     * For a particular mapping/alignment, it may offer low uniqueness in two sense:
     *     <ul>
     *         <li>
     *             low REFERENCE UNIQUENESS: meaning the sequence being mapped match multiple locations on the reference;
     *         </li>
     *         <li>
     *             low READ UNIQUENESS: with only a very short part of the read uniquely explained by this particular alignment;
     *         </li>
     *     </ul>
     *     Good alignments offer both high reference uniqueness and read uniqueness, as judged by the requested
     *     {@code mapQThresholdInclusive} and {@code uniqReadLenInclusive}
     *     (yes we are doing a hard-filtering but more advanced model is not our priority right now 2017-11-20).
     * </p>
     * <p>
     * Note that "original" is meant to be possibly different from the returned configuration,
     * but DOES NOT mean the alignments of the contig as given by the aligner,
     * i.e. the configuration should be one of the best given by
     * {@link AlignedContig#pickBestConfigurations(Set, double)}.
     */
    @VisibleForTesting
    static GoodAndBadMappings removeNonUniqueMappings( final GoodAndBadMappings goodAndBadMappings,
                                                       final int mapQThresholdInclusive,
                                                       final int uniqReadLenInclusive ) {
        final List<AlignmentInterval> inputAlignments = goodAndBadMappings.getGoodMappings();
        if ( inputAlignments.size() <= 2 ) // TODO: 6/27/18 does this need to be change to < 2?
            return goodAndBadMappings;

        // two pass, each focusing on removing the alignments of a contig that offers low uniqueness in one sense:
        final List<AlignmentInterval> selectedAlignments = new ArrayList<>(inputAlignments.size());
        final List<AlignmentInterval> lowUniquenessMappings = new ArrayList<>(goodAndBadMappings.getBadMappings());

        // first pass is for removing alignments with low REFERENCE UNIQUENESS, using low mapping quality as the criterion
        removeDueToLowMQ(inputAlignments, mapQThresholdInclusive, selectedAlignments, lowUniquenessMappings);

        // second pass, the slower one, is to remove alignments offering low READ UNIQUENESS,
        // i.e. with only a very short part of the read being uniquely explained by this particular alignment;
        removeDueToShortReadSpan(selectedAlignments, uniqReadLenInclusive, lowUniquenessMappings);


        return new GoodAndBadMappings(selectedAlignments, lowUniquenessMappings, goodAndBadMappings.getMayBeNullGoodMappingToNonCanonicalChromosome());
    }

    // classify alignments with MQ strictly lower than {@code mapQThresholdInclusive} as bad.
    private static void removeDueToLowMQ( final List<AlignmentInterval> inputAlignments, final int mapQThresholdInclusive,
                                          final List<AlignmentInterval> selectedAlignments, final List<AlignmentInterval> lowUniquenessMappings ) {

        for ( final AlignmentInterval alignment : inputAlignments ) {
            if ( alignment.mapQual >= mapQThresholdInclusive )
                selectedAlignments.add(alignment);
            else
                lowUniquenessMappings.add(alignment);
        }
    }

    // classify alignments offering unique read span strictly shorter than {@code uniqReadLenInclusive} as bad.
    private static void removeDueToShortReadSpan( final List<AlignmentInterval> selectedAlignments, final int uniqReadLenInclusive,
                                                  final List<AlignmentInterval> lowUniquenessMappings ) {
        // the steps are:
        //      search bi-directionally until cannot find overlap any more,
        //      subtract the overlap from the distance covered on the contig by the alignment.
        //      This gives unique read region it explains.
        //      If this unique read region is "short": shorter than {@code uniqReadLenInclusive}), drop it.

        // each alignment has an entry of a tuple2, one for max overlap maxFront, one for max overlap maxRear,
        // max overlap maxFront is a tuple2 registering the index and overlap bases count
        final Map<AlignmentInterval, Tuple2<Integer, Integer>> maxOverlapMap = getMaxOverlapPairs(selectedAlignments);
        for ( Iterator<AlignmentInterval> iterator = selectedAlignments.iterator(); iterator.hasNext(); ) {
            final AlignmentInterval alignment = iterator.next();

            final Tuple2<Integer, Integer> maxOverlapFrontAndRear = maxOverlapMap.get(alignment);
            final int maxOverlapFront = Math.max(0, maxOverlapFrontAndRear._1);
            final int maxOverlapRear = Math.max(0, maxOverlapFrontAndRear._2);

            // theoretically this could be negative for an alignment whose maxFront and maxRear sums together bigger than the read span
            // but earlier configuration scoring would make this impossible because such alignments should be filtered out already
            // considering that it brings more penalty than value, i.e. read bases explained (even if the MQ is 60),
            // but even if it is kept, a negative value won't hurt unless a stupid threshold value is passed in
            final int uniqReadSpan = alignment.endInAssembledContig - alignment.startInAssembledContig + 1
                    - maxOverlapFront - maxOverlapRear;
            if ( uniqReadSpan < uniqReadLenInclusive ) {
                lowUniquenessMappings.add(alignment);
                iterator.remove();
            }
        }
    }

    /**
     * Extract the max overlap information, front and back, for each alignment in {@code configuration}.
     * For each alignment, the corresponding tuple2 has the max (front, rear) overlap base counts.
     */
    @VisibleForTesting
    static Map<AlignmentInterval, Tuple2<Integer, Integer>> getMaxOverlapPairs( final List<AlignmentInterval> configuration ) {

        final List<TempMaxOverlapInfo> intermediateResult =
                new ArrayList<>(Collections.nCopies(configuration.size(), new TempMaxOverlapInfo()));

        // We iterate through all alignments except the last one
        // For the last alignment, which naturally doesn't have any maxRear,
        //     the following implementation sets its maxFront during the iteration, if available at all (it may overlap with nothing)
        for ( int i = 0; i < configuration.size() - 1; ++i ) {

            final AlignmentInterval cur = configuration.get(i);
            // For the i-th alignment, we only look at alignments after it (note j starts from i+1) and find max overlap
            int maxOverlapRearBases = -1;
            int maxOverlapRearIndex = -1;
            for ( int j = i + 1; j < configuration.size(); ++j ) { // note j > i
                final int overlap = AlignmentInterval.overlapOnContig(cur, configuration.get(j));
                if ( overlap > maxOverlapRearBases ) {
                    maxOverlapRearBases = overlap;
                    maxOverlapRearIndex = j;
                } else { // following ones, as guaranteed by the ordering of alignments in the contig, cannot overlap
                    break;
                }
            }

            if ( maxOverlapRearBases > 0 ) {
                // for current alignment (i-th), set its max_overlap_rear, which would not change in later iterations and copy old max_overlap_front
                final Tuple2<Integer, Integer> maxRear = new Tuple2<>(maxOverlapRearIndex, maxOverlapRearBases);
                final Tuple2<Integer, Integer> maxFrontToCopy = intermediateResult.get(i).maxFront;
                intermediateResult.set(i, new TempMaxOverlapInfo(maxFrontToCopy, maxRear));

                // then conditionally set the max_overlap_front of the
                // maxOverlapRearIndex-th alignment
                // that maximally overlaps with the current, i.e. i-th, alignment
                final TempMaxOverlapInfo oldValue = intermediateResult.get(maxOverlapRearIndex);// maxOverlapRearIndex cannot be -1 here
                if ( oldValue.maxFront._2 < maxOverlapRearBases )
                    intermediateResult.set(maxOverlapRearIndex, new TempMaxOverlapInfo(new Tuple2<>(i, maxOverlapRearBases), oldValue.maxRear));
            }
        }

        final Map<AlignmentInterval, Tuple2<Integer, Integer>> maxOverlapMap = new HashMap<>(configuration.size());
        for ( int i = 0; i < configuration.size(); ++i ) {
            maxOverlapMap.put(configuration.get(i),
                    new Tuple2<>(intermediateResult.get(i).maxFront._2, intermediateResult.get(i).maxRear._2));
        }

        return maxOverlapMap;
    }

    /**
     * See the funny alignment signature described in ticket 4951 on GATK github
     *
     * @param intervalOne assumed to start no later   than {@code intervalTwo} on the read
     * @param intervalTwo assumed to start no earlier than {@code intervalOne} on the read
     * @return true if the two given intervals can be stitched together
     * @throws IllegalArgumentException if the two intervals are not sorted according to their {@link AlignmentInterval#startInAssembledContig}
     */
    public static boolean simpleChimeraWithStichableAlignments( final AlignmentInterval intervalOne, final AlignmentInterval intervalTwo ) {
        if ( intervalOne.startInAssembledContig > intervalTwo.startInAssembledContig )
            throw new IllegalArgumentException("Assumption that input intervals are sorted by their starts on read is violated.\tFirst: " +
                    intervalOne.toPackedString() + "\tSecond: " + intervalTwo.toPackedString());
        if ( !intervalOne.referenceSpan.getContig().equals(intervalTwo.referenceSpan.getContig()) )
            return false;
        if ( intervalOne.forwardStrand != intervalTwo.forwardStrand )
            return false;
        if ( intervalOne.containsOnRead(intervalTwo) || intervalTwo.containsOnRead(intervalOne) )
            return false;
        if ( intervalOne.containsOnRef(intervalTwo) || intervalTwo.containsOnRef(intervalOne) )
            return false;
        final boolean refOrderSwap = intervalOne.forwardStrand != (intervalOne.referenceSpan.getStart() < intervalTwo.referenceSpan.getStart());
        if ( refOrderSwap )
            return false;
        final int overlapOnContig = AlignmentInterval.overlapOnContig(intervalOne, intervalTwo);
        final int overlapOnRefSpan = AlignmentInterval.overlapOnRefSpan(intervalOne, intervalTwo);
        if ( overlapOnContig == 0 && overlapOnRefSpan == 0 ) {
            final boolean canBeStitched = intervalTwo.referenceSpan.getStart() - intervalOne.referenceSpan.getEnd() == 1
                    && intervalTwo.startInAssembledContig - intervalOne.endInAssembledContig == 1;
            return canBeStitched;
        } else
            return overlapOnContig == overlapOnRefSpan;
    }

    public List<GoodAndBadMappings> pickAndFilterConfigurations( final Set<String> canonicalChromosomes,
                                                                 final double scoreDiffTolerance ) {
        return filterSecondaryConfigurationsByMappingQualityThreshold(
                pickBestConfigurations(canonicalChromosomes, scoreDiffTolerance),
                SECONDARY_CONFIGURATION_MQ_FILTER_THRESHOLD);
    }

    /**
     * Pick the best configurations based on a heuristic scoring scheme implemented in
     * {@link #computeScoreOfConfiguration(List, Set, int)}.
     *
     * <p>
     * This step itself has its own hard-filtering steps done
     * via {@link #heuristicSpeedUpWhenFacingManyMappings(Set, int)}.
     * </p>
     *
     * <p>
     * Before all possible configurations are scored, the alignments are analyzed
     * via {@link #getBetterNonCanonicalMapping(Set, List, int)} to annotate the
     * contig if a good mapping to non-canonical chromosomes exist.
     * See detailed explanation in {@link #getBetterNonCanonicalMapping(Set, List, int)}.
     * </p>
     *
     * @return a 2-D list, where in the case when multiple configurations are equally top-scored, all such configurations are picked up
     */
    @VisibleForTesting
    public List<GoodAndBadMappings> pickBestConfigurations( final Set<String> canonicalChromosomes,
                                                            final double scoreDiffTolerance ) {
        // nothing to score if only one alignment
        if ( alignmentIntervals.size() == 1 ) {
            return Collections.singletonList(
                    new GoodAndBadMappings(Collections.singletonList(alignmentIntervals.get(0)))
            );
        }

        // step 1: get max aligner score of mappings to canonical chromosomes and speed up in case of too many mappings
        final int maxCanonicalChrAlignerScore = alignmentIntervals.stream()
                .filter(alignmentInterval -> canonicalChromosomes.contains(alignmentInterval.referenceSpan.getContig()))
                .mapToInt(ai -> ai.alnScore).max().orElse(0); // possible that no mapping to canonical chromosomes

        final GoodAndBadMappings preFilteredAlignments =
                heuristicSpeedUpWhenFacingManyMappings(canonicalChromosomes, maxCanonicalChrAlignerScore);
        final List<AlignmentInterval> goodMappings = preFilteredAlignments.goodMappings;
        final List<AlignmentInterval> badMappings = preFilteredAlignments.badMappings;

        final int newMaxCanonicalChrAlignerScore = goodMappings.stream()
                .filter(alignmentInterval -> canonicalChromosomes.contains(alignmentInterval.referenceSpan.getContig()))
                .mapToInt(ai -> ai.alnScore).max().orElse(0); // possible that no mapping to canonical chromosomes

        // annotate contig if a good mapping to non-canonical chromosome exists
        final AlignmentInterval goodMappingToNonCanonicalChromosome =
                getBetterNonCanonicalMapping(canonicalChromosomes, goodMappings, newMaxCanonicalChrAlignerScore);
        if ( goodMappingToNonCanonicalChromosome != null ) { // take it out of consideration
            goodMappings.remove(goodMappingToNonCanonicalChromosome);
        }

        // step 2: generate, score, and pick configurations
        return generateScoreAndPickConfigurations(goodMappings, badMappings, goodMappingToNonCanonicalChromosome, canonicalChromosomes,
                newMaxCanonicalChrAlignerScore, scoreDiffTolerance, contigName);
    }

    // if multiple configurations have equally good scores, return all of them
    private static List<GoodAndBadMappings> generateScoreAndPickConfigurations( final List<AlignmentInterval> goodMappings,
                                                                                final List<AlignmentInterval> badMappings,
                                                                                final AlignmentInterval goodMappingToNonCanonicalChromosome,
                                                                                final Set<String> canonicalChromosomes,
                                                                                final int maxCanonicalChrAlignerScore,
                                                                                final double scoreDiffTolerance,
                                                                                final String contigName ) {
        final List<List<AlignmentInterval>> allConfigurations = Sets.powerSet(new HashSet<>(goodMappings))
                .stream().map(ArrayList::new)
                // make sure within each configuration, alignments would be sorted as they would be in a corresponding AlignedContig
                .map(ls -> ls.stream().sorted(getAlignmentIntervalComparator()).collect(Collectors.toList()))
                .collect(Collectors.toList());

        final List<Double> scores = allConfigurations.stream()
                .map(configuration -> computeScoreOfConfiguration(configuration, canonicalChromosomes, maxCanonicalChrAlignerScore))
                .collect(SVUtils.arrayListCollector(allConfigurations.size()));

        // step 3: pick the best-scored configuration(s) (if multiple configurations have equally good scores, return all of them)
        final double maxScore = scores.stream().mapToDouble(Double::doubleValue).max()
                .orElseThrow(() -> new GATKException("Cannot find best-scoring configuration on alignments of contig: " + contigName));

        return IntStream.range(0, allConfigurations.size())
                .filter(i -> {
                    final double s = scores.get(i);
                    // two configurations with would-be-same-scores can differ by a tolerance due to the sin of comparing
                    // could-be-close floating point values
                    // (see http://www.cygnus-software.com/papers/comparingfloats/Comparing%20floating%20point%20numbers.htm)
                    final double tol = Math.max(Math.ulp(s), scoreDiffTolerance);
                    return s >= maxScore || maxScore - s <= tol;
                })
                .mapToObj(p -> {
                    final ArrayList<AlignmentInterval> copy = new ArrayList<>(goodMappings);
                    final List<AlignmentInterval> pickedAlignments = allConfigurations.get(p);
                    copy.removeAll(pickedAlignments); // remove picked, left are bad
                    copy.addAll(badMappings); // add original bad mappings
                    return new GoodAndBadMappings(pickedAlignments, copy, goodMappingToNonCanonicalChromosome);
                })
                .collect(Collectors.toList());
    }

    /**
     * Speed up if number of alignments is too high (>10):
     * <ul>
     *     <li>if mapped to canonical chromosomes, MQ must be >10;</li>
     *     <li>otherwise, must have AS higher than max canonical aligner score</li>
     * </ul>
     *
     * @param canonicalChromosomes        a set of chromosome names that are defined as canonical, e.g. for Human, chr1-chr22, and chrX and chrY
     * @param maxCanonicalChrAlignerScore among the mappings of {@code alignedContig} that map to canonical chromosomes, the maximum of the aligner scores
     * @return {@link GoodAndBadMappings} with bad mappings defined as those failing hard filters mentioned above.
     */
    @VisibleForTesting
    GoodAndBadMappings heuristicSpeedUpWhenFacingManyMappings( final Set<String> canonicalChromosomes,
                                                               final int maxCanonicalChrAlignerScore ) {
        final List<AlignmentInterval> goods;
        final List<AlignmentInterval> bads;
        if ( alignmentIntervals.size() > 10 ) {
            goods = new ArrayList<>();
            bads = new ArrayList<>();
            for ( final AlignmentInterval alignment : alignmentIntervals ) {
                final boolean isGood = (!canonicalChromosomes.contains(alignment.referenceSpan.getContig()) && alignment.alnScore > maxCanonicalChrAlignerScore)
                        || alignment.mapQual > ALIGNMENT_MQ_THRESHOLD_FOR_SPEED_BOOST;
                if ( isGood )
                    goods.add(alignment);
                else
                    bads.add(alignment);
            }
        } else {
            goods = alignmentIntervals;
            bads = Collections.emptyList();
        }
        return new GoodAndBadMappings(goods, bads);
    }

    /**
     * Some non-canonical reference chromosomes represent alternate haplotypes which are similar to
     * sequence represented in the canonical chromosomes but contain structural rearrangements,
     * including deletions or duplications.
     * In other words, the non-canonical chromosome captures an SV--relative to the canonical chromosomes--of
     * relatively high population frequency.
     * <p>
     * The sample under analysis could have the allele of this non-canonical version,
     * and be marked as having an SV on the corresponding location on the canonical chromosome.
     * <p>
     * An assembly contig from this sample may have two equally well scored alignment configurations, where
     * one configuration has split alignments to canonical chromosomes, hence indicating the SV, whereas
     * the other configuration has a single, often very good (or even better) alignment to a non-canonical chromosome.
     * We send down the chimeric alignment configuration for inference, but note down that a non-canonical chromosome
     * in the reference input could have already captured the SV on this sample.
     *
     * @return {@code null} if the non-canonical chromosome mapping doesn't offer a better score,
     * otherwise the non-canonical chromosome mapping
     */
    @VisibleForTesting
    static AlignmentInterval getBetterNonCanonicalMapping( final Set<String> canonicalChromosomes,
                                                           final List<AlignmentInterval> goodMappings,
                                                           final int maxCanonicalChrAlignerScore ) {
        final List<AlignmentInterval> canonicalMappings = new ArrayList<>(goodMappings.size());
        final List<AlignmentInterval> nonCanonicalMapping = new ArrayList<>();
        for ( final AlignmentInterval alignment : goodMappings ) {
            if ( canonicalChromosomes.contains(alignment.referenceSpan.getContig()) )
                canonicalMappings.add(alignment);
            else
                nonCanonicalMapping.add(alignment);
        }
        if ( canonicalMappings.isEmpty() )
            return null;
        if ( nonCanonicalMapping.size() == 1 &&
                (canonicalMappings.size() > 1 || canonicalMappings.get(0).containsGapOfEqualOrLargerSize(GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY)) ) {
            final double canonicalScore = computeScoreOfConfiguration(canonicalMappings, canonicalChromosomes, maxCanonicalChrAlignerScore);
            final double nonCanonicalScore = computeScoreOfConfiguration(nonCanonicalMapping, canonicalChromosomes, maxCanonicalChrAlignerScore);
            return (canonicalScore > nonCanonicalScore) ? null : nonCanonicalMapping.get(0);
        } else {
            return null;
        }
    }

    /**
     * Computing score of given configuration of alignments.
     * No assumption on the ordering of input alignments.
     */
    @VisibleForTesting
    static double computeScoreOfConfiguration( final List<AlignmentInterval> configuration,
                                               final Set<String> canonicalChromosomes,
                                               final int maxCanonicalChrAlignerScore ) {

        final double tigExplainQual = computeTigExplainQualOfOneConfiguration(configuration, canonicalChromosomes, maxCanonicalChrAlignerScore);

        int redundancy = 0;
        for ( int i = 0; i < configuration.size() - 1; ++i ) {
            for ( int j = i + 1; j < configuration.size(); ++j ) {
                final int overlap = AlignmentInterval.overlapOnContig(configuration.get(i), configuration.get(j));
                redundancy += overlap;
            }
        }

        return tigExplainQual - redundancy;
    }

    private static double computeTigExplainQualOfOneConfiguration( final List<AlignmentInterval> configuration,
                                                                   final Set<String> canonicalChromosomes,
                                                                   final int maxCanonicalChrAlignerScore ) {
        double tigExplainedQual = 0;
        for ( final AlignmentInterval alignmentInterval : configuration ) {
            final int len = alignmentInterval.getSizeOnRead();
            final double weight;
            if ( canonicalChromosomes.contains(alignmentInterval.referenceSpan.getContig()) ) {
                weight = alignmentInterval.mapQual / COVERAGE_MQ_NORMALIZATION_CONST;
            } else {
                weight = Math.max(alignmentInterval.mapQual / COVERAGE_MQ_NORMALIZATION_CONST,
                        alignmentInterval.alnScore > maxCanonicalChrAlignerScore ? 1 : 0);
            }
            tigExplainedQual += weight * len;
        }
        return tigExplainedQual;
    }

    /**
     * For contigs with more than 1 best-scored configurations as determined by
     * {@link AlignedContig#pickBestConfigurations(Set, double)},
     * save the contigs that have one and only one configuration that
     * has all mapping quality strictly above the specified {@code threshold}.
     * Example:
     * if a contig has two "optimal" configurations with MQ's {10, 60, 60}, and {60, 60},
     * this function will favor/pick the {60, 60} configuration if the {@code threshold} is 10,
     * hence breaking the degeneracy;
     * on the other hand if the {@code mqThreshold} is passed in as any value below 10, say 0,
     * then this function returns all both original {@code differentConfigurationsForOneContig}s
     */
    @VisibleForTesting
    static List<GoodAndBadMappings> filterSecondaryConfigurationsByMappingQualityThreshold(
            final List<GoodAndBadMappings> differentConfigurationsForOneContig,
            final int mqThreshold ) {

        if ( differentConfigurationsForOneContig.size() == 1 ) {
            return differentConfigurationsForOneContig;
        } else {
            final List<GoodAndBadMappings> configurationsWithMappingAboveMQThreshold =
                    Utils.stream(differentConfigurationsForOneContig)
                            .filter(rep -> rep.getGoodMappings().stream().mapToInt(ai -> ai.mapQual).min().orElse(mqThreshold) > mqThreshold)
                            .collect(Collectors.toList());
            if ( configurationsWithMappingAboveMQThreshold.size() != 1 ) {
                return differentConfigurationsForOneContig;
            } else {
                return configurationsWithMappingAboveMQThreshold;
            }
        }
    }

    @Override
    public String toString() {
        return "(" + contigName
                + ", "
                + alignmentIntervals.stream().map(AlignmentInterval::toPackedString).collect(Collectors.toList())
                + ")";
    }

    @Override
    public boolean equals( final Object o ) {
        if ( this == o ) return true;
        if ( o == null || getClass() != o.getClass() ) return false;

        final AlignedContig that = (AlignedContig) o;

        if ( !contigName.equals(that.contigName) ) return false;
        if ( !Arrays.equals(contigSequence, that.contigSequence) ) return false;
        return alignmentIntervals.equals(that.alignmentIntervals);
    }

    @Override
    public int hashCode() {
        int result = contigName.hashCode();
        result = 31 * result + Arrays.hashCode(contigSequence);
        result = 31 * result + alignmentIntervals.hashCode();
        return result;
    }

    void serialize( final Kryo kryo, final Output output ) {

        output.writeString(contigName);

        output.writeInt(contigSequence.length);
        for ( final byte base : contigSequence ) {
            output.writeByte(base);
        }

        output.writeInt(alignmentIntervals.size());
        alignmentIntervals.forEach(it -> it.serialize(kryo, output));

    }

    static Comparator<AlignmentInterval> getAlignmentIntervalComparator() {
        Comparator<AlignmentInterval> comparePos = Comparator.comparingInt(aln -> aln.startInAssembledContig);
        Comparator<AlignmentInterval> compareRefTig = Comparator.comparing(aln -> aln.referenceSpan.getContig());
        Comparator<AlignmentInterval> compareRefSpanStart = Comparator.comparingInt(aln -> aln.referenceSpan.getStart());
        return comparePos.thenComparing(compareRefTig).thenComparing(compareRefSpanStart);
    }

    /**
     * After configuration scoring and picking, the original alignments can be classified as
     * good and bad mappings:
     * good: the ones used the picked configuration
     * bad: unused alignments in the chosen configuration; these likely contain more noise than information
     * they can be turned into string representation following the format as in {@link AlignmentInterval#toPackedString()}
     * <p>
     * Note that a special case needs attention:
     * if {@link #getMayBeNullGoodMappingToNonCanonicalChromosome()} returns non-null result,
     * it is indicating an equally good--or better--non-chimeric mapping to a non-canonical chromosome exists,
     * but to preserve the SV signal, we keep the chimeric alignments to canonical chromosomes and
     * signal the situation to downstream units.
     */
    @VisibleForTesting
    public static final class GoodAndBadMappings {

        private final List<AlignmentInterval> goodMappings;
        private final List<AlignmentInterval> badMappings;
        private final AlignmentInterval goodMappingToNonCanonicalChromosome;

        public GoodAndBadMappings( @Nonnull final List<AlignmentInterval> goodMappings ) {
            this(goodMappings, Collections.emptyList(), null);
        }

        public GoodAndBadMappings( @Nonnull final List<AlignmentInterval> goodMappings, @Nonnull final List<AlignmentInterval> badMappings,
                                   final AlignmentInterval goodMappingToNonCanonicalChr ) {
            this.goodMappings = goodMappings;
            this.badMappings = badMappings;

            this.goodMappingToNonCanonicalChromosome = goodMappingToNonCanonicalChr;
        }

        public GoodAndBadMappings( @Nonnull final List<AlignmentInterval> goodMappings, @Nonnull final List<AlignmentInterval> badMappings ) {
            this(goodMappings, badMappings, null);
        }

        public List<AlignmentInterval> getGoodMappings() {
            return goodMappings;
        }

        public List<AlignmentInterval> getBadMappings() {
            return badMappings;
        }

        public AlignmentInterval getMayBeNullGoodMappingToNonCanonicalChromosome() {
            return goodMappingToNonCanonicalChromosome;
        }

        @Override
        public boolean equals( final Object o ) {
            if ( this == o ) return true;
            if ( o == null || getClass() != o.getClass() ) return false;

            final GoodAndBadMappings that = (GoodAndBadMappings) o;

            if ( !goodMappings.equals(that.goodMappings) ) return false;
            if ( !badMappings.equals(that.badMappings) ) return false;
            return Objects.equals(goodMappingToNonCanonicalChromosome, that.goodMappingToNonCanonicalChromosome);
        }

        @Override
        public int hashCode() {
            int result = goodMappings.hashCode();
            result = 31 * result + badMappings.hashCode();
            result = 31 * result + (goodMappingToNonCanonicalChromosome != null ? goodMappingToNonCanonicalChromosome.hashCode() : 0);
            return result;
        }

        @Override
        public String toString() {
            final StringBuilder sb = new StringBuilder("GoodAndBadMappings{");
            sb.append("goodMappings=").append(goodMappings);
            sb.append(", badMappings=").append(badMappings);
            if ( goodMappingToNonCanonicalChromosome != null )
                sb.append(", goodMappingToNonCanonicalChromosome=").append(goodMappingToNonCanonicalChromosome);
            sb.append('}');
            return sb.toString();
        }
    }

    /**
     * Each alignment in a specific configuration has an entry,
     * pointing to the alignments that comes before and after it,
     * that overlaps maximally (i.e. no other front or rear alignments have more overlaps)
     * with the current alignment.
     */
    private static final class TempMaxOverlapInfo {
        final Tuple2<Integer, Integer> maxFront; // 1st holds index pointing to another alignment before this, 2nd holds the count of overlapping bases
        final Tuple2<Integer, Integer> maxRear;  // same intention as above, but for alignments after this

        TempMaxOverlapInfo() {
            maxFront = new Tuple2<>(-1, -1);
            maxRear = new Tuple2<>(-1, -1);
        }

        TempMaxOverlapInfo( final Tuple2<Integer, Integer> maxFront, final Tuple2<Integer, Integer> maxRear ) {
            this.maxFront = maxFront;
            this.maxRear = maxRear;
        }
    }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<AlignedContig> {
        @Override
        public void write( final Kryo kryo, final Output output, final AlignedContig alignedContig ) {
            alignedContig.serialize(kryo, output);
        }

        @Override
        public AlignedContig read( final Kryo kryo, final Input input, final Class<AlignedContig> clazz ) {
            return new AlignedContig(kryo, input);
        }
    }
}
