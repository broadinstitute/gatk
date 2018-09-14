package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceMultiSparkSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SvDiscoveryInputMetaData;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.ContigAlignmentsModifier;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;

import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigAlignmentsSparkArgumentCollection.CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD;
import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.*;

/**
 * This deals with the special case where a contig has multiple (> 2) alignments
 * and seemingly has the complete alt haplotype assembled.
 * See criteria in {@link AssemblyContigWithFineTunedAlignments#hasIncompletePictureFromMultipleAlignments()}.
 * For cases where the contig's alignment shows signature that the assembly doesn't
 * paint the full picture we could decide to emit all BND records,
 * but that could be dealt with later.
 */
public final class CpxVariantInterpreter {

    public static final int MIN_READ_SPAN_AFTER_DEOVERLAP = 2;

    public static List<VariantContext> makeInterpretation(final JavaRDD<AssemblyContigWithFineTunedAlignments> assemblyContigs,
                                                          final SvDiscoveryInputMetaData svDiscoveryInputMetaData) {

        final Broadcast<ReferenceMultiSparkSource> referenceBroadcast = svDiscoveryInputMetaData.getReferenceData().getReferenceBroadcast();
        final Broadcast<SAMSequenceDictionary> referenceSequenceDictionaryBroadcast = svDiscoveryInputMetaData.getReferenceData().getReferenceSequenceDictionaryBroadcast();

        // almost every thing happens in this series of maps
        final JavaPairRDD<CpxVariantCanonicalRepresentation, Iterable<CpxVariantInducingAssemblyContig>> interpretationAndAssemblyEvidence =
                assemblyContigs
                        .mapToPair(tig -> getOneVariantFromOneContig(tig, referenceSequenceDictionaryBroadcast.getValue()))
                        .groupByKey(); // two contigs could give the same variant

        return interpretationAndAssemblyEvidence.map(pair -> turnIntoVariantContext(pair, referenceBroadcast.getValue())).collect();
    }

    private static Tuple2<CpxVariantCanonicalRepresentation, CpxVariantInducingAssemblyContig> getOneVariantFromOneContig
            (final AssemblyContigWithFineTunedAlignments contigWithFineTunedAlignments,
             final SAMSequenceDictionary refSequenceDictionary) {

        final AssemblyContigWithFineTunedAlignments furtherProcessedContig =
                furtherPreprocess(contigWithFineTunedAlignments, refSequenceDictionary);
        final CpxVariantInducingAssemblyContig cpxVariantInducingAssemblyContig =
                new CpxVariantInducingAssemblyContig(furtherProcessedContig, refSequenceDictionary);
        return new Tuple2<>(new CpxVariantCanonicalRepresentation(cpxVariantInducingAssemblyContig), cpxVariantInducingAssemblyContig);
    }

    // =================================================================================================================

    /**
     * Essentially, this step is to de-overlap the alignments
     * (see {@link #deOverlapAlignments(List, SAMSequenceDictionary)})
     * because it would be very difficult to extract retracting jumps
     * (basically indicating homology, which is non-essential to event interpretation)
     * keeping track of them, and making sense of the event.
     *
     * @return the input contig with its alignments de-overlapped
     */
    @VisibleForTesting
    static AssemblyContigWithFineTunedAlignments furtherPreprocess(final AssemblyContigWithFineTunedAlignments contigWithFineTunedAlignments,
                                                                   final SAMSequenceDictionary refSequenceDictionary) {

        final List<AlignmentInterval> deOverlappedAlignmentConfiguration =
                deOverlapAlignments(contigWithFineTunedAlignments.getAlignments(), refSequenceDictionary);

        return new AssemblyContigWithFineTunedAlignments(
                new AlignedContig(contigWithFineTunedAlignments.getContigName(), contigWithFineTunedAlignments.getContigSequence(),
                        deOverlappedAlignmentConfiguration),
                contigWithFineTunedAlignments.getInsertionMappings(),
                contigWithFineTunedAlignments.hasEquallyGoodAlnConfigurations(),
                contigWithFineTunedAlignments.getSAtagForGoodMappingToNonCanonicalChromosome());
    }

    // TODO: 3/9/18 though this function is tested, the test case doesn't cover some edge cases, let's come back later when we have more experience
    /**
     * Here we follow the de-overlapping strategy as implemented in
     * {@link CpxVariantInterpreter#removeOverlap(AlignmentInterval, AlignmentInterval, int, SAMSequenceDictionary, boolean, boolean)},
     * with the following exception:
     * the head and tail alignments are not clipped, i.e. they are always kept as is.
     */
    @VisibleForTesting
    static List<AlignmentInterval> deOverlapAlignments(final List<AlignmentInterval> originalAlignments,
                                                       final SAMSequenceDictionary refSequenceDictionary) {
        final List<AlignmentInterval> result = new ArrayList<>(originalAlignments.size());
        final int totalCount = originalAlignments.size();
        AlignmentInterval one = originalAlignments.get(0);
        for (int i = 1; i < totalCount ; ++i) {
            final AlignmentInterval two = originalAlignments.get(i);
            // TODO: 11/5/17 an edge case is possible where the best configuration contains two alignments,
            //       one of which contains a large gap, and since the gap split happens after the configuration scoring,
            //       (that gap split happens after scoring is due to how MQ and AS are used in the scoring step, gap-split alignment cannot use originating alignment's values, but it takes time to recompute)
            //       one of the alignment from the gap split may be contained in the other original alignment, leading to problems;
            //       here we skip the alignment that is BEFORE the child alignment from the gap-split,
            //       IFF that alignment contains the child alignment in terms of their spans on the read/contig
            //       if you are concerned about the first child alignment from the same gapped alignment being skipped,
            //       don't worry, that is impossible because child alignments of the same gapped alignment cannot overlap on the read.
            if (two.alnModType.equals(ContigAlignmentsModifier.AlnModType.FROM_SPLIT_GAPPED_ALIGNMENT)) {
                if ( one.containsOnRead(two) )
                    continue;
            }
            final int overlapOnContig = AlignmentInterval.overlapOnContig(one, two);
            if ( overlapOnContig == 0 ) { // nothing to remove
                // an extreme rare case is possible when, after de-overlapping, the alignment block is only 1bp long (yes, it occurs),
                // this throws off segmentation algo, hence we drop it here
                // more generally, the case could be that a small (say 2bp, or 5bp) alignment block is left, we need a more general strategy
                if (one.getSizeOnRead() >= MIN_READ_SPAN_AFTER_DEOVERLAP)
                    result.add(one);
                one = two;
            } else {
                final Tuple2<AlignmentInterval, AlignmentInterval> deoverlapped =
                        removeOverlap(one, two, overlapOnContig, refSequenceDictionary,
                                one.equals(originalAlignments.get(0)),
                                two.equals(originalAlignments.get(totalCount - 1)));
                if (deoverlapped._1.getSizeOnRead() >= MIN_READ_SPAN_AFTER_DEOVERLAP) // see comment above
                    result.add(deoverlapped._1);
                one = deoverlapped._2;
            }
        }
        result.add(one);

        final AlignmentInterval head = result.get(0);
        final AlignmentInterval tail = result.get(result.size() - 1);
        final int overlap = AlignmentInterval.overlapOnContig(head, tail);
        if ( overlap > 0 ) {
            throw new UnhandledCaseSeen(
                    "After alignment preprocessing, head and tail alignments should not overlap on read; " +
                    "otherwise it would mean some middle alignments are contained in head/tail which should not happen." +
                    originalAlignments.stream().map(AlignmentInterval::toPackedString).collect(Collectors.toList()));
        }

        return result;
    }

    /**
     * Removes overlap between input {@code contig}'s two alignments.
     * @param one alignment one
     * @param two alignment two
     * @param overlapOnRead non-positive values will throw IllegalArgumentException
     * @param dictionary if null, then {@code one} and {@code two} must be mapped to the same chromosome
     * @param firstIsAlignmentHead if {@code one} is the head alignment of the contig
     * @param secondIsAlignmentTail if {@code two} is the tail alignment of the contig
     */
    @VisibleForTesting
    static Tuple2<AlignmentInterval, AlignmentInterval> removeOverlap(final AlignmentInterval one, final AlignmentInterval two,
                                                                      final int overlapOnRead,
                                                                      final SAMSequenceDictionary dictionary,
                                                                      final boolean firstIsAlignmentHead,
                                                                      final boolean secondIsAlignmentTail) {

        if (overlapOnRead <= 0)
            throw new IllegalArgumentException("Overlap on read is non-positive for two alignments: "
                    + one.toPackedString() + "\t" + two.toPackedString());

        if (one.containsOnRead(two) || two.containsOnRead(one))
            throw new IllegalArgumentException("Two input alignments' overlap on read consumes completely one of them.\t"
                    + one.toPackedString() + "\t" + two.toPackedString());

        final boolean oneYieldToTwo;
        if (firstIsAlignmentHead) {
            oneYieldToTwo = false;
        } else if (secondIsAlignmentTail) {
            oneYieldToTwo = true;
        } else {
            oneYieldToTwo = yieldOverlapToAlignmentTwo(one, two, dictionary);
        }

        final AlignmentInterval reconstructedOne, reconstructedTwo;
        if (oneYieldToTwo) {
            reconstructedOne = ContigAlignmentsModifier.clipAlignmentInterval(one, overlapOnRead, true);
            reconstructedTwo = two;
        } else {
            reconstructedOne = one;
            reconstructedTwo = ContigAlignmentsModifier.clipAlignmentInterval(two, overlapOnRead, false);
        }
        return new Tuple2<>(reconstructedOne, reconstructedTwo);
    }

    /**
     * Implementing homology-yielding strategy between two alignments {@code one} and {@code two}.
     *
     * <p>
     *     The strategy aims to follow the left-align convention:
     *     <ul>
     *         <li>
     *             when two alignments are mapped to different chromosomes,
     *             the homologous sequence is yielded to the alignment block
     *             whose reference contig comes later as defined by {@code refSequenceDictionary}
     *         </li>
     *         <li>
     *             when two alignments are mapped to the same chromosome,
     *             <ul>
     *                 <li>
     *                     if the two alignments are of the same orientation,
     *                     the homologous sequence is yielded to {@code two}
     *                     when both are of the '+' strand
     *                 </li>
     *                 <li>
     *                     if the two alignments are of opposite orientations,
     *                     the homologous sequence is yielded so that
     *                     a) when the two alignments' ref span doesn't overlap,
     *                        we follow the left align convention for the resulting strand-switch breakpoints;
     *                     b) when the two alignments' ref span do overlap,
     *                        we makes it so that the inverted duplicated reference span is minimized
     *                        (this avoids over detection of inverted duplications by
     *                        {@link SimpleChimera#isCandidateInvertedDuplication()}}
     *                 </li>
     *             </ul>
     *         </li>
     *     </ul>
     * </p>
     *
     *
     * @return true if {@code one} should yield the homologous sequence to {@code two}.
     * @throws IllegalArgumentException when the two alignments contains one another in terms of their read span, or
     *                                  when the two alignments map to different chromosomes yet {@code refSequenceDictionary} is {@code null}
     */
    @VisibleForTesting
    static boolean yieldOverlapToAlignmentTwo(final AlignmentInterval one, final AlignmentInterval two,
                                              final SAMSequenceDictionary refSequenceDictionary) {

        Utils.validateArg( !(one.containsOnRead(two) || two.containsOnRead(one)),
                "assumption that two alignments don't contain one another on their read span is violated.\n" +
                        one.toPackedString() + "\n" + two.toPackedString());

        final boolean oneYieldToTwo;
        if ( one.referenceSpan.getContig().equals(two.referenceSpan.getContig()) ) {
            // motivation: when the two alignments' ref span doesn't overlap,
            //             this strategy follows the left align convention for the strand-switch breakpoints;
            //             when the two alignments' ref span do overlap,
            //             this strategy makes it so that the inverted duplicated reference span is minimal.
            if (one.forwardStrand != two.forwardStrand) {
                // jumpStart is for "the starting reference location of a jump that linked two alignment intervals", and
                // jumpLandingRefLoc is for "that jump's landing reference location"
                final int jumpStartRefLoc = one.referenceSpan.getEnd(),
                        jumpLandingRefLoc = two.referenceSpan.getStart();
                oneYieldToTwo = ( (jumpStartRefLoc <= jumpLandingRefLoc) == one.forwardStrand );
            } else {
                oneYieldToTwo = one.forwardStrand;
            }
        } else {
            if (refSequenceDictionary == null)
                throw new IllegalArgumentException("despite input alignments mapped to different chromosomes, " +
                        "input reference sequence dictionary is null. " + one.toPackedString() + "\t" + two.toPackedString());

            oneYieldToTwo = IntervalUtils.compareContigs(one.referenceSpan, two.referenceSpan, refSequenceDictionary) > 0;
        }
        return oneYieldToTwo;
    }

    // =================================================================================================================

    @VisibleForTesting
    static VariantContext turnIntoVariantContext(final Tuple2<CpxVariantCanonicalRepresentation, Iterable<CpxVariantInducingAssemblyContig>> pair,
                                                 final ReferenceMultiSparkSource reference)
            throws IOException {

        final CpxVariantCanonicalRepresentation cpxVariantCanonicalRepresentation = pair._1;
        final byte[] refBases = getRefBases(reference, cpxVariantCanonicalRepresentation);

        final VariantContextBuilder rawVariantContextBuilder = cpxVariantCanonicalRepresentation.toVariantContext(refBases);
        final Iterable<CpxVariantInducingAssemblyContig> evidenceContigs = pair._2;

        if (Utils.stream(evidenceContigs).anyMatch(evidenceContig -> evidenceContig.getPreprocessedTig().getAlignments().size()==0)) {
            throw new GATKException("Some contigs were unmapped, yet seem to be used for inference.\n" + cpxVariantCanonicalRepresentation.toString() +
                    Utils.stream(evidenceContigs).map(tig->tig.getPreprocessedTig().getContigName()).collect(Collectors.toList()));
        }

        final List<String> contigNames = new ArrayList<>();
        final List<String> mayContainNoInfoOnNonCanonicalMapping = new ArrayList<>(); // for storing AssemblyContigWithFineTunedAlignments.NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME
        final List<String> minFlankingMQs = new ArrayList<>();
        final List<String> minFlankingLengths = new ArrayList<>();

        Utils.stream(evidenceContigs)
                .sorted(Comparator.comparing(tig -> tig.getPreprocessedTig().getContigName()))
                .forEach(annotatedCpxInducingContig -> {

                    final AssemblyContigWithFineTunedAlignments preprocessedTig =
                            annotatedCpxInducingContig.getPreprocessedTig();

                    contigNames.add(preprocessedTig.getContigName());

                    final String saTagForGoodMappingToNonCanonicalChromosome =
                            preprocessedTig.getSAtagForGoodMappingToNonCanonicalChromosome();
                    if ( ! saTagForGoodMappingToNonCanonicalChromosome
                            .equals(AssemblyContigWithFineTunedAlignments.NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME)) {
                        mayContainNoInfoOnNonCanonicalMapping.add(saTagForGoodMappingToNonCanonicalChromosome);
                    } else {
                        mayContainNoInfoOnNonCanonicalMapping.add(".");
                    }

                    minFlankingMQs.add( String.valueOf(Math.min(preprocessedTig.getHeadAlignment().mapQual, preprocessedTig.getTailAlignment().mapQual)) );
                    final int minAlnLen = Math.min(preprocessedTig.getHeadAlignment().getSizeOnRead(), preprocessedTig.getTailAlignment().getSizeOnRead());
                    minFlankingLengths.add( String.valueOf(minAlnLen) );
                });

        final Map<String, String> attributeMap = new HashMap<>();

        attributeMap.put(TOTAL_MAPPINGS, String.valueOf(contigNames.size()));

        attributeMap.put(CONTIG_NAMES, String.join(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR, contigNames));

        if ( ! mayContainNoInfoOnNonCanonicalMapping.equals( Collections.nCopies(contigNames.size(), ".") ) ) {
            attributeMap.put(CTG_GOOD_NONCANONICAL_MAPPING, String.join(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR, mayContainNoInfoOnNonCanonicalMapping));
        }

        attributeMap.put(HQ_MAPPINGS,       String.valueOf(minFlankingMQs.stream().mapToInt(Integer::valueOf).filter( mq -> mq >= CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD ).count()) );
        attributeMap.put(MAPPING_QUALITIES, String.join(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR, minFlankingMQs));
        attributeMap.put(ALIGN_LENGTHS,     String.join(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR, minFlankingLengths));
        attributeMap.put(MAX_ALIGN_LENGTH,  String.valueOf(minFlankingLengths.stream().mapToInt(Integer::valueOf).max().orElse(0)));

        attributeMap.forEach(rawVariantContextBuilder::attribute);
        return rawVariantContextBuilder.make();
    }

    private static byte[] getRefBases( final ReferenceMultiSparkSource reference, final CpxVariantCanonicalRepresentation cpxVariantCanonicalRepresentation)
            throws IOException {
        final SimpleInterval affectedRefRegion = cpxVariantCanonicalRepresentation.getAffectedRefRegion();
        SimpleInterval refBase = new SimpleInterval(affectedRefRegion.getContig(), affectedRefRegion.getStart(),
                                                                                   affectedRefRegion.getStart());
        return reference
                .getReferenceBases(refBase)
                .getBases();
    }

    // =================================================================================================================

    static final class UnhandledCaseSeen extends GATKException.ShouldNeverReachHereException {
        private static final long serialVersionUID = 0L;
        UnhandledCaseSeen( final String s ) {
            super(s);
        }

        UnhandledCaseSeen( final String s, final Throwable throwable ) {
            super(s, throwable);
        }

        UnhandledCaseSeen( final Throwable throwable) {this("Seeing unhandled case", throwable);}
    }
}
