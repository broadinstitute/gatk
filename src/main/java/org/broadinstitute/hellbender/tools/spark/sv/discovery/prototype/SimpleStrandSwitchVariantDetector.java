package org.broadinstitute.hellbender.tools.spark.sv.discovery.prototype;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.*;
import org.broadinstitute.hellbender.tools.spark.sv.utils.RDDUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVVCFWriter;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SvCigarUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;
import scala.Tuple3;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;


final class SimpleStrandSwitchVariantDetector implements VariantDetectorFromLocalAssemblyContigAlignments {

    @SuppressWarnings("unchecked")
    private static final List<String> EMPTY_INSERTION_MAPPINGS = Collections.EMPTY_LIST;

    static final int MORE_RELAXED_ALIGNMENT_MIN_LENGTH = 30;
    static final int MORE_RELAXED_ALIGNMENT_MIN_MQ = 20;

    @Override
    public void inferSvAndWriteVCF(final JavaRDD<AlignedContig> contigs, final String vcfOutputFileName,
                                   final Broadcast<ReferenceMultiSource> broadcastReference, final String fastaReference,
                                   final Logger toolLogger) {

        contigs.cache();
        toolLogger.info(contigs.count() + " chimeras indicating either 1) simple strand-switch breakpoints, or 2) inverted duplication.");

        // split between suspected inv dup VS strand-switch breakpoints
        // logic flow: first modify alignments (heuristically) if there are overlaps on read between the alignments,
        //             then split the input reads into two classes--those judged by IsLikelyInvertedDuplication are likely invdup and those aren't
        //             finally send the two split reads down different path, one for invdup and one for BND records
        final Tuple2<JavaRDD<AlignedContig>, JavaRDD<AlignedContig>> split =
                RDDUtils.split(contigs.map(SimpleStrandSwitchVariantDetector::removeOverlap),
                        contig -> BreakpointComplications.isLikelyInvertedDuplication(contig.alignmentIntervals.get(0),
                                contig.alignmentIntervals.get(1)), true);

        final JavaRDD<VariantContext> simpleStrandSwitchBkpts =
                dealWithSimpleStrandSwitchBkpts(split._2, broadcastReference, toolLogger);
        SVVCFWriter.writeVCF(null, vcfOutputFileName.replace(".vcf", "_simpleSS.vcf"),
                fastaReference, simpleStrandSwitchBkpts, toolLogger);

        // TODO: 8/23/17 add inv dup code in the next pr
    }

    // =================================================================================================================

    /**
     * Removes overlap from a designated alignment interval, so that the inverted duplicated reference span is minimal.
     * If the two alignment intervals are NOT overlapping, return the original aligned contig.
     */
    private static AlignedContig removeOverlap(final AlignedContig contig) {
        final int overlapOnRead = AlignmentInterval.overlapOnContig(contig.alignmentIntervals.get(0),
                                                                    contig.alignmentIntervals.get(1));
        if (overlapOnRead==0) {
            return contig;
        } else {
            final AlignmentInterval one = contig.alignmentIntervals.get(0),
                                    two = contig.alignmentIntervals.get(1);
            // jumpStart is for "the starting reference location of a jump that linked two alignment intervals", and
            // jumpLandingRefLoc is for "that jump's landing reference location"
            final int jumpStartRefLoc   = one.referenceSpan.getEnd(),
                      jumpLandingRefLoc = two.referenceSpan.getStart();
            final AlignmentInterval reconstructedOne, reconstructedTwo;
            if (jumpStartRefLoc <= jumpLandingRefLoc ^ one.forwardStrand) {
                reconstructedOne = one;
                reconstructedTwo = clipAlignmentInterval(two, overlapOnRead, false);
            } else {
                reconstructedOne = clipAlignmentInterval(one, overlapOnRead, true);
                reconstructedTwo = two;
            }
            return new AlignedContig(contig.contigName, contig.contigSequence,
                                     Arrays.asList(reconstructedOne, reconstructedTwo),
                                     contig.hasEquallyGoodAlnConfigurations);
        }
    }

    /**
     * Given {@code clipLengthOnRead} to be clipped from an aligned contig's {@link AlignmentInterval} {@code input},
     * return a modified {@link AlignmentInterval} with the requested number of bases clipped away and ref span recomputed.
     * @param input                 alignment to be modified (record not modified in place)
     * @param clipLengthOnRead      number of read bases to be clipped away
     * @param clipFrom3PrimeEnd     to clip from the 3' end of the read or not
     */
    private static AlignmentInterval clipAlignmentInterval(final AlignmentInterval input, final int clipLengthOnRead,
                                                           final boolean clipFrom3PrimeEnd) {
        Utils.validateArg(clipLengthOnRead < input.endInAssembledContig - input.startInAssembledContig + 1,
                            "input alignment to be clipped away: " + input.toPackedString() + "\twith clip length: " + clipLengthOnRead);

        final Tuple2<SimpleInterval, Cigar> result = computeNewRefSpanAndCigar(input, clipLengthOnRead, clipFrom3PrimeEnd);
        final int newTigStart, newTigEnd;
        if (clipFrom3PrimeEnd) {
            newTigStart = input.startInAssembledContig;
            newTigEnd   = input.endInAssembledContig - clipLengthOnRead;
        } else {
            newTigStart = input.startInAssembledContig + clipLengthOnRead;
            newTigEnd   = input.endInAssembledContig;
        }
        return new AlignmentInterval(result._1, newTigStart, newTigEnd, result._2, input.forwardStrand, input.mapQual,
                StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.MISSING_NM,
                input.alnScore, input.isFromSplitGapAlignment, input.hasUndergoneOverlapRemoval);
    }

    /**
     * Given {@code clipLengthOnRead} to be clipped from an aligned contig's {@link AlignmentInterval} {@code input},
     * return a modified reference span with the requested number of bases clipped away and new cigar.
     * @param input                 alignment to be modified (record not modified in place)
     * @param clipLengthOnRead      number of read bases to be clipped away
     * @param clipFrom3PrimeEnd     to clip from the 3' end of the read or not
     * @throws IllegalArgumentException if the {@code input} alignment contains {@link CigarOperator#P} or {@link CigarOperator#N} operations
     */
    @VisibleForTesting
    static Tuple2<SimpleInterval, Cigar> computeNewRefSpanAndCigar(final AlignmentInterval input, final int clipLengthOnRead,
                                                                   final boolean clipFrom3PrimeEnd) {
        Utils.validateArg(input.cigarAlong5to3DirectionOfContig.getCigarElements().stream().map(CigarElement::getOperator)
                        .noneMatch(op -> op.equals(CigarOperator.N) || op.isPadding()),
                "Input alignment contains padding or skip operations, which is currently unsupported: " + input.toPackedString());

        final Tuple3<List<CigarElement>, List<CigarElement>, List<CigarElement>> threeSections = splitCigarByLeftAndRightClipping(input);
        final List<CigarElement> leftClippings = threeSections._1();
        final List<CigarElement> unclippedCigarElementsForThisAlignment = threeSections._2();
        final List<CigarElement> rightClippings = threeSections._3();

        int readBasesConsumed = 0, refBasesConsumed = 0;
        final List<CigarElement> cigarElements = new ArrayList<>(unclippedCigarElementsForThisAlignment);
        if (clipFrom3PrimeEnd) Collections.reverse(cigarElements);
        final List<CigarElement> newMiddleSection = new ArrayList<>(unclippedCigarElementsForThisAlignment.size());
        for (int idx = 0; idx < cigarElements.size(); ++idx) {
            final CigarElement ce = cigarElements.get(idx);
            if (ce.getOperator().consumesReadBases()) {
                if (readBasesConsumed + ce.getLength() < clipLengthOnRead) {
                    readBasesConsumed += ce.getLength();
                } else { // enough read bases would be clipped, note that here we can have only either 'M' or 'I'

                    if (!ce.getOperator().isAlignment() && !ce.getOperator().equals(CigarOperator.I))
                        throw new GATKException.ShouldNeverReachHereException("Logic error, should not reach here: operation consumes read bases but is neither alignment nor insertion.\n " +
                                "Original cigar(" + TextCigarCodec.encode(input.cigarAlong5to3DirectionOfContig) + ")\tclipLengthOnRead(" + clipLengthOnRead + ")\t" + (clipFrom3PrimeEnd ? "clipFrom3PrimeEnd" : "clipFrom5PrimeEnd"));

                    // deal with cigar first
                    newMiddleSection.add( new CigarElement(clipLengthOnRead, CigarOperator.S) );
                    // # of bases left, for the current operation, after requested clipping is done,
                    // may be 0 because we probably don't need the entire current operation to have the requested
                    // number of read bases clipped, e.g. we only 15 more bases from the current operation,
                    // but its length is 20, then 5 bases should be "left over"
                    final int leftOverBasesForCurrOp = readBasesConsumed + ce.getLength() - clipLengthOnRead;
                    if (leftOverBasesForCurrOp!=0) {
                        newMiddleSection.add(// again note that here we can have only either 'M' or 'I'
                                new CigarElement(leftOverBasesForCurrOp, ce.getOperator().isAlignment() ? CigarOperator.M
                                                                                                        : CigarOperator.S) );
                    }
                    newMiddleSection.addAll( cigarElements.subList(idx+1, cigarElements.size()) );

                    // then deal with ref span
                    refBasesConsumed += ce.getOperator().isAlignment() ? (clipLengthOnRead - readBasesConsumed)
                                                                       : ce.getLength();

                    break;
                }
            }
            if ( ce.getOperator().consumesReferenceBases() ) { // if reaches here, not enough read bases have been consumed
                refBasesConsumed += ce.getLength();
            }
        }
        if (clipFrom3PrimeEnd) Collections.reverse(newMiddleSection);
        final Cigar newCigar = constructNewCigar(leftClippings, newMiddleSection, rightClippings);
        if (newCigar.getCigarElements().isEmpty())
            throw new GATKException("Logic error: new cigar is empty.\t" + input.toPackedString() + "\tclip length " +
                    clipLengthOnRead + "\tclip from end " + (clipFrom3PrimeEnd? "3":"5"));

        final SimpleInterval newRefSpan;
        if (clipFrom3PrimeEnd == input.forwardStrand) {
            newRefSpan = new SimpleInterval(input.referenceSpan.getContig(), input.referenceSpan.getStart(),
                    input.referenceSpan.getEnd() - refBasesConsumed);
        } else {
            newRefSpan = new SimpleInterval(input.referenceSpan.getContig(), input.referenceSpan.getStart() + refBasesConsumed,
                    input.referenceSpan.getEnd());
        }

        return new Tuple2<>(newRefSpan, newCigar);
    }

    /**
     * Extract, from provided {@code input} alignment, three parts of cigar elements:
     * <ul>
     *     <li> a list of cigar elements leading to the start position of the input alignment on the read</li>
     *     <li> </li>
     *     <li> a list of cigar elements after the end position of the input alignment on the read</li>
     * </ul>
     * For example, an alignment with cigar "10H20S5M15I25M35D45M30S40H" with start pos 31 and end pos 120,
     * the result would be ([10H, 20S], [5M, 15I, 25M, 35D, 45M], [30S, 40H]).
     */
    @VisibleForTesting
    static Tuple3<List<CigarElement>, List<CigarElement>, List<CigarElement>> splitCigarByLeftAndRightClipping(final AlignmentInterval input) {
        final List<CigarElement> cigarElements = input.cigarAlong5to3DirectionOfContig.getCigarElements();
        final List<CigarElement> left = new ArrayList<>(cigarElements.size()),
                                 middle = new ArrayList<>(cigarElements.size()),
                                 right = new ArrayList<>(cigarElements.size());
        // using input's startInAssembledContig and endInAssembledContig in testing conditions because they must be the boundaries of alignment operations
        int readBasesConsumed = 0;
        for(final CigarElement ce : cigarElements) {
            if (readBasesConsumed < input.startInAssembledContig-1) {
                left.add(ce);
            } else if (readBasesConsumed < input.endInAssembledContig) {
                middle.add(ce);
            } else {
                right.add(ce);
            }
            readBasesConsumed += ce.getOperator().consumesReadBases() || ce.getOperator().equals(CigarOperator.H)? ce.getLength() : 0;
        }

        if (middle.isEmpty())
            throw new GATKException("Logic error: cigar elements corresponding to alignment block is empty. " + input.toPackedString());

        return new Tuple3<>(left, middle, right);
    }

    private static Cigar constructNewCigar(final List<CigarElement> left, final List<CigarElement> middle, final List<CigarElement> right) {
        final Cigar cigar = new Cigar(left);
        middle.forEach(cigar::add);
        right.forEach(cigar::add);
        return new Cigar(SvCigarUtils.compactifyNeighboringSoftClippings(cigar.getCigarElements()));
    }
    // =================================================================================================================

    /**
     * @throws IllegalArgumentException if the assumption that the input aligned assembly contig has 2 alignments
     *                                  mapped to the same chr with strand switch is invalid
     */
    private static Tuple2<ChimericAlignment, byte[]> convertAlignmentIntervalsToChimericAlignment
    (final AlignedContig contigWith2AIMappedToSameChrAndStrandSwitch) {
        Utils.validateArg(InternalSvDiscoverFromLocalAssemblyContigAlignmentsSpark.isLikelyInvBreakpointOrInsInv(contigWith2AIMappedToSameChrAndStrandSwitch),
                "assumption that input aligned assembly contig has 2 alignments mapped to the same chr with strand switch is invalid.\n" +
                        InternalSvDiscoverFromLocalAssemblyContigAlignmentsSpark.onErrorStringRepForAlignedContig(contigWith2AIMappedToSameChrAndStrandSwitch));

        final AlignmentInterval intervalOne = contigWith2AIMappedToSameChrAndStrandSwitch.alignmentIntervals.get(0),
                                intervalTwo = contigWith2AIMappedToSameChrAndStrandSwitch.alignmentIntervals.get(1);

        // TODO: 8/28/17 this default empty insertion mapping treatment is temporary and should be fixed later
        return new Tuple2<>(new ChimericAlignment(intervalOne, intervalTwo, EMPTY_INSERTION_MAPPINGS,
                contigWith2AIMappedToSameChrAndStrandSwitch.contigName), contigWith2AIMappedToSameChrAndStrandSwitch.contigSequence);
    }

    /**
     * Roughly similar to {@link ChimericAlignment#nextAlignmentMayBeNovelInsertion(AlignmentInterval, AlignmentInterval, Integer)}:
     *  1) either alignment may have very low mapping quality (a more relaxed mapping quality threshold);
     *  2) either alignment may consume only a "short" part of the contig, or if assuming that the alignment consumes
     *     roughly the same amount of ref bases and read bases, has isAlignment that is too short
     */
    private static boolean splitPairStrongEnoughEvidenceForCA(final AlignmentInterval intervalOne,
                                                              final AlignmentInterval intervalTwo,
                                                              final int mapQThresholdInclusive,
                                                              final int alignmentLengthThresholdInclusive) {

        if (intervalOne.mapQual < mapQThresholdInclusive || intervalTwo.mapQual < mapQThresholdInclusive)
            return false;

        final int overlap = AlignmentInterval.overlapOnContig(intervalOne, intervalTwo);

        final int x = intervalOne.endInAssembledContig - intervalOne.startInAssembledContig + 1,
                  y = intervalTwo.endInAssembledContig - intervalTwo.startInAssembledContig + 1;

        return Math.min(x - overlap, y - overlap) >= alignmentLengthThresholdInclusive;
    }

    // workflow manager for simple strand-switch alignment contigs
    private JavaRDD<VariantContext> dealWithSimpleStrandSwitchBkpts(final JavaRDD<AlignedContig> contigs,
                                                                    final Broadcast<ReferenceMultiSource> broadcastReference,
                                                                    final Logger toolLogger) {

        final JavaPairRDD<ChimericAlignment, byte[]> simpleStrandSwitchBkpts =
                contigs
                        .filter(tig ->
                                splitPairStrongEnoughEvidenceForCA(tig.alignmentIntervals.get(0), tig.alignmentIntervals.get(1),
                                        MORE_RELAXED_ALIGNMENT_MIN_MQ,  MORE_RELAXED_ALIGNMENT_MIN_LENGTH))
                        .mapToPair(SimpleStrandSwitchVariantDetector::convertAlignmentIntervalsToChimericAlignment).cache();

        toolLogger.info(simpleStrandSwitchBkpts.count() + " chimeras indicating simple strand-switch breakpoints.");

        return simpleStrandSwitchBkpts
                .mapToPair(pair -> new Tuple2<>(new NovelAdjacencyReferenceLocations(pair._1, pair._2), pair._1))
                .groupByKey()
                .mapToPair(noveltyAndEvidence -> inferBNDType(noveltyAndEvidence, broadcastReference.getValue()))
                .flatMap(noveltyTypeAndEvidence ->
                        AnnotatedVariantProducer
                                .produceAnnotatedBNDmatesVcFromNovelAdjacency(noveltyTypeAndEvidence._1,
                                        noveltyTypeAndEvidence._2._1, noveltyTypeAndEvidence._2._2, broadcastReference).iterator());
    }

    private static Tuple2<NovelAdjacencyReferenceLocations, Tuple2<List<SvType>, Iterable<ChimericAlignment>>>
    inferBNDType(final Tuple2<NovelAdjacencyReferenceLocations, Iterable<ChimericAlignment>> noveltyAndEvidence,
                 final ReferenceMultiSource reference) {

        final NovelAdjacencyReferenceLocations novelAdjacency = noveltyAndEvidence._1;
        final Iterable<ChimericAlignment> chimericAlignments = noveltyAndEvidence._2;
        final BreakEndVariantType bkpt_1, bkpt_2;
        if (novelAdjacency.strandSwitch == StrandSwitch.FORWARD_TO_REVERSE) {
            bkpt_1 = new BreakEndVariantType.INV55BND(novelAdjacency, true, reference);
            bkpt_2 = new BreakEndVariantType.INV55BND(novelAdjacency, false, reference);
        } else if (novelAdjacency.strandSwitch == StrandSwitch.REVERSE_TO_FORWARD) {
            bkpt_1 = new BreakEndVariantType.INV33BND(novelAdjacency, true, reference);
            bkpt_2 = new BreakEndVariantType.INV33BND(novelAdjacency, false, reference);
        } else {
            throw new GATKException("Wrong type of novel adjacency sent to wrong analysis pathway: no strand-switch being sent to strand-switch path. \n" +
                    Utils.stream(chimericAlignments).map(ChimericAlignment::onErrStringRep).collect(Collectors.toList()));
        }

        return new Tuple2<>(novelAdjacency, new Tuple2<>(Arrays.asList(bkpt_1, bkpt_2), chimericAlignments));
    }

    // TODO: 8/23/17 add inv dup code in the next PR
}
