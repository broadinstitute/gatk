package org.broadinstitute.hellbender.tools.spark.sv.discovery.prototype;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SvCigarUtils;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;
import scala.Tuple3;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public final class ContigAlignmentsModifier {

    /**
     * Removes overlap between input {@code contig}'s two alignments.
     * If the two alignment intervals are NOT overlapping, return the original aligned contig.
     * @param dictionary if null, then {@code one} and {@code two} must be mapped to the same chromosome
     */
    public static List<AlignmentInterval> removeOverlap(final AlignmentInterval one, final AlignmentInterval two,
                                                        final SAMSequenceDictionary dictionary) {
        if (dictionary == null)
            Utils.validateArg(one.referenceSpan.getContig().equals(two.referenceSpan.getContig()),
                    "despite input alignments mapped to different chromosomes, input reference sequence dictionary is null. \n" +
                            one.toPackedString() + "\t" + two.toPackedString());

        final AlignmentInterval reconstructedOne, reconstructedTwo;
        final int overlapOnRead = AlignmentInterval.overlapOnContig(one, two);
        if (overlapOnRead == 0) {
            reconstructedOne = one;
            reconstructedTwo = two;
        } else {
            final boolean oneYieldToTwo;
            if (one.referenceSpan.getContig().equals(two.referenceSpan.getContig())) {
                if (one.forwardStrand != two.forwardStrand) { // so that the inverted duplicated reference span is minimal.
                    // jumpStart is for "the starting reference location of a jump that linked two alignment intervals", and
                    // jumpLandingRefLoc is for "that jump's landing reference location"
                    final int jumpStartRefLoc = one.referenceSpan.getEnd(),
                            jumpLandingRefLoc = two.referenceSpan.getStart();
                    oneYieldToTwo = jumpStartRefLoc <= jumpLandingRefLoc == one.forwardStrand;
                } else {
                    oneYieldToTwo = one.forwardStrand;
                }
            } else {
                oneYieldToTwo = IntervalUtils.compareContigs(one.referenceSpan, two.referenceSpan, dictionary) > 0;
            }

            if (oneYieldToTwo) {
                reconstructedOne = clipAlignmentInterval(one, overlapOnRead, true);
                reconstructedTwo = two;
            } else {
                reconstructedOne = one;
                reconstructedTwo = clipAlignmentInterval(two, overlapOnRead, false);
            }

        }
        return Arrays.asList(reconstructedOne, reconstructedTwo);
    }

    /**
     * Given {@code clipLengthOnRead} to be clipped from an aligned contig's {@link AlignmentInterval} {@code input},
     * return a modified {@link AlignmentInterval} with the requested number of bases clipped away and ref span recomputed.
     * <p>
     *     Note that the returned AlignmentInterval will have its {@code NM} set to {@link AlignmentInterval#NO_NM}
     *     whereas its {@code AS} and mapping quality are simply copied from the input alignment
     *     (recalculating them doesn't payoff at this stage, and makes the alignment filtering step more complicated since they are used there)
     * </p>
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
                AlignmentInterval.NO_NM, input.alnScore, input.alnModType);
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
}
