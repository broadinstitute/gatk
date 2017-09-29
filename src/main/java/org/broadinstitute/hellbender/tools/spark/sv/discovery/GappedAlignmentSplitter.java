package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.SVConstants;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public final class GappedAlignmentSplitter {

    /**
     * Splits a gapped alignment into multiple alignment regions, based on input sensitivity: i.e. if the gap(s) in the input alignment
     * is equal to or larger than the size specified by {@code sensitivity}, two alignment regions will be generated.
     *
     * To fulfill this functionality, needs to accomplish three key tasks correctly:
     * <ul>
     *     <li>generate appropriate cigar,</li>
     *     <li>infer reference coordinates,</li>
     *     <li>infer contig coordinates</li>
     * </ul>
     * As an example, an alignment with CIGAR "397S118M2D26M6I50M7I26M1I8M13D72M398S", when {@code sensitivity} is set to "1", should be split into 5 alignments
     * <ul>
     *     <li>"397S118M594S"</li>
     *     <li>"515S26M568S"</li>
     *     <li>"547S50M512S"</li>
     *     <li>"631S8M470S"</li>
     *     <li>"639S72M398S"</li>
     * </ul>
     * On the other hand, an alignment with CIGAR "10M10D10M60I10M10I10M50D10M", when {@code sensitivity} is set to "50", should be split into 3 alignments
     * <ul>
     *     <li>"10M10D10M100S"</li>
     *     <li>"80S10M10I10M10S"</li>
     *     <li>"110S10M"</li>
     * </ul>
     * And when an alignment has hard clipping adjacent to soft clippings, e.g. "1H2S3M5I10M20D6M7S8H", it should be split into alignments with CIGAR resembling the original CIGAR as much as possible:
     * <ul>
     *     <li>"1H2S3M28S8H"</li>
     *     <li>"1H10S10M13S8H"</li>
     *     <li>"1H20S6M7S8H"</li>
     * </ul>
     *
     * @return an iterable of size >= 1. if size==1, the returned iterable contains only the input (i.e. either no gap or hasn't reached sensitivity)
     */
    @VisibleForTesting
    public static Iterable<AlignedAssembly.AlignmentInterval> split(final AlignedAssembly.AlignmentInterval oneRegion,
                                                                    final int sensitivity,
                                                                    final int unclippedContigLen) {

        final List<CigarElement> cigarElements = checkCigarAndConvertTerminalInsertionToSoftClip(oneRegion.cigarAlong5to3DirectionOfContig);
        if (cigarElements.size() == 1) return new ArrayList<>( Collections.singletonList(oneRegion) );

        final List<AlignedAssembly.AlignmentInterval> result = new ArrayList<>(3); // blunt guess
        final int originalMapQ = oneRegion.mapQual;

        final List<CigarElement> cigarMemoryList = new ArrayList<>();
        final int clippedNBasesFromStart = SVVariantDiscoveryUtils.getNumClippedBases(true, cigarElements);

        final int hardClippingAtBeginning = cigarElements.get(0).getOperator()== CigarOperator.H ? cigarElements.get(0).getLength() : 0;
        final int hardClippingAtEnd = (cigarElements.get(cigarElements.size()-1).getOperator()== CigarOperator.H)? cigarElements.get(cigarElements.size()-1).getLength() : 0;
        final CigarElement hardClippingAtBeginningMaybeNull = hardClippingAtBeginning==0 ? null : new CigarElement(hardClippingAtBeginning, CigarOperator.H);
        int contigIntervalStart = 1 + clippedNBasesFromStart;
        // we are walking along the contig following the cigar, which indicates that we might be walking backwards on the reference if oneRegion.forwardStrand==false
        int refBoundary1stInTheDirectionOfContig = oneRegion.forwardStrand ? oneRegion.referenceInterval.getStart() : oneRegion.referenceInterval.getEnd();
        for (final CigarElement cigarElement : cigarElements) {
            final CigarOperator op = cigarElement.getOperator();
            final int operatorLen = cigarElement.getLength();
            switch (op) {
                case M: case EQ: case X: case S: case H:
                    cigarMemoryList.add(cigarElement);
                    break;
                case I: case D:
                    if (operatorLen < sensitivity) {
                        cigarMemoryList.add(cigarElement);
                        break;
                    }

                    // collapse cigar memory list into a single cigar for ref & contig interval computation
                    final Cigar memoryCigar = new Cigar(cigarMemoryList);
                    final int effectiveReadLen = memoryCigar.getReadLength() + SVVariantDiscoveryUtils.getTotalHardClipping(memoryCigar) - SVVariantDiscoveryUtils.getNumClippedBases(true, memoryCigar);

                    // task 1: infer reference interval taking into account of strand
                    final SimpleInterval referenceInterval;
                    if (oneRegion.forwardStrand) {
                        referenceInterval = new SimpleInterval(oneRegion.referenceInterval.getContig(),
                                refBoundary1stInTheDirectionOfContig,
                                refBoundary1stInTheDirectionOfContig + (memoryCigar.getReferenceLength()-1));
                    } else {
                        referenceInterval = new SimpleInterval(oneRegion.referenceInterval.getContig(),
                                refBoundary1stInTheDirectionOfContig - (memoryCigar.getReferenceLength()-1), // step backward
                                refBoundary1stInTheDirectionOfContig);
                    }

                    // task 2: infer contig interval
                    final int contigIntervalEnd = contigIntervalStart + effectiveReadLen - 1;

                    // task 3: now add trailing cigar element and create the real cigar for the to-be-returned AR
                    cigarMemoryList.add(new CigarElement(unclippedContigLen-contigIntervalEnd-hardClippingAtEnd, CigarOperator.S));
                    if (hardClippingAtEnd != 0) { // be faithful to hard clipping (as the accompanying bases have been hard-clipped away)
                        cigarMemoryList.add(new CigarElement(hardClippingAtEnd, CigarOperator.H));
                    }
                    final Cigar cigarForNewAlignmentInterval = new Cigar(cigarMemoryList);

                    final AlignedAssembly.AlignmentInterval split = new AlignedAssembly.AlignmentInterval(referenceInterval, contigIntervalStart, contigIntervalEnd, cigarForNewAlignmentInterval, oneRegion.forwardStrand, originalMapQ, SVConstants.DiscoveryStepConstants.ARTIFICIAL_MISMATCH);

                    result.add(split);

                    // update cigar memory
                    cigarMemoryList.clear();
                    if (hardClippingAtBeginningMaybeNull != null) {
                        cigarMemoryList.add(hardClippingAtBeginningMaybeNull); // be faithful about hard clippings
                    }
                    cigarMemoryList.add(new CigarElement(contigIntervalEnd - hardClippingAtBeginning + (op.consumesReadBases() ? operatorLen : 0), CigarOperator.S));

                    // update pointers into reference and contig
                    final int refBoundaryAdvance = op.consumesReadBases() ? memoryCigar.getReferenceLength() : memoryCigar.getReferenceLength() + operatorLen;
                    refBoundary1stInTheDirectionOfContig += oneRegion.forwardStrand ? refBoundaryAdvance : -refBoundaryAdvance;
                    contigIntervalStart += op.consumesReadBases() ? effectiveReadLen + operatorLen : effectiveReadLen;

                    break;
                default:
                    throw new GATKException("Alignment CIGAR contains an unexpected N or P element: " + oneRegion.toPackedString()); // TODO: 1/20/17 still not quite sure if this is quite right, it doesn't blow up on NA12878 WGS, but who knows what happens in the future
            }
        }

        if (result.isEmpty()) {
            return new ArrayList<>(Collections.singletonList(oneRegion));
        }

        final SimpleInterval lastReferenceInterval;
        if (oneRegion.forwardStrand) {
            lastReferenceInterval =  new SimpleInterval(oneRegion.referenceInterval.getContig(), refBoundary1stInTheDirectionOfContig, oneRegion.referenceInterval.getEnd());
        } else {
            lastReferenceInterval =  new SimpleInterval(oneRegion.referenceInterval.getContig(), oneRegion.referenceInterval.getStart(), refBoundary1stInTheDirectionOfContig);
        }

        final Cigar lastForwardStrandCigar = new Cigar(cigarMemoryList);
        int clippedNBasesFromEnd = SVVariantDiscoveryUtils.getNumClippedBases(false, cigarElements);
        result.add(new AlignedAssembly.AlignmentInterval(lastReferenceInterval,
                contigIntervalStart, unclippedContigLen-clippedNBasesFromEnd, lastForwardStrandCigar,
                oneRegion.forwardStrand, originalMapQ, SVConstants.DiscoveryStepConstants.ARTIFICIAL_MISMATCH));


        return result;
    }

    /**
     * Checks the input CIGAR for assumption that operator 'D' is not immediately adjacent to clipping operators.
     * Then convert the 'I' CigarElement, if it is at either end (terminal) of the input cigar, to a corresponding 'S' operator.
     * Note that we allow CIGAR of the format '10H10S10I10M', but disallows the format if after the conversion the cigar turns into a giant clip,
     * e.g. '10H10S10I10S10H' is not allowed (if allowed, it becomes a giant clip of '10H30S10H' which is non-sense).
     *
     * @return a pair of number of clipped (hard and soft, including the ones from the converted terminal 'I') bases at the front and back of the
     *         input {@code cigarAlongInput5to3Direction}.
     *
     * @throws IllegalArgumentException when the checks as described above fail.
     */
    @VisibleForTesting
    public static List<CigarElement> checkCigarAndConvertTerminalInsertionToSoftClip(final Cigar cigarAlongInput5to3Direction) {

        if (cigarAlongInput5to3Direction.numCigarElements()<2 ) return cigarAlongInput5to3Direction.getCigarElements();

        final List<CigarElement> cigarElements = new ArrayList<>(cigarAlongInput5to3Direction.getCigarElements());
        SVVariantDiscoveryUtils.validateCigar(cigarElements);

        final List<CigarElement> convertedList = convertInsToSoftClipFromOneEnd(cigarElements, true);
        return convertInsToSoftClipFromOneEnd(convertedList, false);
    }

    /**
     * Actually convert terminal 'I' to 'S' and in case there's an 'S' comes before 'I', compactify the two neighboring 'S' operations into one.
     *
     * @return the converted and compactified list of cigar elements
     */
    @VisibleForTesting
    public static List<CigarElement> convertInsToSoftClipFromOneEnd(final List<CigarElement> cigarElements,
                                                                    final boolean fromStart) {
        final int numHardClippingBasesFromOneEnd = SVVariantDiscoveryUtils.getNumHardClippingBases(fromStart, cigarElements);
        final int numSoftClippingBasesFromOneEnd = SVVariantDiscoveryUtils.getNumSoftClippingBases(fromStart, cigarElements);

        final int indexOfFirstNonClippingOperation;
        if (numHardClippingBasesFromOneEnd==0 && numSoftClippingBasesFromOneEnd==0) { // no clipping
            indexOfFirstNonClippingOperation = fromStart ? 0 : cigarElements.size()-1;
        } else if (numHardClippingBasesFromOneEnd==0 || numSoftClippingBasesFromOneEnd==0) { // one clipping
            indexOfFirstNonClippingOperation = fromStart ? 1 : cigarElements.size()-2;
        } else {
            indexOfFirstNonClippingOperation = fromStart ? 2 : cigarElements.size()-3;
        }

        final CigarElement element = cigarElements.get(indexOfFirstNonClippingOperation);
        if (element.getOperator() == CigarOperator.I) {

            cigarElements.set(indexOfFirstNonClippingOperation, new CigarElement(element.getLength(), CigarOperator.S));

            return compactifyNeighboringSoftClippings(cigarElements);
        } else {
            return cigarElements;
        }
    }

    /**
     * Compactify two neighboring soft clippings, one of which was converted from an insertion operation.
     * @return the compactified list of operations
     * @throws IllegalArgumentException if there's un-handled edge case where two operations neighboring each other have
     *                                  the same operator (other than 'S') but for some reason was not compactified into one
     */
    @VisibleForTesting
    public static List<CigarElement> compactifyNeighboringSoftClippings(final List<CigarElement> cigarElements) {
        final List<CigarElement> result = new ArrayList<>(cigarElements.size());
        for (final CigarElement element : cigarElements) {
            final int idx = result.size()-1;
            if (result.isEmpty() || result.get(idx).getOperator()!=element.getOperator()) {
                result.add(element);
            } else {
                Utils.validateArg(result.get(idx).getOperator()==CigarOperator.S && element.getOperator()==CigarOperator.S,
                        "Seeing new edge case where two neighboring operations are having the same operator: " + cigarElements.toString());
                result.set(idx, new CigarElement(result.get(idx).getLength()+element.getLength(), CigarOperator.S));
            }
        }
        return result;
    }
}
