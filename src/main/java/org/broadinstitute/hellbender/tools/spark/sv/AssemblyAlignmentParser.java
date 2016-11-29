package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.util.SequenceUtil;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import scala.Tuple2;

import java.io.Serializable;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

import static org.broadinstitute.hellbender.tools.spark.sv.ContigsCollection.loadContigsCollectionKeyedByAssemblyId;

/**
 * Implements a parser for parsing alignments of locally assembled contigs, and generating {@link ChimericAlignment}'s
 * from the alignments, if possible.
 */
class AssemblyAlignmentParser implements Serializable {
    private static final long serialVersionUID = 1L;

    static final int ARTIFICIAL_CIGAR_MISMATCH = -1;

    // TODO: 11/23/16 test
    /**
     * Loads the alignment regions and sequence of all locally-assembled contigs from the text file they are in;
     * one record for each contig.
     * @param ctx                       spark context for IO operations
     * @param pathToInputAlignments     path string to alignments of the contigs; format assumed to be consistent/parsable by {@link AlignmentRegion#toString()}
     * @param pathToInputAssemblies     path string to assembled contigs; format assumed to be consistent/parsable by {@link ContigsCollection#loadContigsCollectionKeyedByAssemblyId(JavaSparkContext, String)}
     * @return                          an PairRDD for all assembled contigs with their alignment regions and sequence
     */
    static JavaPairRDD<Iterable<AlignmentRegion>, byte[]> prepAlignmentRegionsForCalling(final JavaSparkContext ctx,
                                                                                         final String pathToInputAlignments,
                                                                                         final String pathToInputAssemblies) {

        final JavaPairRDD<Tuple2<String, String>, Iterable<AlignmentRegion>> alignmentRegionsKeyedByAssemblyAndContigId
                = ctx.textFile(pathToInputAlignments).map(textLine -> AlignmentRegion.fromString(textLine.split(AlignmentRegion.STRING_REP_SEPARATOR,-1)))
                .flatMap(oneRegion -> breakGappedAlignment(oneRegion, SVConstants.GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY))
                .mapToPair(alignmentRegion -> new Tuple2<>(new Tuple2<>(alignmentRegion.assemblyId, alignmentRegion.contigId), alignmentRegion))
                .groupByKey();


        final JavaPairRDD<Tuple2<String, String>, byte[]> contigSequences
                = loadContigsCollectionKeyedByAssemblyId(ctx, pathToInputAssemblies).flatMapToPair(assemblyIdAndContigsCollection -> {
                    final String assemblyId = assemblyIdAndContigsCollection._1;
                    final ContigsCollection contigsCollection = assemblyIdAndContigsCollection._2;
                    return contigsCollection.getContents().stream().map(pair -> new Tuple2<>(new Tuple2<>(assemblyId, pair._1.toString()), pair._2.toString().getBytes())).collect(Collectors.toList());
                });

        return alignmentRegionsKeyedByAssemblyAndContigId
                .join(contigSequences)
                .mapToPair(Tuple2::_2);
    }

    // TODO: 11/26/16 assuming 'I'/'D' cannot be the beginning or last CIGAR operator AND the following operator cannot be a clipping
    /**
     * Break a gapped alignment into multiple alignment regions, based on input sensitivity: i.e. if the gap(s) in the input alignment
     * is equal to or larger than the size specified, two alignment regions will be generated.
     * @return an iterable of size >= 1. if size==1, the returned iterable contains only the input (i.e. either no gap or hasn't reached sensitivity)
     */
    @VisibleForTesting
    static Iterable<AlignmentRegion> breakGappedAlignment(final AlignmentRegion oneRegion, final int sensitivity) {

        final List<CigarElement> cigarElements = oneRegion.forwardStrandCigar.getCigarElements();
        if(cigarElements.size()==1) return new ArrayList<>( Collections.singletonList(oneRegion) );
        final int gapCount = (int) cigarElements.stream().filter(ce -> ce.getOperator().isIndel() && ce.getLength() >= sensitivity).count();
        if (gapCount==0) return new ArrayList<>( Collections.singletonList(oneRegion) );

        // get necessary information for computing
        final String asmId = oneRegion.assemblyId;
        final String contigId = oneRegion.contigId;
        final String chr = oneRegion.referenceInterval.getContig();
        final Cigar originalForwardStrandCigar = oneRegion.forwardStrandCigar;
        final boolean isForwardStrand = oneRegion.forwardStrand;
        final int newMapQual = oneRegion.mapQual;
        final int refStart = oneRegion.referenceInterval.getStart();
        final int refEnd   = oneRegion.referenceInterval.getEnd();
        final int contigLength = oneRegion.assembledContigLength;

        // needs to accomplish three key functionality: 1. generate appropriate cigar, 2. infer reference coordinates, 3. infer contig coordinates, (cannot do yet: 4. mismatches and mapping quality)
        final List<AlignmentRegion> result = new ArrayList<>(gapCount+1);

        final List<CigarElement> cigarMemoryList = new ArrayList<>();
        final int clippedNBases = SVVariantCallerUtils.getNumClippedBases(true, originalForwardStrandCigar);
        final int hardClippingAtEnd = (cigarElements.get(cigarElements.size()-1).getOperator()== CigarOperator.H)? cigarElements.get(cigarElements.size()-1).getLength() : 0;
        int refIntervalStart = refStart, contigIntervalStart = 1 + clippedNBases;
        int gapIndex=0;
        for(int i=0; i<cigarElements.size(); ++i){ // considers only 'M', 'I', 'D', 'S', and 'H'
            final CigarElement current = cigarElements.get(i);
            final CigarOperator op = current.getOperator();
            final int operatorLen = current.getLength();
            switch (op){
                case M: case S: case H:
                    cigarMemoryList.add(current);
                    break;
                case I: case D:
                    if (operatorLen<sensitivity) {
                        cigarMemoryList.add(current);
                        break;
                    }

                    // collapse cigar memory list into a single cigar
                    final List<CigarElement> cigarMemoryListModifiable = new ArrayList<>(cigarMemoryList);
                    final Cigar memoryCigar = new Cigar(cigarMemoryListModifiable);
                    final int effectiveReadLen = memoryCigar.getReadLength() + SVVariantCallerUtils.getTotalHardClipping(memoryCigar) - SVVariantCallerUtils.getNumClippedBases(true, memoryCigar);

                    // infer reference interval
                    final SimpleInterval referenceInterval = new SimpleInterval(chr, refIntervalStart, refIntervalStart + memoryCigar.getReferenceLength()-1);

                    // infer contig interval
                    final int contigIntervalEnd = contigIntervalStart + effectiveReadLen - 1;

                    // now add trailing cigar element and create the real cigar
                    cigarMemoryListModifiable.add(new CigarElement(contigLength-contigIntervalEnd-hardClippingAtEnd, CigarOperator.S));
                    if(hardClippingAtEnd!=0) cigarMemoryListModifiable.add(new CigarElement(hardClippingAtEnd, CigarOperator.H));
                    final Cigar forwardStrandCigar = new Cigar(cigarMemoryListModifiable);

                    final AlignmentRegion split = new AlignmentRegion(asmId, contigId, referenceInterval, forwardStrandCigar, isForwardStrand, newMapQual, ARTIFICIAL_CIGAR_MISMATCH, contigIntervalStart, contigIntervalEnd);

                    result.add(split);

                    // update pointers into reference and contig
                    final boolean isInsertion = op==CigarOperator.I;
                    refIntervalStart     += isInsertion ? memoryCigar.getReferenceLength()  : memoryCigar.getReferenceLength() + operatorLen;
                    contigIntervalStart  += isInsertion ? effectiveReadLen + operatorLen    : effectiveReadLen;
                    cigarMemoryList.clear();    // update cigar memory (note: might not be an efficient frequent operation on linked-list, at least in theory, in practice, might be ok)
                    final int hardClippingAtBeginning = SVVariantCallerUtils.getTotalHardClipping(memoryCigar);
                    if (hardClippingAtBeginning!=0) cigarMemoryList.add(new CigarElement(hardClippingAtBeginning, CigarOperator.H));
                    cigarMemoryList.add(new CigarElement(contigIntervalEnd - hardClippingAtBeginning + (isInsertion ? operatorLen : 0), CigarOperator.S));

                    if (gapCount==++gapIndex) { // last gap
                        final List<CigarElement> trailingCE = new ArrayList<>( cigarElements.subList(i+1, cigarElements.size()) );
                        if (hardClippingAtBeginning!=0){
                            trailingCE.add(0, new CigarElement(hardClippingAtBeginning, CigarOperator.H));
                            trailingCE.add(1, new CigarElement(contigIntervalStart-hardClippingAtBeginning-1, CigarOperator.S));
                        } else {
                            trailingCE.add(0, new CigarElement(contigIntervalStart-1, CigarOperator.S));
                        }
                        final SimpleInterval lastReferenceInterval =  new SimpleInterval(chr, refIntervalStart, refEnd);
                        final Cigar lastForwardStrandCigar = new Cigar(trailingCE);
                        result.add(new AlignmentRegion(asmId, contigId, lastReferenceInterval, lastForwardStrandCigar,
                                isForwardStrand, newMapQual, ARTIFICIAL_CIGAR_MISMATCH,
                                contigIntervalStart, contigLength-SVVariantCallerUtils.getNumClippedBases(false, originalForwardStrandCigar)));
                    }

                    break;
                default:
                    cigarMemoryList.add(current);
                    break;
            }
        }

        return result;
    }

    /**
     * First step in calling variants: parse all alignment records for a single locally-assembled contig and generate
     * chimeric alignments if available.
     * Applies certain filters to the alignment regions:
     *     1) if the alignment region has too low a mapping quality, it is simply skipped
     *     2) if the alignment region is too small, it is skipped
     *     3) if the alignment region passes the above two filters and the next alignment region could be treated as potential inserted sequence, note down the mapping & alignment information of that region and skip it
     *
     * @param input     made of ((assemblyId, contigId), ({alignmentRegions}, sequence)) of a signalling locally-assembled contig
     *
     * @return          the chimeric alignments of this sequence (empty if the sequence does not have any alignments)
     */
    @VisibleForTesting
    static List<ChimericAlignment> getChimericAlignmentsFromAlignmentRegions(Tuple2<Iterable<AlignmentRegion>, byte[]> input) {

        final byte[] sequence = input._2;
        final List<AlignmentRegion> alignmentRegionList = StreamSupport.stream(input._1().spliterator(), false).sorted(Comparator.comparing(a -> a.startInAssembledContig)).collect(Collectors.toList());
        if (alignmentRegionList.size() < 2) {
            return new ArrayList<>();
        }

        final List<ChimericAlignment> results = new ArrayList<>(alignmentRegionList.size() - 1);
        final List<String> insertionAlignmentRegions = new ArrayList<>();

        final Iterator<AlignmentRegion> iterator = alignmentRegionList.iterator();

        // fast forward to the first alignment region with high MapQ
        AlignmentRegion current = iterator.next();
        while (mapQualTooLow(current) && iterator.hasNext()) {
            current = iterator.next();
        }

        while ( iterator.hasNext() ) {
            final AlignmentRegion next = iterator.next();
            if (currentAlignmentIsTooShort(current, next, SVConstants.DEFAULT_MIN_ALIGNMENT_LENGTH)) {
                continue;
            } else if (nextAlignmentMayBeNovelInsertion(current, next, SVConstants.DEFAULT_MIN_ALIGNMENT_LENGTH)) {
                if (iterator.hasNext()) {
                    insertionAlignmentRegions.add(next.toPackedString());
                    continue;
                } else {
                    break;
                }
            }

            results.add(new ChimericAlignment(current, next, getHomology(current, next, sequence), getInsertedSequence(current, next, sequence), insertionAlignmentRegions));

            current = next;
        }

        return results;
    }

    // TODO: 11/22/16 it might also be suitable to consider the reference context this alignment region is mapped to and not simply apply a hard filter (need to think about how to test)
    static boolean mapQualTooLow(final AlignmentRegion next) {
        return next.mapQual < SVConstants.CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD;
    }

    @VisibleForTesting
    static boolean currentAlignmentIsTooShort(final AlignmentRegion current, final AlignmentRegion next, final Integer minAlignLength) {
        return current.referenceInterval.size() - SVVariantCallerUtils.overlapOnContig(current, next) < minAlignLength;
    }

    /**
     * To implement the idea that for two consecutive alignment regions of a contig, the one with higher reference coordinate might be a novel insertion.
     */
    @VisibleForTesting
    static boolean nextAlignmentMayBeNovelInsertion(final AlignmentRegion current, final AlignmentRegion next, final Integer minAlignLength) {
        return mapQualTooLow(next) ||                               // inserted sequence might have low mapping quality
                currentAlignmentIsTooShort(next, current, minAlignLength) ||     // inserted sequence might be very small
                current.referenceInterval.contains(next.referenceInterval) ||   // one might completely contain the other
                next.referenceInterval.contains(current.referenceInterval);
    }

    /**
     * @return Micro-homology sequence using two alignments of the same contig: as indicated by their overlap on the contig itself.
     *          Empty if they don't overlap on the contig.
     */
    @VisibleForTesting
    static String getHomology(final AlignmentRegion current, final AlignmentRegion next, final byte[] sequence) {

        if (current.endInAssembledContig >= next.startInAssembledContig) {
            final byte[] homologyBytes = Arrays.copyOfRange(sequence, next.startInAssembledContig-1, current.endInAssembledContig);
            if (current.referenceInterval.getStart() > next.referenceInterval.getStart()) {
                SequenceUtil.reverseComplement(homologyBytes, 0, homologyBytes.length);
            }
            return new String(homologyBytes);
        } else {
            return "";
        }
    }

    /**
     * @return Inserted sequence using two alignments of the same contig: as indicated by their separation on the the contig itself.
     */
    @VisibleForTesting
    static String getInsertedSequence(final AlignmentRegion current, final AlignmentRegion next, final byte[] sequence) {

        if (current.endInAssembledContig < next.startInAssembledContig-1) {
            final byte[] insertedSequenceBytes = Arrays.copyOfRange(sequence, current.endInAssembledContig, next.startInAssembledContig-1);
            if (current.referenceInterval.getStart() > next.referenceInterval.getStart()) {
                SequenceUtil.reverseComplement(insertedSequenceBytes, 0, insertedSequenceBytes.length);
            }
            return new String(insertedSequenceBytes);
        } else {
            return "";
        }
    }
}
