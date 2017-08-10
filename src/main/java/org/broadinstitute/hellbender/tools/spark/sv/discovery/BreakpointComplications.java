package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import scala.Tuple2;

import java.util.*;

/**
 * A helper struct for annotating complications that make the locations represented by its associated
 * {@link NovelAdjacencyReferenceLocations} a little ambiguous,
 * so that downstream analysis could infer sv type with these complications.
 * To be updated as more types of complications can be processed and handled by
 * {@link BreakpointComplications( ChimericAlignment )}.
 */
@DefaultSerializer(BreakpointComplications.Serializer.class)
public final class BreakpointComplications {

    @SuppressWarnings("unchecked")
    private static final List<String> DEFAULT_CIGAR_STRINGS_FOR_DUP_SEQ_ON_CTG = new ArrayList<>(Collections.EMPTY_LIST);

    /**
     * '+' strand representations of micro-homology, inserted sequence and duplicated sequence on the reference.
     */
    private String homologyForwardStrandRep = "";
    private String insertedSequenceForwardStrandRep = "";

    private boolean hasDuplicationAnnotation = false;

    private SimpleInterval dupSeqRepeatUnitRefSpan = null;
    private int dupSeqRepeatNumOnRef = 0;
    private int dupSeqRepeatNumOnCtg = 0;
    private List<String> cigarStringsForDupSeqOnCtg = null;
    private boolean dupAnnotIsFromOptimization = false;

    /**
     * @return Intended for use in debugging and exception message only.
     */
    @Override
    public String toString() {
        String toPrint = "homology: " + homologyForwardStrandRep + "\tinserted sequence: " + insertedSequenceForwardStrandRep;

        if (hasDuplicationAnnotation()) {
            toPrint += String.format("\ttandem duplication repeat unit ref span: %s\tref repeat num: %d\tctg repeat num: %d\tcigarStringsForDupSeqOnCtg: %s\ttandupAnnotationIsFromSimpleOptimization: %s",
                    dupSeqRepeatUnitRefSpan, dupSeqRepeatNumOnRef, dupSeqRepeatNumOnCtg, cigarStringsForDupSeqOnCtg, isDupAnnotIsFromOptimization() ? "true" : "false");
        }
        return toPrint;
    }

    boolean hasDuplicationAnnotation() {
        return hasDuplicationAnnotation;
    }

    String getHomologyForwardStrandRep() {
        return homologyForwardStrandRep;
    }

    String getInsertedSequenceForwardStrandRep() {
        return insertedSequenceForwardStrandRep;
    }

    SimpleInterval getDupSeqRepeatUnitRefSpan() {
        return dupSeqRepeatUnitRefSpan;
    }

    int getDupSeqRepeatNumOnRef() {
        return dupSeqRepeatNumOnRef;
    }

    int getDupSeqRepeatNumOnCtg() {
        return dupSeqRepeatNumOnCtg;
    }

    List<String> getCigarStringsForDupSeqOnCtg() {
        return cigarStringsForDupSeqOnCtg;
    }

    boolean isDupAnnotIsFromOptimization() {
        return dupAnnotIsFromOptimization;
    }

    @VisibleForTesting
    BreakpointComplications() {

    }

    // For test purposes only to cover the deficiency in the default initializations
    @VisibleForTesting
    BreakpointComplications(final String homologyForwardStrandRep, final String insertedSequenceForwardStrandRep,
                            final boolean hasDuplicationAnnotation, final SimpleInterval dupSeqRepeatUnitRefSpan,
                            final int dupSeqRepeatNumOnRef, final int dupSeqRepeatNumOnCtg,
                            final List<String> cigarStringsForDupSeqOnCtg, final boolean dupAnnotIsFromOptimization) {
        this.homologyForwardStrandRep = homologyForwardStrandRep;
        this.insertedSequenceForwardStrandRep = insertedSequenceForwardStrandRep;
        this.hasDuplicationAnnotation = hasDuplicationAnnotation;
        this.dupSeqRepeatUnitRefSpan = dupSeqRepeatUnitRefSpan;
        this.dupSeqRepeatNumOnRef = dupSeqRepeatNumOnRef;
        this.dupSeqRepeatNumOnCtg = dupSeqRepeatNumOnCtg;
        this.cigarStringsForDupSeqOnCtg = cigarStringsForDupSeqOnCtg;
        this.dupAnnotIsFromOptimization = dupAnnotIsFromOptimization;
    }

    /**
     * Given an {@link ChimericAlignment} representing two reference intervals rearranged as two intervals on the locally-assembled contig,
     * identify potential complications such as homology and duplication on the reference and/or on the contig.
     */
    BreakpointComplications(final ChimericAlignment chimericAlignment, final byte[] contigSeq) {
        final Tuple2<SimpleInterval, SimpleInterval> referenceSpans = chimericAlignment.getCoordSortedReferenceSpans();
        final SimpleInterval leftReferenceSpan  = referenceSpans._1;
        final SimpleInterval rightReferenceSpan = referenceSpans._2;

        final AlignmentInterval firstAlignmentInterval  = chimericAlignment.regionWithLowerCoordOnContig;
        final AlignmentInterval secondAlignmentInterval = chimericAlignment.regionWithHigherCoordOnContig;

        // a segment with lower coordinate on the locally-assembled contig could map to a higher reference coordinate region
        // under two basic types of SV's: inversion (strand switch necessary) and translocation (no strand switch necessary)
        final boolean isNotSimpleTranslocation =
                ChimericAlignment.isNotSimpleTranslocation(chimericAlignment.regionWithLowerCoordOnContig,
                        chimericAlignment.regionWithHigherCoordOnContig, chimericAlignment.strandSwitch,
                        ChimericAlignment.involvesRefPositionSwitch(firstAlignmentInterval, secondAlignmentInterval));

        // TODO: 12/5/16 simple translocation, don't tackle yet
        if (chimericAlignment.strandSwitch!= StrandSwitch.NO_SWITCH) { // the case involves an inversion
            // TODO: 12/5/16 duplication detection to be done for inversion alleles
            initForSimpleInversion(firstAlignmentInterval, secondAlignmentInterval, contigSeq);
        } else if (isNotSimpleTranslocation) {
            initForInsDel(chimericAlignment, leftReferenceSpan, rightReferenceSpan, contigSeq);
        }
    }

    private void initForSimpleInversion(final AlignmentInterval firstContigRegion, final AlignmentInterval secondContigRegion,
                                        final byte[] contigSeq) {
        homologyForwardStrandRep = getHomology(firstContigRegion, secondContigRegion, contigSeq);
        insertedSequenceForwardStrandRep = getInsertedSequence(firstContigRegion, secondContigRegion, contigSeq);
        dupSeqRepeatUnitRefSpan = null;
        dupSeqRepeatNumOnRef = dupSeqRepeatNumOnCtg = 0;
        cigarStringsForDupSeqOnCtg = null;
        dupAnnotIsFromOptimization = false;
        hasDuplicationAnnotation = false;
    }

    /**
     * todo : see ticket #3529
     * @return true iff the two AI of the {@code longRead} overlaps on reference is more than half of the two AI's minimal read span.
     */
    @VisibleForTesting
    public static boolean isLikelyInvertedDuplication(final AlignmentInterval one, final AlignmentInterval two) {
        return 2 * AlignmentInterval.overlapOnRefSpan(one, two) >
                Math.min(one.endInAssembledContig - one.startInAssembledContig,
                        two.endInAssembledContig - two.startInAssembledContig) + 1;
    }

    private void initForInsDel(final ChimericAlignment chimericAlignment, final SimpleInterval leftReferenceInterval,
                               final SimpleInterval rightReferenceInterval, final byte[] contigSeq) {

        final AlignmentInterval firstContigRegion  = chimericAlignment.regionWithLowerCoordOnContig;
        final AlignmentInterval secondContigRegion = chimericAlignment.regionWithHigherCoordOnContig;

        final int r1e = leftReferenceInterval.getEnd(),
                r2b = rightReferenceInterval.getStart(),
                c1e = firstContigRegion.endInAssembledContig,
                c2b = secondContigRegion.startInAssembledContig;

        final int distBetweenAlignRegionsOnRef = r2b - r1e - 1, // distance-1 between the two regions on reference, denoted as d1 in the comments below
                distBetweenAlignRegionsOnCtg = c2b - c1e - 1; // distance-1 between the two regions on contig, denoted as d2 in the comments below

        if ( distBetweenAlignRegionsOnRef > 0 ) {        // Deletion:
            resolveComplicationForSimpleDel(firstContigRegion, secondContigRegion, distBetweenAlignRegionsOnCtg, contigSeq);
        } else if (distBetweenAlignRegionsOnRef == 0 && distBetweenAlignRegionsOnCtg > 0) { // Insertion: simple insertion, inserted sequence is the sequence [c1e+1, c2b-1] on the contig
            insertedSequenceForwardStrandRep = getInsertedSequence(firstContigRegion, secondContigRegion, contigSeq);
        } else if (distBetweenAlignRegionsOnRef == 0 && distBetweenAlignRegionsOnCtg < 0) { // Tandem repeat contraction: reference has two copies but one copy was deleted on the contig; duplicated sequence on reference are [r1e-|d2|+1, r1e] and [r2b, r2b+|d2|-1]
            resolveComplicationForSimpleTandupContraction(leftReferenceInterval, firstContigRegion, secondContigRegion, r1e, c1e, c2b, contigSeq);
        } else if (distBetweenAlignRegionsOnRef < 0 && distBetweenAlignRegionsOnCtg >= 0) { // Tandem repeat expansion:   reference bases [r1e-|d1|+1, r1e] to contig bases [c1e-|d1|+1, c1e] and [c2b, c2b+|d1|-1] with optional inserted sequence [c1e+1, c2b-1] in between the two intervals on contig
            resolveComplicationForSimpleTandupExpansion(leftReferenceInterval, firstContigRegion, secondContigRegion, r1e, r2b, distBetweenAlignRegionsOnCtg, contigSeq);
        } else if (distBetweenAlignRegionsOnRef < 0 && distBetweenAlignRegionsOnCtg < 0) {  // most complicated case, see below
            // Deletion:  duplication with repeat number N1 on reference, N2 on contig, such that N1 <= 2*N2 (and N2<N1);
            // Insertion: duplication with repeat number N1 on reference, N2 on contig, such that N2 <= 2*N1 (and N1<N2);
            // in both cases, the equal sign on the right can be taken only when there's pseudo-homology between starting bases of the duplicated sequence and starting bases of the right flanking region
            // the reference system with a shorter overlap (i.e. with less-negative distance between regions) has a higher repeat number
            resolveComplicationForComplexTandup(firstContigRegion, secondContigRegion, r1e, distBetweenAlignRegionsOnRef, distBetweenAlignRegionsOnCtg, contigSeq);
        } else if (distBetweenAlignRegionsOnRef == 0 && distBetweenAlignRegionsOnCtg == 0) {// SNP & indel
            throw new GATKException("Detected badly parsed chimeric alignment for identifying SV breakpoints; no rearrangement found: " + chimericAlignment.onErrStringRep());
        }

        if ( insertedSequenceForwardStrandRep.isEmpty() ){
            if ( dupSeqRepeatNumOnCtg != dupSeqRepeatNumOnRef && null == dupSeqRepeatUnitRefSpan )
                throw new GATKException("An identified breakpoint pair seem to suggest insertion but the inserted sequence is empty: " + chimericAlignment.onErrStringRep());
        }
    }

    private void resolveComplicationForSimpleTandupExpansion(final SimpleInterval leftReferenceInterval,
                                                             final AlignmentInterval firstContigRegion,
                                                             final AlignmentInterval secondContigRegion,
                                                             final int r1e, final int r2b,
                                                             final int distBetweenAlignRegionsOnCtg, final byte[] contigSeq) {
        // note this does not incorporate the duplicated reference sequence
        insertedSequenceForwardStrandRep = distBetweenAlignRegionsOnCtg == 0 ? "" : getInsertedSequence(firstContigRegion, secondContigRegion, contigSeq);
        hasDuplicationAnnotation  = true;
        dupSeqRepeatUnitRefSpan   = new SimpleInterval(leftReferenceInterval.getContig(), r2b, r1e);
        dupSeqRepeatNumOnRef      = 1;
        dupSeqRepeatNumOnCtg      = 2;
        cigarStringsForDupSeqOnCtg = new ArrayList<>(2);
        if (firstContigRegion.forwardStrand) {
            cigarStringsForDupSeqOnCtg.add( TextCigarCodec.encode(extractCigarForTandup(firstContigRegion, r1e, r2b)) );
            cigarStringsForDupSeqOnCtg.add( TextCigarCodec.encode(extractCigarForTandup(secondContigRegion, r1e, r2b)) );
        } else {
            cigarStringsForDupSeqOnCtg.add( TextCigarCodec.encode(CigarUtils.invertCigar(extractCigarForTandup(firstContigRegion, r1e, r2b))) );
            cigarStringsForDupSeqOnCtg.add( TextCigarCodec.encode(CigarUtils.invertCigar(extractCigarForTandup(secondContigRegion, r1e, r2b))) );
        }
    }

    private void resolveComplicationForSimpleTandupContraction(final SimpleInterval leftReferenceInterval,
                                                               final AlignmentInterval firstContigRegion,
                                                               final AlignmentInterval secondContigRegion,
                                                               final int r1e, final int c1e, final int c2b,
                                                               final byte[] contigSeq) {
        homologyForwardStrandRep = getHomology(firstContigRegion, secondContigRegion, contigSeq);
        hasDuplicationAnnotation = true;
        dupSeqRepeatUnitRefSpan  = new SimpleInterval(leftReferenceInterval.getContig(), r1e - ( c1e - c2b ), r1e);
        dupSeqRepeatNumOnRef     = 2;
        dupSeqRepeatNumOnCtg     = 1;
        cigarStringsForDupSeqOnCtg = DEFAULT_CIGAR_STRINGS_FOR_DUP_SEQ_ON_CTG;
    }

    private void resolveComplicationForSimpleDel(final AlignmentInterval firstContigRegion,
                                                 final AlignmentInterval secondContigRegion,
                                                 final int distBetweenAlignRegionsOnCtg, final byte[] contigSeq) {
        if (distBetweenAlignRegionsOnCtg>=0) {
            // either: a clean deletion, deleted sequence is [r1e+1, r2b-1] on the reference
            // or    : deletion with scar, i.e. large non-conserved substitution, reference bases [r1e+1, r2b-1] is substituted with contig bases [c1e+1, c2b-1]
            insertedSequenceForwardStrandRep = getInsertedSequence(firstContigRegion, secondContigRegion, contigSeq);
        } else {
            // a sequence of bases of length d1+HOM is deleted, and there's homology (which could be dup, but cannot tell): leftFlank+HOM+[r1e+1, r2b-1]+HOM+rightFlank -> leftFlank+HOM+rightFlank
            homologyForwardStrandRep    = getHomology(firstContigRegion, secondContigRegion, contigSeq);
        }
    }

    private void resolveComplicationForComplexTandup(final AlignmentInterval firstContigRegion,
                                                     final AlignmentInterval secondContigRegion,
                                                     final int r1e, final int distBetweenAlignRegionsOnRef,
                                                     final int distBetweenAlignRegionsOnCtg, final byte[] contigSeq) {

        final TandemRepeatStructure duplicationComplication
                = new TandemRepeatStructure(distBetweenAlignRegionsOnRef, distBetweenAlignRegionsOnCtg);

        final boolean isExpansion     = distBetweenAlignRegionsOnRef<distBetweenAlignRegionsOnCtg;

        final int repeatUnitSpanStart = r1e - duplicationComplication.pseudoHomologyLen
                                            - duplicationComplication.repeatedSeqLen * duplicationComplication.lowerRepeatNumberEstimate
                                            + 1;
        final int repeatUnitSpanEnd   = repeatUnitSpanStart + duplicationComplication.repeatedSeqLen - 1;
        homologyForwardStrandRep      = getHomology(firstContigRegion, secondContigRegion, contigSeq);
        hasDuplicationAnnotation      = true;
        cigarStringsForDupSeqOnCtg    = DEFAULT_CIGAR_STRINGS_FOR_DUP_SEQ_ON_CTG;
        dupSeqRepeatUnitRefSpan       = new SimpleInterval(firstContigRegion.referenceSpan.getContig(), repeatUnitSpanStart, repeatUnitSpanEnd);
        dupSeqRepeatNumOnRef          = isExpansion ? duplicationComplication.lowerRepeatNumberEstimate
                                                    : duplicationComplication.higherRepeatNumberEstimate;
        dupSeqRepeatNumOnCtg          = isExpansion ? duplicationComplication.higherRepeatNumberEstimate
                                                    : duplicationComplication.lowerRepeatNumberEstimate;
        dupAnnotIsFromOptimization    = true;
    }


    /**
     * Given a {@link AlignmentInterval} from a pair of ARs that forms a {@link ChimericAlignment} signalling a tandem duplication,
     * extract a CIGAR from the {@link AlignmentInterval#cigarAlong5to3DirectionOfContig}
     * that corresponds to the alignment between the suspected repeated sequence on reference between
     * [{@code alignmentIntervalTwoReferenceIntervalSpanBegin}, {@code alignmentIntervalOneReferenceIntervalSpanEnd}],
     * and the sequence in {@link AlignmentInterval#referenceSpan}.
     */
    @VisibleForTesting
    static Cigar extractCigarForTandup(final AlignmentInterval contigRegion,
                                       final int alignmentIntervalOneReferenceIntervalSpanEnd,
                                       final int alignmentIntervalTwoReferenceIntervalSpanBegin) {

        final List<CigarElement> elementList = contigRegion.cigarAlong5to3DirectionOfContig.getCigarElements();
        final List<CigarElement> result = new ArrayList<>(elementList.size());
        final int refStart = contigRegion.referenceSpan.getStart(),
                refEnd = contigRegion.referenceSpan.getEnd();
        final boolean isForwardStrand = contigRegion.forwardStrand;
        boolean initiatedCollection = false;
        int refPos = isForwardStrand ? refStart : refEnd;
        for(final CigarElement cigarElement : elementList) {
            final CigarOperator operator = cigarElement.getOperator();
            if ( !operator.isClipping() ) {
                final int opLen = cigarElement.getLength();
                refPos += operator.consumesReferenceBases() ? (isForwardStrand ? opLen : -opLen) : 0;
                final int offsetIntoRepeatRegion = isForwardStrand ? refPos - alignmentIntervalTwoReferenceIntervalSpanBegin
                                                                   : alignmentIntervalOneReferenceIntervalSpanEnd - refPos;
                final int overshootOutOfRepeatRegion = isForwardStrand ? refPos - alignmentIntervalOneReferenceIntervalSpanEnd - 1
                                                                       : alignmentIntervalTwoReferenceIntervalSpanBegin - refPos - 1;

                if ( offsetIntoRepeatRegion > 0 ) {
                    if ( overshootOutOfRepeatRegion <= 0 ) {
                        result.add( initiatedCollection ? cigarElement : new CigarElement(offsetIntoRepeatRegion, operator));
                        initiatedCollection = true;
                    } else {
                        result.add(new CigarElement(opLen-overshootOutOfRepeatRegion, operator));
                        break;
                    }
                }
            }
        }

        return new Cigar(result);
    }

    /**
     * @return Micro-homology sequence using two alignments of the same contig: as indicated by their overlap on the contig itself.
     *          Empty if they don't overlap on the contig.
     */
    @VisibleForTesting
    static String getHomology(final AlignmentInterval current, final AlignmentInterval next, final byte[] contigSequence) {

        if (current.endInAssembledContig >= next.startInAssembledContig) {
            final byte[] homologyBytes = Arrays.copyOfRange(contigSequence,
                    next.startInAssembledContig-1, current.endInAssembledContig);
            if (current.referenceSpan.getStart() > next.referenceSpan.getStart()) {
                SequenceUtil.reverseComplement(homologyBytes, 0, homologyBytes.length);
            }
            return new String(homologyBytes);
        } else {
            return "";
        }
    }

    /**
     * Note: not suitable for the most complicated case dealt with in {@link BreakpointComplications( ChimericAlignment )}
     * @return Inserted sequence using two alignments of the same contig: as indicated by their separation on the the contig itself.
     */
    @VisibleForTesting
    static String getInsertedSequence(final AlignmentInterval current, final AlignmentInterval next, final byte[] contigSequence) {

        if (current.endInAssembledContig < next.startInAssembledContig - 1) {
            final byte[] insertedSequenceBytes = Arrays.copyOfRange(contigSequence,
                    current.endInAssembledContig, next.startInAssembledContig - 1);
            if (current.referenceSpan.getStart() > next.referenceSpan.getStart()) {
                SequenceUtil.reverseComplement(insertedSequenceBytes, 0, insertedSequenceBytes.length);
            }
            return new String(insertedSequenceBytes);
        } else {
            return "";
        }
    }

    // TODO: 03/03/17 this complicated tandem duplication annotation is not exactly reproducible in the following sense:
    //          1) depending on what the assembler might produce, e.g. different runs producing slightly different sequences
    //          hence affecting alignment,
    //          2) the assembler might decide to output RC sequences between runs hence the mapping would be to '+' or '-' strand
    //       these randomness may give slightly different results by this treatment
    /**
     * This auxiliary structure, when constructed given overlaps of two corresponding regions on reference and contig sequences,
     * attempts to find--naively and slowly--the repeat numbers on the reference and on the contig of tandem repeats,
     * as well as the pseudo-homology between the duplicated sequence and the right flanking region.
     *
     * An example might help:
     * an assembled contig that's actually a repeat expansion from 1 repeat to 2 repeats with pseudo-homology:
     * TGCCAGGTTACATGGCAAAGAGGGTAGATATGGGGAGCTGTGAAGAATGGAGCCAGTAATTAAATTCACTGAAGTCTCCACAGGAGGGCAAGGTGGACAATCTGTCCCATAGGAGGGGGATTCATGAGGGGAGCTGTGAAGAATGGAGCCAGTAATTAAATTCACTGAAGTCTCCACAGGAGGGCAAGGTGGACAATCTGTCCCATAGGAGGGGGATTCAGGAGGGCAGCTGTGGATGGTGCAAATGCCATTTATGCTCCTCTCCACCCATATCC
     * can be aligned to chr18,
     * the 1st alignment chr18:312579-718, 140M135S, which can be broken into the following part
     * 31:  TGCCAGGTTACATGGCAAAGAGGGTAGATAT
     * 109: GGGGAGCTGTGAAGAATGGAGCCAGTAATTAAATTCACTGAAGTCTCCACAGGAGGGCAAGGTGGACAATCTGTCCCATAGGAGGGGGATTCATGAGGGGAGCTGTGAA
     * 135: GAATGGAGCCAGTAATTAAATTCACTGAAGTCTCCACAGGAGGGCAAGGTGGACAATCTGTCCCATAGGAGGGGGATTCAGGAGGGCAGCTGTGGATGGTGCAAATGCCATTTATGCTCCTCTCCACCCATATCC
     * And the arithmetic to get the cigar operation length works this way:
     * 31 + 109 = 140
     * 109 = 96 + 13
     * where 31 is the left flanking region before the repeated unit, which itself is 96 bases long (see below),
     * the number 13 is the length of the pseudo-homology between the starting bases of the repeated sequence and the right flanking region
     * a clearer picture emerges when we look at the 2nd alignment
     * chr18:312610-757, 127S148M, which can be broken into
     * 31: TGCCAGGTTACATGGCAAAGAGGGTAGATAT
     * 96: GGGGAGCTGTGAAGAATGGAGCCAGTAATTAAATTCACTGAAGTCTCCACAGGAGGGCAAGGTGGACAATCTGTCCCATAGGAGGGGGATTCATGA
     * 96: GGGGAGCTGTGAAGAATGGAGCCAGTAATTAAATTCACTGAAGTCTCCACAGGAGGGCAAGGTGGACAATCTGTCCCATAGGAGGGGGATTCAGGA
     * 13: GGGCAGCTGTGGA
     * 39: TGGTGCAAATGCCATTTATGCTCCTCTCCACCCATATCC
     * And the arithmetic works this way:
     * 31 + 96 = 127
     * 96 + 13 + 39 = 148
     */
    private static final class TandemRepeatStructure {

        /**
         * In {@link TandemRepeatStructure} where the naive attempt to resolve number of tandem repeats
         * on the reference and sample is done, we assume the lower number of repeats is no higher than this number.
         */
        private static final int MAX_LOWER_CN = 10;

        final int lowerRepeatNumberEstimate;
        final int higherRepeatNumberEstimate;
        final int repeatedSeqLen;
        final int pseudoHomologyLen;


        @VisibleForTesting
        TandemRepeatStructure(final int distBetweenAlignRegionsOnRef, final int distBetweenAlignRegionsOnCtg) {
            // the reference system with a shorter overlap (i.e. with less-negative distance between regions) has a higher repeat number
            final boolean isExpansion = distBetweenAlignRegionsOnRef < distBetweenAlignRegionsOnCtg;
            final int overlapOnLowerCNSequence, overlapOnHigherCNSequence;
            if (isExpansion) {
                overlapOnLowerCNSequence = Math.abs(distBetweenAlignRegionsOnRef);
                overlapOnHigherCNSequence = Math.abs(distBetweenAlignRegionsOnCtg);
            } else {     // d1 is lower absolute value -> reference has higher copy number of the duplication, i.e. Deletion
                overlapOnLowerCNSequence = Math.abs(distBetweenAlignRegionsOnCtg);
                overlapOnHigherCNSequence = Math.abs(distBetweenAlignRegionsOnRef);
            }

            int higherCnEst=0, lowerCnEst=0, unitLen=0, pseudoHomLen=0;
            double err = Double.MAX_VALUE;
            for(int cn2 = 1; cn2< MAX_LOWER_CN; ++cn2) {
                for(int cn1 = cn2 + 1; cn1 <= 2 * cn2; ++cn1) {
                    final int dupLenUpperBound = (cn1 == 2 * cn2) ? overlapOnLowerCNSequence : overlapOnHigherCNSequence;
                    for (int l = 2; l <= dupLenUpperBound; ++l) {
                        for (int lambda = 0; lambda < l; ++lambda) {
                            final int d1 = (2*cn2 - cn1)*l + lambda;
                            final int d2 = cn2*l + lambda;
                            final double newErr = Math.abs(overlapOnHigherCNSequence-d1) + Math.abs(overlapOnLowerCNSequence-d2);
                            if (newErr < err) {
                                err = newErr;
                                higherCnEst = cn1; lowerCnEst = cn2;
                                unitLen= l; pseudoHomLen = lambda;
                            }
                            if (err < 1){
                                lowerRepeatNumberEstimate = lowerCnEst;
                                higherRepeatNumberEstimate = higherCnEst;
                                repeatedSeqLen = unitLen;
                                pseudoHomologyLen = pseudoHomLen;
                                return;
                            }
                        }
                    }
                }
            }

            lowerRepeatNumberEstimate = lowerCnEst;
            higherRepeatNumberEstimate = higherCnEst;
            repeatedSeqLen = unitLen;
            pseudoHomologyLen = pseudoHomLen;
        }
    }

    protected BreakpointComplications(final Kryo kryo, final Input input) {
        homologyForwardStrandRep = input.readString();
        insertedSequenceForwardStrandRep = input.readString();
        hasDuplicationAnnotation = input.readBoolean();
        if (hasDuplicationAnnotation) {
            final String ctg = input.readString();
            final int start = input.readInt();
            final int end = input.readInt();
            dupSeqRepeatUnitRefSpan = new SimpleInterval(ctg, start, end);
            dupSeqRepeatNumOnRef = input.readInt();
            dupSeqRepeatNumOnCtg = input.readInt();
            final int cigarCounts = input.readInt();
            cigarStringsForDupSeqOnCtg = new ArrayList<>(cigarCounts);
            for(int i = 0; i < cigarCounts; ++i) {
                cigarStringsForDupSeqOnCtg.add(input.readString());
            }
            dupAnnotIsFromOptimization = input.readBoolean();
        } else {
            dupSeqRepeatUnitRefSpan = null;
            dupSeqRepeatNumOnRef = 0;
            dupSeqRepeatNumOnCtg = 0;
            cigarStringsForDupSeqOnCtg = null;
            dupAnnotIsFromOptimization = false;
        }
    }

    protected void serialize(final Kryo kryo, final Output output) {
        output.writeString(homologyForwardStrandRep);
        output.writeString(insertedSequenceForwardStrandRep);
        output.writeBoolean(hasDuplicationAnnotation);
        if (hasDuplicationAnnotation) {
            output.writeString(dupSeqRepeatUnitRefSpan.getContig());
            output.writeInt(dupSeqRepeatUnitRefSpan.getStart());
            output.writeInt(dupSeqRepeatUnitRefSpan.getEnd());
            output.writeInt(dupSeqRepeatNumOnRef);
            output.writeInt(dupSeqRepeatNumOnCtg);
            output.writeInt(cigarStringsForDupSeqOnCtg.size());
            cigarStringsForDupSeqOnCtg.forEach(output::writeString);
            output.writeBoolean(dupAnnotIsFromOptimization);
        }
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        BreakpointComplications that = (BreakpointComplications) o;

        if (hasDuplicationAnnotation != that.hasDuplicationAnnotation) return false;
        if (dupSeqRepeatNumOnRef != that.dupSeqRepeatNumOnRef) return false;
        if (dupSeqRepeatNumOnCtg != that.dupSeqRepeatNumOnCtg) return false;
        if (dupAnnotIsFromOptimization != that.dupAnnotIsFromOptimization) return false;
        if (!homologyForwardStrandRep.equals(that.homologyForwardStrandRep)) return false;
        if (!insertedSequenceForwardStrandRep.equals(that.insertedSequenceForwardStrandRep)) return false;
        if (dupSeqRepeatUnitRefSpan != null ? !dupSeqRepeatUnitRefSpan.equals(that.dupSeqRepeatUnitRefSpan) : that.dupSeqRepeatUnitRefSpan != null)
            return false;
        return cigarStringsForDupSeqOnCtg != null ? cigarStringsForDupSeqOnCtg.equals(that.cigarStringsForDupSeqOnCtg)
                : that.cigarStringsForDupSeqOnCtg == null;
    }

    @Override
    public int hashCode() {
        int result = homologyForwardStrandRep.hashCode();
        result = 31 * result + insertedSequenceForwardStrandRep.hashCode();
        result = 31 * result + (hasDuplicationAnnotation ? 1 : 0);
        result = 31 * result + (dupSeqRepeatUnitRefSpan != null ? dupSeqRepeatUnitRefSpan.hashCode() : 0);
        result = 31 * result + dupSeqRepeatNumOnRef;
        result = 31 * result + dupSeqRepeatNumOnCtg;
        result = 31 * result + (cigarStringsForDupSeqOnCtg != null ? cigarStringsForDupSeqOnCtg.hashCode() : 0);
        result = 31 * result + (dupAnnotIsFromOptimization ? 1 : 0);
        return result;
    }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<BreakpointComplications> {
        @Override
        public void write(final Kryo kryo, final Output output, final BreakpointComplications breakpointComplications) {
            breakpointComplications.serialize(kryo, output);
        }

        @Override
        public BreakpointComplications read(final Kryo kryo, final Input input, final Class<BreakpointComplications> klass ) {
            return new BreakpointComplications(kryo, input);
        }
    }
}
