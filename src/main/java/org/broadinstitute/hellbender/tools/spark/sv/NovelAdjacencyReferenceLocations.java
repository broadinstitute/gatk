package org.broadinstitute.hellbender.tools.spark.sv;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import scala.Tuple2;

import java.io.IOException;
import java.util.Arrays;
import java.util.Objects;

import static org.broadinstitute.hellbender.tools.spark.sv.NovelAdjacencyReferenceLocations.EndConnectionType.*;

/**
 * This class represents a pair of inferred genomic locations on the reference whose novel adjacency is generated
 * due to a simple SV event (in other words, a simple rearrangement between two genomic locations)
 * that is suggested by the input {@link ChimericAlignment},
 * and complications in pinning down the locations to exact base pair resolution.
 */
@DefaultSerializer(NovelAdjacencyReferenceLocations.Serializer.class)
class NovelAdjacencyReferenceLocations {

    /**
     * In {@link #resolveComplications(ChimericAlignment)} where the naive attempt to resolve number of tandem repeats
     * on the reference and sample is done, we assume the lower number of repeats is no higher than this number.
     */
    private static final int MAX_LOWER_CN = 10;

    final SimpleInterval leftJustifiedLeftRefLoc;
    final SimpleInterval leftJustifiedRightRefLoc;

    final EndConnectionType endConnectionType;
    final LocationComplication complication;


    /**
     * Represents the strand of evidence that was used in computing this breakpoint pair.
     */
    enum EndConnectionType {
        FIVE_TO_THREE, THREE_TO_THREE, FIVE_TO_FIVE
    }

    /**
     * A helper struct for annotating complications that make the locations represented by this NovelAdjacencyReferenceLocations
     * a little ambiguous, so that downstream analysis could infer sv type with these complications.
     * To be updated as more types of complications can be processed and handled by
     * {@link #resolveComplications(ChimericAlignment)}.
     */
    @DefaultSerializer(LocationComplication.Serializer.class)
    static final class LocationComplication {

        /**
         * '+' strand representations of micro-homology, inserted sequence and duplicated sequence on the reference.
         */
        final String homologyForwardStrandRep;
        final String insertedSequenceForwardStrandRep;
        final String dupSeqForwardStrandRep;
        final int dupSeqRepeatNumOnRef;
        final int dupSeqRepeatNumOnCtg;

        LocationComplication(final String homologyForwardStrandRep,
                             final String insertedSequenceForwardStrandRep,
                             final String dupSeqForwardStrandRep,
                             final int dupSeqRepeatNumOnRef,
                             final int dupSeqRepeatNumOnCtg) {
            this.homologyForwardStrandRep = homologyForwardStrandRep;
            this.insertedSequenceForwardStrandRep = insertedSequenceForwardStrandRep;
            this.dupSeqForwardStrandRep = dupSeqForwardStrandRep;
            this.dupSeqRepeatNumOnRef = dupSeqRepeatNumOnRef;
            this.dupSeqRepeatNumOnCtg = dupSeqRepeatNumOnCtg;
        }

        protected LocationComplication(final Kryo kryo, final Input input) {
            homologyForwardStrandRep = input.readString();
            insertedSequenceForwardStrandRep = input.readString();
            dupSeqForwardStrandRep = input.readString();
            dupSeqRepeatNumOnRef = input.readInt();
            dupSeqRepeatNumOnCtg = input.readInt();
        }

        @Override
        public String toString() {
            return String.format("homology: %s\tinserted sequence: %s\ttandem duplication sequence: %s\tref repeat num: %d\tctg repeat num: %d",
                    homologyForwardStrandRep, insertedSequenceForwardStrandRep, dupSeqForwardStrandRep, dupSeqRepeatNumOnRef, dupSeqRepeatNumOnCtg);
        }

        protected void serialize(final Kryo kryo, final Output output) {
            output.writeString(homologyForwardStrandRep);
            output.writeString(insertedSequenceForwardStrandRep);
            output.writeString(dupSeqForwardStrandRep);
            output.writeInt(dupSeqRepeatNumOnRef);
            output.writeInt(dupSeqRepeatNumOnCtg);
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            LocationComplication that = (LocationComplication) o;

            if (dupSeqRepeatNumOnRef != that.dupSeqRepeatNumOnRef) return false;
            if (dupSeqRepeatNumOnCtg != that.dupSeqRepeatNumOnCtg) return false;
            if (!homologyForwardStrandRep.equals(that.homologyForwardStrandRep)) return false;
            if (!insertedSequenceForwardStrandRep.equals(that.insertedSequenceForwardStrandRep)) return false;
            return dupSeqForwardStrandRep.equals(that.dupSeqForwardStrandRep);
        }

        @Override
        public int hashCode() {
            int result = homologyForwardStrandRep.hashCode();
            result = 31 * result + insertedSequenceForwardStrandRep.hashCode();
            result = 31 * result + dupSeqForwardStrandRep.hashCode();
            result = 31 * result + dupSeqRepeatNumOnRef;
            result = 31 * result + dupSeqRepeatNumOnCtg;
            return result;
        }

        public static final class Serializer extends com.esotericsoftware.kryo.Serializer<LocationComplication> {
            @Override
            public void write(final Kryo kryo, final Output output, final LocationComplication locationComplication) {
                locationComplication.serialize(kryo, output);
            }

            @Override
            public LocationComplication read(final Kryo kryo, final Input input, final Class<LocationComplication> klass ) {
                return new LocationComplication(kryo, input);
            }
        }
    }

    NovelAdjacencyReferenceLocations(final ChimericAlignment chimericAlignment)
            throws IOException {

        // first get endConnectionType, then get complications, finally use complications to justify breakpoints
        endConnectionType = determineEndConnectionType(chimericAlignment);

        complication = resolveComplications(chimericAlignment);

        final Tuple2<SimpleInterval, SimpleInterval> leftJustifiedBreakpoints = leftJustifyBreakpoints(chimericAlignment, this.complication);
        leftJustifiedLeftRefLoc = leftJustifiedBreakpoints._1();
        leftJustifiedRightRefLoc = leftJustifiedBreakpoints._2();
    }

    // TODO: 12/12/16 again, does not work for translocation
    @VisibleForTesting
    static EndConnectionType determineEndConnectionType(final ChimericAlignment chimericAlignment) {
        if (chimericAlignment.regionWithLowerCoordOnContig.forwardStrand == chimericAlignment.regionWithHigherCoordOnContig.forwardStrand) {
            return FIVE_TO_THREE;
        } else {
            return chimericAlignment.regionWithLowerCoordOnContig.forwardStrand ? FIVE_TO_FIVE : THREE_TO_THREE;
        }
    }

    /**
     * Given an ChimericAlignment representing two reference intervals rearranged as two intervals on the locally-assembled contig,
     * identify potential complications such as homology and duplication on the reference and/or on the contig.
     */
    @VisibleForTesting
    static LocationComplication resolveComplications(final ChimericAlignment chimericAlignment)
            throws IOException {

        final SimpleInterval leftReferenceInterval  = chimericAlignment.getCoordSortedReferenceIntervals()._1,
                             rightReferenceInterval = chimericAlignment.getCoordSortedReferenceIntervals()._2;

        final AlignmentRegion firstContigRegion  = chimericAlignment.regionWithLowerCoordOnContig,
                              secondContigRegion = chimericAlignment.regionWithHigherCoordOnContig;

        final byte[] contigSeq = chimericAlignment.contigSeq;

        final boolean involvesSameChromosomeReferenceIntervalSwitch = ChimericAlignment.involvesRefPositionSwitch(firstContigRegion, secondContigRegion);

        // a segment with lower coordinate on the locally-assembled contig could map to a higher reference coordinate region
        // under two basic types of SV's: inversion (strand switch necessary) and translocation (no strand switch necessary)
        if (chimericAlignment.strandSwitch!= ChimericAlignment.StrandSwitch.NO_SWITCH) { // the case involves an inversion
            // TODO: 12/5/16 duplication detection to be done for inversion alleles
            return new LocationComplication(getHomology(firstContigRegion, secondContigRegion, contigSeq), getInsertedSequence(firstContigRegion, secondContigRegion, contigSeq), "", 0, 0);
        } else if (ChimericAlignment.isNotSimpleTranslocation(chimericAlignment.regionWithLowerCoordOnContig, chimericAlignment.regionWithHigherCoordOnContig, chimericAlignment.strandSwitch, involvesSameChromosomeReferenceIntervalSwitch)) { // no translocation
            return getLocationComplicationForInsDel(chimericAlignment, leftReferenceInterval, rightReferenceInterval, firstContigRegion, secondContigRegion, contigSeq);
        } else { // TODO: 12/5/16 simple translocation, don't tackle yet
            return new LocationComplication("", "", "", 0, 0);
        }
    }

    // TODO: 12/12/16 the most complicated tandem duplication case is not exactly reproducible (e.g. '+' and '-' strand may give slightly different results by this treatment)
    private static LocationComplication getLocationComplicationForInsDel(final ChimericAlignment chimericAlignment,
                                                                         final SimpleInterval leftReferenceInterval, final SimpleInterval rightReferenceInterval,
                                                                         final AlignmentRegion firstContigRegion, final AlignmentRegion secondContigRegion,
                                                                         final byte[] contigSeq) {
        final int r1e = leftReferenceInterval.getEnd(),
                  r2b = rightReferenceInterval.getStart(),
                  c1e = chimericAlignment.regionWithLowerCoordOnContig.endInAssembledContig,
                  c2b = chimericAlignment.regionWithHigherCoordOnContig.startInAssembledContig;

        final int distBetweenAlignRegionsOnRef = r2b - r1e - 1, // distance-1 between the two regions on reference, denoted as d1 in the comments below
                  distBetweenAlignRegionsOnCtg = c2b - c1e - 1; // distance-1 between the two regions on contig, denoted as d2 in the comments below

        String homologyForwardStrandRepresentation="", insertedSeqForwardStrandRepresentation="", dupSeqForwardStrandRepresentation="";
        int dupSeqRepeatNumOnRef=0, dupSeqRepeatNumOnCtg=0;
        if (distBetweenAlignRegionsOnRef>0 && distBetweenAlignRegionsOnCtg==0) {        // DEL: simple deletion, deleted sequence is [r1e+1, r2b-1] on the reference
            homologyForwardStrandRepresentation    = "";
            insertedSeqForwardStrandRepresentation = "";
            dupSeqForwardStrandRepresentation      = "";
            dupSeqRepeatNumOnRef = dupSeqRepeatNumOnCtg = 0;
        } else if (distBetweenAlignRegionsOnRef==0 && distBetweenAlignRegionsOnCtg>0) { // INS: simple insertion, inserted sequence is the sequence [c1e+1, c2b-1] on the contig
            homologyForwardStrandRepresentation    = "";
            insertedSeqForwardStrandRepresentation = getInsertedSequence(firstContigRegion, secondContigRegion, contigSeq);
            dupSeqForwardStrandRepresentation      = "";
            dupSeqRepeatNumOnRef = dupSeqRepeatNumOnCtg = 0;
        } else if (distBetweenAlignRegionsOnRef==0 && distBetweenAlignRegionsOnCtg<0) { // DEL: tandem repeat contraction where reference has two copies but one copy was deleted on the contig; duplicated sequence on reference are [r1e-|d2|+1, r1e] and [r2b, r2b+|d2|-1]
            homologyForwardStrandRepresentation    = getHomology(firstContigRegion, secondContigRegion, contigSeq);
            insertedSeqForwardStrandRepresentation = "";
            dupSeqForwardStrandRepresentation      = homologyForwardStrandRepresentation;
            dupSeqRepeatNumOnRef      = 2;
            dupSeqRepeatNumOnCtg      = 1;
        } else if (distBetweenAlignRegionsOnRef<0 && distBetweenAlignRegionsOnCtg==0) { // INS: tandem repeat expansion of reference bases [r1e-|d1|+1, r1e] to contig bases [c1e-|d1|+1, c1e] and [c2b, c2b+|d1|-1]
            final byte[] insertedSeq = Arrays.copyOfRange(contigSeq, c1e-Math.abs(distBetweenAlignRegionsOnRef), c1e);
            if (!chimericAlignment.isForwardStrandRepresentation) SequenceUtil.reverseComplement(insertedSeq);
            homologyForwardStrandRepresentation    = "";
            insertedSeqForwardStrandRepresentation = "";
            dupSeqForwardStrandRepresentation      = new String(insertedSeq);
            dupSeqRepeatNumOnRef      = 1;
            dupSeqRepeatNumOnCtg      = 2;
        } else if (distBetweenAlignRegionsOnRef<0 && distBetweenAlignRegionsOnCtg>0) {  // INS: reference bases [r1e-|d1|+1, r1e] is similar to contig bases [c1e-|d1|+1, c1e] and [c2b, c2b+|d1|-1] and an inserted sequence [c1e+1, c2b-1] in between the two intervals on contig
            final byte[] insertedSeq = Arrays.copyOfRange(contigSeq, c1e-Math.abs(distBetweenAlignRegionsOnRef), c1e);
            if (!chimericAlignment.isForwardStrandRepresentation) SequenceUtil.reverseComplement(insertedSeq);
            homologyForwardStrandRepresentation    = "";
            insertedSeqForwardStrandRepresentation = getInsertedSequence(firstContigRegion, secondContigRegion, contigSeq); // note this does not incorporate the duplicated reference sequence
            dupSeqForwardStrandRepresentation      = new String(insertedSeq);
            dupSeqRepeatNumOnRef      = 1;
            dupSeqRepeatNumOnCtg      = 2;
        } else if (distBetweenAlignRegionsOnRef>0 && distBetweenAlignRegionsOnCtg>0) {  // DEL: deletion with scar, i.e. large non-conserved substitution, reference bases [r1e+1, r2b-1] is substituted with contig bases [c1e+1, c2b-1]
            homologyForwardStrandRepresentation    = "";
            insertedSeqForwardStrandRepresentation = getInsertedSequence(firstContigRegion, secondContigRegion, contigSeq);
            dupSeqForwardStrandRepresentation      = "";
            dupSeqRepeatNumOnRef = dupSeqRepeatNumOnCtg = 0;
        } else if (distBetweenAlignRegionsOnRef>0 && distBetweenAlignRegionsOnCtg<0) {  // DEL: a sequence of bases of length d1+HOM is deleted, and there's homology (which could be dup, but cannot tell): leftFlank+HOM+[r1e+1, r2b-1]+HOM+rightFlank -> leftFlank+HOM+rightFlank
            homologyForwardStrandRepresentation    = getHomology(firstContigRegion, secondContigRegion, contigSeq);
            insertedSeqForwardStrandRepresentation = "";
            dupSeqForwardStrandRepresentation      = "";
            dupSeqRepeatNumOnRef = dupSeqRepeatNumOnCtg = 0;
        } else if (distBetweenAlignRegionsOnRef<0 && distBetweenAlignRegionsOnCtg<0) {  // most complicated case, see below
            // DEL: duplication with repeat number N1 on reference, N2 on contig, such that N1 <= 2*N2 (and N2<N1);
            // INS: duplication with repeat number N1 on reference, N2 on contig, such that N2 <= 2*N1 (and N1<N2);
            // in both cases, the equal sign on the right can be taken only when there's pseudo-homology between starting bases of the duplicated sequence and starting bases of the right flanking region
            final int overlapOnLowerCNSequence, overlapOnHigherCNSequence;
            // the reference system with a shorter overlap (i.e. with less-negative distance between regions) has a higher repeat number
            final boolean isExpansion = distBetweenAlignRegionsOnRef<distBetweenAlignRegionsOnCtg;
            if (isExpansion) {
                overlapOnLowerCNSequence = Math.abs(distBetweenAlignRegionsOnRef);
                overlapOnHigherCNSequence = Math.abs(distBetweenAlignRegionsOnCtg);
            } else {     // d1 is lower absolute value -> reference has higher copy number of the duplication, i.e. DEL
                overlapOnLowerCNSequence = Math.abs(distBetweenAlignRegionsOnCtg);
                overlapOnHigherCNSequence = Math.abs(distBetweenAlignRegionsOnRef);
            }

            final TandemRepeatComplication duplicationComplication = tandemDuplicationComplication(overlapOnHigherCNSequence, overlapOnLowerCNSequence);

            final byte[] dupSeq;
            if (chimericAlignment.isForwardStrandRepresentation){
                dupSeq = Arrays.copyOfRange(contigSeq, c2b-1, c2b+duplicationComplication.repeatedSeqLen-1);
            } else {
                dupSeq = Arrays.copyOfRange(contigSeq, c2b+duplicationComplication.pseduHomologyLen-1, c2b+duplicationComplication.pseduHomologyLen+duplicationComplication.repeatedSeqLen-1);
                SequenceUtil.reverseComplement(dupSeq);
            }

            homologyForwardStrandRepresentation    = getHomology(firstContigRegion, secondContigRegion, contigSeq);
            insertedSeqForwardStrandRepresentation = ""; //isExpansion ? StringUtils.repeat(dupSeqForwardStrandRepresentation, duplicationComplication[1]-duplicationComplication[0]) : "";
            dupSeqForwardStrandRepresentation      = new String(dupSeq);
            dupSeqRepeatNumOnRef      = isExpansion ? duplicationComplication.lowerRepeatNumberEstimate  : duplicationComplication.higherRepeatNumberEstimate; //duplicationComplication[isExpansion?0:1];
            dupSeqRepeatNumOnCtg      = isExpansion ? duplicationComplication.higherRepeatNumberEstimate : duplicationComplication.lowerRepeatNumberEstimate; //duplicationComplication[isExpansion?1:0];
        } else if (distBetweenAlignRegionsOnRef==0 && distBetweenAlignRegionsOnCtg==0) {// SNP & indel
            throw new GATKException("Detected badly parsed chimeric alignment for identifying SV breakpoints; no rearrangement found: " + chimericAlignment.toString());
        }

        if (insertedSeqForwardStrandRepresentation.isEmpty()){
            if (dupSeqRepeatNumOnCtg != dupSeqRepeatNumOnRef && dupSeqForwardStrandRepresentation.isEmpty())
                throw new GATKException("An identified breakpoint pair seem to suggest insertion but the inserted sequence is empty: " + chimericAlignment.toString());
        }

        return new LocationComplication(homologyForwardStrandRepresentation, insertedSeqForwardStrandRepresentation, dupSeqForwardStrandRepresentation, dupSeqRepeatNumOnRef, dupSeqRepeatNumOnCtg);
    }

    /**
     * @return Micro-homology sequence using two alignments of the same contig: as indicated by their overlap on the contig itself.
     *          Empty if they don't overlap on the contig.
     */
    @VisibleForTesting
    static String getHomology(final AlignmentRegion current, final AlignmentRegion next, final byte[] contigSequence) {

        if (current.endInAssembledContig >= next.startInAssembledContig) {
            final byte[] homologyBytes = Arrays.copyOfRange(contigSequence, next.startInAssembledContig-1, current.endInAssembledContig);
            if (current.referenceInterval.getStart() > next.referenceInterval.getStart()) {
                SequenceUtil.reverseComplement(homologyBytes, 0, homologyBytes.length);
            }
            return new String(homologyBytes);
        } else {
            return "";
        }
    }

    /**
     * Note: not suitable for the most complicated case dealt with in {@link #resolveComplications(ChimericAlignment)}
     * @return Inserted sequence using two alignments of the same contig: as indicated by their separation on the the contig itself.
     */
    @VisibleForTesting
    static String getInsertedSequence(final AlignmentRegion current, final AlignmentRegion next, final byte[] contigSequence) {

        if (current.endInAssembledContig < next.startInAssembledContig-1) {
            final byte[] insertedSequenceBytes = Arrays.copyOfRange(contigSequence, current.endInAssembledContig, next.startInAssembledContig-1);
            if (current.referenceInterval.getStart() > next.referenceInterval.getStart()) {
                SequenceUtil.reverseComplement(insertedSequenceBytes, 0, insertedSequenceBytes.length);
            }
            return new String(insertedSequenceBytes);
        } else {
            return "";
        }
    }

    // TODO: 12/11/16 More robust way and ways to test (currently covered, but NOT DIRECTLY)
    /**
     * This function, when given overlaps of two corresponding regions on reference and contig sequences,
     * attempts to find--naively and slowly--the repeat numbers on the reference and on the contig of tandem repeats,
     * as well as the pseudo-homology between the duplicated sequence and the right flanking region.
     *
     * An example might help:
     * an assembled contig that's actually a repeat expansion from 1 repeat to 2 repeats with pseudo-homology:
     * TGCCAGGTTACATGGCAAAGAGGGTAGATATGGGGAGCTGTGAAGAATGGAGCCAGTAATTAAATTCACTGAAGTCTCCACAGGAGGGCAAGGTGGACAATCTGTCCCATAGGAGGGGGATTCATGAGGGGAGCTGTGAAGAATGGAGCCAGTAATTAAATTCACTGAAGTCTCCACAGGAGGGCAAGGTGGACAATCTGTCCCATAGGAGGGGGATTCAGGAGGGCAGCTGTGGATGGTGCAAATGCCATTTATGCTCCTCTCCACCCATATCC
     * can be aligned to chr18:312XXX:
     * the 1st alignment chr18:312579-718 140M135S, which can be broken into the following part
     * 31:  TGCCAGGTTACATGGCAAAGAGGGTAGATAT
     * 109: GGGGAGCTGTGAAGAATGGAGCCAGTAATTAAATTCACTGAAGTCTCCACAGGAGGGCAAGGTGGACAATCTGTCCCATAGGAGGGGGATTCATGAGGGGAGCTGTGAA
     * 135: GAATGGAGCCAGTAATTAAATTCACTGAAGTCTCCACAGGAGGGCAAGGTGGACAATCTGTCCCATAGGAGGGGGATTCAGGAGGGCAGCTGTGGATGGTGCAAATGCCATTTATGCTCCTCTCCACCCATATCC
     * And the arithmetic works this way:
     * 31 + 109 = 140
     * 109 = 96 + 13
     * where 31 is the left flanking region before the repeated unit, which itself is 96 bases long (see below),
     * the number 13 is the length of the pseudo-homology between the starting bases of the repeated sequence and the right flanking region
     * a clearer picture emerges when we look at the 2nd alignment chr18:312610-757 127S148M, which can be broken into
     * 31: TGCCAGGTTACATGGCAAAGAGGGTAGATAT
     * 96: GGGGAGCTGTGAAGAATGGAGCCAGTAATTAAATTCACTGAAGTCTCCACAGGAGGGCAAGGTGGACAATCTGTCCCATAGGAGGGGGATTCATGA
     * 96: GGGGAGCTGTGAAGAATGGAGCCAGTAATTAAATTCACTGAAGTCTCCACAGGAGGGCAAGGTGGACAATCTGTCCCATAGGAGGGGGATTCAGGA
     * 13: GGGCAGCTGTGGA
     * 39: TGGTGCAAATGCCATTTATGCTCCTCTCCACCCATATCC
     * And the arithmetic works this way:
     * 31 + 96 = 127
     * 96 + 13 + 39 = 148
     */
    @VisibleForTesting
    static TandemRepeatComplication tandemDuplicationComplication(final int overlapOnHigherCNSequence, final int overlapOnLowerCNSequence) {

        int higherCnEst=0, lowerCnEst=0, unitLen=0, pseudoHomLen=0;
        double err = Double.MAX_VALUE;
        for(int cn2=1; cn2<MAX_LOWER_CN; ++cn2) {
            for(int cn1=cn2+1; cn1<=2*cn2; ++cn1) {
                final int dupLenUpperBound = cn1==2*cn2 ? overlapOnLowerCNSequence : overlapOnHigherCNSequence;
                for (int l=2; l<=dupLenUpperBound; ++l) {
                    for (int lambda=0; lambda<l; ++lambda) {
                        final int d1 = (2*cn2 - cn1)*l + lambda;
                        final int d2 = cn2*l + lambda;
                        final double newErr = Math.abs(overlapOnHigherCNSequence-d1) + Math.abs(overlapOnLowerCNSequence-d2);
                        if (newErr < err) {
                            err = newErr;
                            higherCnEst = cn1; lowerCnEst = cn2;
                            unitLen= l; pseudoHomLen = lambda;
                        }
                        if (err<1)
                            return new TandemRepeatComplication(lowerCnEst, higherCnEst, unitLen, pseudoHomLen);
                    }
                }
            }
        }

        return new TandemRepeatComplication(lowerCnEst, higherCnEst, unitLen, pseudoHomLen);
    }

    private static final class TandemRepeatComplication {
        final int lowerRepeatNumberEstimate;
        final int higherRepeatNumberEstimate;
        final int repeatedSeqLen;
        final int pseduHomologyLen;

        TandemRepeatComplication(final int lowerRepeatNumberEstimate,
                                 final int higherRepeatNumberEstimate,
                                 final int repeatedSeqLen,
                                 final int pseduHomologyLen){
            this.lowerRepeatNumberEstimate = lowerRepeatNumberEstimate;
            this.higherRepeatNumberEstimate = higherRepeatNumberEstimate;
            this.repeatedSeqLen = repeatedSeqLen;
            this.pseduHomologyLen = pseduHomologyLen;
        }
    }

    // TODO: 12/11/16 again this does not deal with inter-chromosome translocation yet
    /**
     * Returns the reference coordinates of the left and right breakpoints implied by this chimeric alignment.
     * If there is homologous sequence represented in the alignments, it will be assigned to the side of the breakpoint
     * with higher reference coordinates.
     */
    @VisibleForTesting
    static Tuple2<SimpleInterval, SimpleInterval> leftJustifyBreakpoints(final ChimericAlignment ca, final LocationComplication complication) {

        final int homologyLen = complication.homologyForwardStrandRep.length();

        final String leftBreakpointRefContig, rightBreakpointRefContig;
        final int leftBreakpointCoord, rightBreakpointCoord;
        if (complication.dupSeqRepeatNumOnRef !=0 && complication.dupSeqRepeatNumOnCtg !=0) { // todo : development artifact-- assuming tandem duplication is not co-existing with inversion
            final int dupLen = complication.dupSeqForwardStrandRep.length();
            leftBreakpointRefContig = rightBreakpointRefContig = ca.regionWithLowerCoordOnContig.referenceInterval.getContig();
            final SimpleInterval leftReferenceInterval, rightReferenceInterval;
            if (ca.isForwardStrandRepresentation) {
                leftReferenceInterval = ca.regionWithLowerCoordOnContig.referenceInterval;
                rightReferenceInterval = ca.regionWithHigherCoordOnContig.referenceInterval;
            } else {
                leftReferenceInterval = ca.regionWithHigherCoordOnContig.referenceInterval;
                rightReferenceInterval = ca.regionWithLowerCoordOnContig.referenceInterval;
            }
            if (complication.dupSeqRepeatNumOnCtg > complication.dupSeqRepeatNumOnRef) {
                leftBreakpointCoord = leftReferenceInterval.getEnd() - homologyLen - (complication.dupSeqRepeatNumOnCtg - complication.dupSeqRepeatNumOnRef)*dupLen;
            } else {
                leftBreakpointCoord = leftReferenceInterval.getEnd() - homologyLen;
            }
            rightBreakpointCoord = rightReferenceInterval.getStart() - 1;
        } else { // inversion and simple deletion & insertion
            final SimpleInterval leftReferenceInterval, rightReferenceInterval;
            if (ca.isForwardStrandRepresentation) {
                leftReferenceInterval  = ca.regionWithLowerCoordOnContig.referenceInterval;
                rightReferenceInterval = ca.regionWithHigherCoordOnContig.referenceInterval;
            } else {
                leftReferenceInterval  = ca.regionWithHigherCoordOnContig.referenceInterval;
                rightReferenceInterval = ca.regionWithLowerCoordOnContig.referenceInterval;
            }
            leftBreakpointRefContig  = leftReferenceInterval.getContig();
            rightBreakpointRefContig = rightReferenceInterval.getContig();
            if (ca.strandSwitch == ChimericAlignment.StrandSwitch.NO_SWITCH) {
                leftBreakpointCoord  = leftReferenceInterval.getEnd() - homologyLen;
                rightBreakpointCoord = rightReferenceInterval.getStart() - 1;
            } else if (ca.strandSwitch == ChimericAlignment.StrandSwitch.FORWARD_TO_REVERSE){
                leftBreakpointCoord  = leftReferenceInterval.getEnd() - homologyLen;
                rightBreakpointCoord = rightReferenceInterval.getEnd();
            } else {
                leftBreakpointCoord  = leftReferenceInterval.getStart() - 1;
                rightBreakpointCoord = rightReferenceInterval.getStart() + homologyLen - 1;
            }
        }
        final SimpleInterval leftBreakpoint = new SimpleInterval(leftBreakpointRefContig, leftBreakpointCoord, leftBreakpointCoord);
        final SimpleInterval rightBreakpoint = new SimpleInterval(rightBreakpointRefContig, rightBreakpointCoord, rightBreakpointCoord);
        return new Tuple2<>(leftBreakpoint, rightBreakpoint);
    }

    protected NovelAdjacencyReferenceLocations(final Kryo kryo, final Input input) {
        final String contig1 = input.readString();
        final int start1 = input.readInt();
        final int end1 = input.readInt();
        this.leftJustifiedLeftRefLoc = new SimpleInterval(contig1, start1, end1);
        final String contig2 = input.readString();
        final int start2 = input.readInt();
        final int end2 = input.readInt();
        this.leftJustifiedRightRefLoc = new SimpleInterval(contig2, start2, end2);

        this.endConnectionType = EndConnectionType.values()[input.readInt()];
        this.complication = kryo.readObject(input, LocationComplication.class);
    }

    @Override
    public String toString() {
        return String.format("%s\t%s\t%s\t%s", leftJustifiedLeftRefLoc.toString(), leftJustifiedRightRefLoc.toString(),
                                               endConnectionType.name(), complication.toString());
    }

    @VisibleForTesting
    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        final NovelAdjacencyReferenceLocations that = (NovelAdjacencyReferenceLocations) o;

        if (leftJustifiedLeftRefLoc != null ? !leftJustifiedLeftRefLoc.equals(that.leftJustifiedLeftRefLoc) : that.leftJustifiedLeftRefLoc != null)
            return false;
        if (leftJustifiedRightRefLoc != null ? !leftJustifiedRightRefLoc.equals(that.leftJustifiedRightRefLoc) : that.leftJustifiedRightRefLoc != null)
            return false;

        return endConnectionType.equals(that.endConnectionType) && complication.equals(that.complication);
    }

    @VisibleForTesting
    @Override
    public int hashCode() {
        return Objects.hash(leftJustifiedLeftRefLoc, leftJustifiedRightRefLoc, complication, 2659* endConnectionType.ordinal());
    }

    protected void serialize(final Kryo kryo, final Output output) {
        output.writeString(leftJustifiedLeftRefLoc.getContig());
        output.writeInt(leftJustifiedLeftRefLoc.getStart());
        output.writeInt(leftJustifiedLeftRefLoc.getEnd());
        output.writeString(leftJustifiedRightRefLoc.getContig());
        output.writeInt(leftJustifiedRightRefLoc.getStart());
        output.writeInt(leftJustifiedRightRefLoc.getEnd());
        output.writeInt(endConnectionType.ordinal());
        kryo.writeObject(output, complication);
    }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<NovelAdjacencyReferenceLocations> {
        @Override
        public void write(final Kryo kryo, final Output output, final NovelAdjacencyReferenceLocations novelAdjacencyReferenceLocations) {
            novelAdjacencyReferenceLocations.serialize(kryo, output);
        }

        @Override
        public NovelAdjacencyReferenceLocations read(final Kryo kryo, final Input input, final Class<NovelAdjacencyReferenceLocations> klass ) {
            return new NovelAdjacencyReferenceLocations(kryo, input);
        }
    }
}
