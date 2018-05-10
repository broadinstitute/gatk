package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.*;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.SimpleChimera.DistancesBetweenAlignmentsOnRefAndOnRead;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.Strand;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.CigarUtils;

import java.util.*;
import java.util.stream.Collectors;

/**
 * A helper struct for annotating complications that make the locations represented by its associated
 * {@link NovelAdjacencyAndAltHaplotype} a little ambiguous.
 *
 * Inserted sequence contains portions of the contig that are aligned to neither region,
 * and therefore may be inserted in the sample.
 * For example, a translocation breakpoint with a micro-insertion:
 *
 * Contig:
 * ACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG
 * Alignment regions:
 * |-----1:100-200-------|
 *                          |----2:100-200-----|
 * Inserted sequence:
 *  GA
 *
 * Homology represents ambiguity about the exact location of the breakpoint. For example, in this case one alignment
 * region ends with "AC" and the next begins with AC, so we don't know if the AC truly belongs with the first or
 * second alignment region.
 *
 * Contig:
 * ACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG
 * Alignment regions:
 * |-----1:100-200-------|
 *                    |-----2:100-200----------|
 * Homology:
 *  AC
 *
 * We currently handle types for simple chimera-induced precise variants {@link TypeInferredFromSimpleChimera}:
 * <ul>
 *     <li>
 *         deletions
 *     </li>
 *     <li>
 *         small insertions
 *     </li>
 *     <Li>
 *         replacements (i.e. del and ins at the same location)
 *     </Li>
 *     <li>
 *         BND type, which in turn is turned into several subtypes
 *         <ul>
 *             <li>
 *                 intra-chromosome strand-switch BND's, i.e. inversion breakpoint suspects
 *             </li>
 *             <li>
 *                 intra-chromosome no strand-switch but order swap BND's, i.e. large tandem duplication breakpoint suspects
 *             </li>
 *             <li>
 *                 inter-chromosome BND's, with or without strand switch.
 *             </li>
 *         </ul>
 *     </li>
 * </ul>
 */
public abstract class BreakpointComplications {

    /**
     * '+' strand representations of micro-homology, inserted sequence and duplicated sequence on the reference.
     *
     * If there is homologous sequence represented in the {@link SimpleChimera},
     * it will be assigned to the side of the breakpoint with higher reference coordinates
     * (as judged by a {@link SAMSequenceDictionary}), i.e. we follow left alignment convention.
     */
    protected String homologyForwardStrandRep = "";
    protected String insertedSequenceForwardStrandRep = "";

    public String getHomologyForwardStrandRep() {
        return homologyForwardStrandRep;
    }
    public String getInsertedSequenceForwardStrandRep() {
        return insertedSequenceForwardStrandRep;
    }

    /**
     * @return Intended for use in debugging and exception message only.
     */
    @Override
    public String toString() {
        return "homology: " + homologyForwardStrandRep + "\tinserted sequence: " + insertedSequenceForwardStrandRep;
    }

    /**
     * Intended to be overridden by sub classes when more complications are involved.
     */
    public Map<String, Object> toVariantAttributes() {

        final Map<String, Object> attributeMap = new HashMap<>();

        if ( !getInsertedSequenceForwardStrandRep().isEmpty() ) {
            attributeMap.put(GATKSVVCFConstants.INSERTED_SEQUENCE, getInsertedSequenceForwardStrandRep());
            attributeMap.put(GATKSVVCFConstants.INSERTED_SEQUENCE_LENGTH, getInsertedSequenceForwardStrandRep().length());
        }

        if ( !getHomologyForwardStrandRep().isEmpty() ) {
            attributeMap.put(GATKSVVCFConstants.HOMOLOGY, getHomologyForwardStrandRep());
            attributeMap.put(GATKSVVCFConstants.HOMOLOGY_LENGTH, getHomologyForwardStrandRep().length());
        }

        return attributeMap;
    }

    /**
     * To be overridden as appropriate.
     */
    public boolean hasDuplicationAnnotation() {
        return false;
    }

    // =================================================================================================================
    protected BreakpointComplications() {

    }

    @VisibleForTesting
    BreakpointComplications(final String homologyForwardStrandRep, final String insertedSequenceForwardStrandRep) {
        this.homologyForwardStrandRep = homologyForwardStrandRep;
        this.insertedSequenceForwardStrandRep = insertedSequenceForwardStrandRep;
    }

    protected BreakpointComplications(final Kryo kryo, final Input input) {
        homologyForwardStrandRep = input.readString();
        insertedSequenceForwardStrandRep = input.readString();
    }
    protected void serialize(final Kryo kryo, final Output output) {
        output.writeString(homologyForwardStrandRep);
        output.writeString(insertedSequenceForwardStrandRep);
    }
    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        BreakpointComplications that = (BreakpointComplications) o;

        if (!homologyForwardStrandRep.equals(that.homologyForwardStrandRep)) return false;
        return insertedSequenceForwardStrandRep.equals(that.insertedSequenceForwardStrandRep);
    }
    @Override
    public int hashCode() {
        int result = homologyForwardStrandRep.hashCode();
        result = 31 * result + insertedSequenceForwardStrandRep.hashCode();
        return result;
    }

    //==================================================================================================================

    /**
     * @return Micro-homology sequence using two alignments of the same contig: as indicated by their overlap on the contig itself.
     *          Empty if they don't overlap on the contig.
     */
    @VisibleForTesting
    static String inferHomology(final AlignmentInterval first, final AlignmentInterval second, final byte[] contigSequence,
                                final boolean firstAfterSecond) {

        if (first.endInAssembledContig >= second.startInAssembledContig) {
            final byte[] homologyBytes = Arrays.copyOfRange(contigSequence,
                    second.startInAssembledContig-1, first.endInAssembledContig);
            return new String(reverseComplementIfNecessary(homologyBytes, first, second, firstAfterSecond));
        } else {
            return "";
        }
    }

    /**
     * Note: not suitable for the most complicated case dealt with in {@link BreakpointComplications( SimpleChimera )}
     * @return Inserted sequence using two alignments of the same contig: as indicated by their separation on the the contig itself.
     */
    @VisibleForTesting
    static String inferInsertedSequence(final AlignmentInterval first, final AlignmentInterval second, final byte[] contigSequence,
                                        final boolean firstAfterSecond) {

        if (first.endInAssembledContig < second.startInAssembledContig - 1) {
            final byte[] insertedSequenceBytes = Arrays.copyOfRange(contigSequence,
                    first.endInAssembledContig, second.startInAssembledContig - 1);
            return new String(reverseComplementIfNecessary(insertedSequenceBytes, first, second, firstAfterSecond));
        } else {
            return "";
        }
    }

    private static byte[] reverseComplementIfNecessary(final byte[] seq,
                                                       final AlignmentInterval first, final AlignmentInterval second,
                                                       final boolean firstAfterSecond) {
        if ( first.forwardStrand == second.forwardStrand ) { // reference order switch
            if (!first.forwardStrand) {
                SequenceUtil.reverseComplement(seq, 0, seq.length);
            }
        } else {
            if (firstAfterSecond == first.forwardStrand) {
                SequenceUtil.reverseComplement(seq, 0, seq.length);
            }
        }
        return seq;
    }

    //==================================================================================================================

    /**
     * For simple deletion, insertion, and replacement (dep and ins at the same time).
     *
     * NOTE THAT this doesn't work for small duplications (expansion or contraction);
     * for those, check {@link SmallDuplicationBreakpointComplications}.
     */
    @DefaultSerializer(SimpleInsDelOrReplacementBreakpointComplications.Serializer.class)
    public static final class SimpleInsDelOrReplacementBreakpointComplications extends BreakpointComplications {

        @VisibleForTesting
        public SimpleInsDelOrReplacementBreakpointComplications(final String homologyForwardStrandRep, final String insertedSequenceForwardStrandRep) {
            super(homologyForwardStrandRep, insertedSequenceForwardStrandRep);
        }

        SimpleInsDelOrReplacementBreakpointComplications(final SimpleChimera simpleChimera, final byte[] contigSeq, final boolean firstAfterSecond) {

            final DistancesBetweenAlignmentsOnRefAndOnRead distances = simpleChimera.getDistancesBetweenAlignmentsOnRefAndOnRead();
            final int distBetweenAlignRegionsOnRef = distances.gapBetweenAlignRegionsOnRef, // distance-1 between the two regions on reference, denoted as d1 in the comments below
                      distBetweenAlignRegionsOnCtg = distances.gapBetweenAlignRegionsOnCtg; // distance-1 between the two regions on contig, denoted as d2 in the comments below

            if ( distBetweenAlignRegionsOnRef > 0 ) { // some bases deleted: here it could be a simple deletion or replacement
                if (distBetweenAlignRegionsOnCtg>=0) {
                    // either: a clean deletion, deleted sequence is [r1e+1, r2b-1] on the reference
                    // or    : deletion with scar, i.e. large non-conserved substitution, reference bases [r1e+1, r2b-1] is substituted with contig bases [c1e+1, c2b-1]
                    insertedSequenceForwardStrandRep =
                            inferInsertedSequence(simpleChimera.regionWithLowerCoordOnContig,
                                    simpleChimera.regionWithHigherCoordOnContig, contigSeq, firstAfterSecond);
                } else {
                    // a sequence of bases of length d1+HOM is deleted, and there's homology (which could be dup, but cannot tell): leftFlank+HOM+[r1e+1, r2b-1]+HOM+rightFlank -> leftFlank+HOM+rightFlank
                    homologyForwardStrandRep =
                            inferHomology(simpleChimera.regionWithLowerCoordOnContig,
                                    simpleChimera.regionWithHigherCoordOnContig, contigSeq, firstAfterSecond);
                }
            } else if (distBetweenAlignRegionsOnRef == 0 && distBetweenAlignRegionsOnCtg > 0) { // Insertion: simple insertion, inserted sequence is the sequence [c1e+1, c2b-1] on the contig
                insertedSequenceForwardStrandRep =
                        inferInsertedSequence(simpleChimera.regionWithLowerCoordOnContig,
                                simpleChimera.regionWithHigherCoordOnContig, contigSeq, firstAfterSecond);
            } else {
                throw new GATKException.ShouldNeverReachHereException(
                        "Inferring breakpoint complications with the wrong unit: " +
                        "using simple ins-del unit for simple chimera:\n" + simpleChimera.toString());
            }

            if ( distBetweenAlignRegionsOnCtg > 0 && insertedSequenceForwardStrandRep.isEmpty() ){
                    throw new GATKException.ShouldNeverReachHereException(
                            "An identified breakpoint pair seem to suggest insertion " +
                            "but the inserted sequence is empty: " + simpleChimera.toString());
            }
        }

        private SimpleInsDelOrReplacementBreakpointComplications(final Kryo kryo, final Input input) {
            super(kryo, input);
        }
        protected void serialize(final Kryo kryo, final Output output) {
            super.serialize(kryo, output);
        }
        public static final class Serializer extends com.esotericsoftware.kryo.Serializer<SimpleInsDelOrReplacementBreakpointComplications> {
            @Override
            public void write(final Kryo kryo, final Output output, final SimpleInsDelOrReplacementBreakpointComplications breakpointComplications) {
                breakpointComplications.serialize(kryo, output);
            }

            @Override
            public SimpleInsDelOrReplacementBreakpointComplications read(final Kryo kryo, final Input input, final Class<SimpleInsDelOrReplacementBreakpointComplications> klass ) {
                return new SimpleInsDelOrReplacementBreakpointComplications(kryo, input);
            }
        }
        @Override
        public int hashCode() {
            return super.hashCode();
        }
        @Override
        public boolean equals(Object obj) {
            return super.equals(obj);
        }
    }

    /**
     * For duplications small enough that we seemingly have assembled across the whole event.
     * NOTE THAT this does not handle case when one of the duplication copy is inverted;
     * for that see {@link IntraChrStrandSwitchBreakpointComplications} (possibly partially).
     */
    public abstract static class SmallDuplicationBreakpointComplications extends BreakpointComplications {

        protected SimpleInterval dupSeqRepeatUnitRefSpan = null;
        protected int dupSeqRepeatNumOnRef = 0;
        protected int dupSeqRepeatNumOnCtg = 0;
        protected List<Strand> dupSeqStrandOnRef = null;
        protected List<Strand> dupSeqStrandOnCtg = null;

        @VisibleForTesting
        SmallDuplicationBreakpointComplications(final String homologyForwardStrandRep, final String insertedSequenceForwardStrandRep,
                                                final SimpleInterval dupSeqRepeatUnitRefSpan,
                                                final int dupSeqRepeatNumOnRef, final int dupSeqRepeatNumOnCtg,
                                                final List<Strand> dupSeqStrandOnRef, final List<Strand> dupSeqStrandOnCtg) {
            super(homologyForwardStrandRep, insertedSequenceForwardStrandRep);
            this.dupSeqRepeatUnitRefSpan = dupSeqRepeatUnitRefSpan;
            this.dupSeqRepeatNumOnRef = dupSeqRepeatNumOnRef;
            this.dupSeqRepeatNumOnCtg = dupSeqRepeatNumOnCtg;
            this.dupSeqStrandOnRef = dupSeqStrandOnRef;
            this.dupSeqStrandOnCtg = dupSeqStrandOnCtg;
        }

        SmallDuplicationBreakpointComplications() {

        }

        SmallDuplicationBreakpointComplications(final Kryo kryo, final Input input) {
            super(kryo, input);
            final String ctg = input.readString();
            final int start = input.readInt();
            final int end = input.readInt();
            dupSeqRepeatUnitRefSpan = new SimpleInterval(ctg, start, end);
            dupSeqRepeatNumOnRef = input.readInt();
            dupSeqRepeatNumOnCtg = input.readInt();
            dupSeqStrandOnRef = new ArrayList<>(dupSeqRepeatNumOnRef);
            for (int i=0; i<dupSeqRepeatNumOnRef; ++i) {
                dupSeqStrandOnRef.add(Strand.values()[input.readInt()]);
            }
            dupSeqStrandOnCtg = new ArrayList<>(dupSeqRepeatNumOnCtg);
            for (int i=0; i<dupSeqRepeatNumOnCtg; ++i) {
                dupSeqStrandOnCtg.add(Strand.values()[input.readInt()]);
            }
        }

        public final SimpleInterval getDupSeqRepeatUnitRefSpan() {
            return dupSeqRepeatUnitRefSpan;
        }
        public final int getDupSeqRepeatNumOnRef() {
            return dupSeqRepeatNumOnRef;
        }
        public final int getDupSeqRepeatNumOnCtg() {
            return dupSeqRepeatNumOnCtg;
        }
        public final List<Strand> getDupSeqOrientationsOnCtg() {
            return dupSeqStrandOnCtg;
        }

        @Override
        public Map<String, Object> toVariantAttributes() {
            final Map<String, Object> parentAttributesToBeFilled = super.toVariantAttributes();

            parentAttributesToBeFilled.put(GATKSVVCFConstants.DUP_REPEAT_UNIT_REF_SPAN,
                    getDupSeqRepeatUnitRefSpan().toString());
            parentAttributesToBeFilled.put(GATKSVVCFConstants.DUPLICATION_NUMBERS,
                    getDupSeqRepeatNumOnRef() + VCFConstants.INFO_FIELD_ARRAY_SEPARATOR + getDupSeqRepeatNumOnCtg());
            parentAttributesToBeFilled.put(GATKSVVCFConstants.DUP_ORIENTATIONS,
                    getDupSeqOrientationsOnCtg().stream().map(Strand::toString).collect(Collectors.joining()));
            parentAttributesToBeFilled.put(dupSeqRepeatNumOnRef < dupSeqRepeatNumOnCtg ? GATKSVVCFConstants.DUP_TAN_EXPANSION_STRING : GATKSVVCFConstants.DUP_TAN_CONTRACTION_STRING,
                    "");

            return parentAttributesToBeFilled;
        }

        @Override
        public String toString() {
            String toPrint = super.toString();
            toPrint += String.format("\ttandem duplication repeat unit ref span: %s\t"+
                            "ref repeat num: %d\t"+
                            "ctg repeat num: %d\t"+
                            "dupSeqStrandOnRef: %s\t" +
                            "dupSeqStrandOnCtg: %s\t",
                    dupSeqRepeatUnitRefSpan,
                    dupSeqRepeatNumOnRef,
                    dupSeqRepeatNumOnCtg,
                    dupSeqStrandOnRef.stream().map(Strand::toString).collect(SVUtils.arrayListCollector(dupSeqStrandOnRef.size())).toString(),
                    dupSeqStrandOnCtg.stream().map(Strand::toString).collect(SVUtils.arrayListCollector(dupSeqStrandOnCtg.size())).toString());
            return toPrint;
        }

        @Override
        public final boolean hasDuplicationAnnotation() {
            return true;
        }

        public final boolean isDupContraction() {
            return dupSeqRepeatNumOnRef > dupSeqRepeatNumOnCtg;
        }

        protected void serialize(final Kryo kryo, final Output output) {
            super.serialize(kryo, output);
            output.writeString(dupSeqRepeatUnitRefSpan.getContig());
            output.writeInt(dupSeqRepeatUnitRefSpan.getStart());
            output.writeInt(dupSeqRepeatUnitRefSpan.getEnd());
            output.writeInt(dupSeqRepeatNumOnRef);
            output.writeInt(dupSeqRepeatNumOnCtg);
            dupSeqStrandOnRef.forEach(s -> output.writeInt(s.ordinal()));
            dupSeqStrandOnCtg.forEach(s -> output.writeInt(s.ordinal()));
        }

        @Override
        public boolean equals(final Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            if (!super.equals(o)) return false;

            final SmallDuplicationBreakpointComplications that = (SmallDuplicationBreakpointComplications) o;

            if (dupSeqRepeatNumOnRef != that.dupSeqRepeatNumOnRef) return false;
            if (dupSeqRepeatNumOnCtg != that.dupSeqRepeatNumOnCtg) return false;
            if (!dupSeqRepeatUnitRefSpan.equals(that.dupSeqRepeatUnitRefSpan)) return false;
            if (!dupSeqStrandOnRef.equals(that.dupSeqStrandOnRef)) return false;
            return dupSeqStrandOnCtg.equals(that.dupSeqStrandOnCtg);
        }

        @Override
        public int hashCode() {
            int result = super.hashCode();
            result = 31 * result + dupSeqRepeatUnitRefSpan.hashCode();
            result = 31 * result + dupSeqRepeatNumOnRef;
            result = 31 * result + dupSeqRepeatNumOnCtg;
            for (final Strand strand : dupSeqStrandOnRef) {
                result = 31 * result + strand.ordinal();
            }
            for (final Strand strand : dupSeqStrandOnCtg) {
                result = 31 * result + strand.ordinal();
            }
            return result;
        }
    }

    @DefaultSerializer(SmallDuplicationWithPreciseDupRangeBreakpointComplications.Serializer.class)
    public static final class SmallDuplicationWithPreciseDupRangeBreakpointComplications extends SmallDuplicationBreakpointComplications {
        public static final List<String> DEFAULT_CIGAR_STRINGS_FOR_DUP_SEQ_ON_CTG = Collections.emptyList();

        private List<String> cigarStringsForDupSeqOnCtg = null;
        final List<String> getCigarStringsForDupSeqOnCtgForwardStrandRep() {
            return cigarStringsForDupSeqOnCtg;
        }

        @VisibleForTesting
        SmallDuplicationWithPreciseDupRangeBreakpointComplications(final String homologyForwardStrandRep, final String insertedSequenceForwardStrandRep,
                                                                   final SimpleInterval dupSeqRepeatUnitRefSpan,
                                                                   final int dupSeqRepeatNumOnRef, final int dupSeqRepeatNumOnCtg,
                                                                   final List<Strand> dupSeqStrandOnRef, final List<Strand> dupSeqStrandOnCtg, final List<String> cigarStringsForDupSeqOnCtg) {
            super(homologyForwardStrandRep, insertedSequenceForwardStrandRep, dupSeqRepeatUnitRefSpan, dupSeqRepeatNumOnRef, dupSeqRepeatNumOnCtg, dupSeqStrandOnRef, dupSeqStrandOnCtg);
            this.cigarStringsForDupSeqOnCtg = cigarStringsForDupSeqOnCtg;
        }
        private SmallDuplicationWithPreciseDupRangeBreakpointComplications(final Kryo kryo, final Input input) {
            super(kryo, input);

            final int cigarCounts = input.readInt();
            cigarStringsForDupSeqOnCtg = new ArrayList<>(cigarCounts);
            for(int i = 0; i < cigarCounts; ++i) {
                cigarStringsForDupSeqOnCtg.add(input.readString());
            }
        }

        @Override
        public final String toString() {
            return super.toString() +
                    "\tprecise cigarStringsForDupSeqOnCtg:\t" + (cigarStringsForDupSeqOnCtg.isEmpty() ? "" : cigarStringsForDupSeqOnCtg);
        }

        @Override
        public Map<String, Object> toVariantAttributes() {
            final Map<String, Object> parentAttributesToBeFilled = super.toVariantAttributes();
            if ( !getCigarStringsForDupSeqOnCtgForwardStrandRep().isEmpty() ) {
                parentAttributesToBeFilled.put(GATKSVVCFConstants.DUP_SEQ_CIGARS,
                        StringUtils.join(getCigarStringsForDupSeqOnCtgForwardStrandRep(), VCFConstants.INFO_FIELD_ARRAY_SEPARATOR));
            }

            return parentAttributesToBeFilled;
        }

        SmallDuplicationWithPreciseDupRangeBreakpointComplications(final SimpleChimera simpleChimera, final byte[] contigSeq, final boolean firstAfterSecond) {
            final DistancesBetweenAlignmentsOnRefAndOnRead distances = simpleChimera.getDistancesBetweenAlignmentsOnRefAndOnRead();

            if (distances.gapBetweenAlignRegionsOnRef > 0)
                throw new GATKException.ShouldNeverReachHereException(
                        "Simple chimera being sent down the wrong path, where the signature indicates " +
                        "a simple deletion but complication being resolve for small duplication. \n" + simpleChimera.toString());

            if (distances.gapBetweenAlignRegionsOnRef < 0) {
                if (distances.gapBetweenAlignRegionsOnCtg >= 0) {
                    // Tandem repeat expansion:
                    // reference bases [r1e-|d1|+1, r1e] to contig bases [c1e-|d1|+1, c1e] and [c2b, c2b+|d1|-1]
                    // with optional inserted sequence [c1e+1, c2b-1] in between the two intervals on contig
                    resolveComplicationForSimpleTandupExpansion(simpleChimera, distances, contigSeq, firstAfterSecond);
                } else {
                    throw new GATKException.ShouldNeverReachHereException(
                            "Simple chimera being sent down the wrong path, where the signature indicates " +
                            "complex duplication but complication being resolved for simple small duplication. \n" + simpleChimera.toString());
                }
            } else { // here distances.gapBetweenAlignRegionsOnRef == 0
                if (distances.gapBetweenAlignRegionsOnCtg < 0) {
                    // Tandem repeat contraction: reference has two copies but one copy was deleted on the contig;
                    // duplicated sequence on reference are [r1e-|d2|+1, r1e] and [r2b, r2b+|d2|-1]
                    resolveComplicationForSimpleTandupContraction(simpleChimera, distances, contigSeq, firstAfterSecond);
                } else if (distances.gapBetweenAlignRegionsOnCtg > 0) {
                    throw new GATKException.ShouldNeverReachHereException(
                            "Simple chimera being sent down the wrong path, where the signature indicates " +
                            "an insertion but complication being resolve for small duplication. \n" +
                            simpleChimera.toString());
                } else { // gapBetweenAlignRegionsOnCtg == 0
                    throw new GATKException.ShouldNeverReachHereException(
                            "Detected badly parsed chimeric alignment for identifying SV breakpoints; " +
                            "no rearrangement found: " + simpleChimera.toString());
                }
            }

            if ( insertedSequenceForwardStrandRep.isEmpty() ){
                if ( dupSeqRepeatNumOnCtg != dupSeqRepeatNumOnRef && null == dupSeqRepeatUnitRefSpan )
                    throw new GATKException.ShouldNeverReachHereException(
                            "An identified breakpoint pair seem to suggest insertion " +
                            "but the inserted sequence is empty: " + simpleChimera.toString());
            }
        }

        private void resolveComplicationForSimpleTandupExpansion(final SimpleChimera simpleChimera,
                                                                 final DistancesBetweenAlignmentsOnRefAndOnRead distances,
                                                                 final byte[] contigSeq, final boolean firstAfterSecond) {

            final AlignmentInterval firstContigRegion = simpleChimera.regionWithLowerCoordOnContig;
            final AlignmentInterval secondContigRegion = simpleChimera.regionWithHigherCoordOnContig;
            final SimpleInterval leftReferenceInterval;
            if (simpleChimera.isForwardStrandRepresentation)
                leftReferenceInterval = simpleChimera.regionWithLowerCoordOnContig.referenceSpan;
            else
                leftReferenceInterval = simpleChimera.regionWithHigherCoordOnContig.referenceSpan;

            // note this does not incorporate the duplicated reference sequence
            insertedSequenceForwardStrandRep = distances.gapBetweenAlignRegionsOnCtg == 0 ? "" : inferInsertedSequence(firstContigRegion, secondContigRegion, contigSeq, firstAfterSecond);
            dupSeqRepeatUnitRefSpan   = new SimpleInterval(leftReferenceInterval.getContig(), distances.rightAlnRefStart, distances.leftAlnRefEnd);
            dupSeqRepeatNumOnRef      = 1;
            dupSeqRepeatNumOnCtg      = 2;
            dupSeqStrandOnRef         = Collections.singletonList(Strand.POSITIVE);
            dupSeqStrandOnCtg         = Arrays.asList(Strand.POSITIVE, Strand.POSITIVE);
            cigarStringsForDupSeqOnCtg = new ArrayList<>(2);
            if (firstContigRegion.forwardStrand) {
                cigarStringsForDupSeqOnCtg.add( TextCigarCodec.encode(extractCigarForTandupExpansion(firstContigRegion, distances.leftAlnRefEnd, distances.rightAlnRefStart)) );
                cigarStringsForDupSeqOnCtg.add( TextCigarCodec.encode(extractCigarForTandupExpansion(secondContigRegion, distances.leftAlnRefEnd, distances.rightAlnRefStart)) );
            } else {
                cigarStringsForDupSeqOnCtg.add( TextCigarCodec.encode(CigarUtils.invertCigar(extractCigarForTandupExpansion(firstContigRegion, distances.leftAlnRefEnd, distances.rightAlnRefStart))) );
                cigarStringsForDupSeqOnCtg.add( TextCigarCodec.encode(CigarUtils.invertCigar(extractCigarForTandupExpansion(secondContigRegion, distances.leftAlnRefEnd, distances.rightAlnRefStart))) );
                Collections.reverse(cigarStringsForDupSeqOnCtg);
            }
        }

        /**
         * Given a {@link AlignmentInterval} from a pair of ARs that forms a {@link SimpleChimera} signalling a tandem duplication,
         * extract a CIGAR from the {@link AlignmentInterval#cigarAlong5to3DirectionOfContig}
         * that corresponds to the alignment between the suspected repeated sequence on reference between
         * [{@code alignmentIntervalTwoReferenceIntervalSpanBegin}, {@code alignmentIntervalOneReferenceIntervalSpanEnd}],
         * and the sequence in {@link AlignmentInterval#referenceSpan}.
         */
        @VisibleForTesting
        static Cigar extractCigarForTandupExpansion(final AlignmentInterval contigRegion,
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

        private void resolveComplicationForSimpleTandupContraction(final SimpleChimera simpleChimera,
                                                                   final DistancesBetweenAlignmentsOnRefAndOnRead distances,
                                                                   final byte[] contigSeq, final boolean firstAfterSecond) {
            final SimpleInterval leftReferenceInterval;
            if (simpleChimera.isForwardStrandRepresentation)
                leftReferenceInterval = simpleChimera.regionWithLowerCoordOnContig.referenceSpan;
            else
                leftReferenceInterval = simpleChimera.regionWithHigherCoordOnContig.referenceSpan;
            homologyForwardStrandRep = inferHomology(simpleChimera.regionWithLowerCoordOnContig, simpleChimera.regionWithHigherCoordOnContig, contigSeq, firstAfterSecond);
            dupSeqRepeatUnitRefSpan  = new SimpleInterval(leftReferenceInterval.getContig(), distances.leftAlnRefEnd - ( distances.firstAlnCtgEnd - distances.secondAlnCtgStart ), distances.leftAlnRefEnd);
            dupSeqRepeatNumOnRef     = 2;
            dupSeqRepeatNumOnCtg     = 1;
            dupSeqStrandOnRef        = Arrays.asList(Strand.POSITIVE, Strand.POSITIVE);
            dupSeqStrandOnCtg        = Collections.singletonList(Strand.POSITIVE);
            cigarStringsForDupSeqOnCtg = DEFAULT_CIGAR_STRINGS_FOR_DUP_SEQ_ON_CTG;
        }


        protected void serialize(final Kryo kryo, final Output output) {
            super.serialize(kryo, output);
            output.writeInt(cigarStringsForDupSeqOnCtg.size());
            cigarStringsForDupSeqOnCtg.forEach(output::writeString);
        }
        public static final class Serializer extends com.esotericsoftware.kryo.Serializer<SmallDuplicationWithPreciseDupRangeBreakpointComplications> {
            @Override
            public void write(final Kryo kryo, final Output output, final SmallDuplicationWithPreciseDupRangeBreakpointComplications breakpointComplications) {
                breakpointComplications.serialize(kryo, output);
            }

            @Override
            public SmallDuplicationWithPreciseDupRangeBreakpointComplications read(final Kryo kryo, final Input input, final Class<SmallDuplicationWithPreciseDupRangeBreakpointComplications> klass ) {
                return new SmallDuplicationWithPreciseDupRangeBreakpointComplications(kryo, input);
            }
        }

        @Override
        public boolean equals(final Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            if (!super.equals(o)) return false;

            final SmallDuplicationWithPreciseDupRangeBreakpointComplications that = (SmallDuplicationWithPreciseDupRangeBreakpointComplications) o;

            return cigarStringsForDupSeqOnCtg.equals(that.cigarStringsForDupSeqOnCtg);
        }

        @Override
        public int hashCode() {
            int result = super.hashCode();
            result = 31 * result + cigarStringsForDupSeqOnCtg.hashCode();
            return result;
        }
    }

    /**
     * This is for dealing with case when the duplicated range could NOT be inferred exactly,
     * but only from a simple optimization scheme.
     */
    @DefaultSerializer(SmallDuplicationWithImpreciseDupRangeBreakpointComplications.Serializer.class)
    public static final class SmallDuplicationWithImpreciseDupRangeBreakpointComplications extends SmallDuplicationBreakpointComplications {
        private SimpleInterval impreciseDupAffectedRefRange = null;

        SimpleInterval getImpreciseDupAffectedRefRange() {
            return impreciseDupAffectedRefRange;
        }

        @VisibleForTesting
        SmallDuplicationWithImpreciseDupRangeBreakpointComplications(final String homologyForwardStrandRep, final String insertedSequenceForwardStrandRep,
                                                                     final SimpleInterval dupSeqRepeatUnitRefSpan,
                                                                     final int dupSeqRepeatNumOnRef, final int dupSeqRepeatNumOnCtg,
                                                                     final List<Strand> dupSeqStrandOnRef, final List<Strand> dupSeqStrandOnCtg,
                                                                     final SimpleInterval impreciseDupAffectedRefRange) {
            super(homologyForwardStrandRep, insertedSequenceForwardStrandRep, dupSeqRepeatUnitRefSpan, dupSeqRepeatNumOnRef, dupSeqRepeatNumOnCtg, dupSeqStrandOnRef, dupSeqStrandOnCtg);
            this.impreciseDupAffectedRefRange = impreciseDupAffectedRefRange;
        }

        private SmallDuplicationWithImpreciseDupRangeBreakpointComplications(final Kryo kryo, final Input input) {
            super(kryo, input);
            final String chr = input.readString();
            final int start = input.readInt();
            final int end = input.readInt();
            impreciseDupAffectedRefRange = new SimpleInterval(chr, start, end);
        }

        @Override
        public final String toString() {
            return super.toString() + "\timprecise:" ;
        }

        @Override
        public final Map<String, Object> toVariantAttributes() {
            final Map<String, Object> parentAttributesToBeFilled = super.toVariantAttributes();
            parentAttributesToBeFilled.put(GATKSVVCFConstants.DUP_ANNOTATIONS_IMPRECISE, "");
            parentAttributesToBeFilled.put(GATKSVVCFConstants.DUP_IMPRECISE_AFFECTED_RANGE, impreciseDupAffectedRefRange.toString());
            return parentAttributesToBeFilled;
        }

        SmallDuplicationWithImpreciseDupRangeBreakpointComplications(final SimpleChimera simpleChimera, final byte[] contigSeq, final boolean firstAfterSecond) {
            final DistancesBetweenAlignmentsOnRefAndOnRead distances = simpleChimera.getDistancesBetweenAlignmentsOnRefAndOnRead();

            if (distances.gapBetweenAlignRegionsOnRef > 0)
                throw new GATKException.ShouldNeverReachHereException(
                        "Simple chimera being sent down the wrong path, where the signature indicates " +
                        "a simple deletion but complication being resolve for small duplication. \n" + simpleChimera.toString());

            final boolean b = distances.gapBetweenAlignRegionsOnRef < 0 &&  distances.gapBetweenAlignRegionsOnCtg < 0;
            if ( !b) {
                throw new GATKException.ShouldNeverReachHereException(
                        "Simple chimera being sent down the wrong path, where the signature indicates " +
                        "simple duplication but complication being resolved for complex small duplication. \n" + simpleChimera.toString());
            }

            // most complicated case, see below
            // Deletion:  duplication with repeat number N1 on reference, N2 on contig, such that N1 <= 2*N2 (and N2<N1);
            // Insertion: duplication with repeat number N1 on reference, N2 on contig, such that N2 <= 2*N1 (and N1<N2);
            // in both cases, the equal sign on the right can be taken only when there's pseudo-homology between starting bases of the duplicated sequence and starting bases of the right flanking region
            // the reference system with a shorter overlap (i.e. with less-negative distance between regions) has a higher repeat number

            final TandemRepeatStructure duplicationComplication =
                    new TandemRepeatStructure(distances.gapBetweenAlignRegionsOnRef, distances.gapBetweenAlignRegionsOnCtg);

            final boolean isExpansion     = distances.gapBetweenAlignRegionsOnRef < distances.gapBetweenAlignRegionsOnCtg;

            final int repeatUnitSpanStart = distances.leftAlnRefEnd/*r1e*/ - duplicationComplication.pseudoHomologyLen
                    - duplicationComplication.repeatedSeqLen * duplicationComplication.lowerRepeatNumberEstimate
                    + 1;
            final int repeatUnitSpanEnd   = repeatUnitSpanStart + duplicationComplication.repeatedSeqLen - 1;
            homologyForwardStrandRep      = inferHomology(simpleChimera.regionWithLowerCoordOnContig, simpleChimera.regionWithHigherCoordOnContig, contigSeq, firstAfterSecond);

            dupSeqRepeatUnitRefSpan       = new SimpleInterval(simpleChimera.regionWithLowerCoordOnContig.referenceSpan.getContig(), repeatUnitSpanStart, repeatUnitSpanEnd);
            dupSeqRepeatNumOnRef          = isExpansion ? duplicationComplication.lowerRepeatNumberEstimate
                                                        : duplicationComplication.higherRepeatNumberEstimate;
            dupSeqRepeatNumOnCtg          = isExpansion ? duplicationComplication.higherRepeatNumberEstimate
                                                        : duplicationComplication.lowerRepeatNumberEstimate;
            dupSeqStrandOnRef             = new ArrayList<>(Collections.nCopies(dupSeqRepeatNumOnRef, Strand.POSITIVE));
            dupSeqStrandOnCtg             = new ArrayList<>(Collections.nCopies(dupSeqRepeatNumOnCtg, Strand.POSITIVE));

            impreciseDupAffectedRefRange = computeAffectedRefRegion(simpleChimera, isExpansion, distances, duplicationComplication);

            if ( insertedSequenceForwardStrandRep.isEmpty() ){
                if ( dupSeqRepeatNumOnCtg != dupSeqRepeatNumOnRef && null == dupSeqRepeatUnitRefSpan )
                    throw new GATKException.ShouldNeverReachHereException(
                            "An identified breakpoint pair seem to suggest insertion " +
                            "but the inserted sequence is empty: " + simpleChimera.toString());
            }
        }

        private static SimpleInterval computeAffectedRefRegion(final SimpleChimera simpleChimera, final boolean isExpansion,
                                                               final DistancesBetweenAlignmentsOnRefAndOnRead distances,
                                                               final TandemRepeatStructure duplicationComplication) {
            final SimpleInterval leftRefSpan, rightRefSpan;
            if (simpleChimera.isForwardStrandRepresentation) {
                leftRefSpan = simpleChimera.regionWithLowerCoordOnContig.referenceSpan;
                rightRefSpan = simpleChimera.regionWithHigherCoordOnContig.referenceSpan;
            } else {
                leftRefSpan = simpleChimera.regionWithHigherCoordOnContig.referenceSpan;
                rightRefSpan = simpleChimera.regionWithLowerCoordOnContig.referenceSpan;
            }
            if (isExpansion) {
                return new SimpleInterval(leftRefSpan.getContig(), distances.rightAlnRefStart, distances.leftAlnRefEnd);
            } else {
                return new SimpleInterval(leftRefSpan.getContig(),
                        rightRefSpan.getStart() - duplicationComplication.repeatedSeqLen,
                        leftRefSpan.getEnd() + duplicationComplication.repeatedSeqLen);
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

        protected void serialize(final Kryo kryo, final Output output) {
            super.serialize(kryo, output);
            output.writeString(impreciseDupAffectedRefRange.getContig());
            output.writeInt(impreciseDupAffectedRefRange.getStart());
            output.writeInt(impreciseDupAffectedRefRange.getEnd());
        }
        public static final class Serializer extends com.esotericsoftware.kryo.Serializer<SmallDuplicationWithImpreciseDupRangeBreakpointComplications> {
            @Override
            public void write(final Kryo kryo, final Output output, final SmallDuplicationWithImpreciseDupRangeBreakpointComplications breakpointComplications) {
                breakpointComplications.serialize(kryo, output);
            }

            @Override
            public SmallDuplicationWithImpreciseDupRangeBreakpointComplications read(final Kryo kryo, final Input input, final Class<SmallDuplicationWithImpreciseDupRangeBreakpointComplications> klass ) {
                return new SmallDuplicationWithImpreciseDupRangeBreakpointComplications(kryo, input);
            }
        }

        @Override
        public boolean equals(final Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            if (!super.equals(o)) return false;

            final SmallDuplicationWithImpreciseDupRangeBreakpointComplications that = (SmallDuplicationWithImpreciseDupRangeBreakpointComplications) o;

            return impreciseDupAffectedRefRange.equals(that.impreciseDupAffectedRefRange);
        }

        @Override
        public int hashCode() {
            int result = super.hashCode();
            result = 31 * result + impreciseDupAffectedRefRange.hashCode();
            return result;
        }
    }

    /**
     * For events we likely have assembled across one of its breakpoint but not the whole event,
     * hence only an BND record (or two, if both mates) could be generated from the signalling simple chimera.
     */
    abstract static class BNDTypeBreakpointComplications extends BreakpointComplications {
        protected BNDTypeBreakpointComplications() {}

        @VisibleForTesting
        BNDTypeBreakpointComplications(final String homologyForwardStrandRep, final String insertedSequenceForwardStrandRep) {
            super(homologyForwardStrandRep, insertedSequenceForwardStrandRep);
        }

        protected BNDTypeBreakpointComplications(final SimpleChimera simpleChimera, final byte[] contigSeq, final boolean firstAfterSecond) {
            homologyForwardStrandRep = inferHomology(simpleChimera.regionWithLowerCoordOnContig,
                                                    simpleChimera.regionWithHigherCoordOnContig, contigSeq, firstAfterSecond);
            insertedSequenceForwardStrandRep = inferInsertedSequence(simpleChimera.regionWithLowerCoordOnContig,
                                                                    simpleChimera.regionWithHigherCoordOnContig, contigSeq, firstAfterSecond);
        }

        protected BNDTypeBreakpointComplications(final Kryo kryo, final Input input) {
            super(kryo, input);
        }

        protected void serialize(final Kryo kryo, final Output output) {
            super.serialize(kryo, output);
        }
    }

    /**
     * For novel adjacency between reference locations that are
     * on the same chromosome, and with a strand switch.
     *
     * But the the signaling simple chimera cannot have its alignment overlap on reference.
     * If they do, {@link InvertedDuplicationBreakpointComplications}.
     *
     */
    @DefaultSerializer(IntraChrStrandSwitchBreakpointComplications.Serializer.class)
    public static final class IntraChrStrandSwitchBreakpointComplications extends BNDTypeBreakpointComplications {

        //works for constructing a simple inversion breakpoint, without any duplication annotations for inverted duplications.
        @VisibleForTesting
        IntraChrStrandSwitchBreakpointComplications(final String homologyForwardStrandRep, final String insertedSequenceForwardStrandRep) {
            super(homologyForwardStrandRep, insertedSequenceForwardStrandRep);
        }

        IntraChrStrandSwitchBreakpointComplications(final Kryo kryo, final Input input) {
            super(kryo, input);
        }

        IntraChrStrandSwitchBreakpointComplications(final SimpleChimera simpleChimera, final byte[] contigSeq, final boolean firstAfterSecond) {
            resolveComplicationForSimpleStrandSwitch(simpleChimera, contigSeq, firstAfterSecond);
        }

        void resolveComplicationForSimpleStrandSwitch(final SimpleChimera simpleChimera, final byte[] contigSeq, final boolean firstAfterSecond) {

            final AlignmentInterval firstAlignmentInterval  = simpleChimera.regionWithLowerCoordOnContig;
            final AlignmentInterval secondAlignmentInterval = simpleChimera.regionWithHigherCoordOnContig;

            homologyForwardStrandRep = inferHomology(firstAlignmentInterval, secondAlignmentInterval, contigSeq, firstAfterSecond);
            insertedSequenceForwardStrandRep = inferInsertedSequence(firstAlignmentInterval, secondAlignmentInterval, contigSeq, firstAfterSecond);
        }

        public static final class Serializer extends com.esotericsoftware.kryo.Serializer<IntraChrStrandSwitchBreakpointComplications> {
            @Override
            public void write(final Kryo kryo, final Output output, final IntraChrStrandSwitchBreakpointComplications breakpointComplications) {
                breakpointComplications.serialize(kryo, output);
            }

            @Override
            public IntraChrStrandSwitchBreakpointComplications read(final Kryo kryo, final Input input, final Class<IntraChrStrandSwitchBreakpointComplications> klass ) {
                return new IntraChrStrandSwitchBreakpointComplications(kryo, input);
            }
        }
    }

    /**
     * For this specific complication, we support a what could be defined as incomplete picture,
     * that involves inverted duplication:
     *
     * two overlapping alignments to reference
     *
     * first alignment:  -------------------->
     * second alignment:            <---------------------
     *
     *                              |--------||----------|
     *                                 Seg.1      Seg.2
     *
     * At least Seg.1 is invert duplicated, and
     * Seg.2 is inverted trans-inserted between the two copies (one of which is inverted).
     *
     * Yet because such contigs are currently defined to have an incomplete picture,
     * {@link org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments#hasIncompletePictureFromTwoAlignments(AlignmentInterval, AlignmentInterval)}
     * the annotations are actually never used.
     */
    public static final class InvertedDuplicationBreakpointComplications extends BNDTypeBreakpointComplications {
        static final List<Strand> DEFAULT_INV_DUP_REF_ORIENTATION = Collections.singletonList(Strand.POSITIVE);
        static final List<Strand> DEFAULT_INV_DUP_CTG_ORIENTATIONS_FR = Arrays.asList(Strand.POSITIVE, Strand.NEGATIVE);
        static final List<Strand> DEFAULT_INV_DUP_CTG_ORIENTATIONS_RF = Arrays.asList(Strand.NEGATIVE, Strand.POSITIVE);

        private SimpleInterval dupSeqRepeatUnitRefSpan = null;
        private int dupSeqRepeatNumOnRef = 0;
        private int dupSeqRepeatNumOnCtg = 0;
        private List<Strand> dupSeqStrandOnRef = null;
        private List<Strand> dupSeqStrandOnCtg = null;
        private List<String> cigarStringsForDupSeqOnCtg = null;
        private boolean dupAnnotIsFromOptimization = false;
        private SimpleInterval invertedTransInsertionRefSpan = null; // TODO: 10/2/17 see ticket 3647

        // note the following getters might return null, users should take care test null.
        public SimpleInterval getDupSeqRepeatUnitRefSpan() {
            return dupSeqRepeatUnitRefSpan;
        }
        int getDupSeqRepeatNumOnRef() {
            return dupSeqRepeatNumOnRef;
        }
        int getDupSeqRepeatNumOnCtg() {
            return dupSeqRepeatNumOnCtg;
        }
        List<Strand> getDupSeqOrientationsOnCtg() {
            return dupSeqStrandOnCtg;
        }
        List<String> getCigarStringsForDupSeqOnCtg() {
            return cigarStringsForDupSeqOnCtg;
        }
        SimpleInterval getInvertedTransInsertionRefSpan() {
            return invertedTransInsertionRefSpan;
        }
        boolean isDupAnnotIsFromOptimization() {
            return dupAnnotIsFromOptimization;
        }


        @VisibleForTesting
        InvertedDuplicationBreakpointComplications(final String homologyForwardStrandRep, final String insertedSequenceForwardStrandRep,
                                                   final SimpleInterval dupSeqRepeatUnitRefSpan,
                                                   final int dupSeqRepeatNumOnRef, final int dupSeqRepeatNumOnCtg,
                                                   final List<Strand> dupSeqStrandOnRef, final List<Strand> dupSeqStrandOnCtg,
                                                   final List<String> cigarStringsForDupSeqOnCtg, final boolean dupAnnotIsFromOptimization,
                                                   final SimpleInterval invertedTransInsertionRefSpan) {
            super(homologyForwardStrandRep, insertedSequenceForwardStrandRep);
            this.dupSeqRepeatUnitRefSpan = dupSeqRepeatUnitRefSpan;
            this.dupSeqRepeatNumOnRef = dupSeqRepeatNumOnRef;
            this.dupSeqRepeatNumOnCtg = dupSeqRepeatNumOnCtg;
            this.dupSeqStrandOnRef = dupSeqStrandOnRef;
            this.dupSeqStrandOnCtg = dupSeqStrandOnCtg;
            this.cigarStringsForDupSeqOnCtg = cigarStringsForDupSeqOnCtg;
            this.dupAnnotIsFromOptimization = dupAnnotIsFromOptimization;
            this.invertedTransInsertionRefSpan = invertedTransInsertionRefSpan;
        }

        InvertedDuplicationBreakpointComplications(final SimpleChimera simpleChimera, final byte[] contigSeq, final boolean firstAfterSecond) {
            resolveComplicationForInvDup(simpleChimera, contigSeq, firstAfterSecond);
        }

        InvertedDuplicationBreakpointComplications(final Kryo kryo, final Input input) {
            super(kryo, input);
            String ctg = input.readString();
            int start = input.readInt();
            int end = input.readInt();
            dupSeqRepeatUnitRefSpan = new SimpleInterval(ctg, start, end);
            dupSeqRepeatNumOnRef = input.readInt();
            dupSeqRepeatNumOnCtg = input.readInt();
            dupSeqStrandOnRef = new ArrayList<>(dupSeqRepeatNumOnRef);
            for (int i=0; i<dupSeqRepeatNumOnRef; ++i) {
                dupSeqStrandOnRef.add(Strand.values()[input.readInt()]);
            }
            dupSeqStrandOnCtg = new ArrayList<>(dupSeqRepeatNumOnCtg);
            for (int i=0; i<dupSeqRepeatNumOnCtg; ++i) {
                dupSeqStrandOnCtg.add(Strand.values()[input.readInt()]);
            }
            final int cigarCounts = input.readInt();
            cigarStringsForDupSeqOnCtg = new ArrayList<>(cigarCounts);
            for(int i = 0; i < cigarCounts; ++i) {
                cigarStringsForDupSeqOnCtg.add(input.readString());
            }
            dupAnnotIsFromOptimization = input.readBoolean();
            if (input.readBoolean()) {
                ctg = input.readString();
                start = input.readInt();
                end = input.readInt();
                invertedTransInsertionRefSpan = new SimpleInterval(ctg, start, end);
            }
        }

        /**
         * Initialize the fields in this object, assuming the input chimeric alignment is induced by two alignments with
         * "significant" (see {@link SimpleChimera#isCandidateInvertedDuplication()})
         * overlap on their reference spans.
         */
        private void resolveComplicationForInvDup(final SimpleChimera simpleChimera, final byte[] contigSeq, final boolean firstAfterSecond) {

            final AlignmentInterval firstAlignmentInterval  = simpleChimera.regionWithLowerCoordOnContig;
            final AlignmentInterval secondAlignmentInterval = simpleChimera.regionWithHigherCoordOnContig;

            // TODO: 8/8/17 this might be wrong regarding how strand is involved, fix it
            insertedSequenceForwardStrandRep = inferInsertedSequence(firstAlignmentInterval, secondAlignmentInterval, contigSeq, firstAfterSecond);

            dupSeqRepeatNumOnRef = 1;
            dupSeqRepeatNumOnCtg = 2;
            dupSeqStrandOnRef = DEFAULT_INV_DUP_REF_ORIENTATION;

            // jump start and jump landing locations
            final int jumpStartRefLoc = firstAlignmentInterval.forwardStrand ? firstAlignmentInterval.referenceSpan.getEnd()
                    : firstAlignmentInterval.referenceSpan.getStart();
            final int jumpLandingRefLoc = secondAlignmentInterval.forwardStrand ? secondAlignmentInterval.referenceSpan.getStart()
                    : secondAlignmentInterval.referenceSpan.getEnd();

            if (firstAlignmentInterval.forwardStrand) {
                final int alpha = firstAlignmentInterval.referenceSpan.getStart(),
                        omega = secondAlignmentInterval.referenceSpan.getStart();
                dupSeqRepeatUnitRefSpan = new SimpleInterval(firstAlignmentInterval.referenceSpan.getContig(),
                        Math.max(alpha, omega), Math.min(jumpStartRefLoc, jumpLandingRefLoc));
                if ( (alpha <= omega && jumpStartRefLoc < jumpLandingRefLoc) || (alpha > omega && jumpLandingRefLoc < jumpStartRefLoc) ) {
                    invertedTransInsertionRefSpan = new SimpleInterval(firstAlignmentInterval.referenceSpan.getContig(),
                            Math.min(jumpStartRefLoc, jumpLandingRefLoc) + 1, Math.max(jumpStartRefLoc, jumpLandingRefLoc));
                }
                dupSeqStrandOnCtg = DEFAULT_INV_DUP_CTG_ORIENTATIONS_FR;
            } else {
                final int alpha = firstAlignmentInterval.referenceSpan.getEnd(),
                        omega = secondAlignmentInterval.referenceSpan.getEnd();
                dupSeqRepeatUnitRefSpan = new SimpleInterval(firstAlignmentInterval.referenceSpan.getContig(),
                        Math.max(jumpStartRefLoc, jumpLandingRefLoc), Math.min(alpha, omega));
                if ( (alpha >= omega && jumpLandingRefLoc < jumpStartRefLoc) || (alpha < omega && jumpStartRefLoc < jumpLandingRefLoc) ) {
                    invertedTransInsertionRefSpan = new SimpleInterval(firstAlignmentInterval.referenceSpan.getContig(),
                            Math.min(jumpStartRefLoc, jumpLandingRefLoc) + 1, Math.max(jumpStartRefLoc, jumpLandingRefLoc));
                }
                dupSeqStrandOnCtg = DEFAULT_INV_DUP_CTG_ORIENTATIONS_RF;
            }
            // not computing cigars because alt haplotypes will be extracted
            cigarStringsForDupSeqOnCtg = SmallDuplicationWithPreciseDupRangeBreakpointComplications.DEFAULT_CIGAR_STRINGS_FOR_DUP_SEQ_ON_CTG;

            dupAnnotIsFromOptimization = false;
        }

        protected void serialize(final Kryo kryo, final Output output) {
            super.serialize(kryo, output);

            output.writeString(dupSeqRepeatUnitRefSpan.getContig());
            output.writeInt(dupSeqRepeatUnitRefSpan.getStart());
            output.writeInt(dupSeqRepeatUnitRefSpan.getEnd());
            output.writeInt(dupSeqRepeatNumOnRef);
            output.writeInt(dupSeqRepeatNumOnCtg);
            dupSeqStrandOnRef.forEach(s -> output.writeInt(s.ordinal()));
            dupSeqStrandOnCtg.forEach(s -> output.writeInt(s.ordinal()));
            output.writeInt(cigarStringsForDupSeqOnCtg.size());
            cigarStringsForDupSeqOnCtg.forEach(output::writeString);
            output.writeBoolean(dupAnnotIsFromOptimization);
            output.writeBoolean( invertedTransInsertionRefSpan != null );
            if (invertedTransInsertionRefSpan != null) {
                output.writeString(invertedTransInsertionRefSpan.getContig());
                output.writeInt(invertedTransInsertionRefSpan.getStart());
                output.writeInt(invertedTransInsertionRefSpan.getEnd());
            }
        }

        public static final class Serializer extends com.esotericsoftware.kryo.Serializer<InvertedDuplicationBreakpointComplications> {
            @Override
            public void write(final Kryo kryo, final Output output, final InvertedDuplicationBreakpointComplications breakpointComplications) {
                breakpointComplications.serialize(kryo, output);
            }

            @Override
            public InvertedDuplicationBreakpointComplications read(final Kryo kryo, final Input input, final Class<InvertedDuplicationBreakpointComplications> klass ) {
                return new InvertedDuplicationBreakpointComplications(kryo, input);
            }
        }

        @Override
        public final String toString() {
            String toPrint = super.toString();
            toPrint += String.format("\ttandem duplication repeat unit ref span: %s\t"+
                            "ref repeat num: %d\t"+
                            "ctg repeat num: %d\t"+
                            "dupSeqStrandOnRef: %s\t" +
                            "dupSeqStrandOnCtg: %s\t" +
                            "cigarStringsForDupSeqOnCtg: %s\t"+
                            "tandupAnnotationIsFromSimpleOptimization: %s\t" +
                            "invertedTransInsertionRefSpan: %s",
                    dupSeqRepeatUnitRefSpan == null ? "" : dupSeqRepeatUnitRefSpan,
                    dupSeqRepeatNumOnRef, dupSeqRepeatNumOnCtg,
                    dupSeqStrandOnRef == null ? "" : dupSeqStrandOnRef.stream().map(Strand::toString).collect(SVUtils.arrayListCollector(dupSeqStrandOnRef.size())).toString(),
                    dupSeqStrandOnCtg == null ? "" : dupSeqStrandOnCtg.stream().map(Strand::toString).collect(SVUtils.arrayListCollector(dupSeqStrandOnCtg.size())).toString(),
                    cigarStringsForDupSeqOnCtg == null ? "" : cigarStringsForDupSeqOnCtg,
                    isDupAnnotIsFromOptimization() ? "true" : "false",
                    invertedTransInsertionRefSpan == null ? "" : invertedTransInsertionRefSpan);
            return toPrint;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            if (!super.equals(o)) return false;

            InvertedDuplicationBreakpointComplications that = (InvertedDuplicationBreakpointComplications) o;

            if (dupSeqRepeatNumOnRef != that.dupSeqRepeatNumOnRef) return false;
            if (dupSeqRepeatNumOnCtg != that.dupSeqRepeatNumOnCtg) return false;
            if (dupAnnotIsFromOptimization != that.dupAnnotIsFromOptimization) return false;
            if (dupSeqRepeatUnitRefSpan != null ? !dupSeqRepeatUnitRefSpan.equals(that.dupSeqRepeatUnitRefSpan) : that.dupSeqRepeatUnitRefSpan != null)
                return false;
            if (dupSeqStrandOnRef != null ? !dupSeqStrandOnRef.equals(that.dupSeqStrandOnRef) : that.dupSeqStrandOnRef != null)
                return false;
            if (dupSeqStrandOnCtg != null ? !dupSeqStrandOnCtg.equals(that.dupSeqStrandOnCtg) : that.dupSeqStrandOnCtg != null)
                return false;
            if (cigarStringsForDupSeqOnCtg != null ? !cigarStringsForDupSeqOnCtg.equals(that.cigarStringsForDupSeqOnCtg) : that.cigarStringsForDupSeqOnCtg != null)
                return false;
            return invertedTransInsertionRefSpan != null ? invertedTransInsertionRefSpan.equals(that.invertedTransInsertionRefSpan) : that.invertedTransInsertionRefSpan == null;
        }
        @Override
        public int hashCode() {
            int result = super.hashCode();
            result = 31 * result + (dupSeqRepeatUnitRefSpan != null ? dupSeqRepeatUnitRefSpan.hashCode() : 0);
            result = 31 * result + dupSeqRepeatNumOnRef;
            result = 31 * result + dupSeqRepeatNumOnCtg;
            if (dupSeqStrandOnRef != null) {
                for (final Strand strand : dupSeqStrandOnRef) {
                    result = 31 * result + strand.ordinal() ;
                }
            }
            if (dupSeqStrandOnCtg != null) {
                for (final Strand strand : dupSeqStrandOnCtg) {
                    result = 31 * result + strand.ordinal() ;
                }
            }
            result = 31 * result + (cigarStringsForDupSeqOnCtg != null ? cigarStringsForDupSeqOnCtg.hashCode() : 0);
            result = 31 * result + (dupAnnotIsFromOptimization ? 1 : 0);
            result = 31 * result + (invertedTransInsertionRefSpan != null ? invertedTransInsertionRefSpan.hashCode() : 0);
            return result;
        }
    }

    /**
     * For novel adjacency between reference locations that are on the same chromosome,
     * WITHOUT strand switch but with order swap,
     * i.e. a base with higher coordinate on ref has lower coordinate on sample.
     *
     * e.g.
     * first alignment:                                                 -------------------->
     * second alignment: --------------------->
     *
     * novel adjacency   |                                                                  |
     */
    @DefaultSerializer(IntraChrRefOrderSwapBreakpointComplications.Serializer.class)
    static final class IntraChrRefOrderSwapBreakpointComplications extends BNDTypeBreakpointComplications {

        @VisibleForTesting
        IntraChrRefOrderSwapBreakpointComplications(final String homologyForwardStrandRep, final String insertedSequenceForwardStrandRep) {
            super(homologyForwardStrandRep, insertedSequenceForwardStrandRep);
        }

        IntraChrRefOrderSwapBreakpointComplications(final SimpleChimera simpleChimera, final byte[] contigSeq, final boolean firstAfterSecond) {
            super(simpleChimera, contigSeq, firstAfterSecond);
        }

        private IntraChrRefOrderSwapBreakpointComplications(final Kryo kryo, final Input input) {
            super(kryo, input);

        }
        protected void serialize(final Kryo kryo, final Output output) {
            super.serialize(kryo, output);
        }
        public static final class Serializer extends com.esotericsoftware.kryo.Serializer<IntraChrRefOrderSwapBreakpointComplications> {
            @Override
            public void write(final Kryo kryo, final Output output, final IntraChrRefOrderSwapBreakpointComplications breakpointComplications) {
                breakpointComplications.serialize(kryo, output);
            }

            @Override
            public IntraChrRefOrderSwapBreakpointComplications read(final Kryo kryo, final Input input, final Class<IntraChrRefOrderSwapBreakpointComplications> klass ) {
                return new IntraChrRefOrderSwapBreakpointComplications(kryo, input);
            }
        }
    }

    /**
     * For novel adjacency between reference locations that are on the different chromosomes,
     * WITH OR WITHOUT strand switch.
     */
    @DefaultSerializer(InterChromosomeBreakpointComplications.Serializer.class)
    static final class InterChromosomeBreakpointComplications extends BNDTypeBreakpointComplications {

        @VisibleForTesting
        InterChromosomeBreakpointComplications(final String homologyForwardStrandRep, final String insertedSequenceForwardStrandRep) {
            super(homologyForwardStrandRep, insertedSequenceForwardStrandRep);
        }

        InterChromosomeBreakpointComplications(final SimpleChimera simpleChimera, final byte[] contigSeq, final boolean firstAfterSecond) {
            super(simpleChimera, contigSeq, firstAfterSecond);
        }

        private InterChromosomeBreakpointComplications(final Kryo kryo, final Input input) {
            super(kryo, input);
        }
        protected void serialize(final Kryo kryo, final Output output) {
            super.serialize(kryo, output);
        }
        public static final class Serializer extends com.esotericsoftware.kryo.Serializer<InterChromosomeBreakpointComplications> {
            @Override
            public void write(final Kryo kryo, final Output output, final InterChromosomeBreakpointComplications breakpointComplications) {
                breakpointComplications.serialize(kryo, output);
            }

            @Override
            public InterChromosomeBreakpointComplications read(final Kryo kryo, final Input input, final Class<InterChromosomeBreakpointComplications> klass ) {
                return new InterChromosomeBreakpointComplications(kryo, input);
            }
        }
    }
}
