package org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;

/**
 * A wrapper around {@link AlignedContig} to represent mapped assembly contig whose alignments
 * went through {@link AssemblyContigAlignmentsConfigPicker} and may represent SV breakpoints.
 *
 * Example Contig:
 * ACTGACTGCACTGACTGCACTGACTGCACTGACTGCACTGACTGCACTGACTGCACTGACTGCACTGACTGCACTGACTGCACTGACTGCACTGACTGCACTGACTG
 * Alignment regions:
 * |---------1:100-200------------|
 *                                 |----------2:100-200------------------|
 *                                                                       |----------3:100-200-----------------|
 * Assembled breakpoints:
 * 1) links 1:100-200 to 2:100-200
 * 2) links 2:100-200 to 3:100-200
 *
 * A contig could have incomplete picture, see {@link #hasIncompletePicture()}
 */
@DefaultSerializer(AssemblyContigWithFineTunedAlignments.Serializer.class)
public final class AssemblyContigWithFineTunedAlignments {
    private static final AlignedContig.Serializer contigSerializer = new AlignedContig.Serializer();
    public static final List<String> emptyInsertionMappings = Collections.emptyList();

    private final AlignedContig sourceTig;
    // alignments that were given by aligner but thought to be non-good,
    // and better treated as not-so-reliable mappings for inserted sequence;
    // useful for annotation but not essential for reliable event interpretation
    private final List<String> insertionMappings;

    private final boolean hasEquallyGoodAlnConfigurations;

    /**
     * for signalling (i.e. not null) that the alignments went through the special treatment in
     * {@link AssemblyContigAlignmentsConfigPicker#getBetterNonCanonicalMapping(Set, List, int)}
     */
    private final String saTAGForGoodMappingToNonCanonicalChromosome;

    public enum AlignmentSignatureBasicType {
        UNKNOWN, SIMPLE, COMPLEX;
    }

    public enum ReasonForAlignmentClassificationFailure {
        INCOMPLETE, AMBIGUOUS, MIS_ASSEMBLY_OR_MAPPING_SUSPECT;
    }

    public AssemblyContigWithFineTunedAlignments(final AlignedContig contig) {
        this(Utils.nonNull(contig), emptyInsertionMappings, false, NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
    }

    // {@code goodMappingToNonCanonicalChromosome} could be null
    public AssemblyContigWithFineTunedAlignments(final AlignedContig contig,
                                                 final List<String> insertionMappings,
                                                 final boolean hasEquallyGoodAlnConfigurations,
                                                 final AlignmentInterval goodMappingToNonCanonicalChromosome) {
        this(Utils.nonNull(contig), Utils.nonNull(insertionMappings),
                hasEquallyGoodAlnConfigurations,
                goodMappingToNonCanonicalChromosome == null ? NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME
                                                            : goodMappingToNonCanonicalChromosome.toSATagString());
    }

    public AssemblyContigWithFineTunedAlignments(final AlignedContig contig,
                                                 final List<String> insertionMappings,
                                                 final boolean hasEquallyGoodAlnConfigurations,
                                                 final String saTAGForGoodMappingToNonCanonicalChromosome) {
        this.sourceTig = contig;
        this.insertionMappings = insertionMappings;
        this.hasEquallyGoodAlnConfigurations = hasEquallyGoodAlnConfigurations;
        this.saTAGForGoodMappingToNonCanonicalChromosome = saTAGForGoodMappingToNonCanonicalChromosome;
    }

    AssemblyContigWithFineTunedAlignments(final Kryo kryo, final Input input) {
        sourceTig = contigSerializer.read(kryo, input, AlignedContig.class);
        final int insertionMappingsSize = input.readInt();
        insertionMappings = new ArrayList<>(insertionMappingsSize);
        for (int i = 0; i < insertionMappingsSize; ++i) {
            insertionMappings.add(input.readString());
        }

        hasEquallyGoodAlnConfigurations = input.readBoolean();
        saTAGForGoodMappingToNonCanonicalChromosome = input.readString();
    }

    public AlignedContig getSourceContig() {
        return sourceTig;
    }

    public List<AlignmentInterval> getAlignments() {
        return sourceTig.getAlignments();
    }

    public AlignmentInterval getHeadAlignment() {
        return sourceTig.getHeadAlignment();
    }

    public AlignmentInterval getTailAlignment() {
        return sourceTig.getTailAlignment();
    }

    public String getContigName() {
        return sourceTig.getContigName();
    }

    public byte[] getContigSequence() {
        return sourceTig.getContigSequence();
    }

    public List<String> getInsertionMappings() {
        return insertionMappings;
    }

    // after fine tuning, a contig may have no good alignment left, or only 1
    public final boolean isInformative() {
        return sourceTig.getAlignments().size() > 1;
    }

    /**
     * See {@link AssemblyContigAlignmentsConfigPicker#getBetterNonCanonicalMapping(Set, List, int)}.
     */
    public static final String NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME = "NONE";
    public final String getSAtagForGoodMappingToNonCanonicalChromosome() {
        return saTAGForGoodMappingToNonCanonicalChromosome;
    }

    public final boolean hasOnly2GoodAlignments() {
        return sourceTig.hasOnly2Alignments();
    }

    public AlignmentSignatureBasicType getAlignmentSignatureBasicType() {
        if ( hasEquallyGoodAlnConfigurations || hasIncompletePicture() || (! isInformative()) ) {
            return AlignmentSignatureBasicType.UNKNOWN;
        } else if ( hasOnly2GoodAlignments() ) {
            return AlignmentSignatureBasicType.SIMPLE;
        } else {
            return AlignmentSignatureBasicType.COMPLEX;
        }
    }

    public ReasonForAlignmentClassificationFailure getReasonForAlignmentClassificationFailure() {

        if ( hasEquallyGoodAlnConfigurations ) { // ambiguous
            return ReasonForAlignmentClassificationFailure.AMBIGUOUS;
        } else if ( hasIncompletePicture() ) {
            return ReasonForAlignmentClassificationFailure.INCOMPLETE;
        } else if ( ! isInformative() ) { // un-informative contig, mis-assembly suspect
            return ReasonForAlignmentClassificationFailure.MIS_ASSEMBLY_OR_MAPPING_SUSPECT;
        } else {
            throw new UnsupportedOperationException(
                    "operating on contig without a suspicious alignment signature, contig: " + toString());
        }
    }

    public final boolean hasEquallyGoodAlnConfigurations() {
        return hasEquallyGoodAlnConfigurations;
    }

    /**
     * This method exist because extending the assembly contig in either direction
     * could (not always, see ticket #4229) change
     * <ul>
     *     <li>the interpretation and</li>
     *     <li>affected reference region</li>
     * </ul>
     * significantly.
     */
    final boolean hasIncompletePicture() {
        if ( hasOnly2GoodAlignments() )
            return hasIncompletePictureFromTwoAlignments(getHeadAlignment(), getTailAlignment());
        else
            return hasIncompletePictureFromMultipleAlignments();
    }

    //==================================================================================================================

    /**
     * This predicate tests if an assembly contig with two (picked) alignments has the
     * <ul>
     *     <li>full event, or</li>
     *     <li>breakpoint </li>
     * </ul>
     * assembled.
     *
     * If the assembly contig's alignments show any of the following signature,
     * it will be classified as incomplete.
     * it is definitely not giving the whole picture of the alt haplotype,
     * hence without other types of evidence (or linking breakpoints, which itself needs other evidence anyway),
     * human-friendly interpretation for them is unreliable.
     * <ul>
     *     <li>
     *         head and tail alignments contain each other in terms of their reference span,
     *         regardless if strand switch is involved;
     *     </li>
     *     <li>
     *         if the two alignments map to the same chromosome, AND their reference spans overlap, AND in the mean time
     *         <ul>
     *             <li>
     *                 there's reference order switch, or
     *             </li>
     *             <li>
     *                 there's strand switch.
     *             </li>
     *         </ul>
     *     </li>
     * </ul>
     */
    public static boolean hasIncompletePictureFromTwoAlignments(final AlignmentInterval headAlignment,
                                                                final AlignmentInterval tailAlignment) {

        final SimpleInterval referenceSpanOne = headAlignment.referenceSpan;
        final SimpleInterval referenceSpanTwo = tailAlignment.referenceSpan;

        // inter contig mapping will not be treated as incomplete picture for 2-alignment reads
        if ( !referenceSpanOne.getContig().equals(referenceSpanTwo.getContig()))
            return false;

        if ( oneRefSpanContainsTheOther(referenceSpanOne, referenceSpanTwo) )
            return true;

        // Note that this is a little different from the multiple alignment (>2) case:
        // here we allow ref order switch because we can emit BND record, and only one BND is necessary.
        // But for multiple alignments, one BND record is not enough and we need to figure out how to emit records for them.
        if (headAlignment.forwardStrand != tailAlignment.forwardStrand) {
            // TODO: 10/29/17 this obsoletes the inverted duplication call code we have now,
            // but those could be used to figure out how to annotate which known ref regions are invert duplicated
            return referenceSpanOne.overlaps(referenceSpanTwo);
        } else {
            if (headAlignment.forwardStrand) {
                return referenceSpanOne.getStart() > referenceSpanTwo.getStart() &&
                        referenceSpanOne.getStart() <= referenceSpanTwo.getEnd();
            } else {
                return referenceSpanTwo.getStart() > referenceSpanOne.getStart() &&
                        referenceSpanTwo.getStart() <= referenceSpanOne.getEnd();
            }
        }
    }

    //==================================================================================================================

    /**
     * This predicate tests if an assembly contig has the full event (i.e. alt haplotype) assembled
     * Of course, the grand problem of SV is always not getting the big-enough picture
     * but here we have a more workable definition of what is definitely not big-enough:
     *
     * If the assembly contig, with its (picked) alignments, shows any of the following signature,
     * it is definitely not giving the whole picture of the alt haplotype,
     * hence without other types of evidence (or linking breakpoints, which itself needs other evidence anyway),
     * human-friendly interpretation for them is unreliable.
     * <ul>
     *     <li>
     *         head and tail alignment contain each other in terms of their reference span
     *     </li>
     *     <li>
     *         head and tail alignment mapped to different chromosome
     *     </li>
     *     <li>
     *         head and tail alignment mapped to different strands
     *     </li>
     *     <li>
     *         head and tail alignment have reference order switch, indicating we have likely
     *         assembled across a complicated tandem duplication breakpoint
     *     </li>
     *     <li>
     *         valid region, as defined by the region bounded by the outer most boundaries of head and tail alignments,
     *         overlaps with any of the middle alignments but does not completely contain them
     *         (note that middle alignments completely disjoint from the valid region are OK,
     *          they indicate remote parts of reference being duplicated and inserted between head/tail)
     *     </li>
     * </ul>
     * In summary, the assembly contig should "resume the normal flow" as defined by the reference,
     * even though we allow the flow to be interrupted in the middle.
     * This is similar to the case where a particular mapping is unreliable without confident left and right flanking.
     */
    private boolean hasIncompletePictureFromMultipleAlignments() {
        final int goodAlignmentsCount = sourceTig.getAlignments().size();

        final AlignmentInterval head = sourceTig.getHeadAlignment(),
                                tail = sourceTig.getTailAlignment();

        // mapped to different chr or strand
        if ( !head.referenceSpan.getContig().equals(tail.referenceSpan.getContig()) || head.forwardStrand != tail.forwardStrand)
            return true;

        final SimpleInterval referenceSpanHead = head.referenceSpan,
                             referenceSpanTail = tail.referenceSpan;

        // head or tail's ref span contained in the other's
        if (oneRefSpanContainsTheOther(referenceSpanHead, referenceSpanTail))
            return true;

        // middle alignments' ref span should be either 1) disjoint from valid region, or 2) completely contained in valid region
        final SimpleInterval validRegion = new SimpleInterval(referenceSpanHead.getContig(),
                                                              Math.min(referenceSpanHead.getStart(), referenceSpanTail.getStart()),
                                                              Math.max(referenceSpanHead.getEnd(), referenceSpanTail.getEnd()));
        final boolean notCompleteDupRegion =
                sourceTig.getAlignments()
                        .subList(1, goodAlignmentsCount - 1) // take middle (no head no tail) alignments' ref span
                        .stream().map(ai -> ai.referenceSpan)
                        .anyMatch( middle -> middle.overlaps(validRegion) && !validRegion.contains(middle));
        if (notCompleteDupRegion)
            return true;

        // reference order switch, indicating a complicated tandem duplication breakpoint
        if (head.forwardStrand) {
            return referenceSpanHead.getStart() >= referenceSpanTail.getStart();
        } else {
            return referenceSpanHead.getEnd() <= referenceSpanTail.getEnd();
        }
    }

    //==================================================================================================================

    private static boolean oneRefSpanContainsTheOther(final SimpleInterval one, final SimpleInterval two) {
        return one.contains(two) || two.contains(one);
    }

    void serialize(final Kryo kryo, final Output output) {
        contigSerializer.write(kryo, output, sourceTig);
        output.writeInt(insertionMappings.size());
        for (final String mapping : insertionMappings) {
            output.writeString(mapping);
        }
        output.writeBoolean(hasEquallyGoodAlnConfigurations);
       output.writeString(saTAGForGoodMappingToNonCanonicalChromosome);
    }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<AssemblyContigWithFineTunedAlignments> {
        @Override
        public void write(final Kryo kryo, final Output output, final AssemblyContigWithFineTunedAlignments alignedContig) {
            alignedContig.serialize(kryo, output);
        }

        @Override
        public AssemblyContigWithFineTunedAlignments read(final Kryo kryo, final Input input, final Class<AssemblyContigWithFineTunedAlignments> clazz) {
            return new AssemblyContigWithFineTunedAlignments(kryo, input);
        }
    }

    @Override
    public String toString() {
        final StringBuilder sb = new StringBuilder("AssemblyContigWithFineTunedAlignments{");
        sb.append("sourceTig=").append(sourceTig);
        sb.append(", insertionMappings=").append(insertionMappings);
        sb.append(", hasEquallyGoodAlnConfigurations=").append(hasEquallyGoodAlnConfigurations);
        sb.append(", saTAGForGoodMappingToNonCanonicalChromosome='").append(saTAGForGoodMappingToNonCanonicalChromosome).append('\'');
        sb.append('}');
        return sb.toString();
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final AssemblyContigWithFineTunedAlignments that = (AssemblyContigWithFineTunedAlignments) o;

        if (hasEquallyGoodAlnConfigurations != that.hasEquallyGoodAlnConfigurations) return false;
        if (!sourceTig.equals(that.sourceTig)) return false;
        if (!insertionMappings.equals(that.insertionMappings)) return false;
        return saTAGForGoodMappingToNonCanonicalChromosome.equals(that.saTAGForGoodMappingToNonCanonicalChromosome);
    }

    @Override
    public int hashCode() {
        int result = sourceTig.hashCode();
        result = 31 * result + insertionMappings.hashCode();
        result = 31 * result + (hasEquallyGoodAlnConfigurations ? 1 : 0);
        result = 31 * result + saTAGForGoodMappingToNonCanonicalChromosome.hashCode();
        return result;
    }
}
