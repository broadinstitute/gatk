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

@DefaultSerializer(AssemblyContigWithFineTunedAlignments.Serializer.class)
public final class AssemblyContigWithFineTunedAlignments {
    private static final AlignedContig.Serializer contigSerializer = new AlignedContig.Serializer();
    public static final List<String> emptyInsertionMappings = Collections.emptyList();

    private final AlignedContig sourceTig;
    // alignments that were given by aligner but thought to be non-good,
    // and better treated as not-so-reliable mappings for inserted sequence;
    // useful for annotation but not essential for reliable event interpretation
    private final List<String> insertionMappings;

    /**
     * for signalling (i.e. not null) that the alignments went through the special treatment in
     * {@link AssemblyContigAlignmentsConfigPicker#getBetterNonCanonicalMapping(Set, List, int)}
     */
    private final String saTAGForGoodMappingToNonCanonicalChromosome;

    public AssemblyContigWithFineTunedAlignments(final AlignedContig contig) {
        this(Utils.nonNull(contig), emptyInsertionMappings, NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
    }

    // {@code goodMappingToNonCanonicalChromosome} could be null
    public AssemblyContigWithFineTunedAlignments(final AlignedContig contig,
                                                 final List<String> insertionMappings,
                                                 final AlignmentInterval goodMappingToNonCanonicalChromosome) {
        this(Utils.nonNull(contig), Utils.nonNull(insertionMappings),
                goodMappingToNonCanonicalChromosome == null ? NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME
                                                            : goodMappingToNonCanonicalChromosome.toSATagString());
    }

    public AssemblyContigWithFineTunedAlignments(final AlignedContig contig,
                                                 final List<String> insertionMappings,
                                                 final String saTAGForGoodMappingToNonCanonicalChromosome) {
        this.sourceTig = contig;
        this.insertionMappings = insertionMappings;
        this.saTAGForGoodMappingToNonCanonicalChromosome = saTAGForGoodMappingToNonCanonicalChromosome;
    }

    AssemblyContigWithFineTunedAlignments(final Kryo kryo, final Input input) {
        sourceTig = contigSerializer.read(kryo, input, AlignedContig.class);
        final int insertionMappingsSize = input.readInt();
        insertionMappings = new ArrayList<>(insertionMappingsSize);
        for (int i = 0; i < insertionMappingsSize; ++i) {
            insertionMappings.add(input.readString());
        }

        saTAGForGoodMappingToNonCanonicalChromosome = input.readString();
    }

    public AlignedContig getSourceContig() {
        return sourceTig;
    }

    public List<AlignmentInterval> getAlignments() {
        return sourceTig.alignmentIntervals;
    }

    public List<String> getInsertionMappings() {
        return insertionMappings;
    }

    // after fine tuning, a contig may have no good alignment left, or only 1
    public final boolean isInformative() {
        return sourceTig.isInformative();
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

    public final boolean hasIncompletePicture() {
        if ( hasOnly2GoodAlignments() )
            return hasIncompletePictureFromTwoAlignments();
        else
            return hasIncompletePictureFromMultipleAlignments();
    }

    public final boolean hasEquallyGoodAlnConfigurations() {
        return sourceTig.hasEquallyGoodAlnConfigurations;
    }

    //==================================================================================================================

    /**
     * An assembly contig with two good alignments indicates to us that
     * it doesn't have the complete alt haplotype assembled if it shows
     * <ul>
     *     <li>
     *        either alignment's ref span is contained in the other's;
     *     </li>
     *
     *    <li>
     *        if the two alignments map to the same chromosome, and their ref span indicate order switch AND overlap
     *        (this indicate a duplication, inverted or not, but we don't have the full duplicated region assembled)
     *    </li>
     * </ul>
     *
     */
    public boolean hasIncompletePictureFromTwoAlignments() {
        return hasIncompletePictureDueToRefSpanContainment()
                ||
                (firstAndLastAlignmentMappedToSameChr() && hasIncompletePictureDueToOverlappingRefOrderSwitch());
    }

    private boolean hasIncompletePictureDueToRefSpanContainment() {
        return oneRefSpanContainsTheOther(sourceTig.alignmentIntervals.get(0).referenceSpan,
                                          sourceTig.alignmentIntervals.get(1).referenceSpan);
    }

    /**
     * Note that this is a little different from the multiple alignment case (>2):
     * here we allow ref order switch because we can emit BND record, and only one BND is necessary.
     * But for multiple alignments, one BND record is not enough and we need to figure out how to emit records for them.
     */
    private boolean hasIncompletePictureDueToOverlappingRefOrderSwitch() {

        final AlignmentInterval one = sourceTig.alignmentIntervals.get(0);
        final AlignmentInterval two = sourceTig.alignmentIntervals.get(1);
        final SimpleInterval referenceSpanOne = one.referenceSpan;
        final SimpleInterval referenceSpanTwo = two.referenceSpan;

        if (oneRefSpanContainsTheOther(referenceSpanOne, referenceSpanTwo))
            return true;

        if (one.forwardStrand != two.forwardStrand) {
            // TODO: 10/29/17 this obsoletes the inverted duplication call code we have now,
            // but those could be used to figure out how to annotate which known ref regions are invert duplicated
            return referenceSpanOne.overlaps(referenceSpanTwo);
        } else {
            if (one.forwardStrand) {
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
        final int goodAlignmentsCount = sourceTig.alignmentIntervals.size();

        final AlignmentInterval head = sourceTig.alignmentIntervals.get(0),
                                tail = sourceTig.alignmentIntervals.get(goodAlignmentsCount-1);

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
                sourceTig.alignmentIntervals
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

    public boolean firstAndLastAlignmentMappedToSameChr() {

        final String firstMappedChr = this.sourceTig.alignmentIntervals.get(0).referenceSpan.getContig();
        final String lastMappedChr  = this.sourceTig.alignmentIntervals.get(this.sourceTig.alignmentIntervals.size() - 1).referenceSpan.getContig();

        return firstMappedChr.equals(lastMappedChr);
    }

    void serialize(final Kryo kryo, final Output output) {
        contigSerializer.write(kryo, output, sourceTig);
        output.writeInt(insertionMappings.size());
        for (final String mapping : insertionMappings) {
            output.writeString(mapping);
        }
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
        sb.append(", saTAGForGoodMappingToNonCanonicalChromosome='").append(saTAGForGoodMappingToNonCanonicalChromosome).append('\'');
        sb.append('}');
        return sb.toString();
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final AssemblyContigWithFineTunedAlignments that = (AssemblyContigWithFineTunedAlignments) o;

        if (!sourceTig.equals(that.sourceTig)) return false;
        if (!insertionMappings.equals(that.insertionMappings)) return false;
        return saTAGForGoodMappingToNonCanonicalChromosome.equals(that.saTAGForGoodMappingToNonCanonicalChromosome);
    }

    @Override
    public int hashCode() {
        int result = sourceTig.hashCode();
        result = 31 * result + insertionMappings.hashCode();
        result = 31 * result + saTAGForGoodMappingToNonCanonicalChromosome.hashCode();
        return result;
    }
}
