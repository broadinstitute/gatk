package org.broadinstitute.hellbender.tools.spark.sv.discovery.prototype;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.Serializer;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignmentInterval;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

@DefaultSerializer(Serializer.class)
final class AssemblyContigWithFineTunedAlignments {
    private static final AlignedContig.Serializer contigSerializer = new AlignedContig.Serializer();
    public static final List<String> emptyInsertionMappings = Collections.emptyList();

    final AlignedContig contig;
    final List<String> insertionMappings;

    AssemblyContigWithFineTunedAlignments(final AlignedContig contig) {
        this(contig, emptyInsertionMappings);
    }

    AssemblyContigWithFineTunedAlignments(final AlignedContig contig,
                                          final List<String> insertionMappings) {
        this.contig = contig;
        this.insertionMappings = Utils.nonNull(insertionMappings);
    }

    AssemblyContigWithFineTunedAlignments(final Kryo kryo, final Input input) {
        contig = contigSerializer.read(kryo, input, AlignedContig.class);
        final int insertionMappingsSize = input.readInt();
        insertionMappings = new ArrayList<>(insertionMappingsSize);
        if (insertionMappingsSize != 0) {
            for (int i = 0; i < insertionMappingsSize; ++i) {
                insertionMappings.add(input.readString());
            }
        }
    }

    // TODO: 11/21/17 push these predicates to AlignedContig class in a later PR
    // after fine tuning, a contig may have no good alignment left, or only 1
    final boolean isInformative() {
        return contig.alignmentIntervals.size() > 1;
    }

    final boolean hasOnly2GoodAlignments () {
        return contig.alignmentIntervals.size() == 2;
    }

    final boolean hasIncomePicture() {
        if ( hasOnly2GoodAlignments() )
            return hasIncompletePictureFromTwoAlignments();
        else
            return hasIncompletePictureFromMultipleAlignments();
    }

    //==================================================================================================================

    boolean hasIncompletePictureFromTwoAlignments() {
        return hasIncompletePictureDueToRefSpanContainment()
                ||
                (firstAndLastAlignmentMappedToSameChr() && hasIncompletePictureDueToOverlappingRefOrderSwitch());
    }

    boolean hasIncompletePictureDueToRefSpanContainment() {

        final AlignmentInterval one = contig.alignmentIntervals.get(0);
        final AlignmentInterval two = contig.alignmentIntervals.get(1);
        return one.referenceSpan.contains(two.referenceSpan)
                ||
                two.referenceSpan.contains(one.referenceSpan);
    }

    boolean hasIncompletePictureDueToOverlappingRefOrderSwitch() {

        final AlignmentInterval one = contig.alignmentIntervals.get(0);
        final AlignmentInterval two = contig.alignmentIntervals.get(1);
        final SimpleInterval referenceSpanOne = one.referenceSpan;
        final SimpleInterval referenceSpanTwo = two.referenceSpan;

        if (referenceSpanOne.contains(referenceSpanTwo) || referenceSpanTwo.contains(referenceSpanOne))
            return true;

        // TODO: 10/29/17 this obsoletes the inverted duplication call code we have now, but those could be used to figure out how to annotate which known ref regions are invert duplicated
        if (one.forwardStrand != two.forwardStrand) {
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

    boolean firstAndLastAlignmentMappedToSameChr() {
        Utils.validateArg(!hasIncompletePictureDueToRefSpanContainment(),
                "assumption that input contig has alignments whose ref span is not completely covered by the other's is violated \n" +
                        contig.toString());

        final String firstMappedChr = this.contig.alignmentIntervals.get(0).referenceSpan.getContig();
        final String lastMappedChr  = this.contig.alignmentIntervals.get(this.contig.alignmentIntervals.size() - 1).referenceSpan.getContig();

        return firstMappedChr.equals(lastMappedChr);
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
        final int goodAlignmentsCount = contig.alignmentIntervals.size();

        final AlignmentInterval head = contig.alignmentIntervals.get(0),
                                tail = contig.alignmentIntervals.get(goodAlignmentsCount-1);

        // mapped to different chr or strand
        if ( !head.referenceSpan.getContig().equals(tail.referenceSpan.getContig()) || head.forwardStrand != tail.forwardStrand)
            return true;

        final SimpleInterval referenceSpanHead = head.referenceSpan,
                             referenceSpanTail = tail.referenceSpan;

        // head or tail's ref span contained in the other's
        if (referenceSpanHead.contains(referenceSpanTail) || referenceSpanTail.contains(referenceSpanHead))
            return true;

        // middle alignments' ref span should be either 1) disjoint from valid region, or 2) completely contained in valid region
        final SimpleInterval validRegion = new SimpleInterval(referenceSpanHead.getContig(),
                                                              Math.min(referenceSpanHead.getStart(), referenceSpanTail.getStart()),
                                                              Math.max(referenceSpanHead.getEnd(), referenceSpanTail.getEnd()));
        final boolean notCompleteDupRegion =
                contig.alignmentIntervals
                        .subList(1, goodAlignmentsCount - 1) // take middle (no head no tail) alignments' ref span
                        .stream().map(ai -> ai.referenceSpan)
                        .anyMatch( middle -> middle.overlaps(validRegion) && !validRegion.contains(middle));
        if (notCompleteDupRegion)
            return true;

        // reference order switch, indicating a complicated tandem duplication breakpoint
        if (head.forwardStrand) {
            return referenceSpanHead.getStart() >= referenceSpanTail.getStart();
        } else {
            return referenceSpanHead.getEnd() <= referenceSpanHead.getEnd();
        }
    }

    void serialize(final Kryo kryo, final Output output) {
        contigSerializer.write(kryo, output, contig);
        output.writeInt(insertionMappings.size());
        for (final String mapping : insertionMappings) {
            output.writeString(mapping);
        }
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
}
