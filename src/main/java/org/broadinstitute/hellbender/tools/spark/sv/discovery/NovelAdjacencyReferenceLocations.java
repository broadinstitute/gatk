package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.common.annotations.VisibleForTesting;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;

import java.util.Objects;

/**
 * This class represents a pair of inferred genomic locations on the reference whose novel adjacency is generated
 * due to a simple SV event (in other words, a simple rearrangement between two genomic locations)
 * that is suggested by the input {@link AlignedContig},
 * and complications in pinning down the locations to exact base pair resolution.
 */
@DefaultSerializer(NovelAdjacencyReferenceLocations.Serializer.class)
public class NovelAdjacencyReferenceLocations {

    public final SimpleInterval leftJustifiedLeftRefLoc;
    public final SimpleInterval leftJustifiedRightRefLoc;

    public final StrandSwitch strandSwitch;
    public final BreakpointComplications complication;

    public NovelAdjacencyReferenceLocations(final ChimericAlignment chimericAlignment, final byte[] contigSequence) {

        // first get strand switch type, then get complications, finally use complications to justify breakpoints
        strandSwitch = chimericAlignment.strandSwitch;

        complication = new BreakpointComplications(chimericAlignment, contigSequence);

        final Tuple2<SimpleInterval, SimpleInterval> leftJustifiedBreakpoints = leftJustifyBreakpoints(chimericAlignment, complication);
        leftJustifiedLeftRefLoc = leftJustifiedBreakpoints._1();
        leftJustifiedRightRefLoc = leftJustifiedBreakpoints._2();
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

        this.strandSwitch = StrandSwitch.values()[input.readInt()];
        this.complication = kryo.readObject(input, BreakpointComplications.class);
    }

    // TODO: 12/11/16 again this does not deal with inter-chromosome translocation yet
    /**
     * Returns the reference coordinates of the left and right breakpoints implied by this chimeric alignment.
     * If there is homologous sequence represented in the alignments, it will be assigned to the side of the breakpoint
     * with higher reference coordinates.
     */
    @VisibleForTesting
    static Tuple2<SimpleInterval, SimpleInterval> leftJustifyBreakpoints(final ChimericAlignment ca,
                                                                         final BreakpointComplications complication) {
        final int homologyLen = complication.getHomologyForwardStrandRep().length();

        final Tuple2<SimpleInterval, SimpleInterval> coordSortedReferenceSpans = ca.getCoordSortedReferenceSpans();
        final SimpleInterval leftReferenceSpan  = coordSortedReferenceSpans._1,
                             rightReferenceSpan = coordSortedReferenceSpans._2;

        final String leftBreakpointRefContig, rightBreakpointRefContig;
        final int leftBreakpointCoord, rightBreakpointCoord;
        if (complication.hasDuplicationAnnotation()) {

            leftBreakpointRefContig = rightBreakpointRefContig = leftReferenceSpan.getContig();

            if (complication.indicatesInvDup()) { // inverted duplication
                leftBreakpointCoord  = complication.getDupSeqRepeatUnitRefSpan().getStart() - 1;
                rightBreakpointCoord = complication.getDupSeqRepeatUnitRefSpan().getEnd();
            } else { // duplication not involving inverted repeats
                if (leftReferenceSpan.contains(rightReferenceSpan) || rightReferenceSpan.contains(leftReferenceSpan)) { // one ref span contains the other
                    leftBreakpointCoord  = complication.getDupSeqRepeatUnitRefSpan().getStart();
                    rightBreakpointCoord = complication.getDupSeqRepeatUnitRefSpan().getEnd();
                } else {
                    final int expansionDiffMayBeNegative = complication.getDupSeqRepeatNumOnCtg() - complication.getDupSeqRepeatNumOnRef();
                    if (expansionDiffMayBeNegative > 0) {
                        leftBreakpointCoord = leftReferenceSpan.getEnd() - homologyLen + 1
                                              - expansionDiffMayBeNegative * complication.getDupSeqRepeatUnitRefSpan().size();
                    } else {
                        // no shift by 1 because for deletion record, POS is base BEFORE the deleted sequence
                        leftBreakpointCoord = leftReferenceSpan.getEnd() - homologyLen;
                    }
                    // TODO: 9/26/17 if we decide to treat dup as insertion so POS==END, rightBreakpointCoord = rightReferenceInterval.getStart() - 1;
                    rightBreakpointCoord = complication.getDupSeqRepeatUnitRefSpan().getEnd();
                }
            }
        } else { // inversion and simple deletion & insertion
            leftBreakpointRefContig  = leftReferenceSpan.getContig();
            rightBreakpointRefContig = rightReferenceSpan.getContig();
            if (ca.strandSwitch == StrandSwitch.NO_SWITCH) { // simple ins/del
                leftBreakpointCoord  = leftReferenceSpan.getEnd() - homologyLen;
                rightBreakpointCoord = rightReferenceSpan.getStart() - 1;
            } else if (ca.strandSwitch == StrandSwitch.FORWARD_TO_REVERSE){ // detected inv by assembling left breakpoint
                leftBreakpointCoord  = leftReferenceSpan.getEnd() - homologyLen;
                rightBreakpointCoord = rightReferenceSpan.getEnd();
            } else {                                                        // detected inv by assembling right breakpoint
                leftBreakpointCoord  = leftReferenceSpan.getStart() - 1;
                rightBreakpointCoord = rightReferenceSpan.getStart() + homologyLen - 1;
            }
        }

        Utils.validate(leftBreakpointCoord <= rightBreakpointCoord,
                "Inferred novel adjacency reference locations have left location after right location : " +
                        leftBreakpointRefContig + ":" + leftBreakpointCoord + "-" + rightBreakpointCoord + "\n" +
                        ca.onErrStringRep() + "\t" + complication.toString());

        final SimpleInterval leftBreakpoint = new SimpleInterval(leftBreakpointRefContig, leftBreakpointCoord, leftBreakpointCoord);
        final SimpleInterval rightBreakpoint = new SimpleInterval(rightBreakpointRefContig, rightBreakpointCoord, rightBreakpointCoord);
        return new Tuple2<>(leftBreakpoint, rightBreakpoint);
    }


    @VisibleForTesting
    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        final NovelAdjacencyReferenceLocations that = (NovelAdjacencyReferenceLocations) o;

        if (leftJustifiedLeftRefLoc != null ? !leftJustifiedLeftRefLoc.equals(that.leftJustifiedLeftRefLoc)
                : that.leftJustifiedLeftRefLoc != null)
            return false;
        if (leftJustifiedRightRefLoc != null ? !leftJustifiedRightRefLoc.equals(that.leftJustifiedRightRefLoc)
                : that.leftJustifiedRightRefLoc != null)
            return false;

        return strandSwitch.equals(that.strandSwitch) && complication.equals(that.complication);
    }

    @VisibleForTesting
    @Override
    public int hashCode() {
        return Objects.hash(leftJustifiedLeftRefLoc, leftJustifiedRightRefLoc, complication, 2659*strandSwitch.ordinal());
    }

    protected void serialize(final Kryo kryo, final Output output) {
        output.writeString(leftJustifiedLeftRefLoc.getContig());
        output.writeInt(leftJustifiedLeftRefLoc.getStart());
        output.writeInt(leftJustifiedLeftRefLoc.getEnd());
        output.writeString(leftJustifiedRightRefLoc.getContig());
        output.writeInt(leftJustifiedRightRefLoc.getStart());
        output.writeInt(leftJustifiedRightRefLoc.getEnd());
        output.writeInt(strandSwitch.ordinal());
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

    /**
     * @return Intended for use in debugging and exception message only.
     */
    @Override
    public String toString() {
        return String.format("%s\t%s\t%s\t%s", leftJustifiedLeftRefLoc.toString(), leftJustifiedRightRefLoc.toString(),
                strandSwitch.name(), complication.toString());
    }
}
