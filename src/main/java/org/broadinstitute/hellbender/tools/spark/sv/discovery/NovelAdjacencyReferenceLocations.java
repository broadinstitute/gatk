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

    public NovelAdjacencyReferenceLocations(final ChimericAlignment chimericAlignment, final byte[] contigSequence){

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

        final String leftBreakpointRefContig, rightBreakpointRefContig;
        final int leftBreakpointCoord, rightBreakpointCoord;
        if (complication.hasDuplicationAnnotation()) { // todo : development artifact-- assuming tandem duplication is not co-existing with inversion
            leftBreakpointRefContig = rightBreakpointRefContig = ca.regionWithLowerCoordOnContig.referenceSpan.getContig();
            final SimpleInterval leftReferenceInterval, rightReferenceInterval;
            if (ca.isForwardStrandRepresentation) {
                leftReferenceInterval = ca.regionWithLowerCoordOnContig.referenceSpan;
                rightReferenceInterval = ca.regionWithHigherCoordOnContig.referenceSpan;
            } else {
                leftReferenceInterval = ca.regionWithHigherCoordOnContig.referenceSpan;
                rightReferenceInterval = ca.regionWithLowerCoordOnContig.referenceSpan;
            }
            if (complication.getDupSeqRepeatNumOnCtg() > complication.getDupSeqRepeatNumOnRef()) {
                leftBreakpointCoord = leftReferenceInterval.getEnd() - homologyLen
                        - (complication.getDupSeqRepeatNumOnCtg() - complication.getDupSeqRepeatNumOnRef()) * complication.getDupSeqRepeatUnitRefSpan().size();
            } else {
                leftBreakpointCoord = leftReferenceInterval.getEnd() - homologyLen;
            }
            rightBreakpointCoord = rightReferenceInterval.getStart() - 1;
        } else { // inversion and simple deletion & insertion
            final SimpleInterval leftReferenceInterval, rightReferenceInterval;
            if (ca.isForwardStrandRepresentation) {
                leftReferenceInterval  = ca.regionWithLowerCoordOnContig.referenceSpan;
                rightReferenceInterval = ca.regionWithHigherCoordOnContig.referenceSpan;
            } else {
                leftReferenceInterval  = ca.regionWithHigherCoordOnContig.referenceSpan;
                rightReferenceInterval = ca.regionWithLowerCoordOnContig.referenceSpan;
            }
            leftBreakpointRefContig  = leftReferenceInterval.getContig();
            rightBreakpointRefContig = rightReferenceInterval.getContig();
            if (ca.strandSwitch == StrandSwitch.NO_SWITCH) {
                leftBreakpointCoord  = leftReferenceInterval.getEnd() - homologyLen;
                rightBreakpointCoord = rightReferenceInterval.getStart() - 1;
            } else if (ca.strandSwitch == StrandSwitch.FORWARD_TO_REVERSE){
                leftBreakpointCoord  = leftReferenceInterval.getEnd() - homologyLen;
                rightBreakpointCoord = rightReferenceInterval.getEnd();
            } else {
                leftBreakpointCoord  = leftReferenceInterval.getStart() - 1;
                rightBreakpointCoord = rightReferenceInterval.getStart() + homologyLen - 1;
            }
        }

        Utils.validate(leftBreakpointCoord <= rightBreakpointCoord,
                "Inferred novel adjacency reference locations have left location after right location : " +
                        leftBreakpointCoord + "\t" + rightBreakpointCoord + ca.onErrStringRep() + "\n" + complication.toString());

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
