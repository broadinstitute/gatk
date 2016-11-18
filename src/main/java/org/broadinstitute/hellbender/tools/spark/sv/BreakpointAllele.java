package org.broadinstitute.hellbender.tools.spark.sv;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.common.annotations.VisibleForTesting;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import scala.Tuple2;

import java.util.Objects;

import static org.broadinstitute.hellbender.tools.spark.sv.BreakpointAllele.Strandedness.*;

/**
 * This class represents the allele of an SV breakpoint (a novel adjacency between two genomic locations)
 */
@DefaultSerializer(BreakpointAllele.Serializer.class)
class BreakpointAllele {

    final SimpleInterval leftJustifiedLeftBreakpoint;
    final SimpleInterval leftJustifiedRightBreakpoint;
    final boolean leftBreakpointAlignedForward;
    final boolean rightBreakpointAlignedForward;
    final String insertedSequence;
    final String homology;

    // TODO: 11/3/16 see where {@code insertionMappings} could be used
//    final List<String> insertionMappings;

    public BreakpointAllele(final ChimericAlignment chimericAlignment) {

        final Tuple2<SimpleInterval, SimpleInterval> leftJustifiedBreakpoints = chimericAlignment.getLeftJustifiedBreakpoints();
        this.leftJustifiedLeftBreakpoint = leftJustifiedBreakpoints._1();
        this.leftJustifiedRightBreakpoint = leftJustifiedBreakpoints._2();

        this.leftBreakpointAlignedForward = chimericAlignment.regionWithLowerCoordOnContig.forwardStrand;
        this.rightBreakpointAlignedForward = chimericAlignment.regionWithHigherCoordOnContig.forwardStrand;

        this.insertedSequence = chimericAlignment.insertedSequence;
        this.homology = chimericAlignment.homology;
//        this.insertionMappings = chimericAlignment.insertionMappings;
    }

    @SuppressWarnings("unchecked")
    protected BreakpointAllele(final Kryo kryo, final Input input) {
        final String contig1 = input.readString();
        final int start1 = input.readInt();
        final int end1 = input.readInt();
        this.leftJustifiedLeftBreakpoint = new SimpleInterval(contig1, start1, end1);
        final String contig2 = input.readString();
        final int start2 = input.readInt();
        final int end2 = input.readInt();
        this.leftJustifiedRightBreakpoint = new SimpleInterval(contig2, start2, end2);
        this.leftBreakpointAlignedForward = input.readBoolean();
        this.rightBreakpointAlignedForward = input.readBoolean();
        this.insertedSequence = input.readString();
        this.homology = input.readString();
    }

    Strandedness determineStrandedness() {
        if (leftBreakpointAlignedForward == rightBreakpointAlignedForward) {
            return SAME_STRAND;
        } else {
            return leftBreakpointAlignedForward ? FIVE_TO_THREE : THREE_TO_FIVE;
        }
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        final BreakpointAllele that = (BreakpointAllele) o;

        if (leftJustifiedLeftBreakpoint != null ? !leftJustifiedLeftBreakpoint.equals(that.leftJustifiedLeftBreakpoint) : that.leftJustifiedLeftBreakpoint != null)
            return false;
        if (leftJustifiedRightBreakpoint != null ? !leftJustifiedRightBreakpoint.equals(that.leftJustifiedRightBreakpoint) : that.leftJustifiedRightBreakpoint != null)
            return false;

        if (leftBreakpointAlignedForward != that.leftBreakpointAlignedForward)
            return false;

        if (rightBreakpointAlignedForward != that.rightBreakpointAlignedForward)
            return false;

        if (insertedSequence != null ? !insertedSequence.equals(that.insertedSequence) : that.insertedSequence != null)
            return false;

        if (homology==null) {
            return that.homology==null;
        } else {
            return that.homology!=null && homology.equals(that.homology);
        }
    }

    @Override
    public int hashCode() {
        return Objects.hash(leftJustifiedLeftBreakpoint, leftJustifiedRightBreakpoint, leftBreakpointAlignedForward, rightBreakpointAlignedForward, insertedSequence, homology);
    }

    protected void serialize(final Kryo kryo, final Output output) {
        output.writeString(leftJustifiedLeftBreakpoint.getContig());
        output.writeInt(leftJustifiedLeftBreakpoint.getStart());
        output.writeInt(leftJustifiedLeftBreakpoint.getEnd());
        output.writeString(leftJustifiedRightBreakpoint.getContig());
        output.writeInt(leftJustifiedRightBreakpoint.getStart());
        output.writeInt(leftJustifiedRightBreakpoint.getEnd());
        output.writeBoolean(leftBreakpointAlignedForward);
        output.writeBoolean(rightBreakpointAlignedForward);
        output.writeString(insertedSequence);
        output.writeString(homology);
    }

    enum Strandedness {
        SAME_STRAND, THREE_TO_FIVE, FIVE_TO_THREE
    }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<BreakpointAllele> {
        @Override
        public void write(final Kryo kryo, final Output output, final BreakpointAllele breakpointAllele ) {
            breakpointAllele.serialize(kryo, output);
        }

        @Override
        public BreakpointAllele read(final Kryo kryo, final Input input, final Class<BreakpointAllele> klass ) {
            return new BreakpointAllele(kryo, input);
        }
    }
}
