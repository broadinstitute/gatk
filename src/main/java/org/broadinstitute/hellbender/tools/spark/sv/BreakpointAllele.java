package org.broadinstitute.hellbender.tools.spark.sv;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

import static org.broadinstitute.hellbender.tools.spark.sv.BreakpointAllele.InversionType.INV_3_TO_5;
import static org.broadinstitute.hellbender.tools.spark.sv.BreakpointAllele.InversionType.INV_5_TO_3;

/**
 * This class represents the allele of an SV breakpoint (a novel adjacency between two genomic locations)
 */
@DefaultSerializer(BreakpointAllele.Serializer.class)
class BreakpointAllele {

    final SimpleInterval leftAlignedLeftBreakpoint;
    final SimpleInterval leftAlignedRightBreakpoint;
    final String insertedSequence;
    final String homology;
    private final InversionType inversionType;

    public BreakpointAllele(final SimpleInterval leftAlignedLeftBreakpoint, final SimpleInterval leftAlignedRightBreakpoint, final String insertedSequence, final String homology, final boolean fiveToThree, final boolean threeToFive, final List<String> insertionMappings) {
        this.leftAlignedLeftBreakpoint = leftAlignedLeftBreakpoint;
        this.leftAlignedRightBreakpoint = leftAlignedRightBreakpoint;
        this.insertedSequence = insertedSequence;
        this.homology = homology;
        this.inversionType = getInversionType(fiveToThree, threeToFive);
    }

    @SuppressWarnings("unchecked")
    public BreakpointAllele(final Kryo kryo, final Input input) {
        final String contig1 = input.readString();
        final int start1 = input.readInt();
        final int end1 = input.readInt();
        this.leftAlignedLeftBreakpoint = new SimpleInterval(contig1, start1, end1);
        final String contig2 = input.readString();
        final int start2 = input.readInt();
        final int end2 = input.readInt();
        this.leftAlignedRightBreakpoint = new SimpleInterval(contig2, start2, end2);
        this.insertedSequence = input.readString();
        this.homology = input.readString();
        this.inversionType = InversionType.values()[input.readInt()];
    }

    public boolean isInversion() {
        return leftAlignedLeftBreakpoint.getContig().equals(leftAlignedRightBreakpoint.getContig()) && (getInversionType() == INV_3_TO_5 || getInversionType() == INV_5_TO_3);
    }

    public enum InversionType  {
        INV_3_TO_5, INV_5_TO_3, INV_NONE
    }

    public InversionType getInversionType() {
        return inversionType;
    }

    private InversionType getInversionType(final boolean fiveToThree, final boolean threeToFive){
        if(!fiveToThree && threeToFive){
            return InversionType.INV_3_TO_5;
        }else if(fiveToThree && !threeToFive){
            return InversionType.INV_5_TO_3;
        }
        return InversionType.INV_NONE;
    }


    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        final BreakpointAllele that = (BreakpointAllele) o;

        if (leftAlignedLeftBreakpoint != null ? !leftAlignedLeftBreakpoint.equals(that.leftAlignedLeftBreakpoint) : that.leftAlignedLeftBreakpoint != null)
            return false;
        if (leftAlignedRightBreakpoint != null ? !leftAlignedRightBreakpoint.equals(that.leftAlignedRightBreakpoint) : that.leftAlignedRightBreakpoint != null)
            return false;
        if (insertedSequence != null ? !insertedSequence.equals(that.insertedSequence) : that.insertedSequence != null)
            return false;
        if (homology != null ? !homology.equals(that.homology) : that.homology != null) return false;
        return inversionType == that.inversionType;

    }

    @Override
    public int hashCode() {
        return Objects.hash(leftAlignedLeftBreakpoint, leftAlignedRightBreakpoint, insertedSequence, homology, 2659*inversionType.ordinal());
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

    private void serialize(final Kryo kryo, final Output output) {
        output.writeString(leftAlignedLeftBreakpoint.getContig());
        output.writeInt(leftAlignedLeftBreakpoint.getStart());
        output.writeInt(leftAlignedLeftBreakpoint.getEnd());
        output.writeString(leftAlignedRightBreakpoint.getContig());
        output.writeInt(leftAlignedRightBreakpoint.getStart());
        output.writeInt(leftAlignedRightBreakpoint.getEnd());
        output.writeString(insertedSequence);
        output.writeString(homology);
        output.writeInt(inversionType.ordinal());
    }
}
