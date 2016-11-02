package org.broadinstitute.hellbender.tools.spark.sv;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import scala.Tuple2;

import java.util.List;
import java.util.Objects;

import static org.broadinstitute.hellbender.tools.spark.sv.BreakpointAllele.BreakpointAlleleInversion.InversionType.*;

/**
 * This class represents the allele of an SV breakpoint (a novel adjacency between two genomic locations)
 */
class BreakpointAllele {

    final SimpleInterval leftAlignedLeftBreakpoint;
    final SimpleInterval leftAlignedRightBreakpoint;
    final String insertedSequence;
    final String homology;

    // TODO: 11/3/16 see where {@code insertionMappings} could be used
    protected BreakpointAllele(final SimpleInterval leftAlignedLeftBreakpoint, final SimpleInterval leftAlignedRightBreakpoint,
                               final String homology, final String insertedSequence, final List<String> insertionMappings) {
        this.leftAlignedLeftBreakpoint = leftAlignedLeftBreakpoint;
        this.leftAlignedRightBreakpoint = leftAlignedRightBreakpoint;
        this.insertedSequence = insertedSequence;
        this.homology = homology;
    }

    protected BreakpointAllele(final ChimericAlignment chimericAlignment) {

        final Tuple2<SimpleInterval, SimpleInterval> leftAndRightBreakpointsOnReferenceLeftAlignedForHomology = chimericAlignment.getLeftAndRightBreakpointsOnReferenceLeftAlignedForHomology();
        final SimpleInterval leftAlignedLeftBreakpointOnAssembledContig = leftAndRightBreakpointsOnReferenceLeftAlignedForHomology._1();
        final SimpleInterval leftAlignedRightBreakpointOnAssembledContig = leftAndRightBreakpointsOnReferenceLeftAlignedForHomology._2();

        if (!leftAlignedLeftBreakpointOnAssembledContig.getContig().equals(leftAlignedRightBreakpointOnAssembledContig.getContig())) {
            this.leftAlignedLeftBreakpoint = leftAlignedLeftBreakpointOnAssembledContig;
            this.leftAlignedRightBreakpoint = leftAlignedRightBreakpointOnAssembledContig;
        } else if (leftAlignedLeftBreakpointOnAssembledContig.getStart() < leftAlignedRightBreakpointOnAssembledContig.getStart()) {
            this.leftAlignedLeftBreakpoint = leftAlignedLeftBreakpointOnAssembledContig;
            this.leftAlignedRightBreakpoint = leftAlignedRightBreakpointOnAssembledContig;
        } else {
            this.leftAlignedLeftBreakpoint = leftAlignedRightBreakpointOnAssembledContig;
            this.leftAlignedRightBreakpoint = leftAlignedLeftBreakpointOnAssembledContig;
        }

        this.insertedSequence = chimericAlignment.insertedSequence;
        this.homology = chimericAlignment.homology;
    }

    @SuppressWarnings("unchecked")
    protected BreakpointAllele(final Kryo kryo, final Input input) {
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

        if (homology==null) {
            return that.homology==null;
        } else {
            return that.homology!=null && homology.equals(that.homology);
        }
    }

    @Override
    public int hashCode() {
        return Objects.hash(leftAlignedLeftBreakpoint, leftAlignedRightBreakpoint, insertedSequence, homology);
    }

    protected void serialize(final Kryo kryo, final Output output) {
        output.writeString(leftAlignedLeftBreakpoint.getContig());
        output.writeInt(leftAlignedLeftBreakpoint.getStart());
        output.writeInt(leftAlignedLeftBreakpoint.getEnd());
        output.writeString(leftAlignedRightBreakpoint.getContig());
        output.writeInt(leftAlignedRightBreakpoint.getStart());
        output.writeInt(leftAlignedRightBreakpoint.getEnd());
        output.writeString(insertedSequence);
        output.writeString(homology);
    }

    @DefaultSerializer(BreakpointAlleleInversion.Serializer.class)
    static final class BreakpointAlleleInversion extends BreakpointAllele {

        final InversionType inversionType;

        public BreakpointAlleleInversion(final ChimericAlignment chimericAlignment){
            super(chimericAlignment);

            final boolean isFiveToThreeInversion = chimericAlignment.region1.forwardStrand && !chimericAlignment.region2.forwardStrand;
            final boolean isThreeToFiveInversion = !chimericAlignment.region1.forwardStrand && chimericAlignment.region2.forwardStrand;
            this.inversionType = determineInversionType(isFiveToThreeInversion, isThreeToFiveInversion);
        }

        public BreakpointAlleleInversion(final SimpleInterval leftAlignedLeftBreakpoint, final SimpleInterval leftAlignedRightBreakpoint,
                                         final String homology, final String insertedSequence, final List<String> insertionMappings,
                                         final boolean fiveToThree, final boolean threeToFive) {

            super(leftAlignedLeftBreakpoint, leftAlignedRightBreakpoint, homology, insertedSequence, insertionMappings);

            this.inversionType = determineInversionType(fiveToThree, threeToFive);
        }

        @SuppressWarnings("unchecked")
        protected BreakpointAlleleInversion(final Kryo kryo, final Input input) {
            super(kryo, input);
            this.inversionType = BreakpointAlleleInversion.InversionType.values()[input.readInt()];
        }

        @Override
        protected void serialize(final Kryo kryo, final Output output) {
            super.serialize(kryo, output);
            output.writeInt(inversionType.ordinal());
        }

        public static final class Serializer extends com.esotericsoftware.kryo.Serializer<BreakpointAlleleInversion> {
            @Override
            public void write(final Kryo kryo, final Output output, final BreakpointAlleleInversion breakpointAlleleInversion ) {
                breakpointAlleleInversion.serialize(kryo, output);
            }

            @Override
            public BreakpointAlleleInversion read(final Kryo kryo, final Input input, final Class<BreakpointAlleleInversion> klass ) {
                return new BreakpointAlleleInversion(kryo, input);
            }
        }

        @Override
        public boolean equals(final Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            final BreakpointAlleleInversion that = (BreakpointAlleleInversion) o;

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

        enum InversionType  {
            INV_3_TO_5, INV_5_TO_3, INV_NONE
        }

        InversionType getInversionType() {
            return inversionType;
        }

        // TODO: test
        private InversionType determineInversionType(final boolean fiveToThree, final boolean threeToFive){
//            // TODO: the logic could be simplified to this
//            //       it should not be possible for a single assembled contig to be both 5-3 and 3-5, or we must tweak the algo to merge two classes of contigs pointing to the supposedly-same inversion
//            if (fiveToThree || threeToFive) {
//                return fiveToThree ? INV_5_TO_3 : INV_3_TO_5;
//            } else {
//                return INV_NONE;
//            }
            if(!fiveToThree && threeToFive){
                return INV_3_TO_5;
            }else if(fiveToThree && !threeToFive){
                return INV_5_TO_3;
            }
            return INV_NONE;
        }
    }
}
