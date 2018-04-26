package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.BreakEndVariantType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SvType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.StrandSwitch;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import scala.Tuple2;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * This class represents a pair of inferred genomic locations on the reference whose novel adjacency is generated
 * due to an SV event (in other words, a simple rearrangement between two genomic locations)
 * that is suggested by the input {@link SimpleChimera},
 * and complications as enclosed in {@link BreakpointComplications}
 * in pinning down the locations to exact base pair resolution.
 *
 * <p>
 *     It essentially represents a bi-path "big" bubble between two reference locations.
 *     One path is the "reference path" consisting of the contiguous block of bases that can be extracted from the reference,
 *     if possible (i.e. no contiguous block exists between locations from difference chromosomes).
 *     The other path is encoded with the alt haplotype sequence.
 * </p>
 */
@DefaultSerializer(NovelAdjacencyAndAltHaplotype.Serializer.class)
public class NovelAdjacencyAndAltHaplotype {

    private final SimpleInterval leftJustifiedLeftRefLoc;
    private final SimpleInterval leftJustifiedRightRefLoc;

    private final StrandSwitch strandSwitch;
    private final BreakpointComplications complication;

    private final TypeInferredFromSimpleChimera type;
    private final byte[] altHaplotypeSequence;

    public NovelAdjacencyAndAltHaplotype(final SimpleChimera simpleChimera, final byte[] contigSequence,
                                         final SAMSequenceDictionary referenceDictionary) {

        strandSwitch = simpleChimera.strandSwitch;

        try {

            final BreakpointsInference inferredClass =
                    BreakpointsInference.getInferenceClass(simpleChimera, contigSequence, referenceDictionary);

            final Tuple2<SimpleInterval, SimpleInterval> leftJustifiedBreakpoints = inferredClass.getLeftJustifiedBreakpoints();
            leftJustifiedLeftRefLoc = leftJustifiedBreakpoints._1();
            leftJustifiedRightRefLoc = leftJustifiedBreakpoints._2();

            complication = inferredClass.getComplications();

            type = simpleChimera.inferType();

            altHaplotypeSequence = inferredClass.getInferredAltHaplotypeSequence();

        } catch (final IllegalArgumentException iaex) { // catching IAEX specifically because it is the most likely exception thrown if there's bug, this helps quickly debugging what the problem is
            throw new GATKException("Erred when inferring breakpoint location and event type from chimeric alignment:\n" +
                    simpleChimera.toString(), iaex);
        }
    }

    public boolean hasInsertedSequence() {
        return ! complication.getInsertedSequenceForwardStrandRep().isEmpty();
    }

    public boolean hasDuplicationAnnotation() {
        return complication.indicatesRefSeqDuplicatedOnAlt();
    }

    protected NovelAdjacencyAndAltHaplotype(final Kryo kryo, final Input input) {
        final String contig1 = input.readString();
        final int start1 = input.readInt();
        final int end1 = input.readInt();
        this.leftJustifiedLeftRefLoc = new SimpleInterval(contig1, start1, end1);
        final String contig2 = input.readString();
        final int start2 = input.readInt();
        final int end2 = input.readInt();
        this.leftJustifiedRightRefLoc = new SimpleInterval(contig2, start2, end2);

        this.strandSwitch = StrandSwitch.values()[input.readInt()];

        this.complication = (BreakpointComplications) kryo.readClassAndObject(input);

        this.type = TypeInferredFromSimpleChimera.values()[ input.readInt() ];

        final boolean altSeqIsNull = input.readBoolean();
        if (altSeqIsNull) {
            altHaplotypeSequence = null;
        } else {
            final int arraySize = input.readInt();
            altHaplotypeSequence = new byte[arraySize];
            for (int i = 0 ; i < arraySize; ++i) {
                altHaplotypeSequence[i] = input.readByte();
            }
        }
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        NovelAdjacencyAndAltHaplotype that = (NovelAdjacencyAndAltHaplotype) o;

        if (!leftJustifiedLeftRefLoc.equals(that.leftJustifiedLeftRefLoc)) return false;
        if (!leftJustifiedRightRefLoc.equals(that.leftJustifiedRightRefLoc)) return false;
        if (strandSwitch != that.strandSwitch) return false;
        if (!complication.equals(that.complication)) return false;
        return Arrays.equals(altHaplotypeSequence, that.altHaplotypeSequence);
    }

    @Override
    public int hashCode() {
        int result = leftJustifiedLeftRefLoc.hashCode();
        result = 31 * result + leftJustifiedRightRefLoc.hashCode();
        result = 31 * result + strandSwitch.ordinal();
        result = 31 * result + complication.hashCode();
        result = 31 * result + Arrays.hashCode(altHaplotypeSequence);
        return result;
    }

    protected void serialize(final Kryo kryo, final Output output) {
        output.writeString(leftJustifiedLeftRefLoc.getContig());
        output.writeInt(leftJustifiedLeftRefLoc.getStart());
        output.writeInt(leftJustifiedLeftRefLoc.getEnd());
        output.writeString(leftJustifiedRightRefLoc.getContig());
        output.writeInt(leftJustifiedRightRefLoc.getStart());
        output.writeInt(leftJustifiedRightRefLoc.getEnd());
        output.writeInt(strandSwitch.ordinal());

        kryo.writeClassAndObject(output, complication);

        output.writeInt(type.ordinal());

        if (altHaplotypeSequence==null) {
            output.writeBoolean(true);
        } else {
            output.writeBoolean(false);
            output.writeInt(altHaplotypeSequence.length);
            for (final byte b : altHaplotypeSequence) {
                output.writeByte(b);
            }
        }
    }

    public SimpleInterval getLeftJustifiedLeftRefLoc() {
        return leftJustifiedLeftRefLoc;
    }

    public SimpleInterval getLeftJustifiedRightRefLoc() {
        return leftJustifiedRightRefLoc;
    }

    public StrandSwitch getStrandSwitch() {
        return strandSwitch;
    }

    public BreakpointComplications getComplication() {
        return complication;
    }

    TypeInferredFromSimpleChimera getType() {
        return type;
    }

    public byte[] getAltHaplotypeSequence() {
        return altHaplotypeSequence;
    }

    /**
     * @return the inferred type could be a single entry for simple variants, or a list of two entries with BND mates.
     */
    List<SvType> toSimpleOrBNDTypes(final ReferenceMultiSource reference, final SAMSequenceDictionary referenceDictionary) {

        switch (type) {
            case InterChromosome:
            case IntraChrRefOrderSwap:
                final Tuple2<BreakEndVariantType, BreakEndVariantType> orderedMatesForTranslocSuspect =
                        BreakEndVariantType.TransLocBND.getOrderedMates(this, reference, referenceDictionary);
                return Arrays.asList(orderedMatesForTranslocSuspect._1, orderedMatesForTranslocSuspect._2);
            case IntraChrStrandSwitch:
                if ( complication.indicatesRefSeqDuplicatedOnAlt() ) {
                    final int svLength =
                            ((BreakpointComplications.IntraChrStrandSwitchBreakpointComplications) this.getComplication())
                                    .getDupSeqRepeatUnitRefSpan().size();
                    return Collections.singletonList( new SimpleSVType.DuplicationInverted(this, svLength) );
                } else {
                    final Tuple2<BreakEndVariantType, BreakEndVariantType> orderedMatesForInversionSuspect =
                            BreakEndVariantType.InvSuspectBND.getOrderedMates(this, reference);
                    return Arrays.asList(orderedMatesForInversionSuspect._1, orderedMatesForInversionSuspect._2);
                }
            case SIMPLE_DEL:
            {
                final int svLength = leftJustifiedLeftRefLoc.getEnd() - leftJustifiedRightRefLoc.getStart();
                return Collections.singletonList( new SimpleSVType.Deletion(this, svLength) );
            }
            case RPL:
            {
                final int svLength = leftJustifiedLeftRefLoc.getEnd() - leftJustifiedRightRefLoc.getStart();
                return Collections.singletonList( new SimpleSVType.Deletion(this, svLength) );
            }
            case SIMPLE_INS:
            {
                final int svLength = this.getComplication().getInsertedSequenceForwardStrandRep().length();
                return Collections.singletonList( new SimpleSVType.Insertion(this, svLength) );
            }
            case SMALL_DUP_EXPANSION:
            {
                final int svLength = getLengthForDupTandem(this);
                return Collections.singletonList( new SimpleSVType.DuplicationTandem(this, svLength) );
            }
            case DEL_DUP_CONTRACTION:
            {
                final int svLength = leftJustifiedLeftRefLoc.getEnd() - leftJustifiedRightRefLoc.getStart();
                return Collections.singletonList( new SimpleSVType.Deletion(this, svLength) );
            }
            case SMALL_DUP_CPX:
            {
                if ( ((BreakpointComplications.SmallDuplicationWithImpreciseDupRangeBreakpointComplications)
                        this.getComplication()).isDupContraction() ) {
                    final int svLength = leftJustifiedLeftRefLoc.getEnd() - leftJustifiedRightRefLoc.getStart();
                    return Collections.singletonList( new SimpleSVType.Deletion(this, svLength) );
                } else {
                    final int svLength = getLengthForDupTandem(this);
                    return Collections.singletonList( new SimpleSVType.DuplicationTandem(this, svLength) );
                }
            }
            default:
                throw new GATKException.ShouldNeverReachHereException("Inferred type not recognized");
        }
    }

    /**
     * Clean expansion of repeat 1 -> 2, or
     * complex expansion, or
     * expansion of 1 repeat on ref to 2 repeats on alt with inserted sequence in between the 2 repeats
     */
    public static int getLengthForDupTandem(final NovelAdjacencyAndAltHaplotype novelAdjacencyAndAltHaplotype) {
        final BreakpointComplications.SmallDuplicationBreakpointComplications complication = (BreakpointComplications.SmallDuplicationBreakpointComplications)
                novelAdjacencyAndAltHaplotype.getComplication();
        return complication.getInsertedSequenceForwardStrandRep().length()
                + (complication.getDupSeqRepeatNumOnCtg() - complication.getDupSeqRepeatNumOnRef()) * complication.getDupSeqRepeatUnitRefSpan().size();
    }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<NovelAdjacencyAndAltHaplotype> {
        @Override
        public void write(final Kryo kryo, final Output output, final NovelAdjacencyAndAltHaplotype novelAdjacencyAndAltHaplotype) {
            novelAdjacencyAndAltHaplotype.serialize(kryo, output);
        }

        @Override
        public NovelAdjacencyAndAltHaplotype read(final Kryo kryo, final Input input, final Class<NovelAdjacencyAndAltHaplotype> klass ) {
            return new NovelAdjacencyAndAltHaplotype(kryo, input);
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
