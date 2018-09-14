package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceMultiSparkSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.BreakEndVariantType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SvType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.StrandSwitch;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import scala.Tuple2;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.STRUCTURAL_VARIANT_SIZE_LOWER_BOUND;

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

    @VisibleForTesting
    public NovelAdjacencyAndAltHaplotype(final SimpleInterval leftJustifiedLeftRefLoc, final SimpleInterval leftJustifiedRightRefLoc,
                                         final StrandSwitch strandSwitch, final BreakpointComplications complication,
                                         final TypeInferredFromSimpleChimera type, final byte[] altHaplotypeSequence) {
        this.leftJustifiedLeftRefLoc = leftJustifiedLeftRefLoc;
        this.leftJustifiedRightRefLoc = leftJustifiedRightRefLoc;
        this.strandSwitch = strandSwitch;
        this.complication = complication;
        this.type = type;
        this.altHaplotypeSequence = altHaplotypeSequence;
    }

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

            type = simpleChimera.inferType(referenceDictionary);

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
        return complication.hasDuplicationAnnotation();
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

    public TypeInferredFromSimpleChimera getTypeInferredFromSimpleChimera() {
        return type;
    }

    public byte[] getAltHaplotypeSequence() {
        return altHaplotypeSequence;
    }

    public int getDistanceBetweenNovelAdjacencies() {
        if (leftJustifiedLeftRefLoc.getContig().equals(leftJustifiedRightRefLoc.getContig())) {
            return leftJustifiedRightRefLoc.getEnd() - leftJustifiedLeftRefLoc.getStart();
        } else {
            return -1;
        }
    }

    /**
     * For replacement, one could have a case where the sequence being replaced (i.e. deleted sequence)
     * is not large enough to trigger an structural variant DEL call,
     * whereas we emit an INS call, and here we call it "fat" insertion.
     */
    public boolean isCandidateForFatInsertion() {
        return type.equals(TypeInferredFromSimpleChimera.RPL)
                &&
                leftJustifiedRightRefLoc.getEnd() - leftJustifiedLeftRefLoc.getStart() < STRUCTURAL_VARIANT_SIZE_LOWER_BOUND;
    }

    /**
     * Return the interval being replaced, anchored 1 bp from left,
     * i.e.
     * if [2,30] is being replaced, this would return [1,30]
     */
    public SimpleInterval getIntervalForFatInsertion() {
        if (isCandidateForFatInsertion()) {
            return new SimpleInterval(leftJustifiedLeftRefLoc.getContig(), leftJustifiedLeftRefLoc.getStart(), leftJustifiedRightRefLoc.getEnd());
        } else
            throw new UnsupportedOperationException("trying to get interval from a novel adjacency that is not a fat insertion: " + toString());
    }

    /**
     * @return the inferred type could be
     *          1) a single entry for simple variants, or
     *          2) a list of two entries for "replacement" case where the ref- and alt- path are both >=
     *             {@link StructuralVariationDiscoveryArgumentCollection#STRUCTURAL_VARIANT_SIZE_LOWER_BOUND} bp long, or
     *          3) a list of two entries with BND mates.
     *          It is safe, for now, to assume that ({@link SimpleSVType} and {@link BreakEndVariantType} are never mixed)
     */
    @VisibleForTesting
    public List<SvType> toSimpleOrBNDTypes( final ReferenceMultiSparkSource reference, final SAMSequenceDictionary referenceDictionary) {

        switch (type) {
            case INTER_CHR_STRAND_SWITCH_55:
            case INTER_CHR_STRAND_SWITCH_33:
            case INTER_CHR_NO_SS_WITH_LEFT_MATE_FIRST_IN_PARTNER:
            case INTER_CHR_NO_SS_WITH_LEFT_MATE_SECOND_IN_PARTNER:
            {
                final Tuple2<BreakEndVariantType, BreakEndVariantType> orderedMatesForTranslocSuspect =
                        BreakEndVariantType.InterChromosomeBreakend.getOrderedMates(this, reference);
                return Arrays.asList(orderedMatesForTranslocSuspect._1, orderedMatesForTranslocSuspect._2);
            }
            case INTRA_CHR_REF_ORDER_SWAP:
            {
                final Tuple2<BreakEndVariantType, BreakEndVariantType> orderedMatesForTranslocSuspect =
                        BreakEndVariantType.IntraChromosomeRefOrderSwap.getOrderedMates(this, reference);
                return Arrays.asList(orderedMatesForTranslocSuspect._1, orderedMatesForTranslocSuspect._2);
            }
            case INTRA_CHR_STRAND_SWITCH_55:
            case INTRA_CHR_STRAND_SWITCH_33:
                if ( complication.hasDuplicationAnnotation() ) { // inverted duplication
                    return Collections.singletonList( new SimpleSVType.DuplicationInverted(this, reference) );
                } else {
                    if (strandSwitch.equals(StrandSwitch.FORWARD_TO_REVERSE)) {
                        final Tuple2<BreakEndVariantType, BreakEndVariantType> orderedMatesForInversionSuspect =
                                BreakEndVariantType.IntraChromosomalStrandSwitch55BreakEnd.getOrderedMates(this, reference);
                        return Arrays.asList(orderedMatesForInversionSuspect._1, orderedMatesForInversionSuspect._2);
                    } else {
                        final Tuple2<BreakEndVariantType, BreakEndVariantType> orderedMatesForInversionSuspect =
                                BreakEndVariantType.IntraChromosomalStrandSwitch33BreakEnd.getOrderedMates(this, reference);
                        return Arrays.asList(orderedMatesForInversionSuspect._1, orderedMatesForInversionSuspect._2);
                    }
                }
            case SIMPLE_DEL:
            case DEL_DUP_CONTRACTION:
            {
                return Collections.singletonList( new SimpleSVType.Deletion(this, reference) );
            }
            case RPL:
            {
                if ( isCandidateForFatInsertion() ) {
                    return Collections.singletonList( new SimpleSVType.Insertion(this, reference) );
                } else { // "DEL" record with possibly linked "INS"
                    final SimpleSVType.Deletion deletion = new SimpleSVType.Deletion(this, reference);
                    if ( complication.getInsertedSequenceForwardStrandRep().length() < STRUCTURAL_VARIANT_SIZE_LOWER_BOUND ){ // ins seq not long enough for an INS call
                        return Collections.singletonList( deletion );
                    } else {
                        final SimpleSVType.Insertion insertion = new SimpleSVType.Insertion(this, reference);
                        return Arrays.asList( deletion, insertion );
                    }
                }
            }
            case SIMPLE_INS:
            {
                return Collections.singletonList( new SimpleSVType.Insertion(this, reference) );
            }
            case SMALL_DUP_EXPANSION:
            {
                // check if duplicated region is large enough to make an DUP call, if not, emit annotated INS call
                final BreakpointComplications.SmallDuplicationWithPreciseDupRangeBreakpointComplications duplicationComplication =
                        (BreakpointComplications.SmallDuplicationWithPreciseDupRangeBreakpointComplications) this.getComplication();
                if (duplicationComplication.getDupSeqRepeatUnitRefSpan().size() < STRUCTURAL_VARIANT_SIZE_LOWER_BOUND) {
                    return Collections.singletonList( new SimpleSVType.Insertion(this, reference));
                } else {
                    return Collections.singletonList( new SimpleSVType.DuplicationTandem(this, reference) );
                }
            }
            case SMALL_DUP_CPX:
            {
                final BreakpointComplications.SmallDuplicationWithImpreciseDupRangeBreakpointComplications duplicationComplication =
                        (BreakpointComplications.SmallDuplicationWithImpreciseDupRangeBreakpointComplications)
                        this.getComplication();
                if ( duplicationComplication.isDupContraction() ) {
                    return Collections.singletonList( new SimpleSVType.Deletion(this, reference) );
                } else {
                    if (duplicationComplication.getDupSeqRepeatUnitRefSpan().size() < STRUCTURAL_VARIANT_SIZE_LOWER_BOUND) {
                        return Collections.singletonList( new SimpleSVType.Insertion(this, reference));
                    } else {
                        return Collections.singletonList( new SimpleSVType.DuplicationTandem(this, reference) );
                    }
                }
            }
            default:
                throw new GATKException.ShouldNeverReachHereException("Inferred type not recognized");
        }
    }

    /**
     * <ul>
     *     <li>the new copies' length + inserted sequence length for simple expansion</li>
     *     <li>for complex expansion: the difference between affected reference region's size and the alt haplotype's length</li>
     *     <li>throws UnsupportedOperationException otherwise</li>
     * </ul>
     */
    public int getLengthForDupTandemExpansion() {
        // TODO: 6/11/18 the implementation taken below simply performs: inserted sequence length + copy number elevation * copy unit length,
        // which could be wrong when the expansion is complex, or when the copies are different in length due to small insertion and deletions in those copies
        // a better implementation is to simply calculate the difference in length in altseq and ref region size
        final BreakpointComplications.SmallDuplicationBreakpointComplications dupComplication = (BreakpointComplications.SmallDuplicationBreakpointComplications) getComplication();
        if (dupComplication.isDupContraction())
            throw new UnsupportedOperationException("Trying to extract length from a duplication contraction: " + this.toString());
        return dupComplication.getInsertedSequenceForwardStrandRep().length()
                + (dupComplication.getDupSeqRepeatNumOnCtg() - dupComplication.getDupSeqRepeatNumOnRef()) * dupComplication.getDupSeqRepeatUnitRefSpan().size();
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
