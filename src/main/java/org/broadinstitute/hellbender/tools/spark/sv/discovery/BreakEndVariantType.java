package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Allele;
import org.apache.commons.lang3.ArrayUtils;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.StrandSwitch;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.NovelAdjacencyAndAltHaplotype;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVContext;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import scala.Tuple2;

import java.io.IOException;
import java.util.Collections;
import java.util.Map;

import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.*;

public abstract class BreakEndVariantType extends SvType {

    /**
     * Technically, a BND-formatted variant should have two VCF records, for mates, hence we also have this field.
     */
    protected final boolean isTheUpstreamMate;

    public final boolean isTheUpstreamMate() {
        return isTheUpstreamMate;
    }

    protected static SimpleInterval getMateRefLoc(final NovelAdjacencyAndAltHaplotype narl, final boolean forUpstreamLoc) {
        return forUpstreamLoc ? narl.getLeftJustifiedRightRefLoc() : narl.getLeftJustifiedLeftRefLoc();
    }

    protected static String getIDString(final NovelAdjacencyAndAltHaplotype narl, final boolean forUpstreamLoc) {
        return BREAKEND_STR + INTERVAL_VARIANT_ID_FIELD_SEPARATOR +
                narl.getLeftJustifiedLeftRefLoc().getContig() + INTERVAL_VARIANT_ID_FIELD_SEPARATOR +
                narl.getLeftJustifiedLeftRefLoc().getEnd() + INTERVAL_VARIANT_ID_FIELD_SEPARATOR +
                narl.getLeftJustifiedRightRefLoc().getContig() + INTERVAL_VARIANT_ID_FIELD_SEPARATOR +
                narl.getLeftJustifiedRightRefLoc().getStart() + INTERVAL_VARIANT_ID_FIELD_SEPARATOR +
                (forUpstreamLoc ? "1" : "2");
    }

    @Override
    public final boolean isBreakEndOnly() {
        return true;
    }

    @Override
    public String toString() {
        return BREAKEND_STR;
    }

    private static Allele constructAltAllele(final byte[] bases, final boolean bracketPointsLeft, final SimpleInterval novelAdjRefLoc,
                                             final boolean basesFirst) {
        final String s = bracketPointsLeft ? "]" + novelAdjRefLoc.getContig() + ":" + novelAdjRefLoc.getStart() + "]"
                : "[" + novelAdjRefLoc.getContig() + ":" + novelAdjRefLoc.getStart() + "[";
        return Allele.create( basesFirst ? new String(bases) + s  : s + new String(bases) );
    }

    // see VCF spec 4.2 or 4.3 for BND format ALT allele field for SV
    BreakEndVariantType(final String id, final Map<String, String> typeSpecificExtraAttributes,
                        final byte[] bases, final boolean bracketPointsLeft, final SimpleInterval novelAdjRefLoc,
                        final boolean basesFirst, final boolean isTheUpstreamMate) {
        super(id, constructAltAllele(bases, bracketPointsLeft, novelAdjRefLoc, basesFirst), SVContext.NO_LENGTH, typeSpecificExtraAttributes);
        this.isTheUpstreamMate = isTheUpstreamMate;
    }

    //==================================================================================================================

    /**
     * Breakend variant type for inversion suspects: those with novel adjacency between two reference locations
     * on the same chromosome but the novel adjacency brings them together in a strand-switch fashion.
     * This is to be distinguished from the more general "translocation" breakends, which are novel adjacency between
     * reference locations without strand switch if the reference bases are from the same chromosome.
     */
    public static final class InvSuspectBND extends BreakEndVariantType {
        /**
         * for breakends, there's a concept of partner (see VCF spec) relationship between two reference bases.
         * This indicates if the breakpoint, for an inversion suspect specifically, represents (one of the two) bases
         * that is the left of the pair.
         */
        private static final Map<String, String> INV55_FLAG = Collections.singletonMap(INV55, "");
        private static final Map<String, String> INV33_FLAG = Collections.singletonMap(INV33, "");

        public static final String INV55_BND = BREAKEND_STR + "_" + INV55;
        public static final String INV33_BND = BREAKEND_STR + "_" + INV33;

        private final boolean indicatesInv55;

        private static boolean indicatesInv55(boolean bracketPointsLeft, boolean basesFirst) {
            return (bracketPointsLeft && basesFirst);
        }

        @Override
        public String toString() {
            return indicatesInv55 ? INV55_BND : INV33_BND;
        }

        private InvSuspectBND(final String id, final byte[] bases, final SimpleInterval novelAdjRefLoc,
                              final boolean bracketPointsLeft, final boolean basesFirst, boolean isTheUpstreamMate) {
            super(id, indicatesInv55(bracketPointsLeft, basesFirst) ? INV55_FLAG: INV33_FLAG,
                    bases, bracketPointsLeft, novelAdjRefLoc, basesFirst, isTheUpstreamMate);
            indicatesInv55 = indicatesInv55(bracketPointsLeft, basesFirst);
        }

        public static Tuple2<BreakEndVariantType, BreakEndVariantType> getOrderedMates(final NovelAdjacencyAndAltHaplotype narl,
                                                                                       final ReferenceMultiSource reference) {

            // inversion breakend formatted records have "bracketPointsLeft" "basesFirst" taking the same value (see spec)
            final boolean isInv55Suspect;
            if (narl.getStrandSwitch() == StrandSwitch.FORWARD_TO_REVERSE) { // INV55, leftHalfInPartnerPair
                isInv55Suspect = true;
            } else if (narl.getStrandSwitch() == StrandSwitch.REVERSE_TO_FORWARD){
                isInv55Suspect = false;
            } else {
                throw new GATKException("Wrong type of novel adjacency sent to wrong analysis pathway: " +
                        "no strand-switch being sent to strand-switch path. \n" + narl.toString());
            }
            final BreakEndVariantType upstreamBreakpoint = new BreakEndVariantType.InvSuspectBND(getIDString(narl, true),
                    extractBasesForAltAllele(narl, true, reference),
                    getMateRefLoc(narl, true),
                    isInv55Suspect, isInv55Suspect, true);
            final BreakEndVariantType downstreamBreakpoint = new BreakEndVariantType.InvSuspectBND(getIDString(narl, false),
                    extractBasesForAltAllele(narl, false, reference),
                    getMateRefLoc(narl, false),
                    isInv55Suspect, isInv55Suspect, false);
            return new Tuple2<>(upstreamBreakpoint, downstreamBreakpoint);
        }

        static byte[] extractBasesForAltAllele(final NovelAdjacencyAndAltHaplotype narl, final boolean forUpstreamLoc,
                                               final ReferenceMultiSource reference) {
            try {
                final byte[] ref = reference
                        .getReferenceBases(forUpstreamLoc ? narl.getLeftJustifiedLeftRefLoc() :
                                narl.getLeftJustifiedRightRefLoc())
                        .getBases();
                final String ins = narl.getComplication().getInsertedSequenceForwardStrandRep();
                if (ins.isEmpty()) {
                    return ref;
                } else {
                    return forUpstreamLoc ? ArrayUtils.addAll(ref, ins.getBytes())
                            : ArrayUtils.addAll(ref, SequenceUtil.reverseComplement(ins).getBytes());
                }
            } catch (final IOException ioex) {
                throw new GATKException("Could not read reference for extracting reference bases.", ioex);
            }
        }
    }

    //==================================================================================================================
    /**
     * Generic variant type for suspected "translocations", including what could be a breakpoint for
     * a dispersed duplication or an inter-chromosome novel adjacency.
     */
    public static final class TransLocBND extends BreakEndVariantType {
        private static Map<String, String> emptyMap = Collections.emptyMap();

        public static final String STRANDSWITCHLESS_BND = BREAKEND_STR + "_NOSS";

        @Override
        public String toString() {
            return STRANDSWITCHLESS_BND;
        }

        private TransLocBND(final NovelAdjacencyAndAltHaplotype narl, final boolean forUpstreamLoc,
                            final ReferenceMultiSource reference, final SAMSequenceDictionary referenceDictionary,
                            final boolean basesFirst, final boolean bracketPointsLeft) {
            super(getIDString(narl, forUpstreamLoc), emptyMap,
                    extractBasesForAltAllele(narl, forUpstreamLoc, reference, referenceDictionary), bracketPointsLeft,
                    getMateRefLoc(narl, forUpstreamLoc), basesFirst, forUpstreamLoc);
        }

        public static Tuple2<BreakEndVariantType, BreakEndVariantType> getOrderedMates(final NovelAdjacencyAndAltHaplotype narl,
                                                                                       final ReferenceMultiSource reference,
                                                                                       final SAMSequenceDictionary referenceDictionary) {
            final boolean isSameChr = narl.getLeftJustifiedLeftRefLoc().getContig().equals(narl.getLeftJustifiedRightRefLoc().getContig());
            final BreakEndVariantType bkpt_1, bkpt_2;
            if (isSameChr) {
                bkpt_1 = new BreakEndVariantType.TransLocBND(narl, true, reference, referenceDictionary, false, true);
                bkpt_2 = new BreakEndVariantType.TransLocBND(narl, false, reference, referenceDictionary, true, false);
            } else {
                if (narl.getStrandSwitch() == StrandSwitch.NO_SWITCH) {
                    final boolean isFirstOfPartner = IntervalUtils.compareContigs(narl.getLeftJustifiedLeftRefLoc(),
                            narl.getLeftJustifiedRightRefLoc(), referenceDictionary) < 0;
                    bkpt_1 = new BreakEndVariantType.TransLocBND(narl, true, reference, referenceDictionary,
                            isFirstOfPartner, !isFirstOfPartner);
                    bkpt_2 = new BreakEndVariantType.TransLocBND(narl, false, reference, referenceDictionary,
                            !isFirstOfPartner, isFirstOfPartner);
                } else {
                    bkpt_1 = new BreakEndVariantType.TransLocBND(narl, true, reference, referenceDictionary,
                            narl.getStrandSwitch() == StrandSwitch.FORWARD_TO_REVERSE, narl.getStrandSwitch() == StrandSwitch.FORWARD_TO_REVERSE);
                    bkpt_2 = new BreakEndVariantType.TransLocBND(narl, false, reference, referenceDictionary,
                            narl.getStrandSwitch() != StrandSwitch.FORWARD_TO_REVERSE, narl.getStrandSwitch() != StrandSwitch.FORWARD_TO_REVERSE);
                }
            }
            return new Tuple2<>(bkpt_1, bkpt_2);
        }

        static byte[] extractBasesForAltAllele(final NovelAdjacencyAndAltHaplotype narl, final boolean forUpstreamLoc,
                                               final ReferenceMultiSource reference, final SAMSequenceDictionary referenceDictionary) {
            try {
                final SimpleInterval refLoc = forUpstreamLoc ? narl.getLeftJustifiedLeftRefLoc() : narl.getLeftJustifiedRightRefLoc();
                final byte[] ref = reference.getReferenceBases(refLoc).getBases();
                final String ins = narl.getComplication().getInsertedSequenceForwardStrandRep();
                if (ins.isEmpty()) {
                    return ref;
                } else {
                    if (narl.getLeftJustifiedLeftRefLoc().getContig().equals(narl.getLeftJustifiedRightRefLoc().getContig())
                            || (narl.getStrandSwitch() == StrandSwitch.NO_SWITCH && IntervalUtils.compareContigs(narl.getLeftJustifiedLeftRefLoc(), narl.getLeftJustifiedRightRefLoc(), referenceDictionary) > 0)
                            || narl.getStrandSwitch() == StrandSwitch.REVERSE_TO_FORWARD) {
                        return forUpstreamLoc ? ArrayUtils.addAll(ins.getBytes(), ref) : ArrayUtils.addAll(ref, ins.getBytes());
                    } else {
                        return forUpstreamLoc ? ArrayUtils.addAll(ref, ins.getBytes()) : ArrayUtils.addAll(ins.getBytes(), ref);
                    }
                }
            } catch (final IOException ioex) {
                throw new GATKException("Could not read reference for extracting reference bases.", ioex);
            }
        }
    }
}
