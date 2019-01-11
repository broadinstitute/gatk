package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceMultiSparkSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.StrandSwitch;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.NovelAdjacencyAndAltHaplotype;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.TypeInferredFromSimpleChimera;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import scala.Tuple2;

import java.io.IOException;
import java.util.Collections;
import java.util.Map;

import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.*;

public abstract class BreakEndVariantType extends SvType {

    /**
     * Technically, a BND-formatted variant should have two VCF records, for mates, hence we also have this field.
     * Upstream mate is defined as the location in a mate pair that has a lower coordinate according to
     * the reference sequence dictionary.
     */
    private final boolean isTheUpstreamMate;

    protected BreakEndVariantType(final String variantCHR, final int variantPOS, final String variantId,
                                  final Allele refAllele, final Allele altAllele, final Map<String, Object> extraAttributes,
                                  final boolean isTheUpstreamMate) {
        super(variantCHR, variantPOS, NO_APPLICABLE_END, variantId, refAllele, altAllele, NO_APPLICABLE_LEN, extraAttributes);
        this.isTheUpstreamMate = isTheUpstreamMate;
    }

    public final boolean isTheUpstreamMate() {
        return isTheUpstreamMate;
    }

    @Override
    public final boolean hasApplicableEnd() {
        return false;
    }
    @Override
    public final boolean hasApplicableLength() {
        return false;
    }

    @Override
    public final String toString() {
        return BREAKEND_STR;
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        if (!super.equals(o)) return false;

        final BreakEndVariantType that = (BreakEndVariantType) o;

        return isTheUpstreamMate == that.isTheUpstreamMate;
    }

    @Override
    public int hashCode() {
        int result = super.hashCode();
        result = 31 * result + (isTheUpstreamMate ? 1 : 0);
        return result;
    }

    //==================================================================================================================

    private static String getIDString(final NovelAdjacencyAndAltHaplotype narl, final boolean forUpstreamLoc) {
        // if no strand switch or different contig, "", otherwise append INV55/33
        final String bndtype = narl.getStrandSwitch().equals(StrandSwitch.NO_SWITCH) || !narl.getLeftJustifiedLeftRefLoc().getContig().equals(narl.getLeftJustifiedRightRefLoc().getContig())? ""
                : (narl.getStrandSwitch().equals(StrandSwitch.FORWARD_TO_REVERSE) ? INV55 : INV33);
        String locationPartOfString = makeLocationString(narl.getLeftJustifiedLeftRefLoc().getContig(),
                narl.getLeftJustifiedLeftRefLoc().getStart(), narl.getLeftJustifiedRightRefLoc().getContig(),
                narl.getLeftJustifiedRightRefLoc().getEnd());
        return BREAKEND_STR + INTERVAL_VARIANT_ID_FIELD_SEPARATOR +
                (bndtype.isEmpty() ? "" : bndtype + INTERVAL_VARIANT_ID_FIELD_SEPARATOR) +
               locationPartOfString + INTERVAL_VARIANT_ID_FIELD_SEPARATOR + (forUpstreamLoc ? "1" : "2");
    }

    private static String getRefBaseString(final NovelAdjacencyAndAltHaplotype narl, final boolean forUpstreamLoc,
                                           final ReferenceMultiSparkSource reference) {
        try {
            byte[] refBases = reference.getReferenceBases(forUpstreamLoc ? narl.getLeftJustifiedLeftRefLoc() :
                    narl.getLeftJustifiedRightRefLoc())
                    .getBases();
            return new String(refBases);
        } catch (final IOException ioex) {
            throw new GATKException("Could not read reference for extracting reference bases.", ioex);
        }
    }

    public enum SupportedType {
        INTRA_CHR_STRAND_SWITCH_55,// intra-chromosome strand-switch novel adjacency, alignments left-flanking the novel adjacency
        INTRA_CHR_STRAND_SWITCH_33,// intra-chromosome strand-switch novel adjacency, alignments right-flanking the novel adjacency

        INTRA_CHR_REF_ORDER_SWAP,// intra-chromosome reference-order swap, but NO strand-switch, novel adjacency

        INTER_CHR_STRAND_SWITCH_55,// pair WY in Fig.1 in Section 5.4 of VCF spec ver.4.2
        INTER_CHR_STRAND_SWITCH_33,// pair XZ in Fig.1 in Section 5.4 of VCF spec ver.4.2
        INTER_CHR_NO_SS_WITH_LEFT_MATE_FIRST_IN_PARTNER, // the green pair in Fig. 7 in Section 5.4 of VCF spec ver.4.2
        INTER_CHR_NO_SS_WITH_LEFT_MATE_SECOND_IN_PARTNER; // the red pair in Fig. 7 in Section 5.4 of VCF spec ver.4.2
    }

    /**
     * Breakend variant type for inversion suspects: those with novel adjacency between two reference locations
     * on the same chromosome but the novel adjacency brings them together in a strand-switch fashion.
     * This is to be distinguished from the more general "translocation" breakends, which are novel adjacency between
     * reference locations without strand switch if the reference bases are from the same chromosome.
     *
     * Note that dispersed duplication with some copies inverted could also lead to breakpoints with strand switch.
     */
    abstract private static class IntraChromosomalStrandSwitchBreakEnd extends BreakEndVariantType {
        static final Map<String, Object> INV55_FLAG = Collections.singletonMap(INV55, true);
        static final Map<String, Object> INV33_FLAG = Collections.singletonMap(INV33, true);

        private IntraChromosomalStrandSwitchBreakEnd(final String variantCHR, final int variantPOS, final String variantId,
                                                     final Allele refAllele, final Allele altAllele,
                                                     final Map<String, Object> extraAttributes,
                                                     final boolean isTheUpstreamMate) {
            super(variantCHR, variantPOS, variantId, refAllele, altAllele, extraAttributes, isTheUpstreamMate);
        }

        @VisibleForTesting
        static String extractInsertedSequence(final NovelAdjacencyAndAltHaplotype narl, final boolean forUpstreamLoc) {
            final String ins = narl.getComplication().getInsertedSequenceForwardStrandRep();
            return forUpstreamLoc ? ins : SequenceUtil.reverseComplement(ins);
        }
    }

    public static final class IntraChromosomalStrandSwitch55BreakEnd extends IntraChromosomalStrandSwitchBreakEnd {

        @VisibleForTesting
        public IntraChromosomalStrandSwitch55BreakEnd(final String variantCHR, final int variantPOS, final String variantId,
                                                      final Allele refAllele, final Allele altAllele, final Map<String, Object> extraAttributes,
                                                      final boolean isTheUpstreamMate) {
            super(variantCHR, variantPOS, variantId, refAllele, altAllele, extraAttributes, isTheUpstreamMate);
        }

        private IntraChromosomalStrandSwitch55BreakEnd(final NovelAdjacencyAndAltHaplotype narl,
                                                       final ReferenceMultiSparkSource reference,
                                                       final boolean isTheUpstreamMate) {
            super(isTheUpstreamMate ? narl.getLeftJustifiedLeftRefLoc().getContig() : narl.getLeftJustifiedRightRefLoc().getContig(),
                    isTheUpstreamMate ? narl.getLeftJustifiedLeftRefLoc().getStart() : narl.getLeftJustifiedRightRefLoc().getEnd(),
                    BreakEndVariantType.getIDString(narl, isTheUpstreamMate),
                    Allele.create(BreakEndVariantType.getRefBaseString(narl, isTheUpstreamMate, reference), true),
                    constructAltAllele(BreakEndVariantType.getRefBaseString(narl, isTheUpstreamMate, reference),
                            extractInsertedSequence(narl, isTheUpstreamMate),
                            isTheUpstreamMate ? narl.getLeftJustifiedRightRefLoc(): narl.getLeftJustifiedLeftRefLoc()),
                    INV55_FLAG, isTheUpstreamMate);
        }

        public static Tuple2<BreakEndVariantType, BreakEndVariantType> getOrderedMates(final NovelAdjacencyAndAltHaplotype narl,
                                                                                       final ReferenceMultiSparkSource reference) {
            return new Tuple2<>(new IntraChromosomalStrandSwitch55BreakEnd(narl, reference, true),
                                new IntraChromosomalStrandSwitch55BreakEnd(narl, reference, false));
        }

        private static Allele constructAltAllele(final String refBase, final String insertedSequence, final SimpleInterval novelAdjRefLoc) {
            return Allele.create(refBase + insertedSequence + "]" + novelAdjRefLoc.getContig() + ":" + novelAdjRefLoc.getEnd() + "]");
        }
    }

    public static final class IntraChromosomalStrandSwitch33BreakEnd extends IntraChromosomalStrandSwitchBreakEnd {

        @VisibleForTesting
        public IntraChromosomalStrandSwitch33BreakEnd(final String variantCHR, final int variantPOS, final String variantId,
                                                      final Allele refAllele, final Allele altAllele, final Map<String, Object> extraAttributes,
                                                      final boolean isTheUpstreamMate) {
            super(variantCHR, variantPOS, variantId, refAllele, altAllele, extraAttributes, isTheUpstreamMate);
        }

        private IntraChromosomalStrandSwitch33BreakEnd(final NovelAdjacencyAndAltHaplotype narl, final ReferenceMultiSparkSource reference,
                                                       final boolean isTheUpstreamMate) {
            super(isTheUpstreamMate ? narl.getLeftJustifiedLeftRefLoc().getContig() : narl.getLeftJustifiedRightRefLoc().getContig(),
                    isTheUpstreamMate ? narl.getLeftJustifiedLeftRefLoc().getStart() : narl.getLeftJustifiedRightRefLoc().getEnd(),
                    BreakEndVariantType.getIDString(narl, isTheUpstreamMate),
                    Allele.create(BreakEndVariantType.getRefBaseString(narl, isTheUpstreamMate, reference), true),
                    constructAltAllele(BreakEndVariantType.getRefBaseString(narl, isTheUpstreamMate, reference),
                            extractInsertedSequence(narl, isTheUpstreamMate),
                            isTheUpstreamMate ? narl.getLeftJustifiedRightRefLoc(): narl.getLeftJustifiedLeftRefLoc()),
                    INV33_FLAG, isTheUpstreamMate);
        }

        public static Tuple2<BreakEndVariantType, BreakEndVariantType> getOrderedMates(final NovelAdjacencyAndAltHaplotype narl,
                                                                                       final ReferenceMultiSparkSource reference) {
            return new Tuple2<>(new IntraChromosomalStrandSwitch33BreakEnd(narl, reference, true),
                                new IntraChromosomalStrandSwitch33BreakEnd(narl, reference, false));
        }

        private static Allele constructAltAllele(final String refBase, final String insertedSequence, final SimpleInterval novelAdjRefLoc) {
            return Allele.create("[" + novelAdjRefLoc.getContig() + ":" + novelAdjRefLoc.getEnd() + "[" + insertedSequence + refBase);
        }
    }

    public static final class IntraChromosomeRefOrderSwap extends BreakEndVariantType {

        @VisibleForTesting
        public IntraChromosomeRefOrderSwap(final String variantCHR, final int variantPOS, final String variantId,
                                           final Allele refAllele, final Allele altAllele, final Map<String, Object> extraAttributes,
                                           final boolean isTheUpstreamMate) {
            super(variantCHR, variantPOS, variantId, refAllele, altAllele, extraAttributes, isTheUpstreamMate);
        }

        private IntraChromosomeRefOrderSwap(final NovelAdjacencyAndAltHaplotype narl, final ReferenceMultiSparkSource reference,
                                            final boolean isTheUpstreamMate) {
            super(isTheUpstreamMate ? narl.getLeftJustifiedLeftRefLoc().getContig() : narl.getLeftJustifiedRightRefLoc().getContig(),
                    isTheUpstreamMate ? narl.getLeftJustifiedLeftRefLoc().getStart() : narl.getLeftJustifiedRightRefLoc().getEnd(),
                    BreakEndVariantType.getIDString(narl, isTheUpstreamMate),
                    Allele.create(BreakEndVariantType.getRefBaseString(narl, isTheUpstreamMate, reference), true),
                    constructAltAllele(BreakEndVariantType.getRefBaseString(narl, isTheUpstreamMate, reference),
                            narl.getComplication().getInsertedSequenceForwardStrandRep(),
                            isTheUpstreamMate ? narl.getLeftJustifiedRightRefLoc(): narl.getLeftJustifiedLeftRefLoc(),
                            isTheUpstreamMate),
                    noExtraAttributes, isTheUpstreamMate);
        }

        public static Tuple2<BreakEndVariantType, BreakEndVariantType> getOrderedMates(final NovelAdjacencyAndAltHaplotype narl,
                                                                                       final ReferenceMultiSparkSource reference) {
            return new Tuple2<>(new IntraChromosomeRefOrderSwap(narl, reference, true),
                                new IntraChromosomeRefOrderSwap(narl, reference, false));
        }

        private static Allele constructAltAllele(final String refBase, final String insertedSequence, final SimpleInterval novelAdjRefLoc,
                                                 final boolean forUpstreamLoc) {
            if (forUpstreamLoc) {
                return Allele.create("]" + novelAdjRefLoc.getContig() + ":" + novelAdjRefLoc.getEnd() + "]" + insertedSequence + refBase);
            } else {
                return Allele.create(refBase + insertedSequence + "[" + novelAdjRefLoc.getContig() + ":" + novelAdjRefLoc.getEnd() + "[");
            }
        }
    }

    public static final class InterChromosomeBreakend extends BreakEndVariantType {

        @VisibleForTesting
        public InterChromosomeBreakend(final String variantCHR, final int variantPOS, final String variantId,
                                       final Allele refAllele, final Allele altAllele, final Map<String, Object> extraAttributes,
                                       final boolean isTheUpstreamMate) {
            super(variantCHR, variantPOS, variantId, refAllele, altAllele, extraAttributes, isTheUpstreamMate);
        }

        private InterChromosomeBreakend(final NovelAdjacencyAndAltHaplotype narl, final ReferenceMultiSparkSource reference,
                                        final boolean isTheUpstreamMate) {
            super(isTheUpstreamMate ? narl.getLeftJustifiedLeftRefLoc().getContig() : narl.getLeftJustifiedRightRefLoc().getContig(),
                    isTheUpstreamMate ? narl.getLeftJustifiedLeftRefLoc().getStart() : narl.getLeftJustifiedRightRefLoc().getEnd(),
                    BreakEndVariantType.getIDString(narl, isTheUpstreamMate),
                    Allele.create(BreakEndVariantType.getRefBaseString(narl, isTheUpstreamMate, reference), true),
                    constructAltAllele(narl, reference, isTheUpstreamMate),
                    noExtraAttributes, isTheUpstreamMate);
        }

        public static Tuple2<BreakEndVariantType, BreakEndVariantType> getOrderedMates(final NovelAdjacencyAndAltHaplotype narl,
                                                                                       final ReferenceMultiSparkSource reference) {

            return new Tuple2<>(new InterChromosomeBreakend(narl, reference, true),
                                new InterChromosomeBreakend(narl, reference, false));
        }

        // see VCF spec 4.2 for BND format ALT allele field for SV, in particular the examples shown in Fig.1, Fig.2 and Fig.5 of Section 5.4
        private static Allele constructAltAllele(final NovelAdjacencyAndAltHaplotype narl, final ReferenceMultiSparkSource reference,
                                                 final boolean forUpstreamLoc) {
            final String refBase = BreakEndVariantType.getRefBaseString(narl, forUpstreamLoc, reference);
            final String insertedSequence = extractInsertedSequence(narl, forUpstreamLoc);
            final SimpleInterval novelAdjRefLoc = forUpstreamLoc ? narl.getLeftJustifiedRightRefLoc() : narl.getLeftJustifiedLeftRefLoc();

            // see Fig.5 of Section 5.4 of spec Version 4.2 (the green pairs)
            final boolean upstreamLocIsFirstInPartner =
                    narl.getTypeInferredFromSimpleChimera().equals(TypeInferredFromSimpleChimera.INTER_CHR_NO_SS_WITH_LEFT_MATE_FIRST_IN_PARTNER);
            if (narl.getStrandSwitch().equals(StrandSwitch.NO_SWITCH)) {
                if (forUpstreamLoc == upstreamLocIsFirstInPartner) {
                    return Allele.create(refBase + insertedSequence + "[" + novelAdjRefLoc.getContig() + ":" + novelAdjRefLoc.getEnd() + "[");
                } else {
                    return Allele.create("]" + novelAdjRefLoc.getContig() + ":" + novelAdjRefLoc.getStart() + "]" + insertedSequence + refBase);
                }
            } else if (narl.getStrandSwitch().equals(StrandSwitch.FORWARD_TO_REVERSE)){
                return Allele.create(refBase + insertedSequence + "]" + novelAdjRefLoc.getContig() + ":" + novelAdjRefLoc.getEnd() + "]");
            } else {
                return Allele.create("[" + novelAdjRefLoc.getContig() + ":" + novelAdjRefLoc.getEnd() + "[" + insertedSequence + refBase);
            }
        }

        private static String extractInsertedSequence(final NovelAdjacencyAndAltHaplotype narl, final boolean forUpstreamLoc) {
            final String ins = narl.getComplication().getInsertedSequenceForwardStrandRep();
            if (ins.isEmpty() || narl.getStrandSwitch() == StrandSwitch.NO_SWITCH) {
                return ins;
            } else {
                return forUpstreamLoc == (narl.getStrandSwitch().equals(StrandSwitch.FORWARD_TO_REVERSE) ) ? ins: SequenceUtil.reverseComplement(ins);
            }
        }
    }
}
