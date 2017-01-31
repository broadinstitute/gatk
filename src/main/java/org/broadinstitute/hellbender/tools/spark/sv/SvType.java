package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.variant.variantcontext.Allele;

import java.util.Collections;
import java.util.Map;

import static org.broadinstitute.hellbender.tools.spark.sv.NovelAdjacencyReferenceLocations.EndConnectionType.FIVE_TO_FIVE;

/**
 * Various types of structural variations
 */
abstract class SvType {

    static final String TANDEMUPLICATION_CONTRACTION_ID_START_STRING = "DEL-DUPLICATION-TANDEM-CONTRACTION";
    static final String TANDEMUPLICATION_EXPANSION_ID_START_STRING = "INS-DUPLICATION-TANDEM-EXPANSION";


    protected final String variantId;
    protected final Allele altAllele;
    protected final int svLen;
    protected final Map<String, String> extraAttributes;

    enum TYPES {
        INV, DEL, INS, DUP;
    }

    protected SvType(final String id, final Allele altAllele, final int len, final Map<String, String> typeSpecificExtraAttributes) {
        variantId = id;
        this.altAllele = altAllele;
        svLen = len;
        extraAttributes = typeSpecificExtraAttributes;
    }

    final String getVariantId() {
        return variantId;
    }
    final Allele getAltAllele() {
        return altAllele;
    }
    final int getSVLength() {
        return svLen;
    }
    final Map<String, String> getTypeSpecificAttributes() {
        return extraAttributes;
    }

    static final class Inversion extends SvType {

        @Override
        public String toString() {
            return TYPES.INV.name();
        }

        Inversion(final NovelAdjacencyReferenceLocations novelAdjacencyReferenceLocations) {
            super(getIDString(novelAdjacencyReferenceLocations),
                    Allele.create(SVConstants.CallingStepConstants.VCF_ALT_ALLELE_STRING_INV),
                    novelAdjacencyReferenceLocations.leftJustifiedRightRefLoc.getStart() - novelAdjacencyReferenceLocations.leftJustifiedLeftRefLoc.getEnd(),
                    Collections.singletonMap((novelAdjacencyReferenceLocations.endConnectionType == FIVE_TO_FIVE) ? GATKSVVCFHeaderLines.INV55 : GATKSVVCFHeaderLines.INV33, ""));
        }

        private static String getIDString(final NovelAdjacencyReferenceLocations novelAdjacencyReferenceLocations) {
            final String contig = novelAdjacencyReferenceLocations.leftJustifiedLeftRefLoc.getContig();
            final int start = novelAdjacencyReferenceLocations.leftJustifiedLeftRefLoc.getEnd();
            final int end = novelAdjacencyReferenceLocations.leftJustifiedRightRefLoc.getStart();
            final NovelAdjacencyReferenceLocations.EndConnectionType endConnectionType = novelAdjacencyReferenceLocations.endConnectionType;

            return (endConnectionType == FIVE_TO_FIVE ? GATKSVVCFHeaderLines.INV55 : GATKSVVCFHeaderLines.INV33) + SVConstants.CallingStepConstants.VARIANT_ID_FIELD_SEPARATOR +
                    contig + SVConstants.CallingStepConstants.VARIANT_ID_FIELD_SEPARATOR + start + SVConstants.CallingStepConstants.VARIANT_ID_FIELD_SEPARATOR + end;
        }
    }

    static final class Deletion extends SvType {

        @Override
        public String toString() {
            return TYPES.DEL.name();
        }

        @SuppressWarnings("unchecked")
        Deletion(final NovelAdjacencyReferenceLocations novelAdjacencyReferenceLocations) {
            super(getIDString(novelAdjacencyReferenceLocations),
                    Allele.create(SVConstants.CallingStepConstants.VCF_ALT_ALLELE_STRING_DEL),
                    -(novelAdjacencyReferenceLocations.leftJustifiedRightRefLoc.getStart() - novelAdjacencyReferenceLocations.leftJustifiedLeftRefLoc.getEnd()),
                    novelAdjacencyReferenceLocations.complication.dupSeqForwardStrandRep.isEmpty() ? Collections.EMPTY_MAP : Collections.singletonMap(SVConstants.CallingStepConstants.TANDUP_CONTRACTION_STRING, ""));
        }

        private static String getIDString(final NovelAdjacencyReferenceLocations novelAdjacencyReferenceLocations) {

            return  ((novelAdjacencyReferenceLocations.complication.dupSeqForwardStrandRep.isEmpty()) ? TYPES.DEL.name() : TANDEMUPLICATION_CONTRACTION_ID_START_STRING)
                    + SVConstants.CallingStepConstants.VARIANT_ID_FIELD_SEPARATOR
                    + novelAdjacencyReferenceLocations.leftJustifiedLeftRefLoc.getContig() + SVConstants.CallingStepConstants.VARIANT_ID_FIELD_SEPARATOR
                    + novelAdjacencyReferenceLocations.leftJustifiedLeftRefLoc.getEnd() + SVConstants.CallingStepConstants.VARIANT_ID_FIELD_SEPARATOR
                    + novelAdjacencyReferenceLocations.leftJustifiedRightRefLoc.getStart();
        }
    }

    static final class Insertion extends SvType {

        @Override
        public String toString() {
            return TYPES.INS.name();
        }

        @SuppressWarnings("unchecked")
        Insertion(final NovelAdjacencyReferenceLocations novelAdjacencyReferenceLocations) {
            super(getIDString(novelAdjacencyReferenceLocations),
                    Allele.create(SVConstants.CallingStepConstants.VCF_ALT_ALLELE_STRING_INS),
                    novelAdjacencyReferenceLocations.complication.insertedSequenceForwardStrandRep.length(),
                    Collections.EMPTY_MAP);
        }

        private static String getIDString(final NovelAdjacencyReferenceLocations novelAdjacencyReferenceLocations) {

            return TYPES.INS.name() + SVConstants.CallingStepConstants.VARIANT_ID_FIELD_SEPARATOR
                    + novelAdjacencyReferenceLocations.leftJustifiedLeftRefLoc.getContig() + SVConstants.CallingStepConstants.VARIANT_ID_FIELD_SEPARATOR
                    + novelAdjacencyReferenceLocations.leftJustifiedLeftRefLoc.getEnd() + SVConstants.CallingStepConstants.VARIANT_ID_FIELD_SEPARATOR
                    + novelAdjacencyReferenceLocations.leftJustifiedRightRefLoc.getStart();
        }
    }

    static final class DuplicationTandem extends SvType {

        @Override
        public String toString() {
            return TYPES.DUP.name();
        }

        DuplicationTandem(final NovelAdjacencyReferenceLocations novelAdjacencyReferenceLocations) {
            super(getIDString(novelAdjacencyReferenceLocations),
                    Allele.create(SVConstants.CallingStepConstants.VCF_ALT_ALLELE_STRING_DUP),
                    novelAdjacencyReferenceLocations.complication.insertedSequenceForwardStrandRep.length()
                            + (novelAdjacencyReferenceLocations.complication.dupSeqRepeatNumOnCtg - novelAdjacencyReferenceLocations.complication.dupSeqRepeatNumOnRef)*novelAdjacencyReferenceLocations.complication.dupSeqForwardStrandRep.length(),
                    Collections.singletonMap(SVConstants.CallingStepConstants.TANDUP_EXPANSION_STRING, ""));
        }

        private static String getIDString(final NovelAdjacencyReferenceLocations novelAdjacencyReferenceLocations) {

            return TANDEMUPLICATION_EXPANSION_ID_START_STRING + SVConstants.CallingStepConstants.VARIANT_ID_FIELD_SEPARATOR
                    + novelAdjacencyReferenceLocations.leftJustifiedLeftRefLoc.getContig() + SVConstants.CallingStepConstants.VARIANT_ID_FIELD_SEPARATOR
                    + novelAdjacencyReferenceLocations.leftJustifiedLeftRefLoc.getEnd() + SVConstants.CallingStepConstants.VARIANT_ID_FIELD_SEPARATOR
                    + novelAdjacencyReferenceLocations.leftJustifiedRightRefLoc.getStart();
        }
    }
}
