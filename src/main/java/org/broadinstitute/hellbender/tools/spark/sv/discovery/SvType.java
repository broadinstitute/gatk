package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;

import java.util.Collections;
import java.util.Map;

import static org.broadinstitute.hellbender.tools.spark.sv.discovery.NovelAdjacencyReferenceLocations.EndConnectionType;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.NovelAdjacencyReferenceLocations.EndConnectionType.*;

/**
 * Various types of structural variations
 */
public abstract class SvType {

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

    final String getInternalVariantId() {
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
                    Allele.create(createBracketedSymbAlleleString(GATKSVVCFConstants.SYMB_ALT_ALLELE_INV_IN_HEADER)),
                    novelAdjacencyReferenceLocations.leftJustifiedRightRefLoc.getStart() - novelAdjacencyReferenceLocations.leftJustifiedLeftRefLoc.getEnd(),
                    Collections.singletonMap((novelAdjacencyReferenceLocations.endConnectionType == FIVE_TO_FIVE) ? GATKSVVCFConstants.INV55 : GATKSVVCFConstants.INV33, ""));
        }

        private static String getIDString(final NovelAdjacencyReferenceLocations novelAdjacencyReferenceLocations) {
            final String contig = novelAdjacencyReferenceLocations.leftJustifiedLeftRefLoc.getContig();
            final int start = novelAdjacencyReferenceLocations.leftJustifiedLeftRefLoc.getEnd();
            final int end = novelAdjacencyReferenceLocations.leftJustifiedRightRefLoc.getStart();
            final EndConnectionType endConnectionType = novelAdjacencyReferenceLocations.endConnectionType;

            return (endConnectionType == FIVE_TO_FIVE ? GATKSVVCFConstants.INV55 : GATKSVVCFConstants.INV33) + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR +
                    contig + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR + start + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR + end;
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
                    Allele.create(createBracketedSymbAlleleString(GATKSVVCFConstants.SYMB_ALT_ALLELE_DEL_IN_HEADER)),
                    -(novelAdjacencyReferenceLocations.leftJustifiedRightRefLoc.getStart() - novelAdjacencyReferenceLocations.leftJustifiedLeftRefLoc.getEnd()),
                    novelAdjacencyReferenceLocations.complication.hasDuplicationAnnotation() ? Collections.singletonMap(GATKSVVCFConstants.TANDUP_CONTRACTION_STRING, "") : Collections.EMPTY_MAP);
        }

        private static String getIDString(final NovelAdjacencyReferenceLocations novelAdjacencyReferenceLocations) {

            return  ((novelAdjacencyReferenceLocations.complication.hasDuplicationAnnotation()) ? GATKSVVCFConstants.TANDUP_CONTRACTION_INTERNAL_ID_START_STRING : TYPES.DEL.name())
                    + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + novelAdjacencyReferenceLocations.leftJustifiedLeftRefLoc.getContig() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + novelAdjacencyReferenceLocations.leftJustifiedLeftRefLoc.getEnd() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
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
                    Allele.create(createBracketedSymbAlleleString(GATKSVVCFConstants.SYMB_ALT_ALLELE_INS_IN_HEADER)),
                    novelAdjacencyReferenceLocations.complication.getInsertedSequenceForwardStrandRep().length(),
                    Collections.EMPTY_MAP);
        }

        private static String getIDString(final NovelAdjacencyReferenceLocations novelAdjacencyReferenceLocations) {

            return TYPES.INS.name() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + novelAdjacencyReferenceLocations.leftJustifiedLeftRefLoc.getContig() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + novelAdjacencyReferenceLocations.leftJustifiedLeftRefLoc.getEnd() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
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
                    Allele.create(createBracketedSymbAlleleString(GATKSVVCFConstants.SYMB_ALT_ALLELE_DUP_IN_HEADER)),
                    novelAdjacencyReferenceLocations.complication.getInsertedSequenceForwardStrandRep().length()
                            + (novelAdjacencyReferenceLocations.complication.getDupSeqRepeatNumOnCtg() - novelAdjacencyReferenceLocations.complication.getDupSeqRepeatNumOnRef())*novelAdjacencyReferenceLocations.complication.getDupSeqRepeatUnitRefSpan().size(),
                    Collections.singletonMap(GATKSVVCFConstants.TANDUP_EXPANSION_STRING, ""));
        }

        private static String getIDString(final NovelAdjacencyReferenceLocations novelAdjacencyReferenceLocations) {

            return GATKSVVCFConstants.TANDUP_EXPANSION_INTERNAL_ID_START_STRING + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + novelAdjacencyReferenceLocations.leftJustifiedLeftRefLoc.getContig() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + novelAdjacencyReferenceLocations.leftJustifiedLeftRefLoc.getEnd() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + novelAdjacencyReferenceLocations.leftJustifiedRightRefLoc.getStart();
        }
    }

    public static String createBracketedSymbAlleleString(final String vcfHeaderDefinedSymbAltAllele) {
        return "<" + vcfHeaderDefinedSymbAltAllele + ">";
    }
}
