package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;

import java.util.Collections;
import java.util.Map;

public abstract class SimpleSVType extends SvType {
    public static String createBracketedSymbAlleleString(final String vcfHeaderDefinedSymbAltAllele) {
        return "<" + vcfHeaderDefinedSymbAltAllele + ">";
    }

    protected SimpleSVType(final String id, final Allele altAllele, final int len, final Map<String, String> typeSpecificExtraAttributes) {
        super(id, altAllele, len, typeSpecificExtraAttributes);
    }

    public enum TYPES {
        INV, DEL, INS, DUP;
    }

    public static final class Inversion extends SvType {

        @Override
        public String toString() {
            return TYPES.INV.name();
        }

        Inversion(final NovelAdjacencyReferenceLocations novelAdjacencyReferenceLocations) {
            super(getIDString(novelAdjacencyReferenceLocations),
                    Allele.create(createBracketedSymbAlleleString(GATKSVVCFConstants.SYMB_ALT_ALLELE_INV_IN_HEADER)),
                    novelAdjacencyReferenceLocations.leftJustifiedRightRefLoc.getStart() - novelAdjacencyReferenceLocations.leftJustifiedLeftRefLoc.getEnd(),
                    Collections.singletonMap((novelAdjacencyReferenceLocations.strandSwitch == StrandSwitch.FORWARD_TO_REVERSE) ? GATKSVVCFConstants.INV55 : GATKSVVCFConstants.INV33, ""));
        }

        private static String getIDString(final NovelAdjacencyReferenceLocations novelAdjacencyReferenceLocations) {
            final String contig = novelAdjacencyReferenceLocations.leftJustifiedLeftRefLoc.getContig();
            final int start = novelAdjacencyReferenceLocations.leftJustifiedLeftRefLoc.getEnd();
            final int end = novelAdjacencyReferenceLocations.leftJustifiedRightRefLoc.getStart();
            final StrandSwitch strandSwitch = novelAdjacencyReferenceLocations.strandSwitch;

            return (strandSwitch.equals(StrandSwitch.FORWARD_TO_REVERSE) ? GATKSVVCFConstants.INV55 : GATKSVVCFConstants.INV33) + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR +
                    contig + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR + start + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR + end;
        }
    }

    public static final class Deletion extends SvType {

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

    public static final class Insertion extends SvType {

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

    public static final class DuplicationTandem extends SvType {

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
}
