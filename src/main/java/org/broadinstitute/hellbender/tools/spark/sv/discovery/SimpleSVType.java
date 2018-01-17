package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.StrandSwitch;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.NovelAdjacencyReferenceLocations;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.EvidenceTargetLink;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.ReadMetadata;
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
        INV, DEL, INS, DUP, DUP_INV;
    }

    public static final class Inversion extends SimpleSVType {

        @Override
        public String toString() {
            return TYPES.INV.name();
        }

        public Inversion(final NovelAdjacencyReferenceLocations novelAdjacencyReferenceLocations) {
            super(getIDString(novelAdjacencyReferenceLocations),
                    Allele.create(createBracketedSymbAlleleString(GATKSVVCFConstants.SYMB_ALT_ALLELE_INV)),
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

    public static final class Deletion extends SimpleSVType {

        @Override
        public String toString() {
            return TYPES.DEL.name();
        }

        @SuppressWarnings("unchecked")
        public Deletion(final NovelAdjacencyReferenceLocations novelAdjacencyReferenceLocations) {
            super(getIDString(novelAdjacencyReferenceLocations),
                    Allele.create(createBracketedSymbAlleleString(GATKSVVCFConstants.SYMB_ALT_ALLELE_DEL)),
                    -(novelAdjacencyReferenceLocations.leftJustifiedRightRefLoc.getStart() - novelAdjacencyReferenceLocations.leftJustifiedLeftRefLoc.getEnd()),
                    novelAdjacencyReferenceLocations.complication.hasDuplicationAnnotation() ? Collections.singletonMap(GATKSVVCFConstants.DUP_TAN_CONTRACTION_STRING, "") : Collections.EMPTY_MAP);
        }

        private static String getIDString(final NovelAdjacencyReferenceLocations novelAdjacencyReferenceLocations) {

            return  ((novelAdjacencyReferenceLocations.complication.hasDuplicationAnnotation()) ? GATKSVVCFConstants.DUP_TAN_CONTRACTION_INTERNAL_ID_START_STRING : TYPES.DEL.name())
                    + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + novelAdjacencyReferenceLocations.leftJustifiedLeftRefLoc.getContig() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + novelAdjacencyReferenceLocations.leftJustifiedLeftRefLoc.getEnd() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + novelAdjacencyReferenceLocations.leftJustifiedRightRefLoc.getStart();
        }
    }

    public static final class Insertion extends SimpleSVType {

        @Override
        public String toString() {
            return TYPES.INS.name();
        }

        @SuppressWarnings("unchecked")
        public Insertion(final NovelAdjacencyReferenceLocations novelAdjacencyReferenceLocations) {
            super(getIDString(novelAdjacencyReferenceLocations),
                    Allele.create(createBracketedSymbAlleleString(GATKSVVCFConstants.SYMB_ALT_ALLELE_INS)),
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

    public static final class DuplicationTandem extends SimpleSVType {

        @Override
        public String toString() {
            return TYPES.DUP.name();
        }

        public DuplicationTandem(final NovelAdjacencyReferenceLocations novelAdjacencyReferenceLocations) {
            super(getIDString(novelAdjacencyReferenceLocations),
                    Allele.create(createBracketedSymbAlleleString(GATKSVVCFConstants.SYMB_ALT_ALLELE_DUP)),
                    novelAdjacencyReferenceLocations.complication.getInsertedSequenceForwardStrandRep().length()
                            + (novelAdjacencyReferenceLocations.complication.getDupSeqRepeatNumOnCtg() - novelAdjacencyReferenceLocations.complication.getDupSeqRepeatNumOnRef())*novelAdjacencyReferenceLocations.complication.getDupSeqRepeatUnitRefSpan().size(),
                    Collections.singletonMap(GATKSVVCFConstants.DUP_TAN_EXPANSION_STRING, ""));
        }

        private static String getIDString(final NovelAdjacencyReferenceLocations novelAdjacencyReferenceLocations) {

            return GATKSVVCFConstants.DUP_TAN_EXPANSION_INTERNAL_ID_START_STRING + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + novelAdjacencyReferenceLocations.leftJustifiedLeftRefLoc.getContig() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + novelAdjacencyReferenceLocations.leftJustifiedLeftRefLoc.getEnd() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + novelAdjacencyReferenceLocations.leftJustifiedRightRefLoc.getStart();
        }
    }

    public static final class ImpreciseDeletion extends SimpleSVType {

        @Override
        public String toString() {
            return TYPES.DEL.name();
        }

        @SuppressWarnings("unchecked")
        public ImpreciseDeletion(final EvidenceTargetLink evidenceTargetLink, final ReadMetadata metadata) {

            super(getIDString(evidenceTargetLink, metadata),
                    Allele.create(createBracketedSymbAlleleString(GATKSVVCFConstants.SYMB_ALT_ALLELE_DEL)),
                    (evidenceTargetLink.getPairedStrandedIntervals().getLeft().getInterval().midpoint() -
                            evidenceTargetLink.getPairedStrandedIntervals().getRight().getInterval().midpoint()),
                    Collections.EMPTY_MAP);
        }

        private static String getIDString(final EvidenceTargetLink evidenceTargetLink, final ReadMetadata metadata) {

            return TYPES.DEL.name()
                    + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + GATKSVVCFConstants.IMPRECISE + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + metadata.getContigName(evidenceTargetLink.getPairedStrandedIntervals().getLeft().getInterval().getContig())
                    + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + evidenceTargetLink.getPairedStrandedIntervals().getLeft().getInterval().getStart() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + evidenceTargetLink.getPairedStrandedIntervals().getLeft().getInterval().getEnd() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + evidenceTargetLink.getPairedStrandedIntervals().getRight().getInterval().getStart() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + evidenceTargetLink.getPairedStrandedIntervals().getRight().getInterval().getEnd();
        }
    }

    public static final class DuplicationInverted extends SimpleSVType {

        @Override
        public String toString() {
            return "DUP:INV";
        }

        @SuppressWarnings("unchecked")
        public DuplicationInverted(final NovelAdjacencyReferenceLocations novelAdjacencyReferenceLocations) {
            super(getIDString(novelAdjacencyReferenceLocations),
                    Allele.create(createBracketedSymbAlleleString(GATKSVVCFConstants.SYMB_ALT_ALLELE_INVDUP)),
                    novelAdjacencyReferenceLocations.complication.getDupSeqRepeatUnitRefSpan().size(),
                    Collections.EMPTY_MAP);
        }

        private static String getIDString(final NovelAdjacencyReferenceLocations novelAdjacencyReferenceLocations) {
            return GATKSVVCFConstants.DUP_INV_INTERNAL_ID_START_STRING + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + novelAdjacencyReferenceLocations.leftJustifiedLeftRefLoc.getContig() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + novelAdjacencyReferenceLocations.leftJustifiedLeftRefLoc.getEnd() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + novelAdjacencyReferenceLocations.leftJustifiedRightRefLoc.getStart();
        }
    }

}
