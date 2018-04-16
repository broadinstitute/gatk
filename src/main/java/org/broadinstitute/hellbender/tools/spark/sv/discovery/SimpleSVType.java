package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.StrandSwitch;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.NovelAdjacencyAndAltHaplotype;
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

    @Override
    public final boolean isBreakEndOnly() {
        return false;
    }

    public enum TYPES {
        INV, DEL, INS, DUP, DUP_INV;
    }

    public static final class Inversion extends SimpleSVType {

        @Override
        public String toString() {
            return TYPES.INV.name();
        }

        public Inversion(final NovelAdjacencyAndAltHaplotype novelAdjacencyAndAltHaplotype,
                         final int svLength) {
            super(getIDString(novelAdjacencyAndAltHaplotype),
                    Allele.create(createBracketedSymbAlleleString(GATKSVVCFConstants.SYMB_ALT_ALLELE_INV)),
                    svLength,
                    Collections.singletonMap((novelAdjacencyAndAltHaplotype.getStrandSwitch() == StrandSwitch.FORWARD_TO_REVERSE) ? GATKSVVCFConstants.INV55 : GATKSVVCFConstants.INV33, ""));
        }

        private static String getIDString(final NovelAdjacencyAndAltHaplotype novelAdjacencyAndAltHaplotype) {
            final String contig = novelAdjacencyAndAltHaplotype.getLeftJustifiedLeftRefLoc().getContig();
            final int start = novelAdjacencyAndAltHaplotype.getLeftJustifiedLeftRefLoc().getEnd();
            final int end = novelAdjacencyAndAltHaplotype.getLeftJustifiedRightRefLoc().getStart();
            final StrandSwitch strandSwitch = novelAdjacencyAndAltHaplotype.getStrandSwitch();

            return (strandSwitch.equals(StrandSwitch.FORWARD_TO_REVERSE) ? GATKSVVCFConstants.INV55 : GATKSVVCFConstants.INV33) + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR +
                    contig + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR + start + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR + end;
        }
    }

    public static final class Deletion extends SimpleSVType {

        @Override
        public String toString() {
            return TYPES.DEL.name();
        }

        public Deletion(final NovelAdjacencyAndAltHaplotype novelAdjacencyAndAltHaplotype,
                        final int svLength) {
            super(getIDString(novelAdjacencyAndAltHaplotype),
                    Allele.create(createBracketedSymbAlleleString(GATKSVVCFConstants.SYMB_ALT_ALLELE_DEL)),
                    svLength,
                    novelAdjacencyAndAltHaplotype.hasDuplicationAnnotation() ? Collections.singletonMap(GATKSVVCFConstants.DUP_TAN_CONTRACTION_STRING, "") : noExtraAttributes);
        }

        private static String getIDString(final NovelAdjacencyAndAltHaplotype novelAdjacencyAndAltHaplotype) {

            return  ((novelAdjacencyAndAltHaplotype.hasDuplicationAnnotation()) ? GATKSVVCFConstants.DUP_TAN_CONTRACTION_INTERNAL_ID_START_STRING : TYPES.DEL.name())
                    + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + novelAdjacencyAndAltHaplotype.getLeftJustifiedLeftRefLoc().getContig() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + novelAdjacencyAndAltHaplotype.getLeftJustifiedLeftRefLoc().getEnd() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + novelAdjacencyAndAltHaplotype.getLeftJustifiedRightRefLoc().getStart();
        }
    }

    public static final class Insertion extends SimpleSVType {

        @Override
        public String toString() {
            return TYPES.INS.name();
        }

        public Insertion(final NovelAdjacencyAndAltHaplotype novelAdjacencyAndAltHaplotype,
                         final int svLength) {
            super(getIDString(novelAdjacencyAndAltHaplotype),
                    Allele.create(createBracketedSymbAlleleString(GATKSVVCFConstants.SYMB_ALT_ALLELE_INS)),
                    svLength,
                    noExtraAttributes);
        }

        private static String getIDString(final NovelAdjacencyAndAltHaplotype novelAdjacencyAndAltHaplotype) {

            return TYPES.INS.name() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + novelAdjacencyAndAltHaplotype.getLeftJustifiedLeftRefLoc().getContig() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + novelAdjacencyAndAltHaplotype.getLeftJustifiedLeftRefLoc().getEnd() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + novelAdjacencyAndAltHaplotype.getLeftJustifiedRightRefLoc().getStart();
        }
    }

    public static final class DuplicationTandem extends SimpleSVType {

        @Override
        public String toString() {
            return TYPES.DUP.name();
        }

        public DuplicationTandem(final NovelAdjacencyAndAltHaplotype novelAdjacencyAndAltHaplotype,
                                 final int svLength) {
            super(getIDString(novelAdjacencyAndAltHaplotype),
                    Allele.create(createBracketedSymbAlleleString(GATKSVVCFConstants.SYMB_ALT_ALLELE_DUP)),
                    svLength,
                    Collections.singletonMap(GATKSVVCFConstants.DUP_TAN_EXPANSION_STRING, ""));
        }

        private static String getIDString(final NovelAdjacencyAndAltHaplotype novelAdjacencyAndAltHaplotype) {

            return GATKSVVCFConstants.DUP_TAN_EXPANSION_INTERNAL_ID_START_STRING + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + novelAdjacencyAndAltHaplotype.getLeftJustifiedLeftRefLoc().getContig() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + novelAdjacencyAndAltHaplotype.getLeftJustifiedLeftRefLoc().getEnd() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + novelAdjacencyAndAltHaplotype.getLeftJustifiedRightRefLoc().getStart();
        }
    }

    public static final class ImpreciseDeletion extends SimpleSVType {

        @Override
        public String toString() {
            return TYPES.DEL.name();
        }

        public ImpreciseDeletion(final EvidenceTargetLink evidenceTargetLink, final int svLength, final ReadMetadata metadata) {

            super(getIDString(evidenceTargetLink, metadata),
                    Allele.create(createBracketedSymbAlleleString(GATKSVVCFConstants.SYMB_ALT_ALLELE_DEL)),
                    svLength,
                    noExtraAttributes);
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

        public DuplicationInverted(final NovelAdjacencyAndAltHaplotype novelAdjacencyAndAltHaplotype,
                                   final int svLength) {
            super(getIDString(novelAdjacencyAndAltHaplotype),
                    Allele.create(createBracketedSymbAlleleString(GATKSVVCFConstants.SYMB_ALT_ALLELE_INVDUP)),
                    svLength,
                    noExtraAttributes);
        }

        private static String getIDString(final NovelAdjacencyAndAltHaplotype novelAdjacencyAndAltHaplotype) {
            return GATKSVVCFConstants.DUP_INV_INTERNAL_ID_START_STRING + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + novelAdjacencyAndAltHaplotype.getLeftJustifiedLeftRefLoc().getContig() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + novelAdjacencyAndAltHaplotype.getLeftJustifiedLeftRefLoc().getEnd() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + novelAdjacencyAndAltHaplotype.getLeftJustifiedRightRefLoc().getStart();
        }
    }

}
