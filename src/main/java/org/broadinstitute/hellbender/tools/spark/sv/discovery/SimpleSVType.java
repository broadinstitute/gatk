package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceMultiSparkSource;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.StrandSwitch;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.BreakpointComplications;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.NovelAdjacencyAndAltHaplotype;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.EvidenceTargetLink;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.ReadMetadata;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.Collections;
import java.util.Map;

public abstract class SimpleSVType extends SvType {
    public static String createBracketedSymbAlleleString(final String vcfHeaderDefinedSymbAltAllele) {
        return "<" + vcfHeaderDefinedSymbAltAllele + ">";
    }

    protected SimpleSVType(final String variantCHR, final int variantPOS, final int variantEND, final String variantId,
                           final Allele refAllele, final Allele altAllele, final int svLen, final Map<String, Object> extraAttributes) {
        super(variantCHR, variantPOS, variantEND, variantId, refAllele, altAllele, svLen, extraAttributes);
    }

    @Override
    public final boolean hasApplicableEnd() {
        return true;
    }
    @Override
    public final boolean hasApplicableLength() {
        return true;
    }

    public enum SupportedType {
        INV, DEL, INS, DUP, DUP_INV;
    }

    public static final class Inversion extends SimpleSVType {

        @Override
        public String toString() {
            return SupportedType.INV.name();
        }

        @VisibleForTesting
        public Inversion(final String variantCHR, final int variantPOS, final int variantEND, final String variantId,
                         final Allele refAllele, final Allele altAllele, final int svLen, final Map<String, Object> extraAttributes) {
            super(variantCHR, variantPOS, variantEND, variantId, refAllele, altAllele, svLen, extraAttributes);
        }

        @VisibleForTesting

        // TODO: 6/12/18 note the following implementation sets POS and REF at the anchor base, which is not requested by the VCF spec
        // TODO: 6/12/18 also, this interface lets one call inversion with SVLEN !=0, which is not the same as VCF spec examples
        public Inversion(final NovelAdjacencyAndAltHaplotype novelAdjacencyAndAltHaplotype, final int svLength,
                         final ReferenceMultiSparkSource reference) {
            super(novelAdjacencyAndAltHaplotype.getLeftJustifiedLeftRefLoc().getContig(),
                    novelAdjacencyAndAltHaplotype.getLeftJustifiedLeftRefLoc().getStart(),
                    novelAdjacencyAndAltHaplotype.getLeftJustifiedRightRefLoc().getEnd(),
                    getIDString(novelAdjacencyAndAltHaplotype),
                    Allele.create(extractRefBases(novelAdjacencyAndAltHaplotype.getLeftJustifiedLeftRefLoc(), reference), true),
                    Allele.create(createBracketedSymbAlleleString(GATKSVVCFConstants.SYMB_ALT_ALLELE_INV)),
                    svLength,
                    Collections.singletonMap((novelAdjacencyAndAltHaplotype.getStrandSwitch() == StrandSwitch.FORWARD_TO_REVERSE) ? GATKSVVCFConstants.INV55 : GATKSVVCFConstants.INV33, true));
        }

        private static String getIDString(final NovelAdjacencyAndAltHaplotype novelAdjacencyAndAltHaplotype) {
            final StrandSwitch strandSwitch = novelAdjacencyAndAltHaplotype.getStrandSwitch();

            return (strandSwitch.equals(StrandSwitch.FORWARD_TO_REVERSE) ? GATKSVVCFConstants.INV55 : GATKSVVCFConstants.INV33) + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR +
                    makeLocationString(novelAdjacencyAndAltHaplotype);
        }
    }

    public static final class Deletion extends SimpleSVType {

        @Override
        public String toString() {
            return SupportedType.DEL.name();
        }

        @VisibleForTesting
        public Deletion(final String variantCHR, final int variantPOS, final int variantEND, final String variantId,
                        final Allele refAllele, final Allele altAllele, final int svLen, final Map<String, Object> extraAttributes) {
            super(variantCHR, variantPOS, variantEND, variantId, refAllele, altAllele, svLen, extraAttributes);
        }

        public Deletion(final NovelAdjacencyAndAltHaplotype novelAdjacencyAndAltHaplotype,
                        final ReferenceMultiSparkSource reference) {
            super(novelAdjacencyAndAltHaplotype.getLeftJustifiedLeftRefLoc().getContig(),
                    novelAdjacencyAndAltHaplotype.getLeftJustifiedLeftRefLoc().getStart(),
                    novelAdjacencyAndAltHaplotype.getLeftJustifiedRightRefLoc().getEnd(),
                    getIDString(novelAdjacencyAndAltHaplotype),
                    Allele.create(extractRefBases(novelAdjacencyAndAltHaplotype.getLeftJustifiedLeftRefLoc(), reference), true),
                    Allele.create(createBracketedSymbAlleleString(GATKSVVCFConstants.SYMB_ALT_ALLELE_DEL)),
                    - novelAdjacencyAndAltHaplotype.getDistanceBetweenNovelAdjacencies(),
                    novelAdjacencyAndAltHaplotype.hasDuplicationAnnotation() ? Collections.singletonMap(GATKSVVCFConstants.DUP_TAN_CONTRACTION_STRING, true) : noExtraAttributes);
        }

        private static String getIDString(final NovelAdjacencyAndAltHaplotype novelAdjacencyAndAltHaplotype) {

            return  ((novelAdjacencyAndAltHaplotype.hasDuplicationAnnotation()) ? GATKSVVCFConstants.DUP_TAN_CONTRACTION_INTERNAL_ID_START_STRING : SupportedType.DEL.name())
                    + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + makeLocationString(novelAdjacencyAndAltHaplotype);
        }
    }

    public static final class Insertion extends SimpleSVType {

        @Override
        public String toString() {
            return SupportedType.INS.name();
        }

        @VisibleForTesting
        public Insertion(final String variantCHR, final int variantPOS, final int variantEND, final String variantId,
                         final Allele refAllele, final Allele altAllele, final int svLen, final Map<String, Object> extraAttributes) {
            super(variantCHR, variantPOS, variantEND, variantId, refAllele, altAllele, svLen, extraAttributes);
        }

        public Insertion(final NovelAdjacencyAndAltHaplotype novelAdjacencyAndAltHaplotype,
                         final ReferenceMultiSparkSource reference) {
            super(novelAdjacencyAndAltHaplotype.getLeftJustifiedLeftRefLoc().getContig(),
                    novelAdjacencyAndAltHaplotype.getLeftJustifiedLeftRefLoc().getStart(),
                    getEnd(novelAdjacencyAndAltHaplotype),
                    getIDString(novelAdjacencyAndAltHaplotype),
                    Allele.create(getRefBases(novelAdjacencyAndAltHaplotype, reference), true),
                    Allele.create(createBracketedSymbAlleleString(GATKSVVCFConstants.SYMB_ALT_ALLELE_INS)),
                    getLength(novelAdjacencyAndAltHaplotype),
                    noExtraAttributes);
        }

        // these methods exist to distinguish fat insertion and linked del+ins in an RPL event, as well as duplication events whose duplicated unit is not large enough
        private static int getEnd(final NovelAdjacencyAndAltHaplotype novelAdjacencyAndAltHaplotype) {
            return novelAdjacencyAndAltHaplotype.isCandidateForFatInsertion()
                    ? novelAdjacencyAndAltHaplotype.getLeftJustifiedRightRefLoc().getEnd()
                    : novelAdjacencyAndAltHaplotype.getLeftJustifiedLeftRefLoc().getStart();
        }

        private static byte[] getRefBases(final NovelAdjacencyAndAltHaplotype novelAdjacencyAndAltHaplotype,
                                          final ReferenceMultiSparkSource reference) {
            return extractRefBases(novelAdjacencyAndAltHaplotype.isCandidateForFatInsertion()
                    ? novelAdjacencyAndAltHaplotype.getIntervalForFatInsertion()
                    : novelAdjacencyAndAltHaplotype.getLeftJustifiedLeftRefLoc(), reference);
        }

        private static int getLength(final NovelAdjacencyAndAltHaplotype novelAdjacencyAndAltHaplotype) {
            return novelAdjacencyAndAltHaplotype.getComplication().hasDuplicationAnnotation()
                    ? novelAdjacencyAndAltHaplotype.getLengthForDupTandemExpansion()
                    : novelAdjacencyAndAltHaplotype.getComplication().getInsertedSequenceForwardStrandRep().length();
        }

        private static String getIDString(final NovelAdjacencyAndAltHaplotype novelAdjacencyAndAltHaplotype) {
            if (novelAdjacencyAndAltHaplotype.isCandidateForFatInsertion())
                return SupportedType.INS.name() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                        + makeLocationString(novelAdjacencyAndAltHaplotype);
            else {
                return SupportedType.INS.name() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                        + makeLocationString(novelAdjacencyAndAltHaplotype.getLeftJustifiedLeftRefLoc().getContig(),
                                                novelAdjacencyAndAltHaplotype.getLeftJustifiedLeftRefLoc().getStart(),
                                                novelAdjacencyAndAltHaplotype.getLeftJustifiedLeftRefLoc().getContig(),
                                                novelAdjacencyAndAltHaplotype.getLeftJustifiedLeftRefLoc().getStart());
            }
        }
    }

    public static final class DuplicationTandem extends SimpleSVType {

        @Override
        public String toString() {
            return SupportedType.DUP.name();
        }

        @VisibleForTesting
        public DuplicationTandem(final String variantCHR, final int variantPOS, final int variantEND, final String variantId,
                                 final Allele refAllele, final Allele altAllele, final int svLen, final Map<String, Object> extraAttributes) {
            super(variantCHR, variantPOS, variantEND, variantId, refAllele, altAllele, svLen, extraAttributes);
        }

        // TODO: 6/12/18 the following implementation treats DuplicationTandem as simple insertions, and duplication annotations will be saved in INFO columns
        public DuplicationTandem(final NovelAdjacencyAndAltHaplotype novelAdjacencyAndAltHaplotype,
                                 final ReferenceMultiSparkSource reference) {
            super(novelAdjacencyAndAltHaplotype.getLeftJustifiedLeftRefLoc().getContig(),
                    novelAdjacencyAndAltHaplotype.getLeftJustifiedLeftRefLoc().getStart(),
                    novelAdjacencyAndAltHaplotype.getLeftJustifiedLeftRefLoc().getStart(),
                    getIDString(novelAdjacencyAndAltHaplotype),
                    Allele.create(extractRefBases(novelAdjacencyAndAltHaplotype.getLeftJustifiedLeftRefLoc(), reference), true),
                    Allele.create(createBracketedSymbAlleleString(GATKSVVCFConstants.SYMB_ALT_ALLELE_DUP)),
                    novelAdjacencyAndAltHaplotype.getLengthForDupTandemExpansion(),
                    Collections.singletonMap(GATKSVVCFConstants.DUP_TAN_EXPANSION_STRING, true));
        }

        private static String getIDString(final NovelAdjacencyAndAltHaplotype novelAdjacencyAndAltHaplotype) {

            final SimpleInterval dupSeqRepeatUnitRefSpan = ((BreakpointComplications.SmallDuplicationBreakpointComplications)
                    novelAdjacencyAndAltHaplotype.getComplication()).getDupSeqRepeatUnitRefSpan();

            return GATKSVVCFConstants.DUP_TAN_EXPANSION_INTERNAL_ID_START_STRING + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + makeLocationString(dupSeqRepeatUnitRefSpan.getContig(), dupSeqRepeatUnitRefSpan.getStart(),
                    dupSeqRepeatUnitRefSpan.getContig(), dupSeqRepeatUnitRefSpan.getEnd());
        }
    }

    public static final class ImpreciseDeletion extends SimpleSVType {

        @Override
        public String toString() {
            return SupportedType.DEL.name();
        }

        public ImpreciseDeletion(final EvidenceTargetLink evidenceTargetLink, final int svLength, final ReadMetadata metadata,
                                 final ReferenceMultiSparkSource reference) {

            super(metadata.getContigName(evidenceTargetLink.getPairedStrandedIntervals().getLeft().getInterval().getContig()),
                    evidenceTargetLink.getPairedStrandedIntervals().getLeft().getInterval().midpoint(),
                    evidenceTargetLink.getPairedStrandedIntervals().getRight().getInterval().midpoint(),
                    getIDString(evidenceTargetLink, metadata),
                    Allele.create(getRefBases(evidenceTargetLink, metadata, reference), true),
                    Allele.create(createBracketedSymbAlleleString(GATKSVVCFConstants.SYMB_ALT_ALLELE_DEL)),
                    svLength,
                    Collections.singletonMap(GATKSVVCFConstants.IMPRECISE, true));
        }

        private static byte[] getRefBases(final EvidenceTargetLink evidenceTargetLink, final ReadMetadata metadata,
                                          final ReferenceMultiSparkSource reference) {
            final SVInterval leftInterval = evidenceTargetLink.getPairedStrandedIntervals().getLeft().getInterval();
            return extractRefBases(
                    new SimpleInterval(metadata.getContigName(leftInterval.getContig()), leftInterval.midpoint(), leftInterval.midpoint()),
                    reference
            );
        }

        private static String getIDString(final EvidenceTargetLink evidenceTargetLink, final ReadMetadata metadata) {
            final SVInterval leftInterval = evidenceTargetLink.getPairedStrandedIntervals().getLeft().getInterval();
            final SVInterval rightInterval = evidenceTargetLink.getPairedStrandedIntervals().getRight().getInterval();

            return SupportedType.DEL.name()
                    + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + GATKSVVCFConstants.IMPRECISE + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + metadata.getContigName(leftInterval.getContig())
                    + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + leftInterval.getStart() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + leftInterval.getEnd() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + rightInterval.getStart() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + rightInterval.getEnd();
        }
    }

    public static final class DuplicationInverted extends SimpleSVType {

        @Override
        public String toString() {
            return "DUP:INV";
        }

        public DuplicationInverted(final NovelAdjacencyAndAltHaplotype novelAdjacencyAndAltHaplotype,
                                   final ReferenceMultiSparkSource reference) {
            super(novelAdjacencyAndAltHaplotype.getLeftJustifiedLeftRefLoc().getContig(),
                    novelAdjacencyAndAltHaplotype.getLeftJustifiedLeftRefLoc().getStart(),
                    novelAdjacencyAndAltHaplotype.getLeftJustifiedLeftRefLoc().getStart(),
                    getIDString(novelAdjacencyAndAltHaplotype),
                    Allele.create(extractRefBases(novelAdjacencyAndAltHaplotype.getLeftJustifiedLeftRefLoc(), reference), true),
                    Allele.create(createBracketedSymbAlleleString(GATKSVVCFConstants.SYMB_ALT_ALLELE_INVDUP)),
                    getSVLen(novelAdjacencyAndAltHaplotype),
                    noExtraAttributes);
        }

        private static int getSVLen(final NovelAdjacencyAndAltHaplotype novelAdjacencyAndAltHaplotype) {
            return ((BreakpointComplications.InvertedDuplicationBreakpointComplications) novelAdjacencyAndAltHaplotype.getComplication())
                    .getDupSeqRepeatUnitRefSpan().size();
        }

        private static String getIDString(final NovelAdjacencyAndAltHaplotype novelAdjacencyAndAltHaplotype) {
            return GATKSVVCFConstants.DUP_INV_INTERNAL_ID_START_STRING + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + makeLocationString(novelAdjacencyAndAltHaplotype);
        }

    }

}
