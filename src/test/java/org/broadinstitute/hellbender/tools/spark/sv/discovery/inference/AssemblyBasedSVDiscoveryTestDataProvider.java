package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceMultiSparkSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.BreakEndVariantType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SvType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import scala.Tuple2;

import java.util.*;

import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.*;

/**
 * Holding utility methods used for constructing test data
 * as well as expected results.
 *
 * NO TESTS ARE RUN IN THIS PARTICULAR CLASS AND ITS CHILD CLASSES.
 */
public abstract class AssemblyBasedSVDiscoveryTestDataProvider {

    public abstract static class AssemblyBasedSVDiscoveryTestDataForSimpleChimera {
        // these are used for building
        public final AlignmentInterval firstAlignment;
        public final AlignmentInterval secondAlignment;
        public final String evidenceAssemblyContigName;
        public final byte[] evidenceContigSeq;

        // these are used for testing
        public final boolean expectedFirstContigRegionHasLaterReferenceMapping;
        public final SimpleChimera expectedSimpleChimera;
        public final NovelAdjacencyAndAltHaplotype expectedNovelAdjacencyAndAltSeq;

        public final List<SvType> expectedSvTypes;
        public final List<VariantContext> expectedVariantContexts;

        public final Class<? extends BreakpointsInference> expectedInferencerClass;

        public AssemblyBasedSVDiscoveryTestDataForSimpleChimera(final AlignmentInterval firstAlignment, final AlignmentInterval secondAlignment,
                                                                final String evidenceAssemblyContigName, final byte[] evidenceContigSeq,
                                                                final boolean expectedFirstContigRegionHasLaterReferenceMapping,
                                                                final SimpleChimera expectedSimpleChimera,
                                                                final NovelAdjacencyAndAltHaplotype expectedNovelAdjacencyAndAltSeq,
                                                                final List<SvType> expectedSvTypes,
                                                                final List<VariantContext> expectedVariantContexts,
                                                                final Class<? extends BreakpointsInference> expectedInferencerClass) {
            this.firstAlignment = firstAlignment;
            this.secondAlignment = secondAlignment;
            this.evidenceAssemblyContigName = evidenceAssemblyContigName;
            this.evidenceContigSeq = evidenceContigSeq;
            this.expectedFirstContigRegionHasLaterReferenceMapping = expectedFirstContigRegionHasLaterReferenceMapping;
            this.expectedSimpleChimera = expectedSimpleChimera;
            this.expectedNovelAdjacencyAndAltSeq = expectedNovelAdjacencyAndAltSeq;
            this.expectedSvTypes = expectedSvTypes;
            this.expectedVariantContexts = expectedVariantContexts;
            this.expectedInferencerClass = expectedInferencerClass;
        }

        public abstract SAMSequenceDictionary getAppropriateDictionary();

        public abstract ReferenceMultiSparkSource getAppropriateRef();

        public abstract Class<? extends BreakpointsInference> getAppropriateBreakpointInferencer();
    }

    public abstract List<AssemblyBasedSVDiscoveryTestDataForSimpleChimera> getAllTestData();

    // same event, two representations from opposite strands
    public final List<Tuple2<AssemblyBasedSVDiscoveryTestDataForSimpleChimera, AssemblyBasedSVDiscoveryTestDataForSimpleChimera>> getAllTestDataPaired() {
        final List<AssemblyBasedSVDiscoveryTestDataForSimpleChimera> allTestData = getAllTestData();
        final List<Tuple2<AssemblyBasedSVDiscoveryTestDataForSimpleChimera, AssemblyBasedSVDiscoveryTestDataForSimpleChimera>> testDataForSimpleSVs
                = new ArrayList<>(allTestData.size()/2);
        for (int i = 0; i < allTestData.size() - 1; i += 2) {
            testDataForSimpleSVs.add(new Tuple2<>(allTestData.get(i), allTestData.get(i+1)));
        }
        return Collections.unmodifiableList(testDataForSimpleSVs);
    }

    // data block ======================================================================================================

    private static String makeID(final String prefix, final String chr1, final int start, final String chr2, final int stop,
                                 final String postfix) {
        return prefix + INTERVAL_VARIANT_ID_FIELD_SEPARATOR +
                SvType.makeLocationString(chr1, start, chr2, stop) +
                (postfix.isEmpty() ? "" : INTERVAL_VARIANT_ID_FIELD_SEPARATOR + postfix);
    }

    static final byte[] EMPTY_BYTE_ARRAY = new byte[]{};

    private static final Allele INV_SYMB_ALLELE = Allele.create(SimpleSVType.createBracketedSymbAlleleString(SimpleSVType.SupportedType.INV.name()), false);
    private static final Allele DEL_SYMB_ALLELE = Allele.create(SimpleSVType.createBracketedSymbAlleleString(SimpleSVType.SupportedType.DEL.name()), false);
    private static final Allele INS_SYMB_ALLELE = Allele.create(SimpleSVType.createBracketedSymbAlleleString(SimpleSVType.SupportedType.INS.name()), false);
    private static final Allele DUP_SYMB_ALLELE = Allele.create(SimpleSVType.createBracketedSymbAlleleString(SimpleSVType.SupportedType.DUP.name()), false);

    // utils block =====================================================================================================

    /**
     * Note that {@code delRange} is expected to be VCF spec compatible,
     * e.g. if chr1:101-200 is deleted, then {@code delRange} should be chr1:100-200
     */
    static final VariantContextBuilder makeDeletion(final SimpleInterval delRange, final Allele refAllele, final boolean isFromDupContraction) {

        return new VariantContextBuilder()
                .chr(delRange.getContig()).start(delRange.getStart()).stop(delRange.getEnd())
                .alleles(Arrays.asList(refAllele, DEL_SYMB_ALLELE))
                .id(makeID((isFromDupContraction? DUP_TAN_CONTRACTION_INTERNAL_ID_START_STRING : SimpleSVType.SupportedType.DEL.name()), delRange.getContig(), delRange.getStart(), delRange.getContig(), delRange.getEnd(), ""))
                .attribute(VCFConstants.END_KEY, delRange.getEnd())
                .attribute(SVLEN, - delRange.size() + 1)
                .attribute(SVTYPE, SimpleSVType.SupportedType.DEL.name());
    }

    static final SvType makeDeletionType(final SimpleInterval delRange, final Allele refAllele, final boolean isFromDupContraction) {
        return new SimpleSVType.Deletion(
                delRange.getContig(), delRange.getStart(), delRange.getEnd(),
                makeID((isFromDupContraction? DUP_TAN_CONTRACTION_INTERNAL_ID_START_STRING : SimpleSVType.SupportedType.DEL.name()),
                        delRange.getContig(), delRange.getStart(), delRange.getContig(), delRange.getEnd(), ""),
                refAllele,
                DEL_SYMB_ALLELE,
                - delRange.size()+ 1,
                isFromDupContraction ? Collections.singletonMap(DUP_TAN_CONTRACTION_STRING, true) :Collections.emptyMap());
    }

    static final VariantContextBuilder makeInversion(final SimpleInterval invertedRegion, final Allele refAllele) {
        return new VariantContextBuilder()
                .chr(invertedRegion.getContig()).start(invertedRegion.getStart() - 1).stop(invertedRegion.getEnd())     // TODO: 5/2/18 VCF spec doesn't require left shift by 1 for inversion POS
                .alleles(Arrays.asList(refAllele, INV_SYMB_ALLELE))
                .id(makeID(SimpleSVType.SupportedType.INV.name(), invertedRegion.getContig(), invertedRegion.getStart() - 1, invertedRegion.getContig(), invertedRegion.getEnd(), ""))
                .attribute(VCFConstants.END_KEY, invertedRegion.getEnd())
                .attribute(SVLEN, 0)                                                                 // TODO: 5/2/18 this is following VCF spec,
                .attribute(SVTYPE, SimpleSVType.SupportedType.INV.name());
    }

    static final SvType makeInversionType(final SimpleInterval invRange, final Allele refAllele, final boolean isInv55) {
        return new SimpleSVType.Inversion(
                invRange.getContig(),
                invRange.getStart() - 1,
                invRange.getEnd(),
                makeID((isInv55? INV55 : INV33),
                        invRange.getContig(), invRange.getStart() - 1, invRange.getContig(), invRange.getEnd(), ""),
                refAllele,
                INV_SYMB_ALLELE, invRange.size(),
                Collections.singletonMap((isInv55) ? INV55 : INV33, true));
    }

    static final VariantContextBuilder makeInsertion(final String chr, final int pos, final int end, final int svLen,
                                                     final Allele refAllele) {

        return new VariantContextBuilder().chr(chr).start(pos).stop(end)
                .alleles(Arrays.asList(refAllele, INS_SYMB_ALLELE))
                .id(makeID(SimpleSVType.SupportedType.INS.name(), chr, pos, chr, end, ""))
                .attribute(VCFConstants.END_KEY, end)
                .attribute(SVLEN, svLen)
                .attribute(SVTYPE, SimpleSVType.SupportedType.INS.name());
    }

    static final SvType makeInsertionType(final SimpleInterval insertionPos, final Allele refAllele, final int svLen) {
        return new SimpleSVType.Insertion(
                insertionPos.getContig(),
                insertionPos.getStart(),
                insertionPos.getEnd(),
                makeID(SimpleSVType.SupportedType.INS.name(),
                        insertionPos.getContig(), insertionPos.getStart(), insertionPos.getContig(), insertionPos.getEnd(), ""),
                refAllele,
                INS_SYMB_ALLELE, svLen,
                Collections.emptyMap());
    }

    static final VariantContextBuilder makeTandemDuplication(final SimpleInterval duplicatedRange, final Allele refAllele,
                                                             final int svLen) {
        return new VariantContextBuilder().chr(duplicatedRange.getContig()).start(duplicatedRange.getStart() - 1).stop(duplicatedRange.getStart() - 1) // TODO: 5/24/18 by left shifting 1, we are treating it as insertions
                .id(makeID(DUP_TAN_EXPANSION_INTERNAL_ID_START_STRING,
                        duplicatedRange.getContig(), duplicatedRange.getStart() , duplicatedRange.getContig(), duplicatedRange.getEnd(), ""))
                .alleles(Arrays.asList(refAllele, DUP_SYMB_ALLELE))
                .attribute(VCFConstants.END_KEY, duplicatedRange.getStart() - 1) // TODO: 5/24/18 see todo above
                .attribute(SVLEN, svLen)
                .attribute(SVTYPE, SimpleSVType.SupportedType.DUP.name());
    }

    static final SvType makeTandemDuplicationType(final SimpleInterval duplicatedRange, final Allele refAllele, final int svLen) {
        return new SimpleSVType.DuplicationTandem(
                duplicatedRange.getContig(),
                duplicatedRange.getStart() - 1,
                duplicatedRange.getStart() - 1,
                makeID(DUP_TAN_EXPANSION_INTERNAL_ID_START_STRING,
                        duplicatedRange.getContig(), duplicatedRange.getStart(), duplicatedRange.getContig(), duplicatedRange.getEnd(), ""),
                refAllele,
                DUP_SYMB_ALLELE, svLen,
                Collections.singletonMap(DUP_TAN_EXPANSION_STRING, true));
    }

    static final VariantContextBuilder makeBND(final SimpleInterval upstreamLoc, final SimpleInterval dnstreamLoc,
                                               final Allele refAllele, final String insertedSeq, final String bndSubtypeString,
                                               final boolean forUpstreamLoc, final boolean refBaseFirst, final boolean bracketPointsLeft) {

        final String upstreamAltString = upstreamLoc.getContig() + ":" + upstreamLoc.getStart();
        final String dnstreamAltString = dnstreamLoc.getContig() + ":" + dnstreamLoc.getEnd();

        final Allele altAllele;
        if (refBaseFirst) {
            if (bracketPointsLeft)
                altAllele = Allele.create(refAllele.getBaseString() + insertedSeq + "]" + (forUpstreamLoc ? dnstreamAltString : upstreamAltString) + "]");
            else
                altAllele = Allele.create(refAllele.getBaseString() + insertedSeq + "[" + (forUpstreamLoc ? dnstreamAltString : upstreamAltString) + "[");
        } else {
            if (bracketPointsLeft)
                altAllele = Allele.create("]" + (forUpstreamLoc ? dnstreamAltString : upstreamAltString) + "]" + insertedSeq + refAllele.getBaseString());
            else
                altAllele = Allele.create("[" + (forUpstreamLoc ? dnstreamAltString : upstreamAltString) + "[" + insertedSeq + refAllele.getBaseString());
        }
        if (forUpstreamLoc) {
            return new VariantContextBuilder().chr(upstreamLoc.getContig()).start(upstreamLoc.getStart()).stop(upstreamLoc.getEnd())
                    .alleles(Arrays.asList(refAllele, altAllele))
                    .id(makeID(BREAKEND_STR + (bndSubtypeString.isEmpty() ? "" : INTERVAL_VARIANT_ID_FIELD_SEPARATOR + bndSubtypeString),
                            upstreamLoc.getContig(), upstreamLoc.getStart(), dnstreamLoc.getContig(), dnstreamLoc.getEnd(), "1"))
                    .attribute(SVTYPE, BREAKEND_STR);
        } else {
            return new VariantContextBuilder().chr(dnstreamLoc.getContig()).start(dnstreamLoc.getStart()).stop(dnstreamLoc.getEnd())
                    .alleles(Arrays.asList(refAllele, altAllele))
                    .id(makeID(BREAKEND_STR + (bndSubtypeString.isEmpty() ? "" : INTERVAL_VARIANT_ID_FIELD_SEPARATOR + bndSubtypeString),
                            upstreamLoc.getContig(), upstreamLoc.getStart(), dnstreamLoc.getContig(), dnstreamLoc.getEnd(), "2"))
                    .attribute(SVTYPE, BREAKEND_STR);
        }
    }

    static final SvType makeBNDType(final String variantCHR, final int variantPOS, final String variantId,
                                    final Allele refAllele, final Allele altAllele, final Map<String, Object> extraAttributes,
                                    final boolean isTheUpstreamMate,
                                    final BreakEndVariantType.SupportedType type) {
        switch (type) {
            case INTRA_CHR_STRAND_SWITCH_55:
                return new BreakEndVariantType.IntraChromosomalStrandSwitch55BreakEnd(variantCHR, variantPOS, variantId, refAllele, altAllele, extraAttributes, isTheUpstreamMate);
            case INTRA_CHR_STRAND_SWITCH_33:
                return new BreakEndVariantType.IntraChromosomalStrandSwitch33BreakEnd(variantCHR, variantPOS, variantId, refAllele, altAllele, extraAttributes, isTheUpstreamMate);
            case INTRA_CHR_REF_ORDER_SWAP:
                return new BreakEndVariantType.IntraChromosomeRefOrderSwap(variantCHR, variantPOS, variantId, refAllele, altAllele, extraAttributes, isTheUpstreamMate);
            case INTER_CHR_STRAND_SWITCH_55:
            case INTER_CHR_STRAND_SWITCH_33:
            case INTER_CHR_NO_SS_WITH_LEFT_MATE_FIRST_IN_PARTNER:
            case INTER_CHR_NO_SS_WITH_LEFT_MATE_SECOND_IN_PARTNER:
                return new BreakEndVariantType.InterChromosomeBreakend(variantCHR, variantPOS, variantId, refAllele, altAllele, extraAttributes, isTheUpstreamMate);
            default:
                throw new GATKException("Unrecognized type: " + type.name());
        }
    }
}
