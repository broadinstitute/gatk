package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceMultiSparkSource;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SvType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.TestUtilsForAssemblyBasedSVDiscovery;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.ContigAlignmentsModifier;
import org.broadinstitute.hellbender.tools.spark.sv.utils.Strand;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import static org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments.NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.StrandSwitch.NO_SWITCH;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.SimpleChimera.DistancesBetweenAlignmentsOnRefAndOnRead;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.TypeInferredFromSimpleChimera.*;
import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.*;

/**
 * Holding test data and expected values for
 *
 *   INSERTION
 *   DELETION
 *   SMALL-DUPLICATION-EXPANSION
 *   SMALL-DUPLICATION-CONTRACTION
 *
 * By "small" we mean the size of the event is small enough for us to fully assemble across the event,
 * as opposed to large tandem duplications, segmental duplications and/or dispersed duplications.
 */
public final class AssemblyBasedSVDiscoveryTestDataProviderForSimpleSV extends AssemblyBasedSVDiscoveryTestDataProvider {

    public static final class TestDataForSimpleSV extends AssemblyBasedSVDiscoveryTestDataProvider.AssemblyBasedSVDiscoveryTestDataForSimpleChimera {
        public final DistancesBetweenAlignmentsOnRefAndOnRead expectedDistances;


        TestDataForSimpleSV(final AlignmentInterval firstAlignment, final AlignmentInterval secondAlignment, final String evidenceAssemblyContigName, final byte[] evidenceContigSeq,
                            final boolean expectedFirstContigRegionHasLaterReferenceMapping,
                            final SimpleChimera expectedSimpleChimera,
                            final NovelAdjacencyAndAltHaplotype expectedNovelAdjacencyAndAltSeq,
                            final List<SvType> expectedSvTypes,
                            final List<VariantContext> expectedVariantContexts,
                            final DistancesBetweenAlignmentsOnRefAndOnRead expectedDistances,
                            final Class<? extends BreakpointsInference> expectedInferencerClass) {
            super(firstAlignment, secondAlignment, evidenceAssemblyContigName, evidenceContigSeq, expectedFirstContigRegionHasLaterReferenceMapping, expectedSimpleChimera, expectedNovelAdjacencyAndAltSeq, expectedSvTypes, expectedVariantContexts, expectedInferencerClass);
            this.expectedDistances = expectedDistances;
        }

        @Override
        public SAMSequenceDictionary getAppropriateDictionary() {
            return TestUtilsForAssemblyBasedSVDiscovery.b37_seqDict;
        }

        @Override
        public ReferenceMultiSparkSource getAppropriateRef() {
            return TestUtilsForAssemblyBasedSVDiscovery.b37_reference;
        }

        @Override
        public Class<? extends BreakpointsInference> getAppropriateBreakpointInferencer() {
            return expectedInferencerClass;
        }
    }

    final TestDataForSimpleSV forSimpleDeletion_plus;
    private final TestDataForSimpleSV forSimpleDeletion_minus;
    private final TestDataForSimpleSV forDeletionWithHomology_plus;
    final TestDataForSimpleSV forDeletionWithHomology_minus;
    private final TestDataForSimpleSV forSimpleInsertion_plus;
    final TestDataForSimpleSV forSimpleInsertion_minus;
    final TestDataForSimpleSV forLongRangeSubstitution_fudgedDel_plus;
    private final TestDataForSimpleSV forLongRangeSubstitution_fudgedDel_minus;
    private final TestDataForSimpleSV forLongRangeSubstitution_fatIns_plus;
    private final TestDataForSimpleSV forLongRangeSubstitution_fatIns_minus;
    private final TestDataForSimpleSV forLongRangeSubstitution_DelAndIns_plus;
    private final TestDataForSimpleSV forLongRangeSubstitution_DelAndIns_minus;
    final TestDataForSimpleSV forSimpleTanDupContraction_plus;
    private final TestDataForSimpleSV forSimpleTanDupContraction_minus;
    private final TestDataForSimpleSV forSimpleTanDupExpansion_ins_plus;
    final TestDataForSimpleSV forSimpleTanDupExpansion_ins_minus;
    private final TestDataForSimpleSV forSimpleTanDupExpansion_dup_plus;
    private final TestDataForSimpleSV forSimpleTanDupExpansion_dup_minus;
    private final TestDataForSimpleSV forSimpleTanDupExpansionWithNovelIns_ins_plus;
    private final TestDataForSimpleSV forSimpleTanDupExpansionWithNovelIns_ins_minus;
    final TestDataForSimpleSV forSimpleTanDupExpansionWithNovelIns_dup_plus;
    private final TestDataForSimpleSV forSimpleTanDupExpansionWithNovelIns_dup_minus;
    private final TestDataForSimpleSV forComplexTanDup_1to2_pseudoHom_plus;
    final TestDataForSimpleSV forComplexTanDup_1to2_pseudoHom_minus;
    final TestDataForSimpleSV forComplexTanDup_2to1_pseudoHom_plus;
    private final TestDataForSimpleSV forComplexTanDup_2to1_pseudoHom_minus;
    private final TestDataForSimpleSV forComplexTanDup_3to2_noPseudoHom_plus;
    final TestDataForSimpleSV forComplexTanDup_3to2_noPseudoHom_minus;
    final TestDataForSimpleSV forComplexTanDup_2to3_noPseudoHom_plus;
    private final TestDataForSimpleSV forComplexTanDup_2to3_noPseudoHom_minus;
    private final TestDataForSimpleSV forComplexTanDup_1to2_short_pseudoHom_plus;
    private final TestDataForSimpleSV forComplexTanDup_1to2_short_pseudoHom_minus;
    private final TestDataForSimpleSV forComplexTanDup_2to3_short_noPseudoHom_plus;
    private final TestDataForSimpleSV forComplexTanDup_2to3_short_noPseudoHom_minus;

    public AssemblyBasedSVDiscoveryTestDataProviderForSimpleSV() {

        final List<TestDataForSimpleSV> simpleDeletion = forSimpleDeletion();
        forSimpleDeletion_plus = simpleDeletion.get(0);
        forSimpleDeletion_minus = simpleDeletion.get(1);

        final List<TestDataForSimpleSV> deletionWithHomology = forDeletionWithHomology();
        forDeletionWithHomology_plus = deletionWithHomology.get(0);
        forDeletionWithHomology_minus = deletionWithHomology.get(1);

        final List<TestDataForSimpleSV> simpleInsertion = forSimpleInsertion();
        forSimpleInsertion_plus = simpleInsertion.get(0);
        forSimpleInsertion_minus = simpleInsertion.get(1);

        final List<TestDataForSimpleSV> longRangeSubstitution = forLongRangeSubstitution();
        forLongRangeSubstitution_fudgedDel_plus = longRangeSubstitution.get(0);
        forLongRangeSubstitution_fudgedDel_minus = longRangeSubstitution.get(1);
        forLongRangeSubstitution_fatIns_plus = longRangeSubstitution.get(2);
        forLongRangeSubstitution_fatIns_minus = longRangeSubstitution.get(3);
        forLongRangeSubstitution_DelAndIns_plus = longRangeSubstitution.get(4);
        forLongRangeSubstitution_DelAndIns_minus = longRangeSubstitution.get(5);

        final List<TestDataForSimpleSV> simpleTandemDuplicationContraction = forSimpleTandemDuplicationContraction();
        forSimpleTanDupContraction_plus = simpleTandemDuplicationContraction.get(0);
        forSimpleTanDupContraction_minus = simpleTandemDuplicationContraction.get(1);

        final List<TestDataForSimpleSV> simpleTandemDuplicationExpansion = forSimpleTandemDuplicationExpansion();
        forSimpleTanDupExpansion_ins_plus = simpleTandemDuplicationExpansion.get(0);
        forSimpleTanDupExpansion_ins_minus = simpleTandemDuplicationExpansion.get(1);
        forSimpleTanDupExpansion_dup_plus = simpleTandemDuplicationExpansion.get(2);
        forSimpleTanDupExpansion_dup_minus = simpleTandemDuplicationExpansion.get(3);

        final List<TestDataForSimpleSV> simpleTandemDuplicationExpansionWithNovelInsertion = forSimpleTandemDuplicationExpansionWithNovelInsertion();
        forSimpleTanDupExpansionWithNovelIns_ins_plus = simpleTandemDuplicationExpansionWithNovelInsertion.get(0);
        forSimpleTanDupExpansionWithNovelIns_ins_minus = simpleTandemDuplicationExpansionWithNovelInsertion.get(1);
        forSimpleTanDupExpansionWithNovelIns_dup_plus = simpleTandemDuplicationExpansionWithNovelInsertion.get(2);
        forSimpleTanDupExpansionWithNovelIns_dup_minus = simpleTandemDuplicationExpansionWithNovelInsertion.get(3);

        final List<TestDataForSimpleSV> complexTandemDuplication = forComplexTandemDuplication();
        forComplexTanDup_1to2_pseudoHom_plus = complexTandemDuplication.get(0);
        forComplexTanDup_1to2_pseudoHom_minus = complexTandemDuplication.get(1);
        forComplexTanDup_2to1_pseudoHom_plus = complexTandemDuplication.get(2);
        forComplexTanDup_2to1_pseudoHom_minus = complexTandemDuplication.get(3);
        forComplexTanDup_3to2_noPseudoHom_plus = complexTandemDuplication.get(4);
        forComplexTanDup_3to2_noPseudoHom_minus = complexTandemDuplication.get(5);
        forComplexTanDup_2to3_noPseudoHom_plus = complexTandemDuplication.get(6);
        forComplexTanDup_2to3_noPseudoHom_minus = complexTandemDuplication.get(7);

        final List<TestDataForSimpleSV> shortComplexTandemDuplication = forComplexTandemDuplicationIns();
        forComplexTanDup_1to2_short_pseudoHom_plus = shortComplexTandemDuplication.get(0);
        forComplexTanDup_1to2_short_pseudoHom_minus = shortComplexTandemDuplication.get(1);
        forComplexTanDup_2to3_short_noPseudoHom_plus = shortComplexTandemDuplication.get(2);
        forComplexTanDup_2to3_short_noPseudoHom_minus = shortComplexTandemDuplication.get(3);
    }

    @Override
    public List<AssemblyBasedSVDiscoveryTestDataForSimpleChimera> getAllTestData() {
        final List<AssemblyBasedSVDiscoveryTestDataForSimpleChimera> testDataForSimpleSVs = Arrays.asList(
                forSimpleDeletion_plus,
                forSimpleDeletion_minus,
                forSimpleInsertion_plus,
                forSimpleInsertion_minus,
                forLongRangeSubstitution_fudgedDel_plus,
                forLongRangeSubstitution_fudgedDel_minus,
                forLongRangeSubstitution_fatIns_plus,
                forLongRangeSubstitution_fatIns_minus,
                forLongRangeSubstitution_DelAndIns_plus,
                forLongRangeSubstitution_DelAndIns_minus,
                forDeletionWithHomology_plus,
                forDeletionWithHomology_minus,
                forSimpleTanDupContraction_plus,
                forSimpleTanDupContraction_minus,
                forSimpleTanDupExpansion_ins_plus,
                forSimpleTanDupExpansion_ins_minus,
                forSimpleTanDupExpansion_dup_plus,
                forSimpleTanDupExpansion_dup_minus,
                forSimpleTanDupExpansionWithNovelIns_ins_plus,
                forSimpleTanDupExpansionWithNovelIns_ins_minus,
                forSimpleTanDupExpansionWithNovelIns_dup_plus,
                forSimpleTanDupExpansionWithNovelIns_dup_minus,
                forComplexTanDup_1to2_pseudoHom_plus,
                forComplexTanDup_1to2_pseudoHom_minus,
                forComplexTanDup_2to1_pseudoHom_plus,
                forComplexTanDup_2to1_pseudoHom_minus,
                forComplexTanDup_3to2_noPseudoHom_plus,
                forComplexTanDup_3to2_noPseudoHom_minus,
                forComplexTanDup_2to3_noPseudoHom_plus,
                forComplexTanDup_2to3_noPseudoHom_minus,
                forComplexTanDup_1to2_short_pseudoHom_plus,
                forComplexTanDup_1to2_short_pseudoHom_minus,
                forComplexTanDup_2to3_short_noPseudoHom_plus,
                forComplexTanDup_2to3_short_noPseudoHom_minus);
        return Collections.unmodifiableList(testDataForSimpleSVs);
    }


    private static VariantContextBuilder addStandardAttributes(final VariantContextBuilder inputBuilder,
                                                               final int alignLength, final String contigName,
                                                               final String expectedVariantsType, final int end,
                                                               final int svLen, final String altSeq,
                                                               final String homSeq, final String insSeq) {
        VariantContextBuilder builder = inputBuilder
                .attribute(TOTAL_MAPPINGS, 1).attribute(HQ_MAPPINGS, 1)
                .attribute(ALIGN_LENGTHS, alignLength).attribute(MAX_ALIGN_LENGTH, alignLength)
                .attribute(CONTIG_NAMES, contigName)
                .attribute(VCFConstants.END_KEY, end)
                .attribute(SVTYPE, expectedVariantsType)
                .attribute(SVLEN, svLen)
                .attribute(MAPPING_QUALITIES, 60);
        if (!homSeq.isEmpty()) {
            builder.attribute(HOMOLOGY, homSeq).attribute(HOMOLOGY_LENGTH, homSeq.length());
        }
        if (!insSeq.isEmpty()) {
            builder.attribute(INSERTED_SEQUENCE, insSeq).attribute(INSERTED_SEQUENCE_LENGTH, insSeq.length());
        }
        if (!altSeq.isEmpty()) {
            builder.attribute(SEQ_ALT_HAPLOTYPE, altSeq);
        }
        return builder;
    }

    /**
     * 40-'A' + 10-'C'+10-'T' + 40-'G' where the segment 10-'C'+10-'T' is deleted (forward strand representation description).
     *
     * Return a list of two entries for positive and reverse strand representations.
     */
    private static List<TestDataForSimpleSV>
    forSimpleDeletion() {

        final List<TestDataForSimpleSV> result = new ArrayList<>();
        // simple deletion '+' strand representation
        final String leftRefFlank = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence('A', 40);
        final String rightRefFlank = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence('G', 40);
        byte[] contigSeq = (leftRefFlank + rightRefFlank).getBytes();
        String contigName = "simple_del_+";

        final SimpleInterval expectedLeftBreakpoint = new SimpleInterval("21:17000040-17000040");
        final SimpleInterval expectedRightBreakpoint = new SimpleInterval("21:17000060-17000060");
        final BreakpointComplications expectedBreakpointComplications = new BreakpointComplications.SimpleInsDelOrReplacementBreakpointComplications("", "");
        final byte[] expectedAltSeq = EMPTY_BYTE_ARRAY;
        final NovelAdjacencyAndAltHaplotype expectedNovelAdjacencyAndAltHaplotype = new NovelAdjacencyAndAltHaplotype(expectedLeftBreakpoint, expectedRightBreakpoint, NO_SWITCH, expectedBreakpointComplications, SIMPLE_DEL, expectedAltSeq);
        AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("21", 17000001, 17000040), 1 ,40, TextCigarCodec.decode("40M40S"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("21", 17000061, 17000100), 41 ,80, TextCigarCodec.decode("40S40M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        SimpleChimera expectedSimpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, true, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
        DistancesBetweenAlignmentsOnRefAndOnRead expectedDistances = new DistancesBetweenAlignmentsOnRefAndOnRead(20, 0, 17000040, 17000061, 40, 41);
        final List<SvType> expectedSVTypes = Collections.singletonList(makeDeletionType(new SimpleInterval("21:17000040-17000060"), Allele.create("G", true),false));
        final List<VariantContext> expectedVariants = Collections.singletonList(
                addStandardAttributes(makeDeletion(new SimpleInterval("21:17000040-17000060"), Allele.create("G", true), false),
                        40, contigName, SimpleSVType.SupportedType.DEL.name(), 17000060, -20, "", "", "").make());
        result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeq, false, expectedSimpleChimera, expectedNovelAdjacencyAndAltHaplotype, expectedSVTypes, expectedVariants, expectedDistances, BreakpointsInference.SimpleInsertionDeletionBreakpointsInference.class));

        // simple deletion '-' strand representation
        SequenceUtil.reverseComplement(leftRefFlank);
        SequenceUtil.reverseComplement(rightRefFlank);
        contigSeq = (rightRefFlank + leftRefFlank).getBytes();
        contigName = "simple_del_-";
        region1 = new AlignmentInterval(new SimpleInterval("21", 17000061, 17000100), 1 ,40, TextCigarCodec.decode("40M40S"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        region2 = new AlignmentInterval(new SimpleInterval("21", 17000001, 17000040), 41 ,80, TextCigarCodec.decode("40S40M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        expectedSimpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, false, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
        expectedDistances = new DistancesBetweenAlignmentsOnRefAndOnRead(20, 0, 17000040, 17000061, 40, 41);
        result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeq, true, expectedSimpleChimera, expectedNovelAdjacencyAndAltHaplotype, expectedSVTypes, expectedVariants, expectedDistances, BreakpointsInference.SimpleInsertionDeletionBreakpointsInference.class));

        return result;
    }

    /**
     * 40-'C' + 'ATCG' + 34 bases of unique sequence + 'ATCG' + 40-'T' is shrunk to 40-'C' + 'ATCG' + 40-'T' (forward strand representation)
     * Return a list of two entries for positive and reverse strand representations.
     */
    private static List<TestDataForSimpleSV>
    forDeletionWithHomology() {

        final List<TestDataForSimpleSV> result = new ArrayList<>();

        // simple deletion with homology '+' strand representation
        final String leftRefFlank = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence('C', 40);
        final String rightRefFlank = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence('T', 40);
        final String homology = "ATCG";
        byte[] contigSeq = (leftRefFlank + homology + rightRefFlank).getBytes();
        String contigName = "simple_del_with_hom_+";

        final SimpleInterval expectedLeftBreakpoint = new SimpleInterval("21:17000040-17000040");
        final SimpleInterval expectedRightBreakpoint = new SimpleInterval("21:17000078-17000078");
        final BreakpointComplications expectedBreakpointComplications = new BreakpointComplications.SimpleInsDelOrReplacementBreakpointComplications(homology, "");
        final byte[] expectedAltSeq = EMPTY_BYTE_ARRAY;
        final NovelAdjacencyAndAltHaplotype expectedNovelAdjacencyAndAltHaplotype = new NovelAdjacencyAndAltHaplotype(expectedLeftBreakpoint, expectedRightBreakpoint, NO_SWITCH, expectedBreakpointComplications, SIMPLE_DEL, expectedAltSeq);
        AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("21", 17000001, 17000044), 1 ,44, TextCigarCodec.decode("44M40S"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("21", 17000079, 17000122), 41 ,84, TextCigarCodec.decode("40S44M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        SimpleChimera expectedSimpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, true, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
        DistancesBetweenAlignmentsOnRefAndOnRead expectedDistances = new DistancesBetweenAlignmentsOnRefAndOnRead(34, -4, 17000044, 17000079, 44, 41);
        final List<SvType> expectedSVTypes = Collections.singletonList(makeDeletionType(new SimpleInterval("21:17000040-17000078"), Allele.create("G", true),false));
        final List<VariantContext> expectedVariants = Collections.singletonList(
                addStandardAttributes(makeDeletion(new SimpleInterval("21:17000040-17000078"), Allele.create("G", true), false),
                        40, contigName, SimpleSVType.SupportedType.DEL.name(), 17000078, -38, "", homology, "").make());
        result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeq, false, expectedSimpleChimera, expectedNovelAdjacencyAndAltHaplotype, expectedSVTypes, expectedVariants, expectedDistances, BreakpointsInference.SimpleInsertionDeletionBreakpointsInference.class));

        // simple deletion with homology '-' strand representation
        contigSeq = (SequenceUtil.reverseComplement(rightRefFlank) + SequenceUtil.reverseComplement(homology) + SequenceUtil.reverseComplement(leftRefFlank)).getBytes();
        contigName = "simple_del_with_hom_-";
        region1 = new AlignmentInterval(new SimpleInterval("21", 17000079, 17000122), 1 ,44, TextCigarCodec.decode("44M40S"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        region2 = new AlignmentInterval(new SimpleInterval("21", 17000001, 17000044), 41 ,84, TextCigarCodec.decode("40S44M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        expectedSimpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, false, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
        expectedDistances = new DistancesBetweenAlignmentsOnRefAndOnRead(34, -4, 17000044, 17000079, 44, 41);
        result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeq, true, expectedSimpleChimera, expectedNovelAdjacencyAndAltHaplotype, expectedSVTypes, expectedVariants, expectedDistances,  BreakpointsInference.SimpleInsertionDeletionBreakpointsInference.class));

        return result;
    }

    /**
     * 100-'A' + 100-'T' and a 50 bases of 'C' is inserted at the A->T junction point (forward strand description)
     * Return a list of two entries for positive and reverse strand representations.
     */
    private static List<TestDataForSimpleSV>
    forSimpleInsertion() {
        final List<TestDataForSimpleSV> result = new ArrayList<>();

        // simple insertion '+' strand representation
        final String leftRefFlank = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence('A', 100);
        final String insertedSeq  = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence('C', 50);
        final String rightRefFlank = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence('T', 100);
        byte[] contigSeq = (leftRefFlank + insertedSeq + rightRefFlank).getBytes();
        String contigName = "simple_ins_+";

        final SimpleInterval expectedLeftBreakpoint = new SimpleInterval("21:17000100-17000100");
        final SimpleInterval expectedRightBreakpoint = new SimpleInterval("21:17000100-17000100");
        final BreakpointComplications expectedBreakpointComplications = new BreakpointComplications.SimpleInsDelOrReplacementBreakpointComplications("", insertedSeq);
        final byte[] expectedAltSeq = insertedSeq.getBytes();
        final NovelAdjacencyAndAltHaplotype expectedNovelAdjacencyAndAltHaplotype = new NovelAdjacencyAndAltHaplotype(expectedLeftBreakpoint, expectedRightBreakpoint, NO_SWITCH, expectedBreakpointComplications, SIMPLE_INS, expectedAltSeq);
        AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("21", 17000001, 17000100), 1 ,100, TextCigarCodec.decode("100M100S"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("21", 17000101, 17000200), 151 ,250, TextCigarCodec.decode("100S100M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        SimpleChimera expectedSimpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, true, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
        DistancesBetweenAlignmentsOnRefAndOnRead expectedDistances = new DistancesBetweenAlignmentsOnRefAndOnRead(0, 50, 17000100, 17000101, 100, 151);
        final List<SvType> expectedSVTypes = Collections.singletonList(makeInsertionType(new SimpleInterval("21:17000100-17000100"), Allele.create("G", true),50));
        final List<VariantContext> expectedVariants = Collections.singletonList(
                addStandardAttributes(makeInsertion("21", 17000100, 17000100, 50, Allele.create("G", true)),
                        100, contigName, SimpleSVType.SupportedType.INS.name(), 17000100, 50, insertedSeq, "", insertedSeq).make());
        result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeq, false, expectedSimpleChimera, expectedNovelAdjacencyAndAltHaplotype, expectedSVTypes, expectedVariants, expectedDistances, BreakpointsInference.SimpleInsertionDeletionBreakpointsInference.class));

        // simple insertion '-' strand representation
        contigSeq = (SequenceUtil.reverseComplement(rightRefFlank) + SequenceUtil.reverseComplement(insertedSeq) + SequenceUtil.reverseComplement(leftRefFlank)).getBytes();
        contigName = "simple_ins_-";
        region1 = new AlignmentInterval(new SimpleInterval("21", 17000101, 17000200), 1 ,100, TextCigarCodec.decode("100M100S"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        region2 = new AlignmentInterval(new SimpleInterval("21", 17000001, 17000100), 151 ,250, TextCigarCodec.decode("100S100M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        expectedDistances = new DistancesBetweenAlignmentsOnRefAndOnRead(0, 50, 17000100, 17000101, 100, 151);
        expectedSimpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, false, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
        result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeq, true, expectedSimpleChimera, expectedNovelAdjacencyAndAltHaplotype, expectedSVTypes, expectedVariants, expectedDistances, BreakpointsInference.SimpleInsertionDeletionBreakpointsInference.class));

        return result;
    }

    /**
     * fudged deletion case:
     * 100-'A' + 100-'G' where the middle 30-'A'+30-'G' is substituted with 10-'C' (forward strand representation)
     * fat insertion case:
     * 50-'A' + 50-'G' where the middle 10-'A'+10-'G' is substituted with 60-'C' (forward strand representation)
     * Two linked expectedVariantss case:
     * 100-'A' + 100-'G' where the middle 30-'A'+30-'G' is substituted with 55-'C' (forward strand representation)
     */
    private static List<TestDataForSimpleSV>
    forLongRangeSubstitution() {

        final List<TestDataForSimpleSV> result = new ArrayList<>();

        {//fudged deletion case
            // '+' strand representation
            final String leftRefFlank = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence('A', 70);
            final String rightRefFlank = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence('G', 70);
            final String substitution = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence('C', 10);
            byte[] contigSeq = (leftRefFlank + substitution + rightRefFlank).getBytes();
            String contigName = "rpl_fudged_del_+";
            final String insSeqString = substitution;

            final SimpleInterval expectedLeftBreakpoint = new SimpleInterval("21:17000070-17000070");
            final SimpleInterval expectedRightBreakpoint = new SimpleInterval("21:17000130-17000130");
            final BreakpointComplications expectedBreakpointComplications = new BreakpointComplications.SimpleInsDelOrReplacementBreakpointComplications("", insSeqString);
            final byte[] expectedAltSeq = insSeqString.getBytes();
            final NovelAdjacencyAndAltHaplotype expectedNovelAdjacencyAndAltHaplotype = new NovelAdjacencyAndAltHaplotype(expectedLeftBreakpoint, expectedRightBreakpoint, NO_SWITCH, expectedBreakpointComplications, RPL, expectedAltSeq);
            AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("21", 17000001, 17000070), 1, 70, TextCigarCodec.decode("70M80S"), true, 60, 0, 70, ContigAlignmentsModifier.AlnModType.NONE);
            AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("21", 17000131, 17000200), 81, 150, TextCigarCodec.decode("80S70M"), true, 60, 0, 70, ContigAlignmentsModifier.AlnModType.NONE);
            SimpleChimera expectedSimpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, true, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
            DistancesBetweenAlignmentsOnRefAndOnRead expectedDistances = new DistancesBetweenAlignmentsOnRefAndOnRead(60, 10, 17000070, 17000131, 70, 81);
            final List<SvType> expectedSVTypes = Collections.singletonList(makeDeletionType(new SimpleInterval("21:17000070-17000130"), Allele.create("T", true), false));
            final List<VariantContext> expectedVariants = Collections.singletonList(
                    addStandardAttributes(makeDeletion(new SimpleInterval("21:17000070-17000130"), Allele.create("T", true), false),
                            70, contigName, SimpleSVType.SupportedType.DEL.name(), 17000130, -60, insSeqString, "", insSeqString).make());
            result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeq, false, expectedSimpleChimera, expectedNovelAdjacencyAndAltHaplotype, expectedSVTypes, expectedVariants, expectedDistances, BreakpointsInference.SimpleInsertionDeletionBreakpointsInference.class));

            // '-' strand representation
            contigSeq = (SequenceUtil.reverseComplement(rightRefFlank) + SequenceUtil.reverseComplement(substitution) + SequenceUtil.reverseComplement(leftRefFlank)).getBytes();
            contigName = "rpl_fudged_del_-";
            region1 = new AlignmentInterval(new SimpleInterval("21", 17000131, 17000200), 1, 70, TextCigarCodec.decode("70M80S"), false, 60, 0, 70, ContigAlignmentsModifier.AlnModType.NONE);
            region2 = new AlignmentInterval(new SimpleInterval("21", 17000001, 17000070), 81, 150, TextCigarCodec.decode("80S70M"), false, 60, 0, 70, ContigAlignmentsModifier.AlnModType.NONE);
            expectedSimpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, false, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
            expectedDistances = new DistancesBetweenAlignmentsOnRefAndOnRead(60, 10, 17000070, 17000131, 70, 81);
            result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeq, true, expectedSimpleChimera, expectedNovelAdjacencyAndAltHaplotype, expectedSVTypes, expectedVariants, expectedDistances, BreakpointsInference.SimpleInsertionDeletionBreakpointsInference.class));
        }

        {//fat insertion case
            // '+' strand representation
            final String leftRefFlank = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence('A', 40);
            final String rightRefFlank = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence('G', 40);
            final String substitution = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence('C', 60);
            byte[] contigSeq = (leftRefFlank + substitution + rightRefFlank).getBytes();
            String contigName = "rpl_fat_ins_+";
            final String insSeqString = substitution;

            final SimpleInterval expectedLeftBreakpoint = new SimpleInterval("21:17000040-17000040");
            final SimpleInterval expectedRightBreakpoint = new SimpleInterval("21:17000060-17000060");
            final BreakpointComplications expectedBreakpointComplications = new BreakpointComplications.SimpleInsDelOrReplacementBreakpointComplications("", insSeqString);
            final byte[] expectedAltSeq = insSeqString.getBytes();
            final NovelAdjacencyAndAltHaplotype expectedNovelAdjacencyAndAltHaplotype = new NovelAdjacencyAndAltHaplotype(expectedLeftBreakpoint, expectedRightBreakpoint, NO_SWITCH, expectedBreakpointComplications, RPL, expectedAltSeq);
            AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("21", 17000001, 17000040), 1, 40, TextCigarCodec.decode("40M100S"), true, 60, 0, 60, ContigAlignmentsModifier.AlnModType.NONE);
            AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("21", 17000061, 17000100), 101, 140, TextCigarCodec.decode("100S40M"), true, 60, 0, 60, ContigAlignmentsModifier.AlnModType.NONE);
            SimpleChimera expectedSimpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, true, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
            DistancesBetweenAlignmentsOnRefAndOnRead expectedDistances = new DistancesBetweenAlignmentsOnRefAndOnRead(20, 60, 17000040, 17000061, 40, 101);
            final List<SvType> expectedSVTypes = Collections.singletonList(makeInsertionType(new SimpleInterval("21:17000040-17000060"), Allele.create("G", true), 60));
            final List<VariantContext> expectedVariants = Collections.singletonList(
                    addStandardAttributes(makeInsertion("21", 17000040, 17000060, 60, Allele.create("G", true)),
                            40, contigName, SimpleSVType.SupportedType.INS.name(), 17000060, 60, insSeqString, "", insSeqString).make());
            result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeq, false, expectedSimpleChimera, expectedNovelAdjacencyAndAltHaplotype, expectedSVTypes, expectedVariants, expectedDistances, BreakpointsInference.SimpleInsertionDeletionBreakpointsInference.class));

            // '-' strand representation
            contigSeq = (SequenceUtil.reverseComplement(rightRefFlank) + SequenceUtil.reverseComplement(substitution) + SequenceUtil.reverseComplement(leftRefFlank)).getBytes();
            contigName = "rpl_fat_ins_-";
            region1 = new AlignmentInterval(new SimpleInterval("21", 17000061, 17000100), 1, 40, TextCigarCodec.decode("40M100S"), false, 60, 0, 60, ContigAlignmentsModifier.AlnModType.NONE);
            region2 = new AlignmentInterval(new SimpleInterval("21", 17000001, 17000040), 101, 140, TextCigarCodec.decode("100S40M"), false, 60, 0, 60, ContigAlignmentsModifier.AlnModType.NONE);
            expectedSimpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, false, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
            expectedDistances = new DistancesBetweenAlignmentsOnRefAndOnRead(20, 60, 17000040, 17000061, 40, 101);
            result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeq, true, expectedSimpleChimera, expectedNovelAdjacencyAndAltHaplotype, expectedSVTypes, expectedVariants, expectedDistances, BreakpointsInference.SimpleInsertionDeletionBreakpointsInference.class));
        }

        {//two linked expectedVariantss case
            // '+' strand representation
            final String leftRefFlank = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence('A', 70);
            final String rightRefFlank = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence('G', 70);
            final String substitution = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence('C', 55);
            byte[] contigSeq = (leftRefFlank + substitution + rightRefFlank).getBytes();
            String contigName = "rpl_linked_del_ins_+";
            final String insSeqString = substitution;

            final SimpleInterval expectedLeftBreakpoint = new SimpleInterval("21:17000070-17000070");
            final SimpleInterval expectedRightBreakpoint = new SimpleInterval("21:17000130-17000130");
            final BreakpointComplications expectedBreakpointComplications = new BreakpointComplications.SimpleInsDelOrReplacementBreakpointComplications("", insSeqString);
            final byte[] expectedAltSeq = insSeqString.getBytes();
            final NovelAdjacencyAndAltHaplotype expectedNovelAdjacencyAndAltHaplotype = new NovelAdjacencyAndAltHaplotype(expectedLeftBreakpoint, expectedRightBreakpoint, NO_SWITCH, expectedBreakpointComplications, RPL, expectedAltSeq);
            AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("21", 17000001, 17000070), 1, 70, TextCigarCodec.decode("70M125S"), true, 60, 0, 70, ContigAlignmentsModifier.AlnModType.NONE);
            AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("21", 17000131, 17000200), 126, 195, TextCigarCodec.decode("125S70M"), true, 60, 0, 70, ContigAlignmentsModifier.AlnModType.NONE);
            SimpleChimera expectedSimpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, true, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
            DistancesBetweenAlignmentsOnRefAndOnRead expectedDistances = new DistancesBetweenAlignmentsOnRefAndOnRead(60, 55, 17000070, 17000131, 70, 126);
            final List<SvType> expectedSVTypes = Arrays.asList(
                    makeDeletionType(new SimpleInterval("21:17000070-17000130"), Allele.create("T", true), false),
                    makeInsertionType(new SimpleInterval("21:17000070-17000070"), Allele.create("T", true),55)
            );
            final List<VariantContext> expectedVariants = Arrays.asList(
                    addStandardAttributes(makeDeletion(new SimpleInterval("21:17000070-17000130"), Allele.create("T", true), false), 70, contigName, SimpleSVType.SupportedType.DEL.name(), 17000130, -60, "", "", "")
                            .attribute(LINK, SimpleSVType.SupportedType.INS.name() + INTERVAL_VARIANT_ID_FIELD_SEPARATOR + SvType.makeLocationString("21", 17000070, "21", 17000070)).make(),
                    addStandardAttributes(makeInsertion("21", 17000070, 17000070, 55, Allele.create("T", true)), 70, contigName, SimpleSVType.SupportedType.INS.name(), 17000070, 55, insSeqString, "", insSeqString)
                            .attribute(LINK, SimpleSVType.SupportedType.DEL.name() + INTERVAL_VARIANT_ID_FIELD_SEPARATOR + SvType.makeLocationString("21", 17000070, "21", 17000130)).make()
            );
            result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeq, false, expectedSimpleChimera, expectedNovelAdjacencyAndAltHaplotype, expectedSVTypes, expectedVariants, expectedDistances, BreakpointsInference.SimpleInsertionDeletionBreakpointsInference.class));

            // '-' strand representation
            contigSeq = (SequenceUtil.reverseComplement(rightRefFlank) + SequenceUtil.reverseComplement(substitution) + SequenceUtil.reverseComplement(leftRefFlank)).getBytes();
            contigName = "rpl_linked_del_ins_-";
            region1 = new AlignmentInterval(new SimpleInterval("21", 17000131, 17000200), 1, 70, TextCigarCodec.decode("70M125S"), false, 60, 0, 70, ContigAlignmentsModifier.AlnModType.NONE);
            region2 = new AlignmentInterval(new SimpleInterval("21", 17000001, 17000070), 126, 195, TextCigarCodec.decode("125S70M"), false, 60, 0, 70, ContigAlignmentsModifier.AlnModType.NONE);
            expectedSimpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, false, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
            expectedDistances = new DistancesBetweenAlignmentsOnRefAndOnRead(60, 55, 17000070, 17000131, 70, 126);
            result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeq, true, expectedSimpleChimera, expectedNovelAdjacencyAndAltHaplotype, expectedSVTypes, expectedVariants, expectedDistances, BreakpointsInference.SimpleInsertionDeletionBreakpointsInference.class));
        }

        return result;
    }

    /**
     * 40-'A' + 20-'C' + 40-'G' is shrunk to 40-'A' + 10-'C' + 40-'G' (forward strand representation)
     * Return a list of two entries for positive and reverse strand representations.
     */
    private static List<TestDataForSimpleSV>
    forSimpleTandemDuplicationContraction() {

        final List<TestDataForSimpleSV> result = new ArrayList<>();

        // simple tandem duplication contraction '+' strand representation
        final String leftRefFlank = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence('A', 40);
        final String rightRefFlank = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence('G', 40);
        final String doubleDup = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence('C', 20);
        byte[] contigSeq = (leftRefFlank + doubleDup.substring(0, 10) + rightRefFlank).getBytes();
        String contigName = "simple_del_dup_contraction_+";

        final SimpleInterval expectedLeftBreakpoint = new SimpleInterval("21:17000040-17000040");
        final SimpleInterval expectedRightBreakpoint = new SimpleInterval("21:17000050-17000050");
        final BreakpointComplications expectedBreakpointComplications = new BreakpointComplications.SmallDuplicationWithPreciseDupRangeBreakpointComplications("CCCCCCCCCC", "",
                new SimpleInterval("21:17000041-17000050"), 2, 1, Arrays.asList(Strand.POSITIVE, Strand.POSITIVE), Collections.singletonList(Strand.POSITIVE), Collections.emptyList());
        byte[] expectedAltSeq = EMPTY_BYTE_ARRAY;
        final NovelAdjacencyAndAltHaplotype expectedNovelAdjacencyAndAltHaplotype = new NovelAdjacencyAndAltHaplotype(expectedLeftBreakpoint, expectedRightBreakpoint, NO_SWITCH, expectedBreakpointComplications, DEL_DUP_CONTRACTION, expectedAltSeq);
        AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("21", 17000001, 17000050), 1 ,50, TextCigarCodec.decode("50M40S"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("21", 17000051, 17000100), 41 ,90, TextCigarCodec.decode("40S50M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        SimpleChimera expectedSimpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, true, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
        DistancesBetweenAlignmentsOnRefAndOnRead expectedDistances = new DistancesBetweenAlignmentsOnRefAndOnRead(0, -10, 17000050, 17000051, 50, 41);
        final List<SvType> expectedSVTypes = Collections.singletonList(makeDeletionType(new SimpleInterval("21:17000040-17000050"), Allele.create("G", true), true));
        final List<VariantContext> expectedVariants = Collections.singletonList(
                addStandardAttributes(makeDeletion(new SimpleInterval("21:17000040-17000050"), Allele.create("G", true), true), 50 - 10, contigName, SimpleSVType.SupportedType.DEL.name(), 17000050, -10, "", "CCCCCCCCCC", "")
                        .attribute(DUP_REPEAT_UNIT_REF_SPAN, "21:17000041-17000050").attribute(DUPLICATION_NUMBERS, "2,1").attribute(DUP_ORIENTATIONS, "+").attribute(DUP_TAN_CONTRACTION_STRING, "").make()
        );
        result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeq, false, expectedSimpleChimera, expectedNovelAdjacencyAndAltHaplotype, expectedSVTypes, expectedVariants, expectedDistances, BreakpointsInference.SmallDuplicationWithPreciseDupRangeBreakpointsInference.class));

        // simple tandem duplication contraction '-' strand representation
        contigSeq = (SequenceUtil.reverseComplement(rightRefFlank) + SequenceUtil.reverseComplement(doubleDup).substring(0, 10) + SequenceUtil.reverseComplement(leftRefFlank)).getBytes();
        contigName = "simple_del_dup_contraction_-";
        region1 = new AlignmentInterval(new SimpleInterval("21", 17000051, 17000100), 1 ,50, TextCigarCodec.decode("50M40S"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        region2 = new AlignmentInterval(new SimpleInterval("21", 17000001, 17000050), 41 ,90, TextCigarCodec.decode("40S50M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        expectedSimpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, false, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
        expectedDistances = new DistancesBetweenAlignmentsOnRefAndOnRead(0, -10, 17000050, 17000051, 50, 41);
        result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeq, true, expectedSimpleChimera, expectedNovelAdjacencyAndAltHaplotype, expectedSVTypes, expectedVariants, expectedDistances, BreakpointsInference.SmallDuplicationWithPreciseDupRangeBreakpointsInference.class));

        return result;
    }

    /**
     * case that will be called as insertion
     * 40-'A' + 10-'C' + 40-'G' is expanded to 40-'A' + 20-'C' + 40-'G' (forward strand representation)
     *
     * case that will be called as duplication
     * 40-'A' + 55-'C' + 40-'G' is expanded to 40-'A' + 110-'C' + 40-'G' (forward strand representation)
     */
    private static List<TestDataForSimpleSV>
    forSimpleTandemDuplicationExpansion() {

        final List<TestDataForSimpleSV> result = new ArrayList<>();

        {// insertion case
            // '+' strand representation
            final String leftRefFlank = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence('A', 40);
            final String rightRefFlank = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence('G', 40);
            final String doubleDup = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence('C', 20);
            byte[] contigSeq = (leftRefFlank + doubleDup + rightRefFlank).getBytes();
            String contigName = "simple_dup_exp_too_small_+";

            final SimpleInterval expectedLeftBreakpoint = new SimpleInterval("21:17000040-17000040");
            final SimpleInterval expectedRightBreakpoint = new SimpleInterval("21:17000040-17000040");
            final BreakpointComplications expectedBreakpointComplications = new BreakpointComplications.SmallDuplicationWithPreciseDupRangeBreakpointComplications("", "",
                    new SimpleInterval("21:17000041-17000050"), 1, 2, Collections.singletonList(Strand.POSITIVE), Arrays.asList(Strand.POSITIVE, Strand.POSITIVE), Arrays.asList("10M", "10M"));
            final byte[] expectedAltSeq = doubleDup.getBytes();
            final NovelAdjacencyAndAltHaplotype expectedNovelAdjacencyAndAltHaplotype = new NovelAdjacencyAndAltHaplotype(expectedLeftBreakpoint, expectedRightBreakpoint, NO_SWITCH, expectedBreakpointComplications, SMALL_DUP_EXPANSION, expectedAltSeq);
            AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("21", 17000001, 17000050), 1 ,50, TextCigarCodec.decode("50M50S"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
            AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("21", 17000041, 17000090), 51 ,100, TextCigarCodec.decode("50S50M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
            SimpleChimera expectedSimpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, true, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
            DistancesBetweenAlignmentsOnRefAndOnRead expectedDistances = new DistancesBetweenAlignmentsOnRefAndOnRead(-10, 0, 17000050, 17000041, 50, 51);
            final List<SvType> expectedSVTypes = Collections.singletonList(makeInsertionType(new SimpleInterval("21:17000040-17000040"), Allele.create("G", true), 10));
            final List<VariantContext> expectedVariants = Collections.singletonList(
                    addStandardAttributes(makeInsertion("21", 17000040, 17000040, 10, Allele.create("G", true)), 50, contigName, SimpleSVType.SupportedType.INS.name(), 17000040, 10, StringUtil.bytesToString(expectedAltSeq), "", "")
                            .attribute(DUP_REPEAT_UNIT_REF_SPAN, "21:17000041-17000050").attribute(DUP_ORIENTATIONS, "++").attribute(DUP_SEQ_CIGARS, "10M,10M").attribute(DUPLICATION_NUMBERS, "1,2").attribute(DUP_TAN_EXPANSION_STRING, "").make()
            );
            result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeq, false, expectedSimpleChimera, expectedNovelAdjacencyAndAltHaplotype, expectedSVTypes, expectedVariants, expectedDistances, BreakpointsInference.SmallDuplicationWithPreciseDupRangeBreakpointsInference.class));

            // '-' strand representation
            contigSeq = (SequenceUtil.reverseComplement(rightRefFlank) + SequenceUtil.reverseComplement(doubleDup) + SequenceUtil.reverseComplement(leftRefFlank)).getBytes();
            contigName = "simple_dup_exp_too_small_-";
            region1 = new AlignmentInterval(new SimpleInterval("21", 17000041, 17000090), 1 ,50, TextCigarCodec.decode("50M50S"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
            region2 = new AlignmentInterval(new SimpleInterval("21", 17000001, 17000050), 51 ,100, TextCigarCodec.decode("50S50M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
            expectedSimpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, false, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
            expectedDistances = new DistancesBetweenAlignmentsOnRefAndOnRead(-10, 0, 17000050, 17000041, 50, 51);
            result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeq, true, expectedSimpleChimera, expectedNovelAdjacencyAndAltHaplotype, expectedSVTypes, expectedVariants, expectedDistances, BreakpointsInference.SmallDuplicationWithPreciseDupRangeBreakpointsInference.class));
        }

        {// duplication case
            // '+' strand representation
            final String leftRefFlank = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence('A', 40);
            final String rightRefFlank = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence('G', 40);
            final String doubleDup = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence('C', 110);
            byte[] contigSeq = (leftRefFlank + doubleDup + rightRefFlank).getBytes();
            String contigName = "simple_dup_exp_+";

            final SimpleInterval expectedLeftBreakpoint = new SimpleInterval("21:17000040-17000040");
            final SimpleInterval expectedRightBreakpoint = new SimpleInterval("21:17000040-17000040");
            final BreakpointComplications expectedBreakpointComplications = new BreakpointComplications.SmallDuplicationWithPreciseDupRangeBreakpointComplications("", "",
                    new SimpleInterval("21:17000041-17000095"), 1, 2, Collections.singletonList(Strand.POSITIVE), Arrays.asList(Strand.POSITIVE, Strand.POSITIVE), Arrays.asList("55M", "55M"));
            final byte[] expectedAltSeq = doubleDup.getBytes();
            final NovelAdjacencyAndAltHaplotype expectedNovelAdjacencyAndAltHaplotype = new NovelAdjacencyAndAltHaplotype(expectedLeftBreakpoint, expectedRightBreakpoint, NO_SWITCH, expectedBreakpointComplications, SMALL_DUP_EXPANSION, expectedAltSeq);
            AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("21", 17000001, 17000095), 1 ,95, TextCigarCodec.decode("95M95S"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
            AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("21", 17000041, 17000135), 96 ,190, TextCigarCodec.decode("95S95M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
            SimpleChimera expectedSimpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, true, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
            DistancesBetweenAlignmentsOnRefAndOnRead expectedDistances = new DistancesBetweenAlignmentsOnRefAndOnRead(-55, 0, 17000095, 17000041, 95, 96);
            final List<SvType> expectedSVTypes = Collections.singletonList(makeTandemDuplicationType(new SimpleInterval("21:17000041-17000095"), Allele.create("G", true), 55));
            final List<VariantContext> expectedVariants = Collections.singletonList(
                    addStandardAttributes(makeTandemDuplication(new SimpleInterval("21:17000041-17000095"), Allele.create("G", true), 55), 95, contigName, SimpleSVType.SupportedType.DUP.name(), 17000040, 55, doubleDup, "", "")
                            .attribute(DUP_REPEAT_UNIT_REF_SPAN, "21:17000041-17000095").attribute(DUP_ORIENTATIONS, "++").attribute(DUP_SEQ_CIGARS, "55M,55M").attribute(DUPLICATION_NUMBERS, "1,2").attribute(DUP_TAN_EXPANSION_STRING, "").make()
            );
            result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeq, false, expectedSimpleChimera, expectedNovelAdjacencyAndAltHaplotype, expectedSVTypes, expectedVariants, expectedDistances, BreakpointsInference.SmallDuplicationWithPreciseDupRangeBreakpointsInference.class));

            // '-' strand representation
            contigSeq = (SequenceUtil.reverseComplement(rightRefFlank) + SequenceUtil.reverseComplement(doubleDup) + SequenceUtil.reverseComplement(leftRefFlank)).getBytes();
            contigName = "simple_dup_exp_-";
            region1 = new AlignmentInterval(new SimpleInterval("21", 17000041, 17000135), 1 ,95, TextCigarCodec.decode("95M95S"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
            region2 = new AlignmentInterval(new SimpleInterval("21", 17000001, 17000095), 96 ,190, TextCigarCodec.decode("95S95M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
            expectedSimpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, false, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
            expectedDistances = new DistancesBetweenAlignmentsOnRefAndOnRead(-55, 0, 17000095, 17000041, 95, 96);
            result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeq, true, expectedSimpleChimera, expectedNovelAdjacencyAndAltHaplotype, expectedSVTypes, expectedVariants, expectedDistances, BreakpointsInference.SmallDuplicationWithPreciseDupRangeBreakpointsInference.class));
        }

        return result;
    }

    /**
     * Real event, which will be output as INS (but the event was actually from a hg38 sample, but doesn't matter)
     * repeat:     chr21:26849022-26849037
     * repeat sequence: CCGGGAAATGCTTTTT
     * insertedSequenceForwardStrandRep: TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGAGACGGAGTCTCGCTGTCGCCCAGGCTGGAGTGCAGTGGCGCAATCTCGGCTCACTGCAGGCTCCGCCCCCTGGGGTTCACGCCATTCTCCTGCCTCAGCCTCCCGAGTAGCTGGGACTACAGGCGCCCGCCACCTCGCCCGGCTAATTTTTTGTATTTTTAGTAGAGACGGGGTTTCACCGTGTTAGCCAGGATGGTCTCGATCTCCTGACCTCGTGATCCGCCCGCCTCGGCCTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACCGCGCCCGGC
     *
     * Real event, which will be output as DUP
     * leftFlank:  chr21:25297101-25297163
     * repeat:     chr21:25297164-25297252
     * rightFlank: chr21:25297253-25297300
     * GTTAGTAGATATTCTAGCTGACTCAGTTCAGTGTTGCTATGATTAAACAAGAGTGAGTTCCCT
     * AAAAGTAAATGTTATAAGAAATCTTAAGTATTATTTTCTTATGTTTCTAGCCTAATAAAGTGCTTTTATTAAAGCACTTTATTTAAAGG
     * CATTATTGATATTTCATTATGTTCAACAGATGGAGTTAATGTGAATGT
     *
     * insertedSequenceForwardStrandRep: CTCTCTCTCT
     */
    private static List<TestDataForSimpleSV>
    forSimpleTandemDuplicationExpansionWithNovelInsertion() {

        final List<TestDataForSimpleSV> result = new ArrayList<>();

        {// expectedBreakpointComplications slightly different due to small difference in assembly contig sequence
            AlignmentInterval region1 = TestUtilsForAssemblyBasedSVDiscovery.fromSAMRecordString("asm029081:tig00000\t0\t21\t26847644\t60\t1394M1675S\t*\t0\t0\tTATGGTGACAACAGTTACAGCCTCAGAGTGGCTTTGAGATGGAGATTTGAGATGGAGATTTTGGACACTCATAATCTCCATCTCAAGGAGAGGTGGCCCAGAGGCTGAATGAAGATGAGTGAAGGTAGATGTGATTCCCTTAAATTGGGGCAAAAAGGGACAAAAACAGCAGAAAATCTGTATCTTTAAAGACATATGTAATGTATTTCAGTCTATCAACTTCTCTACATAAACTTTAGCTTTTAAAAATATGTTAACGTAAGTTTGACCTTTAGTGTGTTTCTACCTGCAGGGTATTCTTATTGGAGGTTTGTTTAAAAGCATACATTTCTGATCTTGAATGGGTTACTACAAATCCATTATAATTGTTTCATATTTCATGTTGCAGATACAAGTAGGGTTGAAAAAACAGTGAGTTAAAGGCAAAAGGATGGCCGGGAACATGGCTTTTTTATTCTCTGGGTTTCTATCCAGATTTCTGTTCTTTTGCATAATGACTCCAATCTGTTGTGCACCTGTAGTTCTGGGAAATGATTCTTTTTTAATCGCTTCAACAGAGACATGGATGTTGGAGTTGCCAACTACTAAGCTGAAAAACTCCATCTATGCTCAGAAGAACATTTAATCCACTTACTTTTTCTCTTTTATTTAAAGATTAGCACTCATCAGGCATTTGTGGTAATATGCAAATATATACATAGGACATATATGTATATTTATAAGCAAAATGTGAATTGGAAAAACATTTGAATGTAGAAACAAGACCACAGGAGTAAATTTGTACAAGGCACTAGTAAAAGTGACATGTAATATGGGGTTCTTGTAGTGAGTTTCATAATCCAATTTTTGCTCCTTGATTTGAATGGGCACCCAAAATAACACATGCTATCCTAATCCCTACTCCCCATATTTTGGGTTTTATTTTTATAGAATACATATGGGCTTATATAAACATTAATCTCAACATGTTCTAATTTACATATGTAAGCTAATTTTTATTTCTAGAGATAACAGAACAAAACTCAAAACATTTGACATAAAATTATTGGAACAATTAACAGTTTGACCTATTAAACACATTATTGTCCTCTATGAACAGAGGGACTGTCTGAAAAAAAGAACAAGTTGTCTGCATTTTAAAGTGAGAGATAAGCATCAAGGTGTCAATTTCTATTTACACCTTATGTGTTCTTATTTGTTTCACTGATTCATATGTTATAGACACAATATTCTATTCACAATTTTCACGACGTCTATACCAAAGTAAGTATTCAACAAGTAGCCATGAAATGAGGAAATCTGGTAATATACATGAGCTATTAGAATTGTTTTAATGTAAACATTGTCTAGAGAAACAACTAATGTGCATATTTCATAACCGGGAAATGCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGAGACGGAGTCTCGCTGTCGCCCAGGCTGGAGTGCAGTGGCGCAATCTCGGCTCACTGCAGGCTCCGCCCCCTGGGGTTCACGCCATTCTCCTGCCTCAGCCTCCCGAGTAGCTGGGACTACAGGCGCCCGCCACCTCGCCCGGCTAATTTTTTGTATTTTTAGTAGAGACGGGGTTTCACCGTGTTAGCCAGGATGGTCTCGATCTCCTGACCTCGTGATCCGCCCGCCTCGGCCTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACCGCGCCCGGCCCGGGAAATGCTTTTTATTCACATTTTAGTACCATAAGATTGACTAATTAGAAATAGGGAAGCTGTTTGGTATTAGTCCTACTTTTGGGAACATATGGTCTAAAGTAATATTGGACAAATGGATATTTTAATTGATCACAAATGAGAAAGTAGTTAGAAAACTTCTAGTTTAAACAGGTTATATACCCAGAAGTATTGCAAATATTGGAGACAGAAAAATTATTGTAGCTTGCATTTAGACTCAAAATTGATATTCCCTCCACATCATGCCCTCACAGAACTCAGAGGAGTATGATCTGCCAGATCGACTCCCTTTGTGTAACTGATGCAGAAATGAGACTCAAAAGGGTTAAGTAGTTTTCCCAAGATGTCATAGCGAGAGTCAAACTGGGGGCTTAGAAATGACTCTGACTCATAACTTTTAATCAACTGTTTTGACATTTTAACCTATCTAATTGTGTAGGAGGTAATTATATTGTCAGACTTTGGAATGATGTTGTTTCCAGTAAAGTTTTGTTTTAATTATAAATAGGAATTTTCCAGCAATAAAAAATTTCCACCTTAAAAGATTCTCAGACTTTAGTACATCTTTCTCCAAACACAAGGTGGCGATGGTCTACAACAAATGATGTGCGACTTGGTGTTTTTTTTTTTTTTTGTTCTTTCCTTTCCTTTTTTATTCTTAATAGTTCAAGTTAAGAATTTGCAAAAGTTTCACATCTTCTCAATCATGTTTAATAAATTCTAATTAAATATTCTCCTACCTCCTAGTATTATGGAAAATATTTTAAAAATATTACAATGTTAAATGAATTTATTCTTGAGGGCATAATAAAATGCGTTTTTAAATCAACTACTTTTTAATTATGTGTTTGTATTACCATAAACAAAAATCCAATTAAACTTTAAAGAAAGAAAACTGCCTCTGACAAAATAATACTGTGGACCGCTTTTATTCATTACATTTGAGAACTTCTTGTCATTCAAATGAAAAGATTAAGTACATTTGCAATCCACTAAAACAGATTAAAAACTCATTCATTTATTCAATAGATATTAAGTACATACAGTATGTTTAGTATACATTAATACTTGACAATCAATACTGGTTAACTGGTTTCCCTGGTTTAGAAATTTTCCTTAGCAACAACGTAAGGCTTAAAATGAAAAAAGAAAAGTGAGAAAATGTTCTACCACCAGGTGGTGACAAAAGATAAAATTTAAAATCGCTCTTAATGAGCACATACTTCATGTAATTCTTGAATACTGCAAATATAAGTGACTTCCGAATGTCATGTGAATTTAAAATCATATTCTAGGAATATTTTATTAATTAAAGCAAATTAATATTAACATATTATCTCT\t*\tSA:Z:21,26849022,+,1704S657M2I706M,60,2;chr10,97348533,+,1388S317M1364S,0,0;\tMD:Z:1204A189\tRG:Z:GATKSVContigAlignments\tNM:i:1\tAS:i:1389\tXS:i:0", true);
            AlignmentInterval region2 = TestUtilsForAssemblyBasedSVDiscovery.fromSAMRecordString("asm029081:tig00000\t2048\t21\t26849022\t60\t1704H657M2I706M\t*\t0\t0\tCCGGGAAATGCTTTTTATTCACATTTTAGTACCATAAGATTGACTAATTAGAAATAGGGAAGCTGTTTGGTATTAGTCCTACTTTTGGGAACATATGGTCTAAAGTAATATTGGACAAATGGATATTTTAATTGATCACAAATGAGAAAGTAGTTAGAAAACTTCTAGTTTAAACAGGTTATATACCCAGAAGTATTGCAAATATTGGAGACAGAAAAATTATTGTAGCTTGCATTTAGACTCAAAATTGATATTCCCTCCACATCATGCCCTCACAGAACTCAGAGGAGTATGATCTGCCAGATCGACTCCCTTTGTGTAACTGATGCAGAAATGAGACTCAAAAGGGTTAAGTAGTTTTCCCAAGATGTCATAGCGAGAGTCAAACTGGGGGCTTAGAAATGACTCTGACTCATAACTTTTAATCAACTGTTTTGACATTTTAACCTATCTAATTGTGTAGGAGGTAATTATATTGTCAGACTTTGGAATGATGTTGTTTCCAGTAAAGTTTTGTTTTAATTATAAATAGGAATTTTCCAGCAATAAAAAATTTCCACCTTAAAAGATTCTCAGACTTTAGTACATCTTTCTCCAAACACAAGGTGGCGATGGTCTACAACAAATGATGTGCGACTTGGTGTTTTTTTTTTTTTTTGTTCTTTCCTTTCCTTTTTTATTCTTAATAGTTCAAGTTAAGAATTTGCAAAAGTTTCACATCTTCTCAATCATGTTTAATAAATTCTAATTAAATATTCTCCTACCTCCTAGTATTATGGAAAATATTTTAAAAATATTACAATGTTAAATGAATTTATTCTTGAGGGCATAATAAAATGCGTTTTTAAATCAACTACTTTTTAATTATGTGTTTGTATTACCATAAACAAAAATCCAATTAAACTTTAAAGAAAGAAAACTGCCTCTGACAAAATAATACTGTGGACCGCTTTTATTCATTACATTTGAGAACTTCTTGTCATTCAAATGAAAAGATTAAGTACATTTGCAATCCACTAAAACAGATTAAAAACTCATTCATTTATTCAATAGATATTAAGTACATACAGTATGTTTAGTATACATTAATACTTGACAATCAATACTGGTTAACTGGTTTCCCTGGTTTAGAAATTTTCCTTAGCAACAACGTAAGGCTTAAAATGAAAAAAGAAAAGTGAGAAAATGTTCTACCACCAGGTGGTGACAAAAGATAAAATTTAAAATCGCTCTTAATGAGCACATACTTCATGTAATTCTTGAATACTGCAAATATAAGTGACTTCCGAATGTCATGTGAATTTAAAATCATATTCTAGGAATATTTTATTAATTAAAGCAAATTAATATTAACATATTATCTCT\t*\tSA:Z:21,26847644,+,1394M1675S,60,1;chr10,97348533,+,1388S317M1364S,0,0;\tMD:Z:1363\tRG:Z:GATKSVContigAlignments\tNM:i:2\tAS:i:1345\tXS:i:0", true);

            String contigName = "simple_dup_exp_too_small_1_2_with_ins_+";
            byte[] contigSeq = "TATGGTGACAACAGTTACAGCCTCAGAGTGGCTTTGAGATGGAGATTTGAGATGGAGATTTTGGACACTCATAATCTCCATCTCAAGGAGAGGTGGCCCAGAGGCTGAATGAAGATGAGTGAAGGTAGATGTGATTCCCTTAAATTGGGGCAAAAAGGGACAAAAACAGCAGAAAATCTGTATCTTTAAAGACATATGTAATGTATTTCAGTCTATCAACTTCTCTACATAAACTTTAGCTTTTAAAAATATGTTAACGTAAGTTTGACCTTTAGTGTGTTTCTACCTGCAGGGTATTCTTATTGGAGGTTTGTTTAAAAGCATACATTTCTGATCTTGAATGGGTTACTACAAATCCATTATAATTGTTTCATATTTCATGTTGCAGATACAAGTAGGGTTGAAAAAACAGTGAGTTAAAGGCAAAAGGATGGCCGGGAACATGGCTTTTTTATTCTCTGGGTTTCTATCCAGATTTCTGTTCTTTTGCATAATGACTCCAATCTGTTGTGCACCTGTAGTTCTGGGAAATGATTCTTTTTTAATCGCTTCAACAGAGACATGGATGTTGGAGTTGCCAACTACTAAGCTGAAAAACTCCATCTATGCTCAGAAGAACATTTAATCCACTTACTTTTTCTCTTTTATTTAAAGATTAGCACTCATCAGGCATTTGTGGTAATATGCAAATATATACATAGGACATATATGTATATTTATAAGCAAAATGTGAATTGGAAAAACATTTGAATGTAGAAACAAGACCACAGGAGTAAATTTGTACAAGGCACTAGTAAAAGTGACATGTAATATGGGGTTCTTGTAGTGAGTTTCATAATCCAATTTTTGCTCCTTGATTTGAATGGGCACCCAAAATAACACATGCTATCCTAATCCCTACTCCCCATATTTTGGGTTTTATTTTTATAGAATACATATGGGCTTATATAAACATTAATCTCAACATGTTCTAATTTACATATGTAAGCTAATTTTTATTTCTAGAGATAACAGAACAAAACTCAAAACATTTGACATAAAATTATTGGAACAATTAACAGTTTGACCTATTAAACACATTATTGTCCTCTATGAACAGAGGGACTGTCTGAAAAAAAGAACAAGTTGTCTGCATTTTAAAGTGAGAGATAAGCATCAAGGTGTCAATTTCTATTTACACCTTATGTGTTCTTATTTGTTTCACTGATTCATATGTTATAGACACAATATTCTATTCACAATTTTCACGACGTCTATACCAAAGTAAGTATTCAACAAGTAGCCATGAAATGAGGAAATCTGGTAATATACATGAGCTATTAGAATTGTTTTAATGTAAACATTGTCTAGAGAAACAACTAATGTGCATATTTCATAACCGGGAAATGCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGAGACGGAGTCTCGCTGTCGCCCAGGCTGGAGTGCAGTGGCGCAATCTCGGCTCACTGCAGGCTCCGCCCCCTGGGGTTCACGCCATTCTCCTGCCTCAGCCTCCCGAGTAGCTGGGACTACAGGCGCCCGCCACCTCGCCCGGCTAATTTTTTGTATTTTTAGTAGAGACGGGGTTTCACCGTGTTAGCCAGGATGGTCTCGATCTCCTGACCTCGTGATCCGCCCGCCTCGGCCTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACCGCGCCCGGCCCGGGAAATGCTTTTTATTCACATTTTAGTACCATAAGATTGACTAATTAGAAATAGGGAAGCTGTTTGGTATTAGTCCTACTTTTGGGAACATATGGTCTAAAGTAATATTGGACAAATGGATATTTTAATTGATCACAAATGAGAAAGTAGTTAGAAAACTTCTAGTTTAAACAGGTTATATACCCAGAAGTATTGCAAATATTGGAGACAGAAAAATTATTGTAGCTTGCATTTAGACTCAAAATTGATATTCCCTCCACATCATGCCCTCACAGAACTCAGAGGAGTATGATCTGCCAGATCGACTCCCTTTGTGTAACTGATGCAGAAATGAGACTCAAAAGGGTTAAGTAGTTTTCCCAAGATGTCATAGCGAGAGTCAAACTGGGGGCTTAGAAATGACTCTGACTCATAACTTTTAATCAACTGTTTTGACATTTTAACCTATCTAATTGTGTAGGAGGTAATTATATTGTCAGACTTTGGAATGATGTTGTTTCCAGTAAAGTTTTGTTTTAATTATAAATAGGAATTTTCCAGCAATAAAAAATTTCCACCTTAAAAGATTCTCAGACTTTAGTACATCTTTCTCCAAACACAAGGTGGCGATGGTCTACAACAAATGATGTGCGACTTGGTGTTTTTTTTTTTTTTTGTTCTTTCCTTTCCTTTTTTATTCTTAATAGTTCAAGTTAAGAATTTGCAAAAGTTTCACATCTTCTCAATCATGTTTAATAAATTCTAATTAAATATTCTCCTACCTCCTAGTATTATGGAAAATATTTTAAAAATATTACAATGTTAAATGAATTTATTCTTGAGGGCATAATAAAATGCGTTTTTAAATCAACTACTTTTTAATTATGTGTTTGTATTACCATAAACAAAAATCCAATTAAACTTTAAAGAAAGAAAACTGCCTCTGACAAAATAATACTGTGGACCGCTTTTATTCATTACATTTGAGAACTTCTTGTCATTCAAATGAAAAGATTAAGTACATTTGCAATCCACTAAAACAGATTAAAAACTCATTCATTTATTCAATAGATATTAAGTACATACAGTATGTTTAGTATACATTAATACTTGACAATCAATACTGGTTAACTGGTTTCCCTGGTTTAGAAATTTTCCTTAGCAACAACGTAAGGCTTAAAATGAAAAAAGAAAAGTGAGAAAATGTTCTACCACCAGGTGGTGACAAAAGATAAAATTTAAAATCGCTCTTAATGAGCACATACTTCATGTAATTCTTGAATACTGCAAATATAAGTGACTTCCGAATGTCATGTGAATTTAAAATCATATTCTAGGAATATTTTATTAATTAAAGCAAATTAATATTAACATATTATCTCT".getBytes();
            String insertedSeq = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGAGACGGAGTCTCGCTGTCGCCCAGGCTGGAGTGCAGTGGCGCAATCTCGGCTCACTGCAGGCTCCGCCCCCTGGGGTTCACGCCATTCTCCTGCCTCAGCCTCCCGAGTAGCTGGGACTACAGGCGCCCGCCACCTCGCCCGGCTAATTTTTTGTATTTTTAGTAGAGACGGGGTTTCACCGTGTTAGCCAGGATGGTCTCGATCTCCTGACCTCGTGATCCGCCCGCCTCGGCCTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACCGCGCCCGGC";
            final SimpleInterval expectedLeftBreakpoint = new SimpleInterval("21:26849021-26849021");
            final SimpleInterval expectedRightBreakpoint = new SimpleInterval("21:26849021-26849021");
            SimpleChimera expectedSimpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, true, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
            BreakpointComplications expectedBreakpointComplications = new BreakpointComplications.SmallDuplicationWithPreciseDupRangeBreakpointComplications("", insertedSeq,
                    new SimpleInterval("21:26849022-26849037"), 1, 2, Collections.singletonList(Strand.POSITIVE), Arrays.asList(Strand.POSITIVE, Strand.POSITIVE), Arrays.asList("16M", "16M"));
            final String expectedAltSeq = "CCGGGAAATGCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGAGACGGAGTCTCGCTGTCGCCCAGGCTGGAGTGCAGTGGCGCAATCTCGGCTCACTGCAGGCTCCGCCCCCTGGGGTTCACGCCATTCTCCTGCCTCAGCCTCCCGAGTAGCTGGGACTACAGGCGCCCGCCACCTCGCCCGGCTAATTTTTTGTATTTTTAGTAGAGACGGGGTTTCACCGTGTTAGCCAGGATGGTCTCGATCTCCTGACCTCGTGATCCGCCCGCCTCGGCCTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACCGCGCCCGGCCCGGGAAATGCTTTTT";

            final NovelAdjacencyAndAltHaplotype expectedNovelAdjacencyAndAltHaplotype = new NovelAdjacencyAndAltHaplotype(expectedLeftBreakpoint, expectedRightBreakpoint, NO_SWITCH, expectedBreakpointComplications, SMALL_DUP_EXPANSION, expectedAltSeq.getBytes());
            DistancesBetweenAlignmentsOnRefAndOnRead expectedDistances = new DistancesBetweenAlignmentsOnRefAndOnRead(-16, 310, 26849037, 26849022, 1394, 1705);
            final List<SvType> expectedSVTypes = Collections.singletonList(makeInsertionType(new SimpleInterval("21:26849021-26849021"), Allele.create("A", true), insertedSeq.length() + 16));
            final List<VariantContext> expectedVariants = Collections.singletonList(
                    addStandardAttributes(makeInsertion("21", 26849021, 26849021, insertedSeq.length() + 16, Allele.create("A", true)), 1363, contigName, SimpleSVType.SupportedType.INS.name(), 26849021, insertedSeq.length() + 16, expectedAltSeq, "", insertedSeq)
                            .attribute(DUP_REPEAT_UNIT_REF_SPAN, "21:26849022-26849037").attribute(DUP_ORIENTATIONS, "++").attribute(DUP_SEQ_CIGARS, "16M,16M").attribute(DUPLICATION_NUMBERS, "1,2").attribute(DUP_TAN_EXPANSION_STRING, "").make()
            );
            result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeq, false, expectedSimpleChimera, expectedNovelAdjacencyAndAltHaplotype, expectedSVTypes, expectedVariants, expectedDistances, BreakpointsInference.SmallDuplicationWithPreciseDupRangeBreakpointsInference.class));

            contigName = "simple_dup_exp_too_small_1_2_with_ins_-";
            contigSeq = "AGAGATAATATGTTAATATTAATTTGCTTTAATTAATAAAATATTCCTAGAATATGATTTTAAATTCACATGACATTCGGAAGTCACTTATATTTGCAGTATTCAAGAATTACATGAAGTATGTGCTCATTAAGAGCGATTTTAAATTTTATCTTTTGTCACCACCTGGTGGTAGAACATTTTCTCACTTTTCTTTTTTCATTTTAAGCCTTACGTTGTTGCTAAGGAAAATTTCTAAACCAGGGAAACCAGTTAACCAGTATTGATTGTCAAGTATTAATGTATACTAAACATACTGTATGTACTTAATATCTATTGAATAAATGAATGAGTTTTTAATCTGTTTTAGTGGATTGCAAATGTACTTAATCTTTTCATTTGAATGACAAGAAGTTCTCAAATGTAATGAATAAAAGCGGTCCACAGTATTATTTTGTCAGAGGCAGTTTTCTTTCTTTAAAGTTTAATTGGATTTTTGTTTATGGTAATACAAACACATAATTAAAAAGTAGTTGATTTAAAAACGCATTTTATTATGCCCTCAAGAATAAATTCATTTAACATTGTAATATTTTTAAAATATTTTCCATAATACTAGGAGGTAGGAGAATATTTAATTAGAATTTATTAAACATGATTGAGAAGATGTGAAACTTTTGCAAATTCTTAACTTGAACTATTAAGAATAAAAAAGGAAAGGAAAGAACAAAAAAAAAAAAAAAACACCAAGTCGCACATCATTTGTTGTAGACCATCGCCACCTTGTGTTTGGAGAAAGATGTACTAAAGTCTGAGAATCTTTTAAGGTGGAAATTTTTTATTGCTGGAAAATTCCTATTTATAATTAAAACAAAACTTTACTGGAAACAACATCATTCCAAAGTCTGACAATATAATTACCTCCTACACAATTAGATAGGTTAAAATGTCAAAACAGTTGATTAAAAGTTATGAGTCAGAGTCATTTCTAAGCCCCCAGTTTGACTCTCGCTATGACATCTTGGGAAAACTACTTAACCCTTTTGAGTCTCATTTCTGCATCAGTTACACAAAGGGAGTCGATCTGGCAGATCATACTCCTCTGAGTTCTGTGAGGGCATGATGTGGAGGGAATATCAATTTTGAGTCTAAATGCAAGCTACAATAATTTTTCTGTCTCCAATATTTGCAATACTTCTGGGTATATAACCTGTTTAAACTAGAAGTTTTCTAACTACTTTCTCATTTGTGATCAATTAAAATATCCATTTGTCCAATATTACTTTAGACCATATGTTCCCAAAAGTAGGACTAATACCAAACAGCTTCCCTATTTCTAATTAGTCAATCTTATGGTACTAAAATGTGAATAAAAAGCATTTCCCGGGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACGAGGTCAGGAGATCGAGACCATCCTGGCTAACACGGTGAAACCCCGTCTCTACTAAAAATACAAAAAATTAGCCGGGCGAGGTGGCGGGCGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGCAGGAGAATGGCGTGAACCCCAGGGGGCGGAGCCTGCAGTGAGCCGAGATTGCGCCACTGCACTCCAGCCTGGGCGACAGCGAGACTCCGTCTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGCATTTCCCGGTTATGAAATATGCACATTAGTTGTTTCTCTAGACAATGTTTACATTAAAACAATTCTAATAGCTCATGTATATTACCAGATTTCCTCATTTCATGGCTACTTGTTGAATACTTACTTTGGTATAGACGTCGTGAAAATTGTGAATAGAATATTGTGTCTATAACATATGAATCAGTGAAACAAATAAGAACACATAAGGTGTAAATAGAAATTGACACCTTGATGCTTATCTCTCACTTTAAAATGCAGACAACTTGTTCTTTTTTTCAGACAGTCCCTCTGTTCATAGAGGACAATAATGTGTTTAATAGGTCAAACTGTTAATTGTTCCAATAATTTTATGTCAAATGTTTTGAGTTTTGTTCTGTTATCTCTAGAAATAAAAATTAGCTTACATATGTAAATTAGAACATGTTGAGATTAATGTTTATATAAGCCCATATGTATTCTATAAAAATAAAACCCAAAATATGGGGAGTAGGGATTAGGATAGCATGTGTTATTTTGGGTGCCCATTCAAATCAAGGAGCAAAAATTGGATTATGAAACTCACTACAAGAACCCCATATTACATGTCACTTTTACTAGTGCCTTGTACAAATTTACTCCTGTGGTCTTGTTTCTACATTCAAATGTTTTTCCAATTCACATTTTGCTTATAAATATACATATATGTCCTATGTATATATTTGCATATTACCACAAATGCCTGATGAGTGCTAATCTTTAAATAAAAGAGAAAAAGTAAGTGGATTAAATGTTCTTCTGAGCATAGATGGAGTTTTTCAGCTTAGTAGTTGGCAACTCCAACATCCATGTCTCTGTTGAAGCGATTAAAAAAGAATCATTTCCCAGAACTACAGGTGCACAACAGATTGGAGTCATTATGCAAAAGAACAGAAATCTGGATAGAAACCCAGAGAATAAAAAAGCCATGTTCCCGGCCATCCTTTTGCCTTTAACTCACTGTTTTTTCAACCCTACTTGTATCTGCAACATGAAATATGAAACAATTATAATGGATTTGTAGTAACCCATTCAAGATCAGAAATGTATGCTTTTAAACAAACCTCCAATAAGAATACCCTGCAGGTAGAAACACACTAAAGGTCAAACTTACGTTAACATATTTTTAAAAGCTAAAGTTTATGTAGAGAAGTTGATAGACTGAAATACATTACATATGTCTTTAAAGATACAGATTTTCTGCTGTTTTTGTCCCTTTTTGCCCCAATTTAAGGGAATCACATCTACCTTCACTCATCTTCATTCAGCCTCTGGGCCACCTCTCCTTGAGATGGAGATTATGAGTGTCCAAAATCTCCATCTCAAATCTCCATCTCAAAGCCACTCTGAGGCTGTAACTGTTGTCACCATA".getBytes();
            insertedSeq = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGAGACGGAGTCTCGCTGTCGCCCAGGCTGGAGTGCAGTGGCGCAATCTCGGCTCACTGCAGGCTCCGCCCCCTGGGGTTCACGCCATTCTCCTGCCTCAGCCTCCCGAGTAGCTGGGACTACAGGCGCCCGCCACCTCGCCCGGCTAATTTTTTGTATTTTTAGTAGAGACGGGGTTTCACCGTGTTAGCCAGGATGGTCTCGATCTCCTGACCTCGTGATCCGCCCGCCTCGGCCTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACCGCGCCCGGC";
            expectedBreakpointComplications = new BreakpointComplications.SmallDuplicationWithPreciseDupRangeBreakpointComplications("", insertedSeq,
                    new SimpleInterval("21:26849022-26849037"), 1, 2, Collections.singletonList(Strand.POSITIVE), Arrays.asList(Strand.POSITIVE, Strand.POSITIVE), Arrays.asList("16M", "16M"));
            final NovelAdjacencyAndAltHaplotype expectedNovelAdjacencyAndAltHaplotypeDetectedFromReverseStrand = new NovelAdjacencyAndAltHaplotype(expectedLeftBreakpoint, expectedRightBreakpoint, NO_SWITCH, expectedBreakpointComplications, SMALL_DUP_EXPANSION, expectedAltSeq.getBytes());
            expectedDistances = new DistancesBetweenAlignmentsOnRefAndOnRead(-16, 310, 26849037, 26849022, 1366, 1677);
            region1 = TestUtilsForAssemblyBasedSVDiscovery.fromSAMRecordString("asm000001:tig00001\t2064\t21\t26849022\t60\t1704H657M3I706M\t*\t0\t0\tCCGGGAAATGCTTTTTATTCACATTTTAGTACCATAAGATTGACTAATTAGAAATAGGGAAGCTGTTTGGTATTAGTCCTACTTTTGGGAACATATGGTCTAAAGTAATATTGGACAAATGGATATTTTAATTGATCACAAATGAGAAAGTAGTTAGAAAACTTCTAGTTTAAACAGGTTATATACCCAGAAGTATTGCAAATATTGGAGACAGAAAAATTATTGTAGCTTGCATTTAGACTCAAAATTGATATTCCCTCCACATCATGCCCTCACAGAACTCAGAGGAGTATGATCTGCCAGATCGACTCCCTTTGTGTAACTGATGCAGAAATGAGACTCAAAAGGGTTAAGTAGTTTTCCCAAGATGTCATAGCGAGAGTCAAACTGGGGGCTTAGAAATGACTCTGACTCATAACTTTTAATCAACTGTTTTGACATTTTAACCTATCTAATTGTGTAGGAGGTAATTATATTGTCAGACTTTGGAATGATGTTGTTTCCAGTAAAGTTTTGTTTTAATTATAAATAGGAATTTTCCAGCAATAAAAAATTTCCACCTTAAAAGATTCTCAGACTTTAGTACATCTTTCTCCAAACACAAGGTGGCGATGGTCTACAACAAATGATGTGCGACTTGGTGTTTTTTTTTTTTTTTTGTTCTTTCCTTTCCTTTTTTATTCTTAATAGTTCAAGTTAAGAATTTGCAAAAGTTTCACATCTTCTCAATCATGTTTAATAAATTCTAATTAAATATTCTCCTACCTCCTAGTATTATGGAAAATATTTTAAAAATATTACAATGTTAAATGAATTTATTCTTGAGGGCATAATAAAATGCGTTTTTAAATCAACTACTTTTTAATTATGTGTTTGTATTACCATAAACAAAAATCCAATTAAACTTTAAAGAAAGAAAACTGCCTCTGACAAAATAATACTGTGGACCGCTTTTATTCATTACATTTGAGAACTTCTTGTCATTCAAATGAAAAGATTAAGTACATTTGCAATCCACTAAAACAGATTAAAAACTCATTCATTTATTCAATAGATATTAAGTACATACAGTATGTTTAGTATACATTAATACTTGACAATCAATACTGGTTAACTGGTTTCCCTGGTTTAGAAATTTTCCTTAGCAACAACGTAAGGCTTAAAATGAAAAAAGAAAAGTGAGAAAATGTTCTACCACCAGGTGGTGACAAAAGATAAAATTTAAAATCGCTCTTAATGAGCACATACTTCATGTAATTCTTGAATACTGCAAATATAAGTGACTTCCGAATGTCATGTGAATTTAAAATCATATTCTAGGAATATTTTATTAATTAAAGCAAATTAATATTAACATATTATCTCT\t*\tSA:Z:21,26847644,-,1394M1676S,60,1;chr10,97348533,-,1388S317M1365S,0,0;\tMD:Z:1363\tRG:Z:GATKSVContigAlignments\tNM:i:3\tAS:i:1344\tXS:i:0", true);
            region2 = TestUtilsForAssemblyBasedSVDiscovery.fromSAMRecordString("asm000001:tig00001\t16\t21\t26847644\t60\t1394M1676S\t*\t0\t0\tTATGGTGACAACAGTTACAGCCTCAGAGTGGCTTTGAGATGGAGATTTGAGATGGAGATTTTGGACACTCATAATCTCCATCTCAAGGAGAGGTGGCCCAGAGGCTGAATGAAGATGAGTGAAGGTAGATGTGATTCCCTTAAATTGGGGCAAAAAGGGACAAAAACAGCAGAAAATCTGTATCTTTAAAGACATATGTAATGTATTTCAGTCTATCAACTTCTCTACATAAACTTTAGCTTTTAAAAATATGTTAACGTAAGTTTGACCTTTAGTGTGTTTCTACCTGCAGGGTATTCTTATTGGAGGTTTGTTTAAAAGCATACATTTCTGATCTTGAATGGGTTACTACAAATCCATTATAATTGTTTCATATTTCATGTTGCAGATACAAGTAGGGTTGAAAAAACAGTGAGTTAAAGGCAAAAGGATGGCCGGGAACATGGCTTTTTTATTCTCTGGGTTTCTATCCAGATTTCTGTTCTTTTGCATAATGACTCCAATCTGTTGTGCACCTGTAGTTCTGGGAAATGATTCTTTTTTAATCGCTTCAACAGAGACATGGATGTTGGAGTTGCCAACTACTAAGCTGAAAAACTCCATCTATGCTCAGAAGAACATTTAATCCACTTACTTTTTCTCTTTTATTTAAAGATTAGCACTCATCAGGCATTTGTGGTAATATGCAAATATATACATAGGACATATATGTATATTTATAAGCAAAATGTGAATTGGAAAAACATTTGAATGTAGAAACAAGACCACAGGAGTAAATTTGTACAAGGCACTAGTAAAAGTGACATGTAATATGGGGTTCTTGTAGTGAGTTTCATAATCCAATTTTTGCTCCTTGATTTGAATGGGCACCCAAAATAACACATGCTATCCTAATCCCTACTCCCCATATTTTGGGTTTTATTTTTATAGAATACATATGGGCTTATATAAACATTAATCTCAACATGTTCTAATTTACATATGTAAGCTAATTTTTATTTCTAGAGATAACAGAACAAAACTCAAAACATTTGACATAAAATTATTGGAACAATTAACAGTTTGACCTATTAAACACATTATTGTCCTCTATGAACAGAGGGACTGTCTGAAAAAAAGAACAAGTTGTCTGCATTTTAAAGTGAGAGATAAGCATCAAGGTGTCAATTTCTATTTACACCTTATGTGTTCTTATTTGTTTCACTGATTCATATGTTATAGACACAATATTCTATTCACAATTTTCACGACGTCTATACCAAAGTAAGTATTCAACAAGTAGCCATGAAATGAGGAAATCTGGTAATATACATGAGCTATTAGAATTGTTTTAATGTAAACATTGTCTAGAGAAACAACTAATGTGCATATTTCATAACCGGGAAATGCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGAGACGGAGTCTCGCTGTCGCCCAGGCTGGAGTGCAGTGGCGCAATCTCGGCTCACTGCAGGCTCCGCCCCCTGGGGTTCACGCCATTCTCCTGCCTCAGCCTCCCGAGTAGCTGGGACTACAGGCGCCCGCCACCTCGCCCGGCTAATTTTTTGTATTTTTAGTAGAGACGGGGTTTCACCGTGTTAGCCAGGATGGTCTCGATCTCCTGACCTCGTGATCCGCCCGCCTCGGCCTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACCGCGCCCGGCCCGGGAAATGCTTTTTATTCACATTTTAGTACCATAAGATTGACTAATTAGAAATAGGGAAGCTGTTTGGTATTAGTCCTACTTTTGGGAACATATGGTCTAAAGTAATATTGGACAAATGGATATTTTAATTGATCACAAATGAGAAAGTAGTTAGAAAACTTCTAGTTTAAACAGGTTATATACCCAGAAGTATTGCAAATATTGGAGACAGAAAAATTATTGTAGCTTGCATTTAGACTCAAAATTGATATTCCCTCCACATCATGCCCTCACAGAACTCAGAGGAGTATGATCTGCCAGATCGACTCCCTTTGTGTAACTGATGCAGAAATGAGACTCAAAAGGGTTAAGTAGTTTTCCCAAGATGTCATAGCGAGAGTCAAACTGGGGGCTTAGAAATGACTCTGACTCATAACTTTTAATCAACTGTTTTGACATTTTAACCTATCTAATTGTGTAGGAGGTAATTATATTGTCAGACTTTGGAATGATGTTGTTTCCAGTAAAGTTTTGTTTTAATTATAAATAGGAATTTTCCAGCAATAAAAAATTTCCACCTTAAAAGATTCTCAGACTTTAGTACATCTTTCTCCAAACACAAGGTGGCGATGGTCTACAACAAATGATGTGCGACTTGGTGTTTTTTTTTTTTTTTTGTTCTTTCCTTTCCTTTTTTATTCTTAATAGTTCAAGTTAAGAATTTGCAAAAGTTTCACATCTTCTCAATCATGTTTAATAAATTCTAATTAAATATTCTCCTACCTCCTAGTATTATGGAAAATATTTTAAAAATATTACAATGTTAAATGAATTTATTCTTGAGGGCATAATAAAATGCGTTTTTAAATCAACTACTTTTTAATTATGTGTTTGTATTACCATAAACAAAAATCCAATTAAACTTTAAAGAAAGAAAACTGCCTCTGACAAAATAATACTGTGGACCGCTTTTATTCATTACATTTGAGAACTTCTTGTCATTCAAATGAAAAGATTAAGTACATTTGCAATCCACTAAAACAGATTAAAAACTCATTCATTTATTCAATAGATATTAAGTACATACAGTATGTTTAGTATACATTAATACTTGACAATCAATACTGGTTAACTGGTTTCCCTGGTTTAGAAATTTTCCTTAGCAACAACGTAAGGCTTAAAATGAAAAAAGAAAAGTGAGAAAATGTTCTACCACCAGGTGGTGACAAAAGATAAAATTTAAAATCGCTCTTAATGAGCACATACTTCATGTAATTCTTGAATACTGCAAATATAAGTGACTTCCGAATGTCATGTGAATTTAAAATCATATTCTAGGAATATTTTATTAATTAAAGCAAATTAATATTAACATATTATCTCT\t*\tSA:Z:21,26849022,-,1704S657M3I706M,60,3;chr10,97348533,-,1388S317M1365S,0,0;\tMD:Z:1204A189\tRG:Z:GATKSVContigAlignments\tNM:i:1\tAS:i:1384\tXS:i:0", true);
            expectedSimpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, false, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
            result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeq, true, expectedSimpleChimera, expectedNovelAdjacencyAndAltHaplotypeDetectedFromReverseStrand, expectedSVTypes, expectedVariants, expectedDistances, BreakpointsInference.SmallDuplicationWithPreciseDupRangeBreakpointsInference.class));
        }

        {
            // simple tandem duplication expansion with novel insertion '+' strand representation
            final String leftRefFlank = "GTTAGTAGATATTCTAGCTGACTCAGTTCAGTGTTGCTATGATTAAACAAGAGTGAGTTCCCT";                     //63
            final String rightRefFlank = "CATTATTGATATTTCATTATGTTCAACAGATGGAGTTAATGTGAATGT";                                   //48
            final String insertedSeq = "CTCTCTCTCT";                                                                           //10
            final String dup = "AAAAGTAAATGTTATAAGAAATCTTAAGTATTATTTTCTTATGTTTCTAGCCTAATAAAGTGCTTTTATTAAAGCACTTTATTTAAAGG";    //89
            byte[] contigSeq = (leftRefFlank + dup + insertedSeq + dup + rightRefFlank).getBytes();
            String contigName = "simple_dup_exp_1_2_with_ins_+";

            final SimpleInterval expectedLeftBreakpoint = new SimpleInterval("21", 25297163, 25297163);
            final SimpleInterval expectedRightBreakpoint = new SimpleInterval("21", 25297163, 25297163);
            final BreakpointComplications expectedBreakpointComplications = new BreakpointComplications.SmallDuplicationWithPreciseDupRangeBreakpointComplications("", new String(insertedSeq),
                    new SimpleInterval("21", 25297164,25297252), 1, 2, Collections.singletonList(Strand.POSITIVE), Arrays.asList(Strand.POSITIVE, Strand.POSITIVE), Arrays.asList("89M", "89M"));
            final byte[] expectedAltSeq = Arrays.copyOfRange(contigSeq, leftRefFlank.length(), contigSeq.length - rightRefFlank.length());
            final NovelAdjacencyAndAltHaplotype expectedNovelAdjacencyAndAltHaplotype = new NovelAdjacencyAndAltHaplotype(expectedLeftBreakpoint, expectedRightBreakpoint, NO_SWITCH, expectedBreakpointComplications, SMALL_DUP_EXPANSION, expectedAltSeq);
            AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("21", 25297101, 25297252), 1 ,152, TextCigarCodec.decode("152M147S"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
            AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("21", 25297164, 25297300), 163 ,299, TextCigarCodec.decode("162S137M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
            SimpleChimera expectedSimpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, true, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
            DistancesBetweenAlignmentsOnRefAndOnRead expectedDistances = new DistancesBetweenAlignmentsOnRefAndOnRead(-89, 10, 25297252, 25297164, 152, 163);
            final List<SvType> expectedSVTypes = Collections.singletonList(makeTandemDuplicationType(new SimpleInterval("21", 25297164,25297252), Allele.create("T", true), 99));
            final List<VariantContext> expectedVariants = Collections.singletonList(
                    addStandardAttributes(makeTandemDuplication(new SimpleInterval("21", 25297164,25297252), Allele.create("T", true), 99), 137, contigName, SimpleSVType.SupportedType.DUP.name(), 25297163, 99, StringUtil.bytesToString(expectedAltSeq), "", insertedSeq)
                            .attribute(DUP_REPEAT_UNIT_REF_SPAN, "21:25297164-25297252").attribute(DUP_ORIENTATIONS, "++").attribute(DUP_SEQ_CIGARS, "89M,89M").attribute(DUPLICATION_NUMBERS, "1,2").attribute(DUP_TAN_EXPANSION_STRING, "").make()
            );
            result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeq, false, expectedSimpleChimera, expectedNovelAdjacencyAndAltHaplotype, expectedSVTypes, expectedVariants, expectedDistances, BreakpointsInference.SmallDuplicationWithPreciseDupRangeBreakpointsInference.class));

            // simple tandem duplication expansion with novel insertion '-' strand representation
            contigSeq = (SequenceUtil.reverseComplement(rightRefFlank) + SequenceUtil.reverseComplement(dup) + SequenceUtil.reverseComplement(insertedSeq) + SequenceUtil.reverseComplement(dup) + SequenceUtil.reverseComplement(leftRefFlank)).getBytes();
            contigName = "simple_dup_exp_1_2_with_ins_-";

            region1 = new AlignmentInterval(new SimpleInterval("21", 25297164, 25297300), 1 ,137, TextCigarCodec.decode("137M162S"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
            region2 = new AlignmentInterval(new SimpleInterval("21", 25297101, 25297252), 148 ,299, TextCigarCodec.decode("147S152M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
            expectedSimpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, false, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
            expectedDistances = new DistancesBetweenAlignmentsOnRefAndOnRead(-89, 10, 25297252, 25297164, 137, 148);
            result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeq, true, expectedSimpleChimera, expectedNovelAdjacencyAndAltHaplotype, expectedSVTypes, expectedVariants, expectedDistances, BreakpointsInference.SmallDuplicationWithPreciseDupRangeBreakpointsInference.class));
        }

        return result;
    }

    /**
     * These test data was based on a real observation on a locally-assembled contig
     * "TGCCAGGTTACATGGCAAAGAGGGTAGATATGGGGAGCTGTGAAGAATGGAGCCAGTAATTAAATTCACTGAAGTCTCCACAGGAGGGCAAGGTGGACAATCTGTCCCATAGGAGGGGGATTCATGAGGGGAGCTGTGAAGAATGGAGCCAGTAATTAAATTCACTGAAGTCTCCACAGGAGGGCAAGGTGGACAATCTGTCCCATAGGAGGGGGATTCAGGAGGGCAGCTGTGGATGGTGCAAATGCCATTTATGCTCCTCTCCACCCATATCC"
     * with two alignment records chr18:312579-312718 140M135S
     *                            chr18:312610-312757 127S148M
     * for a tandem repeat expansion event from 1 copy to 2 copies with also a pseudo-homology

     * Return a list of eight entries for positive and reverse strand representations for:
     * 1. expansion from 1 unit to 2 units with pseudo-homology
     * 2. contraction from 2 units to 1 unit with pseudo-homology
     * 3. contraction from 3 units to 2 units without pseudo-homology
     * 4. expansion from 2 units to 3 units without pseudo-homology
     */
    private static List<TestDataForSimpleSV>
    forComplexTandemDuplication() {

        final List<TestDataForSimpleSV> result = new ArrayList<>();
        final String leftRefFlank       = "TGCCAGGTTACATGGCAAAGAGGGTAGATAT";                                                                    // 31
        final String rightRefFlank      = "TGGTGCAAATGCCATTTATGCTCCTCTCCACCCATATCC";                                                            // 39
        final String firstRepeat        = "GGGGAGCTGTGAAGAATGGAGCCAGTAATTAAATTCACTGAAGTCTCCACAGGAGGGCAAGGTGGACAATCTGTCCCATAGGAGGGGGATTCATGA";   // 96
        final String secondRepeat       = "GGGGAGCTGTGAAGAATGGAGCCAGTAATTAAATTCACTGAAGTCTCCACAGGAGGGCAAGGTGGACAATCTGTCCCATAGGAGGGGGATTCAGGA";   // 96
        final String pseudoHomology     = "GGGCAGCTGTGGA";                                                                                      // 13


        // first test (the original observed event, but assigned to a different chromosome): expansion from 1 unit to 2 units with pseudo-homology
        String contigName = "cpx_dup_exp_1_2_pseudo_+";
        final byte[] fakeRefSeqForComplexExpansionWithPseudoHomology = String.format("%s%s%s%s", leftRefFlank, firstRepeat, pseudoHomology, rightRefFlank).getBytes();
        final byte[] contigSeqForComplexExpansionWithPseudoHomology = String.format("%s%s%s%s%s", leftRefFlank, firstRepeat, secondRepeat, pseudoHomology, rightRefFlank).getBytes();
        SimpleInterval expectedLeftBreakpoint = new SimpleInterval("20", 312609, 312609);
        SimpleInterval expectedRightBreakpoint = new SimpleInterval("20", 312609, 312609);
        BreakpointComplications expectedBreakpointComplications = new BreakpointComplications.SmallDuplicationWithImpreciseDupRangeBreakpointComplications(pseudoHomology, "",
                new SimpleInterval("20", 312610, 312705), 1, 2, Collections.singletonList(Strand.POSITIVE), Arrays.asList(Strand.POSITIVE, Strand.POSITIVE),
                new SimpleInterval("20", 312610, 312718));
        byte[] expectedAltSeq = (firstRepeat+secondRepeat+pseudoHomology).getBytes();
        NovelAdjacencyAndAltHaplotype expectedNovelAdjacencyAndAltHaplotype = new NovelAdjacencyAndAltHaplotype(expectedLeftBreakpoint, expectedRightBreakpoint, NO_SWITCH, expectedBreakpointComplications, SMALL_DUP_CPX, expectedAltSeq);
        AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("20", 312579, 312718), 1 ,140, TextCigarCodec.decode("140M135S"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("20", 312610, 312757), 128 ,275, TextCigarCodec.decode("127S148M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        SimpleChimera expectedSimpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, true, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
        DistancesBetweenAlignmentsOnRefAndOnRead expectedDistances = new DistancesBetweenAlignmentsOnRefAndOnRead(-109, -13, 312718, 312610, 140, 128);
        List<SvType> expectedSVTypes = Collections.singletonList(makeTandemDuplicationType(new SimpleInterval("20:312610-312705"), Allele.create("T", true), 96));
        List<VariantContext> expectedVariants = Collections.singletonList(
                addStandardAttributes(makeTandemDuplication(new SimpleInterval("20:312610-312705"), Allele.create("T", true), 96), 140 - pseudoHomology.length(), contigName, SimpleSVType.SupportedType.DUP.name(), 312609, 96, StringUtil.bytesToString(expectedAltSeq), pseudoHomology, "")
                        .attribute(DUP_REPEAT_UNIT_REF_SPAN, "20:312610-312705").attribute(DUP_IMPRECISE_AFFECTED_RANGE, "20:312610-312718").attribute(DUP_ORIENTATIONS, "++").attribute(DUPLICATION_NUMBERS, "1,2").attribute(DUP_TAN_EXPANSION_STRING, "").attribute(DUP_ANNOTATIONS_IMPRECISE, "").make()
        );
        result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeqForComplexExpansionWithPseudoHomology, false, expectedSimpleChimera, expectedNovelAdjacencyAndAltHaplotype, expectedSVTypes, expectedVariants, expectedDistances, BreakpointsInference.SmallDuplicationWithImpreciseDupRangeBreakpointsInference.class));

        contigName = "cpx_dup_exp_1_2_pseudo_-";
        final byte[] contigSeqForComplexExpansionWithPseudoHomology_reverseStrand = Arrays.copyOf(contigSeqForComplexExpansionWithPseudoHomology, contigSeqForComplexExpansionWithPseudoHomology.length);
        SequenceUtil.reverseComplement(contigSeqForComplexExpansionWithPseudoHomology_reverseStrand);
        expectedNovelAdjacencyAndAltHaplotype = new NovelAdjacencyAndAltHaplotype(expectedLeftBreakpoint, expectedRightBreakpoint, NO_SWITCH, expectedBreakpointComplications, SMALL_DUP_CPX, expectedAltSeq);
        region1 = new AlignmentInterval(new SimpleInterval("20", 312610, 312757), 1 ,148, TextCigarCodec.decode("148M127S"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        region2 = new AlignmentInterval(new SimpleInterval("20", 312579, 312718), 136 ,275, TextCigarCodec.decode("135S140M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        expectedSimpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, false, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
        expectedDistances = new DistancesBetweenAlignmentsOnRefAndOnRead(-109, -13, 312718, 312610, 148, 136);
        result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeqForComplexExpansionWithPseudoHomology_reverseStrand, true, expectedSimpleChimera, expectedNovelAdjacencyAndAltHaplotype, expectedSVTypes, expectedVariants, expectedDistances, BreakpointsInference.SmallDuplicationWithImpreciseDupRangeBreakpointsInference.class));

        // second test: contraction from 2 units to 1 unit with pseudo-homology
        contigName = "cpx_dup_contract_2_1_pseudo_+";
        final byte[] contigSeqForComplexContractionWithPseudoHomology = fakeRefSeqForComplexExpansionWithPseudoHomology;
        expectedLeftBreakpoint = new SimpleInterval("20", 312609, 312609);
        expectedRightBreakpoint = new SimpleInterval("20", 312705, 312705);
        expectedBreakpointComplications = new BreakpointComplications.SmallDuplicationWithImpreciseDupRangeBreakpointComplications(firstRepeat+pseudoHomology, "",
                new SimpleInterval("20", 312610, 312705), 2, 1, Arrays.asList(Strand.POSITIVE, Strand.POSITIVE), Collections.singletonList(Strand.POSITIVE),
                new SimpleInterval("20", 312610, 312814));
        expectedAltSeq= (firstRepeat+pseudoHomology).getBytes();
        expectedNovelAdjacencyAndAltHaplotype = new NovelAdjacencyAndAltHaplotype(expectedLeftBreakpoint, expectedRightBreakpoint, NO_SWITCH, expectedBreakpointComplications, SMALL_DUP_CPX, expectedAltSeq);
        region1 = new AlignmentInterval(new SimpleInterval("20", 312579, 312718), 1, 140, TextCigarCodec.decode("140M39S"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        region2 = new AlignmentInterval(new SimpleInterval("20", 312706, 312853), 32, 179, TextCigarCodec.decode("31S148M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        expectedSimpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, true, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
        expectedDistances = new DistancesBetweenAlignmentsOnRefAndOnRead(-13, -109, 312718, 312706, 140, 32);
        expectedSVTypes = Collections.singletonList(makeDeletionType(new SimpleInterval("20:312609-312705"), Allele.create("T", true), true));
        expectedVariants = Collections.singletonList(
                addStandardAttributes(makeDeletion(new SimpleInterval("20:312609-312705"), Allele.create("T", true), true), 140 - expectedBreakpointComplications.getHomologyForwardStrandRep().length(), contigName, SimpleSVType.SupportedType.DEL.name(), 312705, -96, StringUtil.bytesToString(expectedAltSeq), firstRepeat+pseudoHomology, "")
                        .attribute(DUP_REPEAT_UNIT_REF_SPAN, "20:312610-312705").attribute(DUP_IMPRECISE_AFFECTED_RANGE, "20:312610-312814").attribute(DUPLICATION_NUMBERS, "2,1").attribute(DUP_ORIENTATIONS, "+").attribute(DUP_TAN_CONTRACTION_STRING, "").attribute(DUP_ANNOTATIONS_IMPRECISE, "").make()
        );
        result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeqForComplexContractionWithPseudoHomology, false, expectedSimpleChimera, expectedNovelAdjacencyAndAltHaplotype, expectedSVTypes, expectedVariants, expectedDistances, BreakpointsInference.SmallDuplicationWithImpreciseDupRangeBreakpointsInference.class));

        contigName = "cpx_dup_contract_2_1_pseudo_-";
        final byte[] contigSeqForComplexContractionWithPseudoHomology_reverseStrand = Arrays.copyOf(contigSeqForComplexContractionWithPseudoHomology, contigSeqForComplexContractionWithPseudoHomology.length);
        SequenceUtil.reverseComplement(contigSeqForComplexContractionWithPseudoHomology_reverseStrand);
        expectedNovelAdjacencyAndAltHaplotype =  new NovelAdjacencyAndAltHaplotype(expectedLeftBreakpoint, expectedRightBreakpoint, NO_SWITCH, expectedBreakpointComplications, SMALL_DUP_CPX, expectedAltSeq);
        region1 = new AlignmentInterval(new SimpleInterval("20", 312706, 312853), 1, 148, TextCigarCodec.decode("148M31S"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        region2 = new AlignmentInterval(new SimpleInterval("20", 312579, 312718), 40, 179, TextCigarCodec.decode("39S140M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        expectedSimpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, false, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
        expectedDistances = new DistancesBetweenAlignmentsOnRefAndOnRead(-13, -109, 312718, 312706, 148, 40);
        result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeqForComplexContractionWithPseudoHomology_reverseStrand, true, expectedSimpleChimera, expectedNovelAdjacencyAndAltHaplotype, expectedSVTypes, expectedVariants, expectedDistances, BreakpointsInference.SmallDuplicationWithImpreciseDupRangeBreakpointsInference.class));

        // third test: contraction from 3 units to 2 units without pseudo-homology
        contigName = "cpx_dup_contract_3_2_+";
        final byte[] fakeRefSeqForComplexContractionNoPseudoHomology = String.format("%s%s%s%s%s", leftRefFlank, firstRepeat, secondRepeat, firstRepeat, rightRefFlank).getBytes();
        final byte[] contigSeqForComplexContractionNoPseudoHomology = String.format("%s%s%s%s", leftRefFlank, firstRepeat, secondRepeat, rightRefFlank).getBytes();
        expectedLeftBreakpoint = new SimpleInterval("20", 312609, 312609);
        expectedRightBreakpoint = new SimpleInterval("20", 312705, 312705);
        expectedBreakpointComplications = new BreakpointComplications.SmallDuplicationWithImpreciseDupRangeBreakpointComplications(firstRepeat+secondRepeat, "",
                new SimpleInterval("20", 312610, 312705), 3, 2, Arrays.asList(Strand.POSITIVE, Strand.POSITIVE, Strand.POSITIVE), Arrays.asList(Strand.POSITIVE, Strand.POSITIVE),
                new SimpleInterval("20", 312610, 312897));
        expectedAltSeq= (firstRepeat+secondRepeat).getBytes();
        expectedNovelAdjacencyAndAltHaplotype = new NovelAdjacencyAndAltHaplotype(expectedLeftBreakpoint, expectedRightBreakpoint, NO_SWITCH, expectedBreakpointComplications, SMALL_DUP_CPX, expectedAltSeq);
        region1 = new AlignmentInterval(new SimpleInterval("20", 312579, 312801), 1, 223, TextCigarCodec.decode("223M39S"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        region2 = new AlignmentInterval(new SimpleInterval("20", 312706, 312936), 32, 262, TextCigarCodec.decode("31S231M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        expectedSimpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, true, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
        expectedDistances = new DistancesBetweenAlignmentsOnRefAndOnRead(-96, -192, 312801, 312706, 223, 32);
        expectedSVTypes = Collections.singletonList(makeDeletionType(new SimpleInterval("20:312609-312705"), Allele.create("T", true), true));
        expectedVariants = Collections.singletonList(
                addStandardAttributes(makeDeletion(new SimpleInterval("20:312609-312705"), Allele.create("T", true), true), 223 - expectedBreakpointComplications.getHomologyForwardStrandRep().length(), contigName, SimpleSVType.SupportedType.DEL.name(), 312705, -96, StringUtil.bytesToString(expectedAltSeq), firstRepeat+secondRepeat, "")
                        .attribute(DUP_REPEAT_UNIT_REF_SPAN, "20:312610-312705").attribute(DUP_IMPRECISE_AFFECTED_RANGE, "20:312610-312897").attribute(DUPLICATION_NUMBERS, "3,2").attribute(DUP_ORIENTATIONS, "++").attribute(DUP_TAN_CONTRACTION_STRING, "").attribute(DUP_ANNOTATIONS_IMPRECISE, "").make()
        );
        result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeqForComplexContractionNoPseudoHomology, false, expectedSimpleChimera, expectedNovelAdjacencyAndAltHaplotype, expectedSVTypes, expectedVariants, expectedDistances, BreakpointsInference.SmallDuplicationWithImpreciseDupRangeBreakpointsInference.class));

        contigName = "cpx_dup_contract_3_2_-";
        final byte[] contigSeqForComplexContractionNoPseudoHomology_reverseStrand = Arrays.copyOf(contigSeqForComplexContractionNoPseudoHomology, contigSeqForComplexContractionNoPseudoHomology.length);
        SequenceUtil.reverseComplement(contigSeqForComplexContractionNoPseudoHomology_reverseStrand);
        expectedNovelAdjacencyAndAltHaplotype = new NovelAdjacencyAndAltHaplotype(expectedLeftBreakpoint, expectedRightBreakpoint, NO_SWITCH, expectedBreakpointComplications, SMALL_DUP_CPX, expectedAltSeq);
        region1 = new AlignmentInterval(new SimpleInterval("20", 312706, 312936), 1, 231, TextCigarCodec.decode("231M31S"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        region2 = new AlignmentInterval(new SimpleInterval("20", 312579, 312801), 40, 262, TextCigarCodec.decode("39S223M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        expectedSimpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, false, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
        expectedDistances = new DistancesBetweenAlignmentsOnRefAndOnRead(-96, -192, 312801, 312706, 231, 40);
        result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeqForComplexContractionNoPseudoHomology_reverseStrand, true, expectedSimpleChimera, expectedNovelAdjacencyAndAltHaplotype, expectedSVTypes, expectedVariants, expectedDistances, BreakpointsInference.SmallDuplicationWithImpreciseDupRangeBreakpointsInference.class));

        // fourth test: expansion from 2 units to 3 units without pseudo-homology
        contigName = "cpx_dup_exp_2_3_+";
        final byte[] contigSeqForComplexExpansionNoPseudoHomology = fakeRefSeqForComplexContractionNoPseudoHomology;
        expectedLeftBreakpoint = new SimpleInterval("20", 312609, 312609);
        expectedRightBreakpoint = new SimpleInterval("20", 312609, 312609);
        expectedBreakpointComplications = new BreakpointComplications.SmallDuplicationWithImpreciseDupRangeBreakpointComplications(firstRepeat, "",
                new SimpleInterval("20", 312610, 312705), 2, 3, Arrays.asList(Strand.POSITIVE, Strand.POSITIVE), Arrays.asList(Strand.POSITIVE, Strand.POSITIVE, Strand.POSITIVE),
                new SimpleInterval("20", 312610, 312801));
         expectedAltSeq= (firstRepeat+secondRepeat+firstRepeat).getBytes();
        expectedNovelAdjacencyAndAltHaplotype = new NovelAdjacencyAndAltHaplotype(expectedLeftBreakpoint, expectedRightBreakpoint, NO_SWITCH, expectedBreakpointComplications, SMALL_DUP_CPX, expectedAltSeq);
        region1 = new AlignmentInterval(new SimpleInterval("20", 312579, 312801), 1, 223, TextCigarCodec.decode("223M135S"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        region2 = new AlignmentInterval(new SimpleInterval("20", 312610, 312840), 128, 358, TextCigarCodec.decode("127S231M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        expectedSimpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, true, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
        expectedDistances = new DistancesBetweenAlignmentsOnRefAndOnRead(-192, -96, 312801, 312610, 223, 128);
        expectedSVTypes = Collections.singletonList(makeTandemDuplicationType(new SimpleInterval("20:312610-312705"), Allele.create("T", true), 96));
        expectedVariants = Collections.singletonList(
                addStandardAttributes(makeTandemDuplication(new SimpleInterval("20:312610-312705"), Allele.create("T", true), 96), 223 - 96, contigName, SimpleSVType.SupportedType.DUP.name(), 312609, 96, StringUtil.bytesToString(expectedAltSeq), firstRepeat, "")
                        .attribute(DUP_REPEAT_UNIT_REF_SPAN, "20:312610-312705").attribute(DUP_IMPRECISE_AFFECTED_RANGE, "20:312610-312801").attribute(DUP_ORIENTATIONS, "+++").attribute(DUPLICATION_NUMBERS, "2,3").attribute(DUP_TAN_EXPANSION_STRING, "").attribute(DUP_ANNOTATIONS_IMPRECISE, "").make()
        );
        result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeqForComplexExpansionNoPseudoHomology, false, expectedSimpleChimera, expectedNovelAdjacencyAndAltHaplotype, expectedSVTypes, expectedVariants, expectedDistances, BreakpointsInference.SmallDuplicationWithImpreciseDupRangeBreakpointsInference.class));

        contigName = "cpx_dup_exp_2_3_-";
        final byte[] contigSeqForComplexExpansionNoPseudoHomology_reverseStrand = Arrays.copyOf(contigSeqForComplexExpansionNoPseudoHomology, contigSeqForComplexExpansionNoPseudoHomology.length);
        SequenceUtil.reverseComplement(contigSeqForComplexExpansionNoPseudoHomology_reverseStrand);
        expectedNovelAdjacencyAndAltHaplotype = new NovelAdjacencyAndAltHaplotype(expectedLeftBreakpoint, expectedRightBreakpoint, NO_SWITCH, expectedBreakpointComplications, SMALL_DUP_CPX, expectedAltSeq);
        region1 = new AlignmentInterval(new SimpleInterval("20", 312610, 312840), 1, 231, TextCigarCodec.decode("231M127S"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        region2 = new AlignmentInterval(new SimpleInterval("20", 312579, 312801), 136, 358, TextCigarCodec.decode("135S223M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        expectedSimpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, false, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
        expectedDistances = new DistancesBetweenAlignmentsOnRefAndOnRead(-192, -96, 312801, 312610, 231, 136);
        result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeqForComplexExpansionNoPseudoHomology_reverseStrand, true, expectedSimpleChimera, expectedNovelAdjacencyAndAltHaplotype, expectedSVTypes, expectedVariants, expectedDistances, BreakpointsInference.SmallDuplicationWithImpreciseDupRangeBreakpointsInference.class));

        return result;
    }

    /**
     * See {@link #forComplexTandemDuplication()} .
     * Here we are simply making the repeat sequence shorter than 50 bases so that the end expectedSVTypes will be INS instead of DUP
     */
    private static List<TestDataForSimpleSV>
    forComplexTandemDuplicationIns() {

        final List<TestDataForSimpleSV> result = new ArrayList<>();
        final String leftRefFlank       = "TGCCAGGTTACATGGCAAAGAGGGTAGATAT";              // 31
        final String rightRefFlank      = "TGGTGCAAATGCCATTTATGCTCCTCTCCACCCATATCC";      // 39
        final String firstRepeat        = "GGGGAGCTGTGAAGAATGGAGCCAGTAATTAAATTCACTGAA";   // 42
        final String secondRepeat       = "GGGGAGCTGTGAAGAATGGAGCCAGTAATTAAATTCACTGAA";   // 42
        final String pseudoHomology     = "GGGCAGCTGTGGA";                                // 13


        // first test : expansion from 1 unit to 2 units with pseudo-homology
        String contigName = "cpx_dup_exp_small_1_2_pseudo_+";
        final byte[] fakeRefSeqForComplexExpansionWithPseudoHomology = String.format("%s%s%s%s", leftRefFlank, firstRepeat, pseudoHomology, rightRefFlank).getBytes();
        final byte[] contigSeqForComplexExpansionWithPseudoHomology = String.format("%s%s%s%s%s", leftRefFlank, firstRepeat, secondRepeat, pseudoHomology, rightRefFlank).getBytes();
        SimpleInterval expectedLeftBreakpoint = new SimpleInterval("20", 312609, 312609);
        SimpleInterval expectedRightBreakpoint = new SimpleInterval("20", 312609, 312609);
        BreakpointComplications expectedBreakpointComplications = new BreakpointComplications.SmallDuplicationWithImpreciseDupRangeBreakpointComplications(pseudoHomology, "",
                new SimpleInterval("20", 312610, 312651), 1, 2, Collections.singletonList(Strand.POSITIVE), Arrays.asList(Strand.POSITIVE, Strand.POSITIVE),
                new SimpleInterval("20", 312610, 312664));
        byte[] expectedAltSeq = String.format("%s%s%s", firstRepeat, secondRepeat, pseudoHomology).getBytes();
        NovelAdjacencyAndAltHaplotype expectedNovelAdjacencyAndAltHaplotype = new NovelAdjacencyAndAltHaplotype(expectedLeftBreakpoint, expectedRightBreakpoint, NO_SWITCH, expectedBreakpointComplications, SMALL_DUP_CPX, expectedAltSeq);
        AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("20", 312579, 312664), 1 ,86, TextCigarCodec.decode("86M81S"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("20", 312610, 312703), 74 ,167, TextCigarCodec.decode("73S94M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        SimpleChimera expectedSimpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, true, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
        DistancesBetweenAlignmentsOnRefAndOnRead expectedDistances = new DistancesBetweenAlignmentsOnRefAndOnRead(-55, -13, 312664, 312610, 86, 74);
        List<SvType> expectedSVTypes = Collections.singletonList(makeInsertionType(new SimpleInterval("20:312609-312609"), Allele.create("T", true), 42));
        List<VariantContext> expectedVariants = Collections.singletonList(
                addStandardAttributes(makeInsertion("20", 312609, 312609, 42, Allele.create("T", true)), 86 - pseudoHomology.length(), contigName, SimpleSVType.SupportedType.INS.name(), 312609, 42, StringUtil.bytesToString(expectedAltSeq), "", "")
                        .attribute(HOMOLOGY, pseudoHomology).attribute(HOMOLOGY_LENGTH, pseudoHomology.length())
                        .attribute(DUP_REPEAT_UNIT_REF_SPAN, "20:312610-312651").attribute(DUP_IMPRECISE_AFFECTED_RANGE, "20:312610-312664").attribute(DUP_ORIENTATIONS, "++").attribute(DUPLICATION_NUMBERS, "1,2").attribute(DUP_TAN_EXPANSION_STRING, "").attribute(DUP_ANNOTATIONS_IMPRECISE, "").make()
        );
        result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeqForComplexExpansionWithPseudoHomology, false, expectedSimpleChimera, expectedNovelAdjacencyAndAltHaplotype, expectedSVTypes, expectedVariants, expectedDistances, BreakpointsInference.SmallDuplicationWithImpreciseDupRangeBreakpointsInference.class));

        contigName = "cpx_dup_exp_small_1_2_pseudo_-";
        final byte[] contigSeqForComplexExpansionWithPseudoHomology_reverseStrand = Arrays.copyOf(contigSeqForComplexExpansionWithPseudoHomology, contigSeqForComplexExpansionWithPseudoHomology.length);
        SequenceUtil.reverseComplement(contigSeqForComplexExpansionWithPseudoHomology_reverseStrand);
        expectedNovelAdjacencyAndAltHaplotype = new NovelAdjacencyAndAltHaplotype(expectedLeftBreakpoint, expectedRightBreakpoint, NO_SWITCH, expectedBreakpointComplications, SMALL_DUP_CPX, expectedAltSeq);
        region1 = new AlignmentInterval(new SimpleInterval("20", 312610, 312703), 1 ,94, TextCigarCodec.decode("94M73S"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        region2 = new AlignmentInterval(new SimpleInterval("20", 312579, 312664), 82 ,167, TextCigarCodec.decode("81S86M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        expectedSimpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, false, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
        expectedDistances = new DistancesBetweenAlignmentsOnRefAndOnRead(-55, -13, 312664, 312610, 94, 82);
        result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeqForComplexExpansionWithPseudoHomology_reverseStrand, true, expectedSimpleChimera, expectedNovelAdjacencyAndAltHaplotype, expectedSVTypes, expectedVariants, expectedDistances, BreakpointsInference.SmallDuplicationWithImpreciseDupRangeBreakpointsInference.class));

        // second test: expansion from 2 units to 3 units without pseudo-homology
        contigName = "cpx_dup_exp_too_small_2_3_+";
        final byte[] contigSeqForComplexExpansionNoPseudoHomology = String.format("%s%s%s%s%s", leftRefFlank, firstRepeat, secondRepeat, firstRepeat, rightRefFlank).getBytes();
        expectedLeftBreakpoint = new SimpleInterval("20", 312609, 312609);
        expectedRightBreakpoint = new SimpleInterval("20", 312609, 312609);
        expectedBreakpointComplications = new BreakpointComplications.SmallDuplicationWithImpreciseDupRangeBreakpointComplications(firstRepeat, "",
                new SimpleInterval("20", 312610, 312651), 2, 3, Arrays.asList(Strand.POSITIVE, Strand.POSITIVE), Arrays.asList(Strand.POSITIVE, Strand.POSITIVE, Strand.POSITIVE),
                new SimpleInterval("20", 312610, 312693));
        expectedAltSeq= String.format("%s%s%s", firstRepeat, secondRepeat, firstRepeat).getBytes();
        expectedNovelAdjacencyAndAltHaplotype = new NovelAdjacencyAndAltHaplotype(expectedLeftBreakpoint, expectedRightBreakpoint, NO_SWITCH, expectedBreakpointComplications, SMALL_DUP_CPX, expectedAltSeq);
        region1 = new AlignmentInterval(new SimpleInterval("20", 312579, 312693), 1, 115, TextCigarCodec.decode("115M81S"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        region2 = new AlignmentInterval(new SimpleInterval("20", 312610, 312732), 74, 196, TextCigarCodec.decode("73S123M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        expectedSimpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, true, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
        expectedDistances = new DistancesBetweenAlignmentsOnRefAndOnRead(-84, -42, 312693, 312610, 115, 74);
        expectedSVTypes = Collections.singletonList(makeInsertionType(new SimpleInterval("20:312609-312609"), Allele.create("T", true), 42));
        expectedVariants = Collections.singletonList(
                addStandardAttributes(makeInsertion("20", 312609, 312609, 42, Allele.create("T", true)), 115 - firstRepeat.length(), contigName, SimpleSVType.SupportedType.INS.name(), 312609, 42, StringUtil.bytesToString(expectedAltSeq), "", "")
                        .attribute(HOMOLOGY, firstRepeat).attribute(HOMOLOGY_LENGTH, firstRepeat.length())
                        .attribute(DUP_REPEAT_UNIT_REF_SPAN, "20:312610-312651").attribute(DUP_IMPRECISE_AFFECTED_RANGE, "20:312610-312693").attribute(DUP_ORIENTATIONS, "+++").attribute(DUPLICATION_NUMBERS, "2,3").attribute(DUP_TAN_EXPANSION_STRING, "").attribute(DUP_ANNOTATIONS_IMPRECISE, "").make()
        );
        result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeqForComplexExpansionNoPseudoHomology, false, expectedSimpleChimera, expectedNovelAdjacencyAndAltHaplotype, expectedSVTypes, expectedVariants, expectedDistances, BreakpointsInference.SmallDuplicationWithImpreciseDupRangeBreakpointsInference.class));

        contigName = "cpx_dup_exp_too_small_2_3_-";
        final byte[] contigSeqForComplexExpansionNoPseudoHomology_reverseStrand = Arrays.copyOf(contigSeqForComplexExpansionNoPseudoHomology, contigSeqForComplexExpansionNoPseudoHomology.length);
        SequenceUtil.reverseComplement(contigSeqForComplexExpansionNoPseudoHomology_reverseStrand);
        expectedNovelAdjacencyAndAltHaplotype = new NovelAdjacencyAndAltHaplotype(expectedLeftBreakpoint, expectedRightBreakpoint, NO_SWITCH, expectedBreakpointComplications, SMALL_DUP_CPX, expectedAltSeq);
        region1 = new AlignmentInterval(new SimpleInterval("20", 312610, 312732), 1, 123, TextCigarCodec.decode("123M73S"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        region2 = new AlignmentInterval(new SimpleInterval("20", 312579, 312693), 82, 196, TextCigarCodec.decode("81S115M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        expectedSimpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, false, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
        expectedDistances = new DistancesBetweenAlignmentsOnRefAndOnRead(-84, -42, 312693, 312610, 123, 82);
        result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeqForComplexExpansionNoPseudoHomology_reverseStrand, true, expectedSimpleChimera, expectedNovelAdjacencyAndAltHaplotype, expectedSVTypes, expectedVariants, expectedDistances, BreakpointsInference.SmallDuplicationWithImpreciseDupRangeBreakpointsInference.class));

        return result;
    }
}
