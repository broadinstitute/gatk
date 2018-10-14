package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceMultiSparkSource;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.BreakEndVariantType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SvType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.TestUtilsForAssemblyBasedSVDiscovery;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.ContigAlignmentsModifier;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.StrandSwitch;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import static org.broadinstitute.hellbender.tools.spark.sv.discovery.TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments.NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.TypeInferredFromSimpleChimera.INTRA_CHR_STRAND_SWITCH_33;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.TypeInferredFromSimpleChimera.INTRA_CHR_STRAND_SWITCH_55;
import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.*;

/**
 * Holding test data and expected values for inversion breakpoints
 * NOT FOR REAL SYMBOLIC INVERSIONS WHOSE BREAKPOINTS AND FLANKING COPY NUMBER EVENTS ARE FULLY RESOLVED.
 *
 * Extracted here because inversions are treated differently
 * by {@link SimpleNovelAdjacencyInterpreter}, which outputs BND's,
 * and by {@link ContigChimericAlignmentIterativeInterpreter}, which was historically developed first and outputs symbolic inversion calls.
 */
public class AssemblyBasedSVDiscoveryTestDataProviderForInversionBreakpoints extends AssemblyBasedSVDiscoveryTestDataProvider {

    public static final class TestDataForInversion extends AssemblyBasedSVDiscoveryTestDataForSimpleChimera {

        private TestDataForInversion(final AlignmentInterval firstAlignment, final AlignmentInterval secondAlignment, final String evidenceAssemblyContigName, final byte[] evidenceContigSeq,
                                     final boolean expectedFirstContigRegionHasLaterReferenceMapping,
                                     final SimpleChimera expectedSimpleChimera,
                                     final NovelAdjacencyAndAltHaplotype expectedNovelAdjacencyAndAltSeq,
                                     final List<SvType> expectedSvTypes,
                                     final List<VariantContext> expectedVariantContexts) {
            super(firstAlignment, secondAlignment, evidenceAssemblyContigName, evidenceContigSeq, expectedFirstContigRegionHasLaterReferenceMapping, expectedSimpleChimera, expectedNovelAdjacencyAndAltSeq, expectedSvTypes, expectedVariantContexts,
                    BreakpointsInference.IntraChrStrandSwitchBreakpointInference.class);
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
            return BreakpointsInference.IntraChrStrandSwitchBreakpointInference.class;
        }
    }

    // the chromosome that the long contig1 is supposed to be mapped to is actually chr19, but to make tests runnable, we could only use "20" or "21"
    public static final String chrForLongContig1 = "20";
    public static final String homologyForLongContig1 = "GGAGGGAGGTGGGGGGGGTCAGCCCCCCGCCCAGCCAGCCGCCCCGTCCGGGAGGTGGGGGGTGCCTCTGCCCGGCCGCCCCTACTGGGAAGTGAGGAGCCCCTCTGCCCGGCCACCACCCCGTCTGGGAGGTGTACCCAACAGCTCATTGAGAACGGGCCATGATGACAATGGCAGTTTTGTGGAATAGAAACGGGGGAAAGGTGGGGAAAAGATTGAGAAATCGGATGGTTGCTGTGTCTGTGTAGAAAGAGGTAGACATGGGAGACTTCATTTTGTTCTGTACTAAGAAAAATTCTTCTGCCTTGGGATCCTGTTGACCTATGACCTTACCCCCAACCCTGTGCTCTCTGAAACATGTGCTGTGTCCACTCAGGGTTAAATGGATTAAGGGCGGTGCAAGATGTGCTTTGTTAAACAGATGCTTGAAGGCAGCATGCTCGTTAAGAGTCATCACCACTCCCTAATCTCAAGTACCCGGGGACACAAACACTGCGGAAGGCCGCAGGGTCCTCTGCCTAGGAAAACCAGAGACCTTTGTTCACTTGTTTATCTGCTGACCTTCCCTCCACTGTTGTCCTATGACCCTGCCAAATCCCCCTCTGCGAGAAACACCCAAGAATGATCAACTAAAAAAAAAAAAGAAAAAAGAAA";
    public static final String LONG_CONTIG1 =
            "TTTTTTTTTTTTTTTCTGAGACCGAATCTCGCTCTGTCACCCAGGCTGGAGTGCAGTGGCACGATCTTGGCTTACTGCAAGCTCTGCCTCCTGGGTTCATGCCATTCTCCTGCCTCAGCCCCACCCCCCCACCCCCCCAGGTAGCTG" +
                    "GGACTACAGGTGTCTGCCACCACACCTGGCTAAATTTTTTTGTATTTTTAGTAGAGACGGGGTTTCACCGTGTTAGCCAAGATGGTTTCGCTCTCCTGACCTCGCGATCCGCCCACCTCGGCCTCTCAAAGTGCTGGGATTACAGGCCTGAGCCACTGCGCCC" +
                    "TGCCTGACCCCTTCTTTTAAAACAAATCTTTTGTCTTTGTCTTCATTTCTGCATTCGTCCCCTTCGTTCAATCCTGTAGGAATTGACAGTGATATTGGGGAGTTCCCATCCCTGGATTTGGGATTTCCTCGAGTTTCCAGCCCTGTCCTTGTGGCCAAAAAGT" +
                    "GGCACAGTTTGGTTTTTGAATCCTGTTCGTGACGCCAAAAATATTCTGTGTGGGAAACATTCAAGGGGAGAAGAAAAGACACACACACAATACCTTTAAGGGTAAATAAGCTTTGTGCCACGTAAATGGAAATGCAGATATAATTAGCAAATTATATAATAAG" +
                    "CAAATCAATGTAATAACCAGATTGATATAATAAGCATGTTGATATAATAAGCAAATTGCAATGGGAAGAGGAGAAAGGAAAAGAGATATATATATTTACACTCACCAGACTATGGAGGATTCACCACCAGACTGGGAAGCAACAGCTTGGGCTCCAGAGTCAG" +
                    "CTACTCATTCCTGCACAGATGAGGAGGGTCTAATGAAGCTTCAGCACAATCTGATACCCTAGCTCTTTTGTAACGAGTTGTTTGGCATAAGGCCCAGTCATGAGGGCCATTTGCAACTGGGGTCAAGGAACACAAAACTGTCAACTTGTTTTTGCGATTGTCT" +
                    "ATTGTTTTTCAACAACTAATATATAGAAATAGATTGAAATAGAGATTTCTCTGAAACAGCGCTGAATGAATGCCTCAAGGGCCTCACACAACCTGTTTGAGAACTTGGTGACTACTGAGTTTGTCCACGTTCAATTAAGTTCAAATTTAGTATTTAACTTTTC" +
                    "CTCCACAAATTACCCAGTCTTGATGAATCAGCTCTGTCTAGACATAGGGCAAGGTGAACCCCTTGGGCAGTTACACAACCTCCGCCTTCTGGGTTTAAGCAATTCTCCTGCCTCAGCCTCCGGACTAGCTGGGTCTACAGGTGTGCAGCACCACACCCAGCTA" +
                    "GTTATTTGTACTTTTAGTAGAAATGGGGTTTCACCATGTTGGCCAGGCTGGTCTTGAACTCCTGACCTCAAGTGATCCACCCACCTTGGCCTCCCAAAGTGTTGCGATTACAGGCACTTGCCAGTGAACCTGGCCCTAAATGACTTCTTTCTATCTCCTATAA" +
                    "CAATTTGAAATTACTTAAAGGTGGTTTCAAATTGAAAAAATAAAAAGAATTTGGATAAAAATAAAATGTAAACAGTTTTTAAAAATTACAAGAGATTACAAAATATACATGTAAAACCGGAGTGGTCAAAAATGACAAATTTGATTTATTTATAAGGTTTATT" +
                    "AAAATTAGCTTTAGTATTGATAATACACTATTACCAAAGTAAAAGCTGATTTTCTCTTGAAAAAAATTTTATGTATTATTAATATGACAGCAAAATACTTCTGTTCACCTTTTGAATATATTCAAAAAGAGAGAGAGTAAAAAAAAAGATAGAATTTTCCCAT" +
                    "GATCTGGGTTGGGCCTGGCTTAGCTCAGGGAGGAAGCCCTGCCTGAAAAACGCTGCAGCTTAGGCTGTGACTTTTTCTTCACTCAGTCCTGATCACATTTTCTGTTACTCAGGGCCTGAGGGGGCGGGGGCCTTAAGCATTATCCAATCAGAAACGCTGGGCT" +
                    "GACGCCCCGTCCGGGAGGGAGGTGGGGGGGGTCAGCCCCCCGCCCAGCCAGCCGCCCCGTCCGGGAGGTGGGGGGTGCCTCTGCCCGGCCGCCCCTACTGGGAAGTGAGGAGCCCCTCTGCCCGGCCACCACCCCGTCTGGGAGGTGTACCCAACAGCTCATT" +
                    "GAGAACGGGCCATGATGACAATGGCAGTTTTGTGGAATAGAAACGGGGGAAAGGTGGGGAAAAGATTGAGAAATCGGATGGTTGCTGTGTCTGTGTAGAAAGAGGTAGACATGGGAGACTTCATTTTGTTCTGTACTAAGAAAAATTCTTCTGCCTTGGGATC" +
                    "CTGTTGACCTATGACCTTACCCCCAACCCTGTGCTCTCTGAAACATGTGCTGTGTCCACTCAGGGTTAAATGGATTAAGGGCGGTGCAAGATGTGCTTTGTTAAACAGATGCTTGAAGGCAGCATGCTCGTTAAGAGTCATCACCACTCCCTAATCTCAAGTA" +
                    "CCCGGGGACACAAACACTGCGGAAGGCCGCAGGGTCCTCTGCCTAGGAAAACCAGAGACCTTTGTTCACTTGTTTATCTGCTGACCTTCCCTCCACTGTTGTCCTATGACCCTGCCAAATCCCCCTCTGCGAGAAACACCCAAGAATGATCAACTAAAAAAAA" +
                    "AAAAGAAAAAAGAAAAAAGAAAAAATACACCCTGGCAGAAAACTTCTCGAATGTAAGTAAAGTGGGAGTTTTCAACTTGCCCACCTCCTTACAAAAACATGCAATAAGAATGTGCTAGATATAAAGCTTACAAATGTGAGAAATACGGGAAGTTTTTCAATTT" +
                    "TAACAAATACTTTCAAAGTCATGTGAAAATTCTCACTAGAAAGAAATCCAATAAAGGCAAGTTGAGGCATAAAAAAAAAAAAAAGAAAAAGCTCAGAATACATATAAGATTGACTAAGTCAGCCAGAAAATAATCCCCTAAAAGAAATTTCTCTCTAAACACC" +
                    "CAATGTGCACAGCTACTCTCAGCATGAGAAACATGAGCTTTATGAAGAAAGGTGGCAGATTTTCAGAAGAATTTTATAAAAGTTTCTTTTCCATCTCTGCTGCTCTCTCATCTCCTAGCCATTGAATGGGGGTTCTATATTGAAATACATCTGACAACTTCCA" +
                    "ACAACACTTTTTGATCAAGAAATAGAATTTGACTATGTTCGTATAGTGGAATATATTAGAACTTGTAACACAGCTAACTGAATAGCTATTATGGTGTTTGGGTGGCCACATCACCTGTCTTTATTTGTCCGGTAATAGCAGCATTCCAATTTAAAGAAATAAA" +
                    "AGATACCAAAATTGTGTTTACTTTTAATTATTCCTATTGAATAAAGTAATAAGCATGTCAGACTGATATCTATCATAACAATAAATTTTGTTTGGATATTATATTAGATATAAATATTTAAGTATGAATAATTTTAATGAACTAGTCATAATGTATGTAGCAT" +
                    "TTTTAAAAATTGTAACTATACTTCAGTTAAAACACTTTATATTTCAAAAGCATAAATAACAATATTAAAATAACAATTTAGGTGATTCATTCAAAGTAAGTATTGGGGCTTTATGTTCATACTATTGTAGAAAATACTGCTTATGGCTCATGCCTGTAATCCC" +
                    "AGCACATTGGGAGGCTGAGGTGGGTAGATCACCTGAGGTCAGGAGTTCCTGATCCCATCTTTACTAAAAATACAAAACTTACCCAGGGTGGTTGTGCACACTTGTAATCCCAGCTACTTGGGAGGCTGAGGCAGGAGAATTGCTTGAACAAGGGAGGAAATGG" +
                    "TTGCAGTGAGCCATGATCATGCCACTGAACCCCAGCCTGGGCAAGAGAGTGAGACTGTCTCAAAAAAAAAAAAAACTGTTTAATTTTTATGAATGCAGGTTTTCTGCAAACACTACACATAACTATGCTAATTGTTCTGAAGTAATAAATAGAAAGCAAGGCA" +
                    "CAACTACAGACTCCACTGTTCAGTTTATGCACTGAACTGTTCTTGCTTTTGCAGTGTAAGTATTTCTGCCTGCAAATACTGGATAATTACCTTGGATCATCAGATTTCTATCAAAGGAATTTAGTATCTTTTAGTCTTTATCATTTTGTATTGCTAAATTTAT" +
                    "CTGTGTGTTAAGCTTCTGTGTGCTCTTAAAATGAGGTTTTATCTAAACAAACCTGTGTCTACTTTAAAAGACTAAACATGAAAAAACTAAACTTTTCAGAACCAAAAACAAAGCAATAAATCTGAAGTACTAGATAGTCTGGAGTGAGATTTATTTAGCTTTT" +
                    "TTTTTTTTTTGAGATGGAGTCTCGCTCTGTCACCGAGGCTGGAGTGCAGTGGCACGAACTCGGCTCACTGCAAAAGCTCTGCCTCCCAGCTTCATGCCATTCTCCTACCTCAGCCTCCCAAGTAGCTGGGATTACAGGCAACTGCCACCACGCCCAGCTAATT" +
                    "TTTTTGTATTTTTAGTAGAGACGGGGTTTCACCGTGTTAGCCAGGATGGTCTTGATCTCCTGACCTCATGATCTGCCTGCCTCGGCCTCCCAAAGTGCTGGGATTACAGGCATGCAAAACCGCGCCCTGCCCTTATTTAGCTCTTAATATGATTTACATATAT" +
                    "TTCAAAAAAGCAGAGAAAAATATCTACATATAATCTAAATCCCTTAAGAGAGAATAGCAAAAAAATTTTTGGAACTTTTTAAAGAGTTTTGAACTCTTGGACATCTGAATTTTGCACACTATATGCACTTGAAAGAATGTTAATGGGGAAAAAGCAGAAGAGA" +
                    "AAAAGACGTTATAAAAAAATCCATGAGTGCACAACACCTATGAAACAGAATAGAGAGCCCAGAAATAATGCCTTTCACCTGCAACCATCAGATTTCTGACAAAGCTGACAAGAGGAATGTGGGAAGAATTCTCTCTTTCATAAATGGTGCTGGAATAACTATC" +
                    "TACCACTATGTAGAAGACTGAAGTGGACCCCTTCATTACACCATATAAAAAAATCAACTGAAGATAAATTAAGGACTTAAATGTAAAACTTAAAATTATAAGAAACCCTGCAAGATAACCTAGGAAATAGCATTCTAGACACAGAAACAGGTAAAGACTTCAT" +
                    "GATGAAGCTACCAAAAGCAACTGCAACAGAAGTAAATTGACAAATGGGATGTATTTAAACTTAAGAGCTTCTTCACAGCAAAGGAAACTATCAACAGAGTAAACAGACAAACTAGAGAATAAAAGAATATATTTGTAAATTTTGCCTCTGAAAAAGGTCTAAT" +
                    "ATACAGAATTTATTAGGAACTTAAACAAGTTTACAAGAAAAAAACACACTCATTAAAAAGTATGCAAAAAACATGAACAGATGCCTTTCAATAGAAGATACACATAGGGCTAACAAGCATATGAAAAAAAAAATGCTCATTGCTAATCATTAGAGAAATTAGA" +
                    "AGAGAAATCACAGGCTGGGTGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCAAGGCAGGCAGATTACAAGGTCAGGAGATCAAGATCATCCTGGCTAACATGGTGAAACCCTGTCTCTACTAAAAATACAAAAATCAGCCAGGTGTGGCAGTGT" +
                    "GCACCTGTAGTCCCAGCTACTCAGGAGGCTGAGGCAGGAGAATTGCTTGAATCTGGTAGGCAGAGGTTGCAGTGAGCTGAGATCACACCACTGCACTCCTGCCTGGGCAACAGAGCAAGACTCCGTCTCAAACACACACACACAGACACACACACACACACAC" +
                    "ACACACACACACACACACACGCAGAGAAACCACAATGAGATACCACCTCATACCAGTCAGAATGGCTATTTTTAAAAAGTCAAAAGATAACAGATGCTGACAAAGTTGCAGAGAAAAGGGAATGCTTATACTCCTCTGGTGGGAGTGTAAATTATTTCAACAA" +
                    "CTGTAAAAATCAGTGTGACAATTCCTCACAGAATGAAAAACAGAATTATCATTCGACTGGGAAACCTCATAATTGAGTATATACCCAAAGAAATATGAAATATTATAAAGACACATCCACATGCATGTTCACTGCAGCACTATTCACAATAGCAAAGACACGG" +
                    "ACAGACTAAATGCCTATCAATGGCAGACTGGATCAAGAAAATATGGTATGGTCAGATGCGGTGGCTCATGCCTGTAATTCCAGCCCTTTGGGAGGCTGAGGCAGGTGGATTGCCTGAGCTTAGAAGTTTGAGACCACTCTGGGCAACATGGCAAAATTTTGTC" +
                    "TCCACAGAAGATACAAAAAAAAAAAAAAAAAA";
    private final TestDataForInversion forSimpleInversionWithInsertion_LeftBreakpoint_plus;
    final TestDataForInversion forSimpleInversionWithHomology_RightBreakpoint_minus;
    final TestDataForInversion forSimpleInversionWithHom_leftPlus;
    private final TestDataForInversion forSimpleInversionWithHom_leftMinus;
    private final TestDataForInversion forSimpleInversionWithHom_rightPlus;
    private final TestDataForInversion forSimpleInversionWithHom_rightMinus;

    public AssemblyBasedSVDiscoveryTestDataProviderForInversionBreakpoints() {
        forSimpleInversionWithInsertion_LeftBreakpoint_plus = forSimpleInversionWithInsertion_LeftBreakpoint_plus();
        forSimpleInversionWithHomology_RightBreakpoint_minus = forSimpleInversionWithHomology_RightBreakpoint_minus();
        final List<TestDataForInversion> inversion4 = forSimpleInversionWithHomology();
        forSimpleInversionWithHom_leftPlus = inversion4.get(0);
        forSimpleInversionWithHom_leftMinus = inversion4.get(1);
        forSimpleInversionWithHom_rightPlus = inversion4.get(2);
        forSimpleInversionWithHom_rightMinus = inversion4.get(3);
    }

    public List<AssemblyBasedSVDiscoveryTestDataForSimpleChimera> getAllTestData() {
        final List<TestDataForInversion> testDataForInversions =
                Arrays.asList(
                        forSimpleInversionWithInsertion_LeftBreakpoint_plus,
                        forSimpleInversionWithHomology_RightBreakpoint_minus,
                        forSimpleInversionWithHom_leftPlus,
                        forSimpleInversionWithHom_leftMinus,
                        forSimpleInversionWithHom_rightPlus,
                        forSimpleInversionWithHom_rightMinus);
        return Collections.unmodifiableList(testDataForInversions);
    }

    private static TestDataForInversion
    forSimpleInversionWithInsertion_LeftBreakpoint_plus() {
        // inversion with inserted sequence
        String contigName = "simple_inv_55_with_ins_+";
        final byte[] leftFlank = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence(146, (byte)'A');
        final byte[] rightFlankRC = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence(50, (byte)'C');
        final byte[] contigSeq = new byte[leftFlank.length+1+rightFlankRC.length];
        System.arraycopy(leftFlank, 0, contigSeq, 0, leftFlank.length);
        contigSeq[leftFlank.length] = (byte) 'T';
        System.arraycopy(rightFlankRC, 0, contigSeq, leftFlank.length+1, rightFlankRC.length);

        // reference intervals are changed from real values, which were generated on a different chromosome, to a fake but valid value.
        final AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("21", 17069151, 17069296), 1, 146, TextCigarCodec.decode("146M51S"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        final AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("21", 17069315, 17069364), 148, 197, TextCigarCodec.decode("147S50M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        final AlignedContig alignedContig = new AlignedContig(contigName, contigSeq, Arrays.asList(region1, region2));
        final SimpleChimera expectedSimpleChimera = new SimpleChimera(contigName, region1, region2, StrandSwitch.FORWARD_TO_REVERSE, true, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);

        final SimpleInterval expectedLeftBreakpoint = new SimpleInterval("21", 17069296, 17069296);
        final SimpleInterval expectedRightBreakpoint = new SimpleInterval("21", 17069364, 17069364);
        final BreakpointComplications expectedBreakpointComplications = new BreakpointComplications.IntraChrStrandSwitchBreakpointComplications("", "T");

        final NovelAdjacencyAndAltHaplotype expectedNovelAdjacencyAndAltHaplotype = new NovelAdjacencyAndAltHaplotype(expectedLeftBreakpoint, expectedRightBreakpoint,
                StrandSwitch.FORWARD_TO_REVERSE, expectedBreakpointComplications, INTRA_CHR_STRAND_SWITCH_55,  "T".getBytes());

        final List<SvType> expectedSVTypes = Arrays.asList(
                makeBNDType("21", 17069296, "BND_INV55_21_17069296_17069364_1", Allele.create("A", true), Allele.create("AT]21:17069364]"), Collections.singletonMap(INV55, true), true, BreakEndVariantType.SupportedType.INTRA_CHR_STRAND_SWITCH_55),
                makeBNDType("21", 17069364, "BND_INV55_21_17069296_17069364_2", Allele.create("G", true), Allele.create("GA]21:17069296]"), Collections.singletonMap(INV55, true), false, BreakEndVariantType.SupportedType.INTRA_CHR_STRAND_SWITCH_55)
        );

        final VariantContext expectedFirstMate = makeBND(new SimpleInterval("21:17069296-17069296"), new SimpleInterval("21:17069364-17069364"), Allele.create("A", true), "T", INV55, true, true, true)
                .attribute(INV55, true).attribute(INSERTED_SEQUENCE, "T").attribute(INSERTED_SEQUENCE_LENGTH, 1).attribute(TOTAL_MAPPINGS, 1).attribute(HQ_MAPPINGS, 1).attribute(MAPPING_QUALITIES, 60).attribute(ALIGN_LENGTHS, 50).attribute(MAX_ALIGN_LENGTH, 50).attribute(CONTIG_NAMES, contigName)
                .make();
        final VariantContext expectedSecondMate = makeBND(new SimpleInterval("21:17069296-17069296"), new SimpleInterval("21:17069364-17069364"), Allele.create("G", true), "A", INV55, false, true, true)
                .attribute(INV55, true).attribute(INSERTED_SEQUENCE, "T").attribute(INSERTED_SEQUENCE_LENGTH, 1).attribute(TOTAL_MAPPINGS, 1).attribute(HQ_MAPPINGS, 1).attribute(MAPPING_QUALITIES, 60).attribute(ALIGN_LENGTHS, 50).attribute(MAX_ALIGN_LENGTH, 50).attribute(CONTIG_NAMES, contigName).make();
        final List<VariantContext> expectedVariants = Arrays.asList(
                new VariantContextBuilder(expectedFirstMate).attribute(BND_MATEID_STR, expectedSecondMate.getID()).make(),
                new VariantContextBuilder(expectedSecondMate).attribute(BND_MATEID_STR, expectedFirstMate.getID()).make()
        );
        return new TestDataForInversion(region1, region2, contigName, contigSeq, false, expectedSimpleChimera, expectedNovelAdjacencyAndAltHaplotype, expectedSVTypes, expectedVariants);
    }

    private static TestDataForInversion
    forSimpleInversionWithHomology_RightBreakpoint_minus() {
        // inversion with strange left breakpoint
        final String contigName = "simple_inv_33_with_hom_-";
        final byte[] contigSequence = LONG_CONTIG1.getBytes(); // 6210 length
        final AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("20", 20262839, 20265443), 1, contigSequence.length - 3603, TextCigarCodec.decode("1970M1I611M1I24M3603H"), false, 60, 8, 2561, ContigAlignmentsModifier.AlnModType.NONE);
        final AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("20", 20248787, 20253040), 1954, contigSequence.length, TextCigarCodec.decode("1953S10M3I9M1I246M2D1572M1I798M5D730M1I347M4I535M"), true, 60, 41, 4068, ContigAlignmentsModifier.AlnModType.NONE);
        final AlignedContig alignedContig = new AlignedContig(contigName, contigSequence, Arrays.asList(region1, region2));

        final SimpleChimera expectedSimpleChimera = new SimpleChimera(contigName, region1, region2, StrandSwitch.REVERSE_TO_FORWARD, false, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);

        final SimpleInterval expectedLeftBreakpoint = new SimpleInterval("20", 20248787, 20248787);
        final SimpleInterval expectedRightBreakpoint = new SimpleInterval("20", 20263493, 20263493);
        final BreakpointComplications expectedBreakpointComplications = new BreakpointComplications.IntraChrStrandSwitchBreakpointComplications(homologyForLongContig1, "");

        final NovelAdjacencyAndAltHaplotype expectedNovelAdjacencyAndAltHaplotype = new NovelAdjacencyAndAltHaplotype(expectedLeftBreakpoint, expectedRightBreakpoint,
                StrandSwitch.REVERSE_TO_FORWARD, expectedBreakpointComplications, TypeInferredFromSimpleChimera.INTRA_CHR_STRAND_SWITCH_33, EMPTY_BYTE_ARRAY);

        final List<SvType> expectedSVTypes = Arrays.asList(
                makeBNDType("20", 20248787, "BND_INV33_20_20248787_20263493_1", Allele.create("C", true), Allele.create("[20:20263493[C"), Collections.singletonMap(INV33, true), true, BreakEndVariantType.SupportedType.INTRA_CHR_STRAND_SWITCH_33),
                makeBNDType("20", 20263493, "BND_INV33_20_20248787_20263493_2", Allele.create("G", true), Allele.create("[20:20248787[G"), Collections.singletonMap(INV33, true), false, BreakEndVariantType.SupportedType.INTRA_CHR_STRAND_SWITCH_33)
        );

        final VariantContext expectedFirstMate = makeBND(new SimpleInterval("20:20248787-20248787"), new SimpleInterval("20:20263493-20263493"), Allele.create("C", true), "", INV33, true, false, false)
                .attribute(INV33, true).attribute(HOMOLOGY, homologyForLongContig1).attribute(HOMOLOGY_LENGTH, homologyForLongContig1.length()).attribute(TOTAL_MAPPINGS, 1).attribute(HQ_MAPPINGS, 1).attribute(MAPPING_QUALITIES, 60).attribute(ALIGN_LENGTHS, 1951).attribute(MAX_ALIGN_LENGTH, 1951).attribute(CONTIG_NAMES, contigName).make();
        final VariantContext expectedSecondMate = makeBND(new SimpleInterval("20:20248787-20248787"), new SimpleInterval("20:20263493-20263493"), Allele.create("G", true), "", INV33, false, false, false)
                .attribute(INV33, true).attribute(HOMOLOGY, homologyForLongContig1).attribute(HOMOLOGY_LENGTH, homologyForLongContig1.length()).attribute(TOTAL_MAPPINGS, 1).attribute(HQ_MAPPINGS, 1).attribute(MAPPING_QUALITIES, 60).attribute(ALIGN_LENGTHS, 1951).attribute(MAX_ALIGN_LENGTH, 1951).attribute(CONTIG_NAMES, contigName).make();

        final List<VariantContext> expectedVariants = Arrays.asList(
                new VariantContextBuilder(expectedFirstMate).attribute(BND_MATEID_STR, expectedSecondMate.getID()).make(),
                new VariantContextBuilder(expectedSecondMate).attribute(BND_MATEID_STR, expectedFirstMate.getID()).make()
        );
        return new TestDataForInversion(region1, region2, contigName, contigSequence, true, expectedSimpleChimera, expectedNovelAdjacencyAndAltHaplotype, expectedSVTypes, expectedVariants);
    }

    /**
     * The following four tests are all going to be for the same inversion, testing if implementations are correct for
     * identifying the breakpoints by looking at different representations of evidence.
     * The inversion we are looking at is
     *
     * '+' strand representation: G....100....G|ACACA|C....100....C               A....100....A|TGTGT|T....100....T
     *
     * 100-bases of 'G' is the left flanking before the homology |ACACA| and the region starting with 100-bases of 'C' and
     * ending with 100-bases of 'A' and maybe (homologyForwardStrandRep uncertainty) the homology |TGTGT| is inverted.
     * 100-bases of 'T' is the right flanking region.
     *
     * Returns a list of four Tuple4's with left flanking evidence '+'/'-' strand representation and right flanking side.
     */
    private static List<TestDataForInversion>
    forSimpleInversionWithHomology() {

        final List<TestDataForInversion> result = new ArrayList<>();

        final String leftLeftPlus = makeDummySequence('G', 100);
        final String leftLeftMinus = makeDummySequence('C', 100);
        final String leftRightPlus = makeDummySequence('C', 100);
        final String leftRightMinus = makeDummySequence('G', 100);
        final String rightLeftPlus = makeDummySequence('A', 100);
        final String rightLeftMinus = makeDummySequence('T', 100);
        final String rightRightPlus = makeDummySequence('T', 100);
        final String rightRightMinus = makeDummySequence('A', 100);
        final String leftHomology = "ACACA";
        final String rightHomology = "TGTGT";
        {// left breakpoint '+'/'-' strand representation
            byte[] contigSeq = (leftLeftPlus + leftHomology + rightLeftMinus).getBytes();

            String contigName = "simple_inv_55_with_hom_+";

            AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("20", 200101, 200205), 1, 105, TextCigarCodec.decode("105M100S"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
            AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("20", 200501, 200605), 101, 205, TextCigarCodec.decode("100S105M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
            SimpleChimera expectedSimpleChimera = new SimpleChimera(contigName, region1, region2, StrandSwitch.FORWARD_TO_REVERSE, true, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
            final SimpleInterval expectedLeftBreakpoint = new SimpleInterval("20", 200200,200200);
            final SimpleInterval expectedRightBreakpoint = new SimpleInterval("20", 200605, 200605);
            final BreakpointComplications expectedBreakpointComplications = new BreakpointComplications.IntraChrStrandSwitchBreakpointComplications(leftHomology, "");
            final NovelAdjacencyAndAltHaplotype breakpoints = new NovelAdjacencyAndAltHaplotype(expectedLeftBreakpoint, expectedRightBreakpoint, StrandSwitch.FORWARD_TO_REVERSE, expectedBreakpointComplications, INTRA_CHR_STRAND_SWITCH_55, EMPTY_BYTE_ARRAY);
            final List<SvType> expectedSVTypes = Arrays.asList(
                    makeBNDType("20", 200200, "BND_INV55_20_200200_200605_1", Allele.create("C", true), Allele.create("C]20:200605]"), Collections.singletonMap(INV55, true), true, BreakEndVariantType.SupportedType.INTRA_CHR_STRAND_SWITCH_55),
                    makeBNDType("20", 200605, "BND_INV55_20_200200_200605_2", Allele.create("T", true), Allele.create("T]20:200200]"), Collections.singletonMap(INV55, true), false, BreakEndVariantType.SupportedType.INTRA_CHR_STRAND_SWITCH_55)
            );
            final VariantContext expectedFirstMate = makeBND(new SimpleInterval("20:200200-200200"), new SimpleInterval("20:200605-200605"), Allele.create("C", true), "", INV55, true, true, true)
                    .attribute(HOMOLOGY, leftHomology).attribute(HOMOLOGY_LENGTH, leftHomology.length()).attribute(INV55, true).attribute(TOTAL_MAPPINGS, 1).attribute(HQ_MAPPINGS, 1).attribute(MAPPING_QUALITIES, 60).attribute(ALIGN_LENGTHS, 100).attribute(MAX_ALIGN_LENGTH, 100).attribute(CONTIG_NAMES, contigName).make();
            final VariantContext expectedSecondMate = makeBND(new SimpleInterval("20:200200-200200"), new SimpleInterval("20:200605-200605"), Allele.create("T", true), "", INV55, false, true, true)
                    .attribute(HOMOLOGY, leftHomology).attribute(HOMOLOGY_LENGTH, leftHomology.length()).attribute(INV55, true).attribute(TOTAL_MAPPINGS, 1).attribute(HQ_MAPPINGS, 1).attribute(MAPPING_QUALITIES, 60).attribute(ALIGN_LENGTHS, 100).attribute(MAX_ALIGN_LENGTH, 100).attribute(CONTIG_NAMES, contigName).make();

            final List<VariantContext> expectedVariants = Arrays.asList(
                    new VariantContextBuilder(expectedFirstMate).attribute(BND_MATEID_STR, expectedSecondMate.getID()).make(),
                    new VariantContextBuilder(expectedSecondMate).attribute(BND_MATEID_STR, expectedFirstMate.getID()).make()
            );

            result.add(new TestDataForInversion(region1, region2, contigName, contigSeq, false, expectedSimpleChimera, breakpoints, expectedSVTypes, expectedVariants));

            contigSeq = (rightLeftPlus + rightHomology + leftLeftMinus).getBytes();
            contigName = "simple_inv_55_with_hom_-";
            region1 = new AlignmentInterval(new SimpleInterval("20", 200501, 200605), 1, 105, TextCigarCodec.decode("105M100S"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
            region2 = new AlignmentInterval(new SimpleInterval("20", 200101, 200205), 101, 205, TextCigarCodec.decode("100S105M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
            expectedSimpleChimera = new SimpleChimera(contigName, region1, region2, StrandSwitch.FORWARD_TO_REVERSE, false, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
            result.add(new TestDataForInversion(region1, region2, contigName, contigSeq, true, expectedSimpleChimera, breakpoints, expectedSVTypes, expectedVariants));
        }
        {// right breakpoint '+'/'-' strand representation
            byte[] contigSeq = (leftRightMinus + rightHomology + rightRightPlus).getBytes();

            String contigName = "simple_inv_33_with_hom_+";

            AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("20", 200201, 200305), 1, 105, TextCigarCodec.decode("105M100S"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
            AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("20", 200601, 200705), 101, 205, TextCigarCodec.decode("100S105M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
            SimpleChimera expectedSimpleChimera = new SimpleChimera(contigName, region1, region2, StrandSwitch.REVERSE_TO_FORWARD, true, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
            final SimpleInterval expectedLeftBreakpoint = new SimpleInterval("20", 200201,200201);
            final SimpleInterval expectedRightBreakpoint = new SimpleInterval("20", 200606, 200606);
            final BreakpointComplications expectedBreakpointComplications = new BreakpointComplications.IntraChrStrandSwitchBreakpointComplications(leftHomology, "");
            final NovelAdjacencyAndAltHaplotype breakpoints = new NovelAdjacencyAndAltHaplotype(expectedLeftBreakpoint, expectedRightBreakpoint, StrandSwitch.REVERSE_TO_FORWARD, expectedBreakpointComplications, INTRA_CHR_STRAND_SWITCH_33, EMPTY_BYTE_ARRAY);
            final List<SvType> expectedSVTypes = Arrays.asList(
                    makeBNDType("20", 200201, "BND_INV33_20_200201_200606_1", Allele.create("A", true), Allele.create("[20:200606[A"), Collections.singletonMap(INV33, true), true, BreakEndVariantType.SupportedType.INTRA_CHR_STRAND_SWITCH_33),
                    makeBNDType("20", 200606, "BND_INV33_20_200201_200606_2", Allele.create("G", true), Allele.create("[20:200201[G"), Collections.singletonMap(INV33, true), false, BreakEndVariantType.SupportedType.INTRA_CHR_STRAND_SWITCH_33)
            );
            final VariantContext expectedFirstMate = makeBND(new SimpleInterval("20:200201-200201"), new SimpleInterval("20:200606-200606"), Allele.create("A", true), "", INV33, true, false, false)
                    .attribute(HOMOLOGY, leftHomology).attribute(HOMOLOGY_LENGTH, leftHomology.length()).attribute(INV33, true).attribute(TOTAL_MAPPINGS, 1).attribute(HQ_MAPPINGS, 1).attribute(MAPPING_QUALITIES, 60).attribute(ALIGN_LENGTHS, 100).attribute(MAX_ALIGN_LENGTH, 100).attribute(CONTIG_NAMES, contigName).make();
            final VariantContext expectedSecondMate = makeBND(new SimpleInterval("20:200201-200201"), new SimpleInterval("20:200606-200606"), Allele.create("G", true), "", INV33, false, false, false)
                    .attribute(HOMOLOGY, leftHomology).attribute(HOMOLOGY_LENGTH, leftHomology.length()).attribute(INV33, true).attribute(TOTAL_MAPPINGS, 1).attribute(HQ_MAPPINGS, 1).attribute(MAPPING_QUALITIES, 60).attribute(ALIGN_LENGTHS, 100).attribute(MAX_ALIGN_LENGTH, 100).attribute(CONTIG_NAMES, contigName).make();
            final List<VariantContext> expectedVariants = Arrays.asList(
                    new VariantContextBuilder(expectedFirstMate).attribute(BND_MATEID_STR, expectedSecondMate.getID()).make(),
                    new VariantContextBuilder(expectedSecondMate).attribute(BND_MATEID_STR, expectedFirstMate.getID()).make()
            );
            result.add(new TestDataForInversion(region1, region2, contigName, contigSeq, false, expectedSimpleChimera, breakpoints, expectedSVTypes, expectedVariants));

            contigSeq = (rightRightMinus + leftHomology + leftRightPlus).getBytes();

            contigName = "simple_inv_33_with_hom_-";
            region1 = new AlignmentInterval(new SimpleInterval("20", 200601, 200705), 1, 105, TextCigarCodec.decode("105M100S"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
            region2 = new AlignmentInterval(new SimpleInterval("20", 200201, 200305), 101, 205, TextCigarCodec.decode("100S105M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
            expectedSimpleChimera = new SimpleChimera(contigName, region1, region2, StrandSwitch.REVERSE_TO_FORWARD, false, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
            result.add(new TestDataForInversion(region1, region2, contigName, contigSeq, true, expectedSimpleChimera, breakpoints, expectedSVTypes, expectedVariants));
        }
        return result;
    }
}
