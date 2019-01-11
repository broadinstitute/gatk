package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceMultiSparkSource;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryPipelineSpark;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.StrandSwitch;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.*;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.EvidenceTargetLink;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.ReadMetadata;
import org.broadinstitute.hellbender.tools.spark.sv.integration.SVIntegrationTestDataProvider;
import org.broadinstitute.hellbender.tools.spark.sv.utils.PairedStrandedIntervalTree;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.tools.spark.sv.utils.StrandedInterval;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.mockito.Mockito;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import scala.Tuple2;

import java.util.*;

import static org.broadinstitute.hellbender.tools.spark.sv.discovery.AnnotatedVariantProducer.produceAnnotatedVcFromAssemblyEvidence;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.AnnotatedVariantProducer.produceLinkedAssemblyBasedVariants;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments.NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME;
import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.*;
import static org.mockito.Mockito.when;
public class AnnotatedVariantProducerUnitTest extends GATKBaseTest {

    // -----------------------------------------------------------------------------------------------
    // Evidence summary annotation
    // -----------------------------------------------------------------------------------------------
    @DataProvider
    private Object[][] forAssemblyBasedAnnotation() {
        final List<Object[]> data = new ArrayList<>(20);

        final String testSample = "testSample";
        final JavaSparkContext testSparkContext = SparkContextFactory.getTestSparkContext();
        final Broadcast<ReferenceMultiSparkSource> referenceBroadcast = testSparkContext.broadcast(TestUtilsForAssemblyBasedSVDiscovery.b37_reference);
        final Broadcast<SAMSequenceDictionary> refSeqDictBroadcast = testSparkContext.broadcast(TestUtilsForAssemblyBasedSVDiscovery.b37_seqDict);

        for (final AssemblyBasedSVDiscoveryTestDataProvider.AssemblyBasedSVDiscoveryTestDataForSimpleChimera testData : new AssemblyBasedSVDiscoveryTestDataProviderForSimpleSV().getAllTestData()) {
            data.add(new Object[]{testData.expectedSvTypes, new SimpleNovelAdjacencyAndChimericAlignmentEvidence(testData.expectedNovelAdjacencyAndAltSeq, Collections.singletonList(testData.expectedSimpleChimera)), referenceBroadcast, refSeqDictBroadcast, testSample, LINK,
                                  testData.expectedVariantContexts});
        }
        for (final AssemblyBasedSVDiscoveryTestDataProvider.AssemblyBasedSVDiscoveryTestDataForSimpleChimera testData : new AssemblyBasedSVDiscoveryTestDataProviderForInversionBreakpoints().getAllTestData()) {
            data.add(new Object[]{testData.expectedSvTypes, new SimpleNovelAdjacencyAndChimericAlignmentEvidence(testData.expectedNovelAdjacencyAndAltSeq, Collections.singletonList(testData.expectedSimpleChimera)), referenceBroadcast, refSeqDictBroadcast, testSample, BND_MATEID_STR,
                                  testData.expectedVariantContexts});
        }

        final Broadcast<ReferenceMultiSparkSource> referenceBroadcast_b38 = testSparkContext.broadcast(TestUtilsForAssemblyBasedSVDiscovery.b38_reference_chr20_chr21);
        final Broadcast<SAMSequenceDictionary> refSeqDictBroadcast_b38 = testSparkContext.broadcast(TestUtilsForAssemblyBasedSVDiscovery.b38_seqDict_chr20_chr21);

        for (final AssemblyBasedSVDiscoveryTestDataProvider.AssemblyBasedSVDiscoveryTestDataForSimpleChimera testData : new AssemblyBasedSVDiscoveryTestDataProviderForBreakEndVariants().getAllTestData()) {
            data.add(new Object[]{testData.expectedSvTypes, new SimpleNovelAdjacencyAndChimericAlignmentEvidence(testData.expectedNovelAdjacencyAndAltSeq, Collections.singletonList(testData.expectedSimpleChimera)), referenceBroadcast_b38, refSeqDictBroadcast_b38, testSample, BND_MATEID_STR,
                                  testData.expectedVariantContexts});
        }

        return data.toArray(new Object[data.size()][]);
    }
    @Test(groups = "sv", dataProvider = "forAssemblyBasedAnnotation")
    public void testAssemblyBasedAnnotation(final List<SvType> inferredTypes,
                                            final SimpleNovelAdjacencyAndChimericAlignmentEvidence novelAdjacencyAndAssemblyEvidence,
                                            final Broadcast<ReferenceMultiSparkSource> broadcastReference,
                                            final Broadcast<SAMSequenceDictionary> broadcastSequenceDictionary,
                                            final String sampleId,
                                            final String linkKey,
                                            final List<VariantContext> expectedVariants) {
        if (inferredTypes.size() == 1) {
            VariantContextTestUtils.assertVariantContextsAreEqual(produceAnnotatedVcFromAssemblyEvidence(inferredTypes.get(0), novelAdjacencyAndAssemblyEvidence, broadcastReference, broadcastSequenceDictionary, null, sampleId).make(),
                    expectedVariants.get(0), Arrays.asList(CONTIG_NAMES, HQ_MAPPINGS)); // these two are omitted because: 1) HQ_MAPPINGS equal testing code is wrong (saying "1"!=1), 2)CONTIG_NAMES will be tested in another test
        } else if (inferredTypes.size() == 2){
            final List<VariantContext> variantContexts = produceLinkedAssemblyBasedVariants(new Tuple2<>(inferredTypes.get(0), inferredTypes.get(1)), novelAdjacencyAndAssemblyEvidence, broadcastReference, broadcastSequenceDictionary, null, sampleId, linkKey);
            VariantContextTestUtils.assertVariantContextsAreEqual(variantContexts.get(0),
                    expectedVariants.get(0), Arrays.asList(CONTIG_NAMES, HQ_MAPPINGS));
            VariantContextTestUtils.assertVariantContextsAreEqual(variantContexts.get(1),
                    expectedVariants.get(1), Arrays.asList(CONTIG_NAMES, HQ_MAPPINGS));
        }
    }

    @Test(groups = "sv")
    public void testMiscCases() {
        final JavaSparkContext testSparkContext = SparkContextFactory.getTestSparkContext();
        Broadcast<ReferenceMultiSparkSource> referenceBroadcast = testSparkContext.broadcast(TestUtilsForAssemblyBasedSVDiscovery.b38_reference_chr20_chr21);
        Broadcast<SAMSequenceDictionary> refSeqDictBroadcast = testSparkContext.broadcast(TestUtilsForAssemblyBasedSVDiscovery.b38_seqDict_chr20_chr21);

        // the following works for: multiple evidence contigs, insertion sequence mapping available, good non-canonical chromosome mapping available
        final AlignmentInterval asm018485_tig00004_1 = TestUtilsForAssemblyBasedSVDiscovery.fromSAMRecordString("asm018485:tig00004\t0\tchr20\t28831147\t38\t275M2I114M229S\t*\t0\t0\tTATCTTCACATAAAAACTACACAGTATCATTCTGAGAAACTTGTTTGTGATGTGTGCATTCATCTCACAGATTTGAACCCTTCCATCTTTTGAGCAGTTTGTACACCTTCTTTTTGTAAAATCTACAAGTGGATATATGGAGCGCTTTGAGGCCTATTGTGGAAAAGGAAATACCTTCACATAAAAACTACACAGAAGCATTCTGAGAAACTTCTTTTTGACGTGTGCATTCATCTCACAGAGTTGAACATTTCATGTGATTGAGCAGCTTTGAAACACTCTTTTTGTAAAATCTGCAAGTGGGTATTTGCAGCACTTTGAGGCCTATTTTGGAAAAGGAAATATCTTCCCATAAAAACTACATAGAAACATTCTCAGAAACTTCTTTGTGTTATCTGCATGTATGTAGAGAGTTGAACCTTTCATTTGATTTAGCAGTTTGCAAACACTCTTTTAGTAGAATCTGCAAGTAGATATTTGAAGCCCTTGGGGCCTATTATGGAAAAGGAAATATCTTCACATAAAAACTATGCAAAAGCGTTCTGAGAAACTTCATTGTGATGTGTGCATTCACTTAACAGAGTTGAACCTTTCTTTGAATTAAGCAGTTTTGAAACACT\t*\tSA:Z:chr3,90319741,-,121M499S,0,11;chr20,29212084,+,409S71M140S,38,4;chrUn_KN707904v1_decoy,742,+,620M,60,8;\tMD:Z:69A12T1G22C43T21T24A18G3T6C5T8A1A3C0G5T11T18G13A6G3C5C0A5C0G0G6C12A13C4G6G3C11\tRG:Z:GATKSVContigAlignments\tNM:i:34\tAS:i:211\tXS:i:181", true);
        final AlignmentInterval asm018485_tig00004_2 = TestUtilsForAssemblyBasedSVDiscovery.fromSAMRecordString("asm018485:tig00004\t2048\tchr20\t29212084\t38\t409H71M140H\t*\t0\t0\tAGAGTTGAACCTTTCATTTGATTTAGCAGTTTGCAAACACTCTTTTAGTAGAATCTGCAAGTAGATATTTG\t*SA:Z:chr20,28831147,+,275M2I114M229S,38,34;chr3,90319741,-,121M499S,0,11;chrUn_KN707904v1_decoy,742,+,620M,60,8;\tMD:Z:32T0G12T15G8\tRG:Z:GATKSVContigAlignments\tNM:i:4\tAS:i:51\tXS:i:33", true);
        final AlignmentInterval asm018485_tig00004_nonCanonical = TestUtilsForAssemblyBasedSVDiscovery.fromSAMRecordString("asm018485:tig00004\t2048\tchrUn_KN707904v1_decoy\t742\t60\t620M\t*\t0\t0\tTATCTTCACATAAAAACTACACAGTATCATTCTGAGAAACTTGTTTGTGATGTGTGCATTCATCTCACAGATTTGAACCCTTCCATCTTTTGAGCAGTTTGTACACCTTCTTTTTGTAAAATCTACAAGTGGATATATGGAGCGCTTTGAGGCCTATTGTGGAAAAGGAAATACCTTCACATAAAAACTACACAGAAGCATTCTGAGAAACTTCTTTTTGACGTGTGCATTCATCTCACAGAGTTGAACATTTCATGTGATTGAGCAGCTTTGAAACACTCTTTTTGTAAAATCTGCAAGTGGGTATTTGCAGCACTTTGAGGCCTATTTTGGAAAAGGAAATATCTTCCCATAAAAACTACATAGAAACATTCTCAGAAACTTCTTTGTGTTATCTGCATGTATGTAGAGAGTTGAACCTTTCATTTGATTTAGCAGTTTGCAAACACTCTTTTAGTAGAATCTGCAAGTAGATATTTGAAGCCCTTGGGGCCTATTATGGAAAAGGAAATATCTTCACATAAAAACTATGCAAAAGCGTTCTGAGAAACTTCATTGTGATGTGTGCATTCACTTAACAGAGTTGAACCTTTCTTTGAATTAAGCAGTTTTGAAACACT\t*\tSA:Z:chr20,28831147,+,275M2I114M229S,38,34;chr3,90319741,-,121M499S,0,11;chr20,29212084,+,409S71M140S,38,4;\tMD:Z:107C34T55T22T81A25G52C12G224\tRG:Z:GATKSVContigAlignments\tNM:i:8\tAS:i:580\tXS:i:293", true);
        final AlignmentInterval asm018485_tig00004_insmapping = TestUtilsForAssemblyBasedSVDiscovery.fromSAMRecordString("asm018485:tig00004\t2064\tchr3\t90319741\t0\t121M499H\t*\t0\t0\tAGTGTTTCAAAACTGCTTAATTCAAAGAAAGGTTCAACTCTGTTAAGTGAATGCACACATCACAATGAAGTTTCTCAGAACGCTTTTGCATAGTTTTTATGTGAAGATATTTCCTTTTCCA\t*\tSA:Z:chr20,28831147,+,275M2I114M229S,38,34;chr20,29212084,+,409S71M140S,38,4;chrUn_KN707904v1_decoy,742,+,620M,60,8;\tMD:Z:5A11C3C0A23T8G24T4C2T0C17C13\tRG:Z:GATKSVContigAlignments\tNM:i:11\tAS:i:66\tXS:i:62", true);
        final AlignmentInterval asm028167_tig00007_1 = TestUtilsForAssemblyBasedSVDiscovery.fromSAMRecordString("asm028167:tig00007\t0\tchr20\t28831011\t60\t72M2I339M2I114M104S\t*\t0\t0\tGAGAAACTTCTTTGTGATGTGTGCATTCATCTCACAGAGATGAACCTATCTTTTCATAGAGCAGTTTTGAAACTCTCTTTCTGTAGAATCTGCGACTGGATATTTGGAGCCCTTAGCGGCCTATGGTGGAAACGGAATTATCTTCACATAAAAACTACACAGTATCATTCTGAGAAACTTGTTTGTGATGTGTGCATTCATCTCACAGATTTGAACCCTTCCATCTTTTGAGCAGTTTGTACACCTTCTTTTTGTAAAATCTACAAGTGGATATATGGAGCGCTTTGAGGCCTATTGTGGAAAAGGAAATACCTTCACATAAAAACTACACAGAAGCATTCTGAGAAACTTCTTTTTGACGTGTGCATTCATCTCACAGAGTTGAACATTTCATGTGATTGAGCAGCTTTGAAACACTCTTTTTGTAAAATCTGCAAGTGGGTATTTGCAGCACTTTGAGGCCTATTTTGGAAAAGGAAATATCTTCCCATAAAAACTACATAGAAACATTCTCAGAAACTTCTTTGTGTTATCTGCATGTATGTAGAGAGTTGAACCTTTCATTTGATTTAGCAGTTTGCAAACACTCTTTTAGTAGAATCTGCAAGTAGATATTTGAAGCCCTTGGGGCCT\t*\tSA:Z:chr20,29212084,+,547S71M15S,21,4;chrUn_KN707904v1_decoy,604,+,633M,60,10;\tMD:Z:24G29G36A1A36A74A12T1G22C43T21T24A18G3T6C5T8A1A3C0G5T11T18G13A6G3C5C0A5C0G0G6C12A13C4G6G3C11\tRG:Z:GATKSVContigAlignments\tNM:i:41\tAS:i:304\tXS:i:179", true);
        final AlignmentInterval asm028167_tig00007_2 = TestUtilsForAssemblyBasedSVDiscovery.fromSAMRecordString("asm028167:tig00007\t2048\tchr20\t29212084\t21\t547H71M15H\t*\t0\t0\tAGAGTTGAACCTTTCATTTGATTTAGCAGTTTGCAAACACTCTTTTAGTAGAATCTGCAAGTAGATATTTG\t*SA:Z:chr20,28831011,+,72M2I339M2I114M104S,60,41;chrUn_KN707904v1_decoy,604,+,633M,60,10;\tMD:Z:32T0G12T15G8\tRG:Z:GATKSVContigAlignments\tNM:i:4\tAS:i:51\tXS:i:41", true);
        final AlignmentInterval asm028167_tig00007_nonCanonical = TestUtilsForAssemblyBasedSVDiscovery.fromSAMRecordString("asm028167:tig00007\t2048\tchrUn_KN707904v1_decoy\t604\t60\t633M\t*\t0\t0\tGAGAAACTTCTTTGTGATGTGTGCATTCATCTCACAGAGATGAACCTATCTTTTCATAGAGCAGTTTTGAAACTCTCTTTCTGTAGAATCTGCGACTGGATATTTGGAGCCCTTAGCGGCCTATGGTGGAAACGGAATTATCTTCACATAAAAACTACACAGTATCATTCTGAGAAACTTGTTTGTGATGTGTGCATTCATCTCACAGATTTGAACCCTTCCATCTTTTGAGCAGTTTGTACACCTTCTTTTTGTAAAATCTACAAGTGGATATATGGAGCGCTTTGAGGCCTATTGTGGAAAAGGAAATACCTTCACATAAAAACTACACAGAAGCATTCTGAGAAACTTCTTTTTGACGTGTGCATTCATCTCACAGAGTTGAACATTTCATGTGATTGAGCAGCTTTGAAACACTCTTTTTGTAAAATCTGCAAGTGGGTATTTGCAGCACTTTGAGGCCTATTTTGGAAAAGGAAATATCTTCCCATAAAAACTACATAGAAACATTCTCAGAAACTTCTTTGTGTTATCTGCATGTATGTAGAGAGTTGAACCTTTCATTTGATTTAGCAGTTTGCAAACACTCTTTTAGTAGAATCTGCAAGTAGATATTTGAAGCCCTTGGGGCCT\t*\tSA:Z:chr20,28831011,+,72M2I339M2I114M104S,60,41;chr20,29212084,+,547S71M15S,21,4;\tMD:Z:108A23A112C34T55T22T81A25G52C12G99\tRG:Z:GATKSVContigAlignments\tNM:i:10\tAS:i:583\tXS:i:304", true);

        final SimpleChimera asm018485_tig00004_chimera = new SimpleChimera("asm018485:tig00004", asm018485_tig00004_1, asm018485_tig00004_2, StrandSwitch.NO_SWITCH,
                true, Collections.singletonList(asm018485_tig00004_insmapping.toPackedString()), asm018485_tig00004_nonCanonical.toSATagString());
        final SimpleChimera asm028167_tig00007_chimera = new SimpleChimera("asm028167:tig00007", asm028167_tig00007_1, asm028167_tig00007_2, StrandSwitch.NO_SWITCH,
                true, Collections.emptyList(), asm028167_tig00007_nonCanonical.toSATagString());
        List<SimpleChimera> simpleChimeras = Arrays.asList(asm018485_tig00004_chimera, asm028167_tig00007_chimera);
        SimpleInterval leftBreakpoint = new SimpleInterval("chr20:28831535-28831535");
        SimpleInterval rightBreakpoint = new SimpleInterval("chr20:29212083-29212083");
        BreakpointComplications complications = new BreakpointComplications.SimpleInsDelOrReplacementBreakpointComplications("", "TTATCTGCATGTATGTAG");
        NovelAdjacencyAndAltHaplotype novelAdjacencyAndAltHaplotype = new NovelAdjacencyAndAltHaplotype(leftBreakpoint, rightBreakpoint, StrandSwitch.NO_SWITCH, complications, TypeInferredFromSimpleChimera.RPL, "TTATCTGCATGTATGTAG".getBytes());
        SimpleNovelAdjacencyAndChimericAlignmentEvidence simpleNovelAdjacencyAndChimericAlignmentEvidence = new SimpleNovelAdjacencyAndChimericAlignmentEvidence(novelAdjacencyAndAltHaplotype, simpleChimeras);
        Allele refAllele = Allele.create("G", true);
        final SimpleSVType.Deletion del_chr20_28831535_29212083 = new SimpleSVType.Deletion("chr20", 28831535, 29212083, "DEL_chr20_28831535_29212083", refAllele, Allele.create("<DEL>"), -380548, Collections.emptyMap());
        VariantContext expected = new VariantContextBuilder().chr("chr20").start(28831535).stop(29212083).id("DEL_chr20_28831535_29212083")
                .alleles(Arrays.asList(refAllele, Allele.create("<DEL>")))
                .attribute(VCFConstants.END_KEY, 29212083).attribute(SVLEN, -380548).attribute(SVTYPE, "DEL").attribute(CONTIG_NAMES, "asm018485:tig00004,asm028167:tig00007")
                .attribute(TOTAL_MAPPINGS, 2).attribute(HQ_MAPPINGS, 0).attribute(MAPPING_QUALITIES, "38,21").attribute(ALIGN_LENGTHS, "71,71").attribute(MAX_ALIGN_LENGTH, 71)
                .attribute(SEQ_ALT_HAPLOTYPE, "TTATCTGCATGTATGTAG").attribute(INSERTED_SEQUENCE, "TTATCTGCATGTATGTAG").attribute(INSERTED_SEQUENCE_LENGTH, 18)
                .attribute(INSERTED_SEQUENCE_MAPPINGS, "500_620_chr3:90319741-90319861_-_499H121M_0_11_66_O").attribute(CTG_GOOD_NONCANONICAL_MAPPING, "chrUn_KN707904v1_decoy,742,+,620M,60,8,580,chrUn_KN707904v1_decoy,604,+,633M,60,10,583").make();
        VariantContext actual = AnnotatedVariantProducer.produceAnnotatedVcFromAssemblyEvidence(del_chr20_28831535_29212083, simpleNovelAdjacencyAndChimericAlignmentEvidence, referenceBroadcast, refSeqDictBroadcast, null, "testSample").make();
        VariantContextTestUtils.assertVariantContextsAreEqual(actual, expected, Collections.singletonList(HQ_MAPPINGS));

        // cnv call available
        referenceBroadcast = testSparkContext.broadcast(TestUtilsForAssemblyBasedSVDiscovery.b37_reference);
        refSeqDictBroadcast = testSparkContext.broadcast(TestUtilsForAssemblyBasedSVDiscovery.b37_seqDict);
        final SAMFileHeader samFileHeader = new SAMFileHeader(refSeqDictBroadcast.getValue());
        SAMReadGroupRecord test = new SAMReadGroupRecord("test");
        test.setSample("sample");
        samFileHeader.addReadGroup(test);
        final Broadcast<SVIntervalTree<VariantContext>> cnvCallsBroadcast =
                StructuralVariationDiscoveryPipelineSpark.broadcastCNVCalls(testSparkContext, samFileHeader, SVIntegrationTestDataProvider.EXTERNAL_CNV_CALLS);
        final AlignmentInterval asm000000_tig00006_1 = TestUtilsForAssemblyBasedSVDiscovery.fromSAMRecordString("asm000000:tig00006\t16\t21\t43349675\t60\t448M229S\t*\t0\t0\tACTAGGTGGGTTATAACTTTTATTTAAAACTTTCAGTTCCAGCTGATGGTTATACCATTGGGAGCCTCCATTTACTTAGAAATGAAACTGAAAACAGACAACTAAAGCATGTCCAGGACTCCTGGCTCCACACCATGCCAGGCGACATCACTCAAGTCTCCAAAGATCACCAAGTGTCCAGCTCAGCTCCTGCCCTCATCAGCAAGTTTTCCAAATGAAAGTTACGTTGAAAGCCACAGTTACCATACTGTAACCAGAATTCAGGCAGTGGCTGCTAGCAGAGTATGATGAACAAGAGCAGGTCTGGTATAAAGACAGTGACTTTGCATTCCAAAGCTTAGCTTAGGGGAAGAACAGGCTTCTGCCTTAAGGGTACCCCTTTGCTTTCGGGGCAGAAAGCAGGCACTTTCAAAAGGGGGCTTGGCATGAATGTCATGAAAGGGAGGAACACCACTGTGAACCCGCTGCCCTACACGGCAGTTCTAGGGCTGAACTCACCGAACAGTGTTAACAAAAAGAGGCCTTGCTGTCTTATCATTTTTATTTAACGCACGAACATTAAGCAGTGTCTCACCCTGGACATTTTACAAGAGATTAAGCTGGCTGGATGCCTTTGCAAAAACAGTGCCCTAAAAATGTGTCATGTTTGGCCAAGATGCTCATCCAAGAATGGAAAA\t*\tSA:Z:21,43353486,-,442S235M,60,0;\tMD:Z:387T60\tRG:Z:GATKSVContigAlignments\tNM:i:1\tAS:i:443\tXS:i:0", true);
        final AlignmentInterval asm000000_tig00006_2 = TestUtilsForAssemblyBasedSVDiscovery.fromSAMRecordString("asm000000:tig00006\t2064\t21\t43353486\t60\t442H235M\t*\t0\t0\tGAGGAACACCACTGTGAACCCGCTGCCCTACACGGCAGTTCTAGGGCTGAACTCACCGAACAGTGTTAACAAAAAGAGGCCTTGCTGTCTTATCATTTTTATTTAACGCACGAACATTAAGCAGTGTCTCACCCTGGACATTTTACAAGAGATTAAGCTGGCTGGATGCCTTTGCAAAAACAGTGCCCTAAAAATGTGTCATGTTTGGCCAAGATGCTCATCCAAGAATGGAAAA\t*\tSA:Z:21,43349675,-,448M229S,60,1;\tMD:Z:235\tRG:Z:GATKSVContigAlignments\tNM:i:0\tAS:i:235\tXS:i:0", true);
        final AlignmentInterval asm000001_tig00001_1 = TestUtilsForAssemblyBasedSVDiscovery.fromSAMRecordString("asm000001:tig00001\t16\t21\t43349641\t60\t25M1D456M417S\t*\t0\t0\tGTCTCCCTGGCTTCTGAGATGGGCCTTCCCCCGACTAGGTGGGTTATAACTTTTATTTAAAACTTTCAGTTCCAGCTGATGGTTATACCATTGGGAGCCTCCATTTACTTAGAAATGAAACTGAAAACAGACAACTAAAGCATGTCCAGGACTCCTGGCTCCACACCATGCCAGGCGACATCACTCAAGTCTCCAAAGATCACCAAGTGTCCAGCTCAGCTCCTGCCCTCATCAGCAAGTTTTCCAAATGAAAGTTACGTTGAAAGCCACAGTTACCATACTGTAACCAGAATTCAGGCAGTGGCTGCTAGCAGAGTATGATGAACAAGAGCAGGTCTGGTATAAAGACAGTGACTTTGCATTCCAAAGCTTAGCTTAGGGGAAGAACAGGCTTCTGCCTTAAGGGTACCCCTTTGCTTTCGGGGCAGAAAGCAGGCACTTTCAAAAGGGGGCTTGGCATGAATGTCATGAAAGGGAGGAACACCACTGTGAACCCGCTGCCCTACACGGCAGTTCTAGGGCTGAACTCACCGAACAGTGTTAACAAAAAGAGGCCTTGCTGTCTTATCATTTTTATTTAACGCACGAACATTAAGCAGTGTCTCACCCTGGACATTTTACAAGAGATTAAGCTGGCTGGATGCCTTTGCAAAAACAGTGCCCTAAAAATGTGTCATGTTTGGCCAAGATGCTCATCCAAGAATGGAAAAGGCCATGTACACAATCCAAGCACCCGAGGGTGTTCTACTCCCAACTGACCCTTCCCAGGAGCCCGGGCAGATCCCAACAGGACTTCCTCCTTGTGGGTATGCATAGGATCCAGGCTGGCAAGAGCGACCAGGCTCCTCCTCCCGCACTCACAGCCCCGTGAAAGGGGAGGGGAGGGGAGGGAACCCGT\t*\tSA:Z:21,43353486,-,475S423M,60,0;\tMD:Z:25^T6A388T60\tRG:Z:GATKSVContigAlignments\tNM:i:3\tAS:i:454\tXS:i:0", true);
        final AlignmentInterval asm000001_tig00001_2 = TestUtilsForAssemblyBasedSVDiscovery.fromSAMRecordString("asm000001:tig00001\t2064\t21\t43353486\t60\t475H423M\t*\t0\t0\tGAGGAACACCACTGTGAACCCGCTGCCCTACACGGCAGTTCTAGGGCTGAACTCACCGAACAGTGTTAACAAAAAGAGGCCTTGCTGTCTTATCATTTTTATTTAACGCACGAACATTAAGCAGTGTCTCACCCTGGACATTTTACAAGAGATTAAGCTGGCTGGATGCCTTTGCAAAAACAGTGCCCTAAAAATGTGTCATGTTTGGCCAAGATGCTCATCCAAGAATGGAAAAGGCCATGTACACAATCCAAGCACCCGAGGGTGTTCTACTCCCAACTGACCCTTCCCAGGAGCCCGGGCAGATCCCAACAGGACTTCCTCCTTGTGGGTATGCATAGGATCCAGGCTGGCAAGAGCGACCAGGCTCCTCCTCCCGCACTCACAGCCCCGTGAAAGGGGAGGGGAGGGGAGGGAACCCGT\t*SA:Z:21,43349641,-,25M1D456M417S,60,3;\tMD:Z:423\tRG:Z:GATKSVContigAlignments\tNM:i:0\tAS:i:423\tXS:i:0", true);
        final SimpleChimera asm000000_tig00006_chimera = new SimpleChimera("asm000000:tig00006", asm000000_tig00006_1, asm000000_tig00006_2, StrandSwitch.NO_SWITCH,
                false, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
        final SimpleChimera asm000001_tig00001_chimera = new SimpleChimera("asm000001:tig00001", asm000001_tig00001_1, asm000001_tig00001_2, StrandSwitch.NO_SWITCH,
                false, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
        simpleChimeras = Arrays.asList(asm000001_tig00001_chimera, asm000000_tig00006_chimera);
        leftBreakpoint = new SimpleInterval("21:43350116-43350116");
        rightBreakpoint = new SimpleInterval("21:43353485-43353485");
        complications = new BreakpointComplications.SimpleInsDelOrReplacementBreakpointComplications("GAGGAA", "");
        novelAdjacencyAndAltHaplotype = new NovelAdjacencyAndAltHaplotype(leftBreakpoint, rightBreakpoint, StrandSwitch.NO_SWITCH, complications, TypeInferredFromSimpleChimera.SIMPLE_DEL, new byte[]{});
        simpleNovelAdjacencyAndChimericAlignmentEvidence = new SimpleNovelAdjacencyAndChimericAlignmentEvidence(novelAdjacencyAndAltHaplotype, simpleChimeras);
        final SimpleSVType.Deletion del_21_43350116_43353485 = new SimpleSVType.Deletion("21", 43350116, 43353485, "DEL_21_43350116_43353485", refAllele, Allele.create("<DEL>"), -3369, Collections.emptyMap());
        expected = new VariantContextBuilder().chr("21").start(43350116).stop(43353485).id("DEL_21_43350116_43353485")
                .alleles(Arrays.asList(refAllele, Allele.create("<DEL>"))).attribute(VCFConstants.END_KEY, 43353485)
                .attribute(SVLEN, -3369).attribute(SVTYPE, "DEL").attribute(CONTIG_NAMES, "asm000000:tig00006,asm000001:tig00001").attribute(TOTAL_MAPPINGS, 2)
                .attribute(HQ_MAPPINGS, 2).attribute(MAPPING_QUALITIES, "60,60").attribute(ALIGN_LENGTHS, "229,417").attribute(MAX_ALIGN_LENGTH, 417).attribute(HOMOLOGY, "GAGGAA")
                .attribute(HOMOLOGY_LENGTH, 6).attribute(EXTERNAL_CNV_CALLS, "CNV_21_43350200_43353400:1:80").make();
        actual = AnnotatedVariantProducer.produceAnnotatedVcFromAssemblyEvidence(del_21_43350116_43353485, simpleNovelAdjacencyAndChimericAlignmentEvidence, referenceBroadcast, refSeqDictBroadcast, cnvCallsBroadcast, "sample").make();
        VariantContextTestUtils.assertVariantContextsAreEqual(actual, expected, Collections.singletonList(HQ_MAPPINGS));
    }

    // -----------------------------------------------------------------------------------------------
    // CI test
    // -----------------------------------------------------------------------------------------------
    @DataProvider(name = "CIIntervals")
    private Object[][] getCIIntervalTests() {
        return new Object[][] {
                new Object[] { 200, new SVInterval(1, 190, 225), "-10,25", null},
                new Object[] { 200, new SVInterval(1, 200, 225), "0,25", null},
                new Object[] { 200, new SVInterval(1, 201, 225), null, new IllegalStateException("Interval must contain point")}
        };
    }

    @Test(dataProvider = "CIIntervals")
    public void testProduceCIInterval(final int point, final SVInterval interval, final String expected, final Exception expectedException) {
        if (expectedException == null) {
            Assert.assertEquals(AnnotatedVariantProducer.produceCIInterval(point, interval), expected);
        } else {
            try {
                AnnotatedVariantProducer.produceCIInterval(point, interval);
                Assert.fail("did not throw expected exception " + expectedException);
            } catch (Throwable e) {
                Assert.assertEquals(e.getClass(), expectedException.getClass());
                Assert.assertEquals(e.getMessage(), expectedException.getMessage());
            }
        }
    }


    // -----------------------------------------------------------------------------------------------
    // EvidenceTargetLink-based annotations
    // -----------------------------------------------------------------------------------------------
    @DataProvider(name = "evidenceTargetLinksAndPreciseVariants")
    private Object[][] getEvidenceTargetLinksAndPreciseVariants() {

        final VariantContext unAnnotatedVC = new VariantContextBuilder()
                .id("TESTID")
                .chr("20").start(200).stop(300)
                .alleles("N", SimpleSVType.ImpreciseDeletion.createBracketedSymbAlleleString(SYMB_ALT_ALLELE_DEL))
                .attribute(VCFConstants.END_KEY, 300)
                .attribute(SVTYPE, SimpleSVType.SupportedType.DEL.toString())
                .make();

        final VariantContext annotatedVC = new VariantContextBuilder()
                .id("TESTID")
                .chr("20").start(200).stop(300)
                .alleles("N", SimpleSVType.ImpreciseDeletion.createBracketedSymbAlleleString(SYMB_ALT_ALLELE_DEL))
                .attribute(VCFConstants.END_KEY, 300)
                .attribute(SVTYPE, SimpleSVType.SupportedType.DEL.toString())
                .attribute(READ_PAIR_SUPPORT, 7)
                .attribute(SPLIT_READ_SUPPORT, 5)
                .make();

        List<Object[]> tests = new ArrayList<>();
        tests.add(new Object[] {
                Arrays.asList(
                        new EvidenceTargetLink(
                                new StrandedInterval(new SVInterval(0, 190, 210), true),
                                new StrandedInterval(new SVInterval(0, 310, 320), false),
                                5, 7, new HashSet<>(), new HashSet<>())),
                Arrays.asList( unAnnotatedVC ),
                Arrays.asList( annotatedVC ) }
        );
        tests.add(new Object[] {
                Arrays.asList(
                        new EvidenceTargetLink(
                                new StrandedInterval(new SVInterval(0, 190, 210), true),
                                new StrandedInterval(new SVInterval(0, 310, 320), true),
                                5, 7, new HashSet<>(), new HashSet<>())),
                Arrays.asList( unAnnotatedVC ),
                Arrays.asList( unAnnotatedVC ) }
        );
        tests.add(new Object[] {
                Arrays.asList(
                        new EvidenceTargetLink(
                                new StrandedInterval(new SVInterval(0, 190, 210), true),
                                new StrandedInterval(new SVInterval(0, 310, 320), false),
                                3, 4, new HashSet<>(), new HashSet<>()),
                        new EvidenceTargetLink(
                                new StrandedInterval(new SVInterval(0, 192, 215), true),
                                new StrandedInterval(new SVInterval(0, 299, 303), false),
                                2, 3, new HashSet<>(), new HashSet<>())),
                Arrays.asList( unAnnotatedVC ),
                Arrays.asList( annotatedVC ) }
        );

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "evidenceTargetLinksAndPreciseVariants", groups = "sv")
    public void testProcessEvidenceTargetLinks(final List<EvidenceTargetLink> etls,
                                               final List<VariantContext> inputVariants,
                                               final List<VariantContext> expectedVariants) {

        final Logger localLogger = LogManager.getLogger(AnnotatedVariantProducer.class);
        final StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigAlignmentsSparkArgumentCollection params =
                new StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigAlignmentsSparkArgumentCollection();

        ReadMetadata metadata = Mockito.mock(ReadMetadata.class);
        when(metadata.getMaxMedianFragmentSize()).thenReturn(300);
        when(metadata.getContigName(0)).thenReturn("20");

        PairedStrandedIntervalTree<EvidenceTargetLink> evidenceTree = new PairedStrandedIntervalTree<>();
        etls.forEach(e -> evidenceTree.put(e.getPairedStrandedIntervals(), e));

        final List<VariantContext> processedVariantContexts =
                AnnotatedVariantProducer.annotateBreakpointBasedCallsWithImpreciseEvidenceLinks(inputVariants,
                        evidenceTree, metadata, TestUtilsForAssemblyBasedSVDiscovery.b37_reference, params, localLogger);

        VariantContextTestUtils.assertEqualVariants(processedVariantContexts, expectedVariants);
    }
}