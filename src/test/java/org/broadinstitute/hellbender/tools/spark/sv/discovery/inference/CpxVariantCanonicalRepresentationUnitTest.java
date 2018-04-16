package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SVTestUtils;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.VariantContextTestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.util.*;
import java.util.stream.Collectors;

import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.*;

public class CpxVariantCanonicalRepresentationUnitTest extends GATKBaseTest {

    @DataProvider(name = "forSerialization")
    private Object[][] forSerialization() {
        final List<Object[]> data = new ArrayList<>(20);
        for (final CpxSVInferenceTestUtils.PreprocessedAndAnalysisReadyContigWithExpectedResults x : CpxSVInferenceTestUtils.PREPROCESSED_AND_ANALYSIS_READY_CONTIGS_AND_EXPECTED_RESULTS) {
            data.add(new Object[]{x.expectedCpxVariantCanonicalRepresentation});
        }
        return data.toArray(new Object[data.size()][]);
    }

    @Test(groups = "sv", dataProvider = "forSerialization")
    public void testSerialization(final CpxVariantCanonicalRepresentation cpxVariantCanonicalRepresentation) {
        final ByteArrayOutputStream bos = new ByteArrayOutputStream();
        final Output out = new Output(bos);
        final Kryo kryo = new Kryo();
        kryo.writeClassAndObject(out, cpxVariantCanonicalRepresentation);
        out.flush();

        final ByteArrayInputStream bis = new ByteArrayInputStream(bos.toByteArray());
        final Input in = new Input(bis);
        @SuppressWarnings("unchecked")
        final CpxVariantCanonicalRepresentation roundTrip = (CpxVariantCanonicalRepresentation) kryo.readClassAndObject(in);
        Assert.assertEquals(roundTrip, cpxVariantCanonicalRepresentation);
    }

    // =================================================================================================================

    @Test
    public void testSpecialCtor() {
        final AlignedContig alignedContig = SVTestUtils.fromPrimarySAMRecordString("asm000308:tig00000\t0\tchr1\t14491391\t60\t1276M530S\t*\t0\t0\tTTGCTGCACCCATCAATCCGTCGTCTACATTAGGTATTTCTCCTAATGCTATCCCTCCCCTAGTCCCCTACCCGCCGACAGGTCCCGGTGTGTGATATTCCCCTCCCTGTGTCCATGTTACTCTTTTGATATTACCAGGGACACCTGGATTTCTACTGATTTTAATGAGAATACCTTCTGTATTCACCATTAAATATGATGGTAGTTGCTGGTTTTAGCTCGATATTATTTGTCATGCTGATTAAGTGGACTTGCATTCCTAGCTTTCAAAGAGGTTTTCTTTCCTTTTTAATAAGGAGTGGGTGTTGCATGTTATCTAATACATTGTCAGTGTGTGGCTATTTATAAGTTCTGTGTGATTATACCATTTATCTATATAAATGTCCCCTTTGTTCTATTTAATAATGCTTTTCCCCCTTCTGAATTCCACTTTGAATTTGAATTCTACTTTGTCTGAAATAGATCCTGCCACCCCTGCTTTTAAAAAGAAAAAAATCTTTTTGCTTGTATTATTTAACTTTTTTGCCTATCCCTCCCTTTTTAATCTTTTCATACCATTGCTTTTCAGTGTCTCGAGCAGTAAGACATTTAACAATTATCAGCCCCATGCTTACTTTGTGCCAGACACTGGATTAAACAAAAATGGAAAAAGAGGATAGAATGTGCTGGAAGGGGTACATTCAAACCCAGTCTGAACTGGCCACTGCTGTGAGCAGGTTTGGGGACAGCAGTAGATCCTAGAAGGGCCTGACCAGCTGGGGAAACTGGCCAGGCTGTCCAGAGGTGACAAGAGGATTGTCACCCAGACTTGCCCAAGAAGAGTGAATCTGAGTCTTGGAGAGAACAGGAGTTTGGGTTCTTCTGGGCCCAGATGGCCTCAGGGCTCCCTGGAATTTGGGGACCCCACAGTTGGTCGCCACCATGAATTGAGGAGCCTTGCTTCTCTCCACACTGTCTTTTCCCTGCCTCCTCGTGGCTTCTGCTTCACTCATTCACTCATTCTGTCAGTGAATGATTCTTCAGCACCTGCCCTGCATAGGATGCCATTGTAGGTGCTGGGAAATCAACGGGAAGAAGATGGAAAACGAGACTTCCCTTATGAAGCTTCTGTTCTACAGAGGTGGGCAGACATGGCCAGAAAAAGCACAAGGCCATTTCCAATGGTGGAAGGGCCAGGACTGCTGCCCTTTCTGATAGCTTCTCTTTACACTTAGGAGAAAATTCAGGGCCCCATAATCCCTAGGCCCTACATAACCACACATGCACACACCACAAACCACACACACACACCACACACACCACACACCACACGCTACACACACCATGCACATACCACAAACCACACACACACCACACACACACCACACATTACACACACCACACAGACACCACACACCACACAAACATCACACACCACAGACACATCACACACCACACACACACCACACATACACAGCACACAGCTCACGCATACACAGCACACACATCACACATACACATACCACATACACACCACACACACACCACAAACCACATATACACAGCACACACATCACACAAACACATACCACATACACACCACACACCGTACATACACAGCATACACATTACACATACACACACGACACACCACACACAGACCACACACCACACAGATACAGCACAGAGACACTACACATACCACATACACACTACACACCACACACACACCACACATACACAGCGCACATACACACACCACACACACCACATACAAACCACATACCACATACACACCACACATATACCACACAGACACCATACATA\t*\tSA:Z:chr1,14492666,+,1684S71M51S,60,0,71;chr4,8687087,-,422S46M1338S,60,1,41;\tMD:Z:144T1131\tRG:Z:GATKSVContigAlignments\tNM:i:1\tAS:i:1271\tXS:i:70", true);
        final AssemblyContigWithFineTunedAlignments preprocessedTig = new AssemblyContigWithFineTunedAlignments(alignedContig);
        final CpxVariantInducingAssemblyContig analysisReadyContig = new CpxVariantInducingAssemblyContig(preprocessedTig, CpxSVInferenceTestUtils.bareBoneHg38SAMSeqDict);
        final SimpleInterval manuallyCalculatedAffectedRefRegion = new SimpleInterval("chr1", 14492666, 14492666);
        final byte[] manuallyCalculatedAltSeq = Arrays.copyOfRange(alignedContig.getContigSequence(), 1275, 1685);
        final List<String> manuallyCalculatedAltArrangements = Arrays.asList("1", "UINS-62", "-chr4:8687087-8687132", "UINS-300", "1");
        final List<SimpleInterval> manuallyCalculatedSegments = Collections.singletonList(manuallyCalculatedAffectedRefRegion);

        final CpxVariantCanonicalRepresentation cpxVariantCanonicalRepresentation = new CpxVariantCanonicalRepresentation(analysisReadyContig);

        Assert.assertEquals(cpxVariantCanonicalRepresentation.getAffectedRefRegion(),
                            manuallyCalculatedAffectedRefRegion);
        Assert.assertEquals(cpxVariantCanonicalRepresentation.getReferenceSegments(),
                            manuallyCalculatedSegments);
        Assert.assertEquals(cpxVariantCanonicalRepresentation.getEventDescriptions(),
                            manuallyCalculatedAltArrangements);
        Assert.assertTrue(Arrays.equals(cpxVariantCanonicalRepresentation.getAltSeq(),
                                        manuallyCalculatedAltSeq));

        final byte[] dummyRefSequence = "TCGA".getBytes();
        final VariantContextBuilder variantContextBuilder = new VariantContextBuilder()
                .chr("chr1").start(14492666).stop(14492666)
                .alleles(Arrays.asList(Allele.create(dummyRefSequence, true), Allele.create(SimpleSVType.createBracketedSymbAlleleString(CPX_SV_SYB_ALT_ALLELE_STR), false)))
                .id("CPX_chr1:14492666-14492666")
                .attribute(VCFConstants.END_KEY, 14492666)
                .attribute(SVTYPE, GATKSVVCFConstants.CPX_SV_SYB_ALT_ALLELE_STR)
                .attribute(SVLEN, 409)
                .attribute(SEQ_ALT_HAPLOTYPE, new String(manuallyCalculatedAltSeq));
        variantContextBuilder.attribute(CPX_EVENT_ALT_ARRANGEMENTS,
                String.join(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR, manuallyCalculatedAltArrangements));
        variantContextBuilder.attribute(CPX_SV_REF_SEGMENTS,
                String.join(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR, manuallyCalculatedSegments.stream().map(SimpleInterval::toString).collect(Collectors.toList())));
        final VariantContext manuallyCalculatedVariantContext = variantContextBuilder.make();

        VariantContextTestUtils.assertVariantContextsAreEqual(
                cpxVariantCanonicalRepresentation.toVariantContext(dummyRefSequence).make(),
                manuallyCalculatedVariantContext,
                Collections.emptyList());
    }

    @DataProvider(name = "forGeneralCtor")
    private Object[][] forGeneralCtor() {
        final List<Object[]> data = new ArrayList<>(20);
        for (final CpxSVInferenceTestUtils.PreprocessedAndAnalysisReadyContigWithExpectedResults x : CpxSVInferenceTestUtils.PREPROCESSED_AND_ANALYSIS_READY_CONTIGS_AND_EXPECTED_RESULTS) {
            data.add(new Object[]{x.expectedCpxVariantInducingAssemblyContig, x.expectedCpxVariantCanonicalRepresentation});
        }
        return data.toArray(new Object[data.size()][]);
    }

    @Test(groups = "sv", dataProvider = "forGeneralCtor")
    public void testGeneralCtor(final CpxVariantInducingAssemblyContig cpxVariantInducingAssemblyContig,
                                final CpxVariantCanonicalRepresentation expectedResult) {
        Assert.assertEquals(new CpxVariantCanonicalRepresentation(cpxVariantInducingAssemblyContig), expectedResult);
    }

    @DataProvider(name = "forExtractRefSegments")
    private Object[][] forExtractRefSegments() {
        final List<Object[]> data = new ArrayList<>(20);
        for (final CpxSVInferenceTestUtils.PreprocessedAndAnalysisReadyContigWithExpectedResults x : CpxSVInferenceTestUtils.PREPROCESSED_AND_ANALYSIS_READY_CONTIGS_AND_EXPECTED_RESULTS) {
            data.add(new Object[]{x.expectedCpxVariantInducingAssemblyContig.getBasicInfo(),
                                  x.expectedCpxVariantInducingAssemblyContig.getEventPrimaryChromosomeSegmentingLocations(),
                                  x.expectedCpxVariantCanonicalRepresentation.getReferenceSegments(),
                                  x.expectedCpxVariantInducingAssemblyContig.getTwoBaseBoundaries()
            });
        }
        return data.toArray(new Object[data.size()][]);
    }

    @Test(groups = "sv", dataProvider = "forExtractRefSegments")
    public void testExtractRefSegments(final CpxVariantInducingAssemblyContig.BasicInfo basicInfo,
                                       final List<SimpleInterval> segmentingLocations,
                                       final List<SimpleInterval> expectedResult,
                                       final Set<SimpleInterval> twoBaseBoundaries) {
        Assert.assertEquals(CpxVariantCanonicalRepresentation.extractRefSegments(basicInfo, segmentingLocations, twoBaseBoundaries), expectedResult);
    }

    @DataProvider(name = "forExtractAltArrangements")
    private Object[][] forExtractAltArrangements() {
        final List<Object[]> data = new ArrayList<>(20);
        for (final CpxSVInferenceTestUtils.PreprocessedAndAnalysisReadyContigWithExpectedResults x : CpxSVInferenceTestUtils.PREPROCESSED_AND_ANALYSIS_READY_CONTIGS_AND_EXPECTED_RESULTS) {
            final CpxVariantInducingAssemblyContig cpxVariantInducingAssemblyContig = x.expectedCpxVariantInducingAssemblyContig;
            final CpxVariantInducingAssemblyContig.BasicInfo basicInfo = cpxVariantInducingAssemblyContig.getBasicInfo();
            final List<SimpleInterval> eventPrimaryChromosomeSegmentingLocations = cpxVariantInducingAssemblyContig.getEventPrimaryChromosomeSegmentingLocations();
            final List<SimpleInterval> segments = CpxVariantCanonicalRepresentation.extractRefSegments(basicInfo, eventPrimaryChromosomeSegmentingLocations, cpxVariantInducingAssemblyContig.getTwoBaseBoundaries());
            data.add(new Object[]{basicInfo,
                                  cpxVariantInducingAssemblyContig.getPreprocessedTig().getAlignments(),
                                  cpxVariantInducingAssemblyContig.getJumps(),
                                  segments,
                                  x.expectedCpxVariantCanonicalRepresentation.getEventDescriptions()
            });
        }
        return data.toArray(new Object[data.size()][]);
    }

    @Test(groups = "sv", dataProvider = "forExtractAltArrangements")
    public void testExtractAltArrangements(final CpxVariantInducingAssemblyContig.BasicInfo basicInfo,
                                       final List<AlignmentInterval> contigAlignments,
                                       final List<CpxVariantInducingAssemblyContig.Jump> jumps,
                                       final List<SimpleInterval> segments,
                                       final List<String> expectedResult) {
        Assert.assertEquals(CpxVariantCanonicalRepresentation.extractAltArrangements(basicInfo, contigAlignments, jumps, segments),
                            expectedResult);
    }

    @DataProvider(name = "forExtractAltHaplotypeSeq")
    private Object[][] forExtractAltHaplotypeSeq() {
        final List<Object[]> data = new ArrayList<>(20);
        for (final CpxSVInferenceTestUtils.PreprocessedAndAnalysisReadyContigWithExpectedResults x : CpxSVInferenceTestUtils.PREPROCESSED_AND_ANALYSIS_READY_CONTIGS_AND_EXPECTED_RESULTS) {
            final CpxVariantInducingAssemblyContig cpxVariantInducingAssemblyContig = x.expectedCpxVariantInducingAssemblyContig;
            final CpxVariantInducingAssemblyContig.BasicInfo basicInfo = cpxVariantInducingAssemblyContig.getBasicInfo();
            final List<SimpleInterval> eventPrimaryChromosomeSegmentingLocations = cpxVariantInducingAssemblyContig.getEventPrimaryChromosomeSegmentingLocations();
            final List<SimpleInterval> segments = CpxVariantCanonicalRepresentation.extractRefSegments(basicInfo, eventPrimaryChromosomeSegmentingLocations, cpxVariantInducingAssemblyContig.getTwoBaseBoundaries());

            data.add(new Object[]{cpxVariantInducingAssemblyContig.getPreprocessedTig(),
                                  segments,
                                  basicInfo,
                                  x.expectedCpxVariantCanonicalRepresentation.getAltSeq()
            });
        }
        return data.toArray(new Object[data.size()][]);
    }

    @Test(groups = "sv", dataProvider = "forExtractAltHaplotypeSeq")
    public void testExtractAltHaplotypeSeq(final AssemblyContigWithFineTunedAlignments tigWithInsMappings,
                                           final List<SimpleInterval> segments,
                                           final CpxVariantInducingAssemblyContig.BasicInfo basicInfo,
                                           final byte[] expectedResult) {
        Assert.assertTrue(
                Arrays.equals(CpxVariantCanonicalRepresentation.extractAltHaplotypeSeq(tigWithInsMappings, segments, basicInfo),
                              expectedResult)
        );
    }

    @DataProvider(name = "forGetAffectedReferenceRegion")
    private Object[][] forGetAffectedReferenceRegion() {
        final List<Object[]> data = new ArrayList<>(20);
        for (final CpxSVInferenceTestUtils.PreprocessedAndAnalysisReadyContigWithExpectedResults x : CpxSVInferenceTestUtils.PREPROCESSED_AND_ANALYSIS_READY_CONTIGS_AND_EXPECTED_RESULTS) {

            data.add(new Object[]{x.expectedCpxVariantInducingAssemblyContig.getEventPrimaryChromosomeSegmentingLocations(),
                                  x.expectedCpxVariantCanonicalRepresentation.getAffectedRefRegion()
            });
        }
        return data.toArray(new Object[data.size()][]);
    }

    @Test(groups = "sv", dataProvider = "forGetAffectedReferenceRegion")
    public void testGetAffectedReferenceRegion(final List<SimpleInterval> eventPrimaryChromosomeSegmentingLocations,
                                               final SimpleInterval expectedResult) {
        Assert.assertEquals(CpxVariantCanonicalRepresentation.getAffectedReferenceRegion(eventPrimaryChromosomeSegmentingLocations),
                            expectedResult);
    }

    @DataProvider(name = "forToVariantContext")
    private Object[][] forToVariantContext() {
        final List<Object[]> data = new ArrayList<>(20);
        for (final CpxSVInferenceTestUtils.PreprocessedAndAnalysisReadyContigWithExpectedResults x : CpxSVInferenceTestUtils.PREPROCESSED_AND_ANALYSIS_READY_CONTIGS_AND_EXPECTED_RESULTS) {

            data.add(new Object[]{x.expectedCpxVariantCanonicalRepresentation,
                                  x.assumedReferenceSequence,
                                  x.expectedVariantContext
            });
        }
        return data.toArray(new Object[data.size()][]);
    }

    @Test(groups = "sv", dataProvider = "forToVariantContext")
    public void testToVariantContext(final CpxVariantCanonicalRepresentation tobetested, final byte[] refBases,
                                     final VariantContext expectedResult) {
        VariantContextTestUtils.assertVariantContextsAreEqual(tobetested.toVariantContext(refBases).make(),
                                                              expectedResult,
                                                              Collections.emptyList());
    }
}
