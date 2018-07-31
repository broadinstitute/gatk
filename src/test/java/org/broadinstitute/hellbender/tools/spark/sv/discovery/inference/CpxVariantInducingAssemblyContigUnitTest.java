package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import org.apache.commons.collections4.CollectionUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.TestUtilsForAssemblyBasedSVDiscovery;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.util.*;

import static org.broadinstitute.hellbender.tools.spark.sv.discovery.TestUtilsForAssemblyBasedSVDiscovery.bareBoneHg38SAMSeqDict;

public class CpxVariantInducingAssemblyContigUnitTest extends GATKBaseTest {
    @DataProvider(name = "forBasicInfoCtor")
    private Object[][] forBasicInfoCtor() {
        final List<Object[]> data = new ArrayList<>(20);

        data.add(new Object[]{new AlignedContig("unmapped", "AAAAA".getBytes(), Collections.emptyList()),
                              null, GATKException.class
        });

        for (final TestUtilsCpxVariantInference.PreprocessedAndAnalysisReadyContigWithExpectedResults x : TestUtilsCpxVariantInference.PREPROCESSED_AND_ANALYSIS_READY_CONTIGS_AND_EXPECTED_RESULTS) {

            data.add(new Object[]{x.expectedCpxVariantInducingAssemblyContig.getPreprocessedTig().getSourceContig(),
                                  x.expectedCpxVariantInducingAssemblyContig.getBasicInfo(),
                                  null
            });

        }

        return data.toArray(new Object[data.size()][]);
    }

    @Test(groups = "sv", dataProvider = "forBasicInfoCtor")
    @SuppressWarnings("rawtypes")
    public void testBasicInfoCtor(final AlignedContig contig,
                                  final CpxVariantInducingAssemblyContig.BasicInfo expectedResults,
                                  final Class expectedExceptionClass) {
        try {
            final CpxVariantInducingAssemblyContig.BasicInfo basicInfo = new CpxVariantInducingAssemblyContig.BasicInfo(contig);
            Assert.assertEquals(basicInfo, expectedResults);
        } catch (final Exception e) {
            Assert.assertEquals(e.getClass(), expectedExceptionClass);
        }
    }

    @DataProvider(name = "forJumpCtor")
    private Object[][] forJumpCtor() {
        final List<Object[]> data = new ArrayList<>(20);

        data.add(new Object[]{TestUtilsForAssemblyBasedSVDiscovery.fromSAMRecordString("asm025517:tig00006\t2048\tchr17\t82596440\t60\t109M142H\t*\t0\t0\tAGCCTCAGAGGTGGCGTCAGGACGTGCCTGCCCCACCGGTTGCCCTGGTGCCCTCGTCACGCCCCGGGACCGTGCACACGTGGGGACTGTTTCCAGACGCACTTTCTGC\t*\tSA:Z:chr17,82596480,+,106S145M,60,0;\tMD:Z:72C36\tRG:Z:GATKSVContigAlignments\tNM:i:1\tAS:i:104\tXS:i:19", true),
                              TestUtilsForAssemblyBasedSVDiscovery.fromSAMRecordString("asm025517:tig00006\t0\tchr17\t82596480\t60\t106S145M\t*\t0\t0\tAGCCTCAGAGGTGGCGTCAGGACGTGCCTGCCCCACCGGTTGCCCTGGTGCCCTCGTCACGCCCCGGGACCGTGCACACGTGGGGACTGTTTCCAGACGCACTTTCTGCCCTGGTGCCCTCGTCACGCCCCGGGACCGCGCACACGTGGGGACTGTTTCCAGACGCACTTTCTGCGGGCAGTCTGTGTGGCAGGGCTCCCTGCCCAGCTCCTGCAGCCTCATCAAGTCTCCCACTAAGGAGGTGTCGCTCC\t*\tSA:Z:chr17,82596440,+,109M142S,60,1;\tMD:Z:145\tRG:Z:GATKSVContigAlignments\tNM:i:0\tAS:i:145\tXS:i:20", true),
                              null,
                              IllegalArgumentException.class
        });

        for (final TestUtilsCpxVariantInference.PreprocessedAndAnalysisReadyContigWithExpectedResults x : TestUtilsCpxVariantInference.PREPROCESSED_AND_ANALYSIS_READY_CONTIGS_AND_EXPECTED_RESULTS) {
            final CpxVariantInducingAssemblyContig expectedCpxVariantInducingAssemblyContig = x.expectedCpxVariantInducingAssemblyContig;
            final List<AlignmentInterval> selectedAlignments = expectedCpxVariantInducingAssemblyContig.getPreprocessedTig().getAlignments();
            for (int i = 0; i < selectedAlignments.size() - 1; ++i) {
                data.add(new Object[]{selectedAlignments.get(i), selectedAlignments.get(i + 1),
                                      expectedCpxVariantInducingAssemblyContig.getJumps().get(i),
                                      null});
            }
        }

        return data.toArray(new Object[data.size()][]);
    }

    @Test(groups = "sv", dataProvider = "forJumpCtor")
    @SuppressWarnings("rawtypes")
    public void testJumpCtor(final AlignmentInterval one, final AlignmentInterval two,
                             final CpxVariantInducingAssemblyContig.Jump expectedResults,
                             final Class expectedExceptionClass) {
        try {
            final CpxVariantInducingAssemblyContig.Jump jump = new CpxVariantInducingAssemblyContig.Jump(one, two);
            Assert.assertEquals(jump, expectedResults);
        } catch (final Exception e) {
            Assert.assertEquals(e.getClass(), expectedExceptionClass);
        }
    }

    @DataProvider(name = "forSerialization")
    private Object[][] forSerialization() {
        final List<Object[]> data = new ArrayList<>(20);
        for (final TestUtilsCpxVariantInference.PreprocessedAndAnalysisReadyContigWithExpectedResults x : TestUtilsCpxVariantInference.PREPROCESSED_AND_ANALYSIS_READY_CONTIGS_AND_EXPECTED_RESULTS) {
            data.add(new Object[]{x.expectedCpxVariantInducingAssemblyContig});
        }
        return data.toArray(new Object[data.size()][]);
    }

    @Test(groups = "sv", dataProvider = "forSerialization")
    public void testSerialization(final CpxVariantInducingAssemblyContig cpxVariantInducingAssemblyContig) {
        final ByteArrayOutputStream bos = new ByteArrayOutputStream();
        final Output out = new Output(bos);
        final Kryo kryo = new Kryo();
        kryo.writeClassAndObject(out, cpxVariantInducingAssemblyContig);
        out.flush();

        final ByteArrayInputStream bis = new ByteArrayInputStream(bos.toByteArray());
        final Input in = new Input(bis);
        @SuppressWarnings("unchecked")
        final CpxVariantInducingAssemblyContig roundTrip = (CpxVariantInducingAssemblyContig) kryo.readClassAndObject(in);
        Assert.assertEquals(cpxVariantInducingAssemblyContig, roundTrip);
        Assert.assertEquals(cpxVariantInducingAssemblyContig.hashCode(), roundTrip.hashCode());
    }

    // =================================================================================================================

    @DataProvider(name = "forExtractJumpsOnReference")
    private Object[][] forExtractJumpsOnReference() {
        final List<Object[]> data = new ArrayList<>(20);

        for (final TestUtilsCpxVariantInference.PreprocessedAndAnalysisReadyContigWithExpectedResults x : TestUtilsCpxVariantInference.PREPROCESSED_AND_ANALYSIS_READY_CONTIGS_AND_EXPECTED_RESULTS) {
            final CpxVariantInducingAssemblyContig expectedCpxVariantInducingAssemblyContig = x.expectedCpxVariantInducingAssemblyContig;
            final List<AlignmentInterval> selectedAlignments = expectedCpxVariantInducingAssemblyContig.getPreprocessedTig().getAlignments();
            data.add(new Object[]{selectedAlignments, expectedCpxVariantInducingAssemblyContig.getJumps()});
        }

        return data.toArray(new Object[data.size()][]);
    }

    @Test(groups = "sv", dataProvider = "forExtractJumpsOnReference")
    public void testExtractJumpsOnReference(final List<AlignmentInterval> alignments,
                                            final List<CpxVariantInducingAssemblyContig.Jump> expectedResults) {
        final List<CpxVariantInducingAssemblyContig.Jump> result = CpxVariantInducingAssemblyContig.extractJumpsOnReference(alignments);
        Assert.assertEquals(result, expectedResults);
    }

    @DataProvider(name = "forExtractSegmentingRefLocationsOnEventPrimaryChromosome")
    private Object[][] forExtractSegmentingRefLocationsOnEventPrimaryChromosome() {
        final List<Object[]> data = new ArrayList<>(20);

        for (final TestUtilsCpxVariantInference.PreprocessedAndAnalysisReadyContigWithExpectedResults x : TestUtilsCpxVariantInference.PREPROCESSED_AND_ANALYSIS_READY_CONTIGS_AND_EXPECTED_RESULTS) {
            final CpxVariantInducingAssemblyContig expectedCpxVariantInducingAssemblyContig = x.expectedCpxVariantInducingAssemblyContig;
            final List<CpxVariantInducingAssemblyContig.Jump> jumps = expectedCpxVariantInducingAssemblyContig.getJumps();
            final CpxVariantInducingAssemblyContig.BasicInfo basicInfo = expectedCpxVariantInducingAssemblyContig.getBasicInfo();
            data.add(new Object[]{jumps, basicInfo, expectedCpxVariantInducingAssemblyContig.getEventPrimaryChromosomeSegmentingLocations()});
        }

        return data.toArray(new Object[data.size()][]);
    }

    @Test(groups = "sv", dataProvider = "forExtractSegmentingRefLocationsOnEventPrimaryChromosome")
    public void testExtractSegmentingRefLocationsOnEventPrimaryChromosome(final List<CpxVariantInducingAssemblyContig.Jump> jumps,
                                                                          final CpxVariantInducingAssemblyContig.BasicInfo basicInfo,
                                                                          List<SimpleInterval> expectedResults) {
        final List<SimpleInterval> result =
                CpxVariantInducingAssemblyContig.extractSegmentingRefLocationsOnEventPrimaryChromosome(jumps, basicInfo,
                        bareBoneHg38SAMSeqDict);
        Assert.assertEquals(result, expectedResults);
    }

    // =================================================================================================================
    @DataProvider(name = "forGetAlignmentsBoundedBySegmentingLocations")
    private Object[][] forGetAlignmentsBoundedBySegmentingLocations() {
        final List<Object[]> data = new ArrayList<>(20);
        for (final TestUtilsCpxVariantInference.PreprocessedAndAnalysisReadyContigWithExpectedResults x : TestUtilsCpxVariantInference.PREPROCESSED_AND_ANALYSIS_READY_CONTIGS_AND_EXPECTED_RESULTS) {
            final CpxVariantInducingAssemblyContig expectedCpxVariantInducingAssemblyContig = x.expectedCpxVariantInducingAssemblyContig;
            final Set<SimpleInterval> expectedTwoBaseBoundaries = x.expectedTwoBaseBoundaries;
            data.add(new Object[]{expectedCpxVariantInducingAssemblyContig, expectedTwoBaseBoundaries});
        }
        return data.toArray(new Object[data.size()][]);
    }

    @Test(groups = "sv", dataProvider = "forGetAlignmentsBoundedBySegmentingLocations")
    public void testGetAlignmentsBoundedBySegmentingLocations(final CpxVariantInducingAssemblyContig cpxVariantInducingAssemblyContig,
                                                              final Set<SimpleInterval> expectedAlignmentsBoundedBySegmentingLocations) {
        Assert.assertTrue(CollectionUtils.isEqualCollection(cpxVariantInducingAssemblyContig.getTwoBaseBoundaries(),
                                                            expectedAlignmentsBoundedBySegmentingLocations));
    }
}
