package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.TestUtilsForAssemblyBasedSVDiscovery;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.ContigAlignmentsModifier;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.StrandSwitch;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import scala.Tuple2;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import static org.broadinstitute.hellbender.tools.spark.sv.discovery.TestUtilsForAssemblyBasedSVDiscovery.*;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments.NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.AssemblyBasedSVDiscoveryTestDataProviderForSimpleSV.TestDataForSimpleSV;

public class SimpleChimeraUnitTest extends AssemblyBasedSVDiscoveryBaseTest {

    @DataProvider
    private Object[][] forSplitPairStrongEnoughEvidenceForCA() {
        final List<Object[]> data = new ArrayList<>(20);

        String alignmentOneSAMString = "asm018495:tig00014\t16\tchr20\t29851171\t60\t72S582M\t*\t0\t0\tATCTGTAAGTGGATATTTGGAGCCCTTTGTGGCCTATGGTGGAAAAGGAAATATCTTCAAATAAAAAGTATGCAGAAGCATTCTGAGAAACTTTTTTGTGCTGTGTGCATTCATCTCACAGGGTTGAACCTATCTTATGATTGAGCAGTTTTGAATCACTCTTTATGTAGAATCTGCAAGTGGATATTTGGAGCGTTTTGAGACCTACCGTGGAAAAGCAATTATCTTCAGATAAAAACTACACAGAAGCATTCTGAGAAACTGCTTTATGATGTGTGCATTCATCTCACAGAGTTGAACCTTTCTTTTGATTGAGCAGCTTTGAAACACTCTTTTTGTAGAATCTGCAAGTGGATATTTCGTGTGCTTTGAGTCCTACCATGGAAAAGGAAATATCTTCACATAAAAAATACTCAGAGGAATTCTGAGAAACTTCTTTGTGATGTGTGCATTCAACTCACAGAATTGAACCTATCTTTAGATTGAGCAGTTTAGAATCTCTCTTTTTGCAGTATCTGCAAGTGGATATTTGGAGCCCTTTGCAGCCTGTGGTGGAAAAGAAAATATCTTCACACAAAAACTACTCGGAAGCATTCTTAGAAACTTCTTTTTGATGTGTGCATTCAAATCACAGAGTTGAACCTATATTTTCAT\t*\tSA:Z:chr20,28843254,-,141M513S,60,15;\tMD:Z:5T15C6A553\tRG:Z:GATKSVContigAlignments\tNM:i:3\tAS:i:567\tXS:i:282";
        String alignmentTwoSAMString = "asm018495:tig00014\t2064\tchr20\t28843254\t60\t141M513H\t*\t0\t0\tATCTGTAAGTGGATATTTGGAGCCCTTTGTGGCCTATGGTGGAAAAGGAAATATCTTCAAATAAAAAGTATGCAGAAGCATTCTGAGAAACTTTTTTGTGCTGTGTGCATTCATCTCACAGGGTTGAACCTATCTTATGAT\t*\tSA:Z:chr20,29851171,-,72S582M,60,3;\tMD:Z:5C7A0G52C2C0A2A2A11G2A0C6A4G14A13A6\tRG:Z:GATKSVContigAlignments\tNM:i:15\tAS:i:66\tXS:i:0";

        AlignmentInterval alignmentOne = fromSAMRecordString(alignmentOneSAMString, true);
        AlignmentInterval alignmentTwo = fromSAMRecordString(alignmentTwoSAMString, true);
        data.add(new Object[]{alignmentOne, alignmentTwo, SimpleNovelAdjacencyInterpreter.MORE_RELAXED_ALIGNMENT_MIN_MQ, SimpleNovelAdjacencyInterpreter.MORE_RELAXED_ALIGNMENT_MIN_LENGTH, true});

        alignmentOneSAMString = "asm018495:tig00014\t16\tchr20\t29851171\t19\t72S582M\t*\t0\t0\tATCTGTAAGTGGATATTTGGAGCCCTTTGTGGCCTATGGTGGAAAAGGAAATATCTTCAAATAAAAAGTATGCAGAAGCATTCTGAGAAACTTTTTTGTGCTGTGTGCATTCATCTCACAGGGTTGAACCTATCTTATGATTGAGCAGTTTTGAATCACTCTTTATGTAGAATCTGCAAGTGGATATTTGGAGCGTTTTGAGACCTACCGTGGAAAAGCAATTATCTTCAGATAAAAACTACACAGAAGCATTCTGAGAAACTGCTTTATGATGTGTGCATTCATCTCACAGAGTTGAACCTTTCTTTTGATTGAGCAGCTTTGAAACACTCTTTTTGTAGAATCTGCAAGTGGATATTTCGTGTGCTTTGAGTCCTACCATGGAAAAGGAAATATCTTCACATAAAAAATACTCAGAGGAATTCTGAGAAACTTCTTTGTGATGTGTGCATTCAACTCACAGAATTGAACCTATCTTTAGATTGAGCAGTTTAGAATCTCTCTTTTTGCAGTATCTGCAAGTGGATATTTGGAGCCCTTTGCAGCCTGTGGTGGAAAAGAAAATATCTTCACACAAAAACTACTCGGAAGCATTCTTAGAAACTTCTTTTTGATGTGTGCATTCAAATCACAGAGTTGAACCTATATTTTCAT\t*\tSA:Z:chr20,28843254,-,141M513S,60,15;\tMD:Z:5T15C6A553\tRG:Z:GATKSVContigAlignments\tNM:i:3\tAS:i:567\tXS:i:282";
        alignmentTwoSAMString = "asm018495:tig00014\t2064\tchr20\t28843254\t60\t141M513H\t*\t0\t0\tATCTGTAAGTGGATATTTGGAGCCCTTTGTGGCCTATGGTGGAAAAGGAAATATCTTCAAATAAAAAGTATGCAGAAGCATTCTGAGAAACTTTTTTGTGCTGTGTGCATTCATCTCACAGGGTTGAACCTATCTTATGAT\t*\tSA:Z:chr20,29851171,-,72S582M,60,3;\tMD:Z:5C7A0G52C2C0A2A2A11G2A0C6A4G14A13A6\tRG:Z:GATKSVContigAlignments\tNM:i:15\tAS:i:66\tXS:i:0";
        alignmentOne = fromSAMRecordString(alignmentOneSAMString, true);
        alignmentTwo = fromSAMRecordString(alignmentTwoSAMString, true);
        data.add(new Object[]{alignmentOne, alignmentTwo, SimpleNovelAdjacencyInterpreter.MORE_RELAXED_ALIGNMENT_MIN_MQ, SimpleNovelAdjacencyInterpreter.MORE_RELAXED_ALIGNMENT_MIN_LENGTH, false});

        alignmentOneSAMString = "asm018495:tig00014\t16\tchr20\t29851171\t60\t72S582M\t*\t0\t0\tATCTGTAAGTGGATATTTGGAGCCCTTTGTGGCCTATGGTGGAAAAGGAAATATCTTCAAATAAAAAGTATGCAGAAGCATTCTGAGAAACTTTTTTGTGCTGTGTGCATTCATCTCACAGGGTTGAACCTATCTTATGATTGAGCAGTTTTGAATCACTCTTTATGTAGAATCTGCAAGTGGATATTTGGAGCGTTTTGAGACCTACCGTGGAAAAGCAATTATCTTCAGATAAAAACTACACAGAAGCATTCTGAGAAACTGCTTTATGATGTGTGCATTCATCTCACAGAGTTGAACCTTTCTTTTGATTGAGCAGCTTTGAAACACTCTTTTTGTAGAATCTGCAAGTGGATATTTCGTGTGCTTTGAGTCCTACCATGGAAAAGGAAATATCTTCACATAAAAAATACTCAGAGGAATTCTGAGAAACTTCTTTGTGATGTGTGCATTCAACTCACAGAATTGAACCTATCTTTAGATTGAGCAGTTTAGAATCTCTCTTTTTGCAGTATCTGCAAGTGGATATTTGGAGCCCTTTGCAGCCTGTGGTGGAAAAGAAAATATCTTCACACAAAAACTACTCGGAAGCATTCTTAGAAACTTCTTTTTGATGTGTGCATTCAAATCACAGAGTTGAACCTATATTTTCAT\t*\tSA:Z:chr20,28843254,-,141M513S,60,15;\tMD:Z:5T15C6A553\tRG:Z:GATKSVContigAlignments\tNM:i:3\tAS:i:567\tXS:i:282";
        alignmentTwoSAMString = "asm018495:tig00014\t2064\tchr20\t28843254\t19\t141M513H\t*\t0\t0\tATCTGTAAGTGGATATTTGGAGCCCTTTGTGGCCTATGGTGGAAAAGGAAATATCTTCAAATAAAAAGTATGCAGAAGCATTCTGAGAAACTTTTTTGTGCTGTGTGCATTCATCTCACAGGGTTGAACCTATCTTATGAT\t*\tSA:Z:chr20,29851171,-,72S582M,60,3;\tMD:Z:5C7A0G52C2C0A2A2A11G2A0C6A4G14A13A6\tRG:Z:GATKSVContigAlignments\tNM:i:15\tAS:i:66\tXS:i:0";
        alignmentOne = fromSAMRecordString(alignmentOneSAMString, true);
        alignmentTwo = fromSAMRecordString(alignmentTwoSAMString, true);
        data.add(new Object[]{alignmentOne, alignmentTwo, SimpleNovelAdjacencyInterpreter.MORE_RELAXED_ALIGNMENT_MIN_MQ, SimpleNovelAdjacencyInterpreter.MORE_RELAXED_ALIGNMENT_MIN_LENGTH, false});

        alignmentOne = new AlignmentInterval(new SimpleInterval("chr1:10000-10028"), 1, 29, TextCigarCodec.decode("29M"), true, 60, 0, 29, ContigAlignmentsModifier.AlnModType.NONE);
        alignmentTwo = new AlignmentInterval(new SimpleInterval("chr1:10201-10501"), 30, 330, TextCigarCodec.decode("301M"), true, 60, 6, 295, ContigAlignmentsModifier.AlnModType.NONE);
        data.add(new Object[]{alignmentOne, alignmentTwo, SimpleNovelAdjacencyInterpreter.MORE_RELAXED_ALIGNMENT_MIN_MQ, SimpleNovelAdjacencyInterpreter.MORE_RELAXED_ALIGNMENT_MIN_LENGTH, false});

        alignmentOne = new AlignmentInterval(new SimpleInterval("chr1:10201-10501"), 30, 330, TextCigarCodec.decode("301M"), false, 60, 6, 295, ContigAlignmentsModifier.AlnModType.NONE);
        alignmentTwo = new AlignmentInterval(new SimpleInterval("chr1:10000-10028"), 1, 29, TextCigarCodec.decode("29M"), false, 60, 0, 29, ContigAlignmentsModifier.AlnModType.NONE);
        data.add(new Object[]{alignmentOne, alignmentTwo, SimpleNovelAdjacencyInterpreter.MORE_RELAXED_ALIGNMENT_MIN_MQ, SimpleNovelAdjacencyInterpreter.MORE_RELAXED_ALIGNMENT_MIN_LENGTH, false});


        return data.toArray(new Object[data.size()][]);
    }
    @Test(dataProvider = "forSplitPairStrongEnoughEvidenceForCA", groups = "sv")
    public void testSplitPairStrongEnoughEvidenceForCA(final AlignmentInterval intervalOne,
                                                       final AlignmentInterval intervalTwo,
                                                       final int mapQThresholdInclusive,
                                                       final int alignmentLengthThresholdInclusive,
                                                       final boolean expected) {
        Assert.assertEquals(SimpleChimera.splitPairStrongEnoughEvidenceForCA(intervalOne, intervalTwo, mapQThresholdInclusive, alignmentLengthThresholdInclusive),
                expected);
    }

    @Test(expectedExceptions = GATKException.class)
    public void testTypeInference_expectException() {
        final AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("21", 100001, 100100), 1 ,100, TextCigarCodec.decode("100M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        final AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("21", 100101, 100200), 101 ,200, TextCigarCodec.decode("100M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        final SimpleChimera simpleChimera = new SimpleChimera(region1, region2, Collections.emptyList(), "1",
                AssemblyContigWithFineTunedAlignments.NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME, b37_seqDict);
        simpleChimera.inferType(null);
    }


    private static final class TestData {
        private final AlignmentInterval one;
        private final AlignmentInterval two;
        private final SAMSequenceDictionary refDict;

        private final SimpleChimera.DistancesBetweenAlignmentsOnRefAndOnRead expectedDistances;
        private final StrandSwitch expectedStrandSwitch;
        private final boolean expectedIsForwardStrandRepresentation;
        private final boolean expectedFirstContigRegionHasLaterReferenceMapping;

        private final boolean isSimpleTranslocation;
        private final boolean isInversionOrInversionBreakpoint;

        private final TypeInferredFromSimpleChimera typeInferred;

        TestData(final AlignmentInterval one, final AlignmentInterval two, final SAMSequenceDictionary refDict, final SimpleChimera.DistancesBetweenAlignmentsOnRefAndOnRead expectedDistances, final StrandSwitch expectedStrandSwitch, final boolean expectedIsForwardStrandRepresentation, final boolean expectedFirstContigRegionHasLaterReferenceMapping, final boolean isSimpleTranslocation, final boolean isInversionOrInversionBreakpoint, final TypeInferredFromSimpleChimera typeInferred) {
            this.one = one;
            this.two = two;
            this.refDict = refDict;
            this.expectedDistances = expectedDistances;
            this.expectedStrandSwitch = expectedStrandSwitch;
            this.expectedIsForwardStrandRepresentation = expectedIsForwardStrandRepresentation;
            this.expectedFirstContigRegionHasLaterReferenceMapping = expectedFirstContigRegionHasLaterReferenceMapping;
            this.isSimpleTranslocation = isSimpleTranslocation;
            this.isInversionOrInversionBreakpoint = isInversionOrInversionBreakpoint;
            this.typeInferred = typeInferred;
        }
    }

    private List<TestData> casesForSimpleSymbolicVariants() {
        final List<TestData> result = new ArrayList<>(20);
        for (final AssemblyBasedSVDiscoveryTestDataProvider.AssemblyBasedSVDiscoveryTestDataForSimpleChimera testDataForSimpleSV : forSimpleSV.getAllTestData()) {
            result.add(new TestData(testDataForSimpleSV.firstAlignment, testDataForSimpleSV.secondAlignment, testDataForSimpleSV.getAppropriateDictionary(),
                    ((TestDataForSimpleSV)testDataForSimpleSV).expectedDistances, testDataForSimpleSV.expectedSimpleChimera.strandSwitch,
                    testDataForSimpleSV.expectedSimpleChimera.isForwardStrandRepresentation,
                    testDataForSimpleSV.expectedFirstContigRegionHasLaterReferenceMapping,
                    false, false,
                    testDataForSimpleSV.expectedNovelAdjacencyAndAltSeq.getTypeInferredFromSimpleChimera()));
        }
        return result;
    }

    private List<TestData> casesForInversion() {

        final List<TestData> result = new ArrayList<>(20);

        for (final AssemblyBasedSVDiscoveryTestDataProvider.AssemblyBasedSVDiscoveryTestDataForSimpleChimera testDataForSimpleSV : forInversionBreakpoints.getAllTestData()) {
            result.add(new TestData(testDataForSimpleSV.firstAlignment, testDataForSimpleSV.secondAlignment, testDataForSimpleSV.getAppropriateDictionary(),
                    null, testDataForSimpleSV.expectedSimpleChimera.strandSwitch,
                    testDataForSimpleSV.expectedSimpleChimera.isForwardStrandRepresentation,
                    testDataForSimpleSV.expectedFirstContigRegionHasLaterReferenceMapping, false, false,
                    testDataForSimpleSV.expectedNovelAdjacencyAndAltSeq.getTypeInferredFromSimpleChimera()));
        }

        return result;
    }

    private List<TestData> casesForInvertedDuplication() {

        final List<TestData> result = new ArrayList<>(20);

        AlignmentInterval intervalOne = new AlignmentInterval(new SimpleInterval("chr21:25625477-25625587"), 1, 111, TextCigarCodec.decode("111M212H"), false, 60, 0, 111, ContigAlignmentsModifier.AlnModType.NONE);
        AlignmentInterval intervalTwo = new AlignmentInterval(new SimpleInterval("chr21:25625379-25625595"), 107, 323, TextCigarCodec.decode("106S217M"), true, 60, 0, 127, ContigAlignmentsModifier.AlnModType.NONE);
        result.add(new TestData(intervalOne, intervalTwo, bareBoneHg38SAMSeqDict, null,
                StrandSwitch.REVERSE_TO_FORWARD, false, true,
                false, true, TypeInferredFromSimpleChimera.INTRA_CHR_STRAND_SWITCH_33));

        intervalOne = new AlignmentInterval(new SimpleInterval("chr20", 48513458, 48513545), 1, 88, TextCigarCodec.decode("88M227H"), true, 39, 1, 83, ContigAlignmentsModifier.AlnModType.NONE);
        intervalTwo = new AlignmentInterval(new SimpleInterval("chr20", 48513297, 48513578), 84, 365, TextCigarCodec.decode("83S282M"), false, 60, 0, 282, ContigAlignmentsModifier.AlnModType.NONE);
        result.add(new TestData(intervalOne, intervalTwo, TestUtilsForAssemblyBasedSVDiscovery.bareBoneHg38SAMSeqDict, null,
                StrandSwitch.FORWARD_TO_REVERSE, true, true,
                false, true, TypeInferredFromSimpleChimera.INTRA_CHR_STRAND_SWITCH_55));
        return result;
    }

    private List<TestData> casesForBreakEndVariants() {

        final List<TestData> result = new ArrayList<>(20);
        for (final AssemblyBasedSVDiscoveryTestDataProvider.AssemblyBasedSVDiscoveryTestDataForSimpleChimera testDataForSimpleSV : forBreakEndVariants.getAllTestData()) {
            result.add(new TestData(testDataForSimpleSV.firstAlignment, testDataForSimpleSV.secondAlignment, testDataForSimpleSV.getAppropriateDictionary(),
                    null, testDataForSimpleSV.expectedSimpleChimera.strandSwitch,
                    testDataForSimpleSV.expectedSimpleChimera.isForwardStrandRepresentation,
                    testDataForSimpleSV.expectedFirstContigRegionHasLaterReferenceMapping,
                    true, false,
                    testDataForSimpleSV.expectedNovelAdjacencyAndAltSeq.getTypeInferredFromSimpleChimera()));
        }

        return result;
    }


    @DataProvider(name = "testRepresentationAndSerialization")
    private Object[][] testRepresentationAndSerialization() {
        final List<Object[]> data = new ArrayList<>(50);

        casesForSimpleSymbolicVariants().forEach(obj -> data.add(new Object[]{obj}));
        casesForInversion().forEach(obj -> data.add(new Object[]{obj}));
        casesForInvertedDuplication().forEach(obj -> data.add(new Object[]{obj}));
        casesForBreakEndVariants().forEach(obj -> data.add(new Object[]{obj}));

        return data.toArray(new Object[data.size()][]);
    }

    @Test(dataProvider = "testRepresentationAndSerialization", groups = "sv")
    public void testRepresentationAndSerialization(final TestData testData) {

        final AlignmentInterval alignmentOne = testData.one;
        final AlignmentInterval alignmentTwo = testData.two;
        final SAMSequenceDictionary refDict = testData.refDict;
        final StrandSwitch expectedStrandSwitch = testData.expectedStrandSwitch;
        final boolean expectedIsForwardStrandRepresentation = testData.expectedIsForwardStrandRepresentation;
        final boolean expectedFirstContigRegionHasLaterReferenceMapping = testData.expectedFirstContigRegionHasLaterReferenceMapping;
        final SimpleChimera.DistancesBetweenAlignmentsOnRefAndOnRead expectedDistances = testData.expectedDistances;
        final boolean expectedIsSimpleTranslocation = testData.isSimpleTranslocation;

        Assert.assertEquals(SimpleChimera.determineStrandSwitch(alignmentOne, alignmentTwo), expectedStrandSwitch);
        Assert.assertEquals(SimpleChimera.isForwardStrandRepresentation(alignmentOne, alignmentTwo, expectedStrandSwitch, refDict),
                            expectedIsForwardStrandRepresentation);

        final SimpleChimera simpleChimera = new SimpleChimera(alignmentOne, alignmentTwo, Collections.emptyList(),
                "dummyName", NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME, refDict);
        Tuple2<SimpleInterval, SimpleInterval> coordinateSortedRefSpans = simpleChimera.getCoordinateSortedRefSpans(refDict);
        Assert.assertTrue( IntervalUtils.compareLocatables(coordinateSortedRefSpans._1(), coordinateSortedRefSpans._2(), refDict) < 0 );
        Assert.assertEquals(IntervalUtils.compareLocatables(alignmentOne.referenceSpan, alignmentTwo.referenceSpan, refDict) > 0,
                            expectedFirstContigRegionHasLaterReferenceMapping);
        if (testData.expectedDistances != null) {
            Assert.assertEquals(simpleChimera.getDistancesBetweenAlignmentsOnRefAndOnRead(),
                                expectedDistances);
        }
        Assert.assertEquals(simpleChimera.isCandidateSimpleTranslocation(), expectedIsSimpleTranslocation);
        Assert.assertEquals(simpleChimera.inferType(testData.refDict), testData.typeInferred);

        final ByteArrayOutputStream bos = new ByteArrayOutputStream();
        final Output out = new Output(bos);
        final Kryo kryo = new Kryo();
        kryo.writeClassAndObject(out, simpleChimera);
        out.flush();

        final ByteArrayInputStream bis = new ByteArrayInputStream(bos.toByteArray());
        final Input in = new Input(bis);
        @SuppressWarnings("unchecked")
        final SimpleChimera roundTrip = (SimpleChimera) kryo.readClassAndObject(in);
        Assert.assertEquals(roundTrip, simpleChimera);
        Assert.assertEquals(roundTrip.hashCode(), simpleChimera.hashCode());
    }
}