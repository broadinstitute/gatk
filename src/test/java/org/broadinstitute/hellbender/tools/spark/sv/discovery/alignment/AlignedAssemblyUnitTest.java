package org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment;


import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.TestUtilsForAssemblyBasedSVDiscovery;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.AlignedAssemblyOrExcuse;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class AlignedAssemblyUnitTest extends GATKBaseTest {
    private static final String dummyRefName = "1";
    private static final List<String> refNames = Collections.singletonList(dummyRefName);


    /**
     * @see AlignmentIntervalUnitTest#createInputsAndExpectedResults_BwaMemAlignmentConstruction()
     * 4 contigs of 3 assemblies pointing to the same inversion event.
     *  contig 1 belongs to assembly 1
     *  contig 2 and 3 belong to assembly 2
     *  contig 4 belongs to assembly 3.
     */
    @DataProvider(name = "AlignedAssemblySerializationTest")
    private Object[][] createInputsAndExpectedResults_Serialization() {

        final int[] alignmentStartsOnRef_0Based = {96, 196, 195, 95, 101, 201, 101, 201};
        final int[] alignmentStartsOnTig_0BasedInclusive = {0, 4, 0, 5, 0, 6, 0, 7};
        final int[] alignmentEndsOnTig_0BasedExclusive = {4, 8, 5, 10, 6, 12, 7, 14};
        final int[] seqLen = {8, 8, 10, 10, 12, 12, 14, 14};
        final int[] mapQual = {0, 1, 10, 20, 30, 40, 50, 60};
        final int[] mismatches = {0, 1, 1, 0, 2, 3, 3, 2};
        final boolean[] strandedness = {true, false, true, false, false, true, false, true};
        final String[] cigarStrings = {"4M4S", "4M4H", "5M5S", "5M5H", "6S6M", "6H6M", "7S7M", "7H7M"}; // each different number represent a different contig's pair of chimeric alignments
        final Cigar[] cigars = Arrays.stream(cigarStrings).map(TextCigarCodec::decode).toArray(Cigar[]::new);

        // these sequence are technically wrong the for the inversion event, but the test purpose is for serialization so it is irrelevant
        final byte[] dummySequenceForContigOne = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence(seqLen[0], (byte)'A');

        final byte[] dummySequenceForContigTwo = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence(seqLen[0], (byte)'T');
        final byte[] dummySequenceForContigThree = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence(seqLen[0], (byte)'C');

        final byte[] dummySequenceForContigFour = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence(seqLen[0], (byte)'G');


        final List<AlignedContig> allContigs = new ArrayList<>();

        for(int pair=0; pair<cigars.length/2; ++pair) {

            final List<AlignmentInterval> alignmentIntervalsForSimpleInversion = new ArrayList<>(8);
            final SimpleInterval referenceIntervalLeft = new SimpleInterval(refNames.get(0), alignmentStartsOnRef_0Based[2*pair]+1, alignmentStartsOnRef_0Based[2*pair]+cigars[2*pair].getReferenceLength());
            final AlignmentInterval alignmentIntervalLeft = new AlignmentInterval(referenceIntervalLeft, alignmentStartsOnTig_0BasedInclusive[2*pair]+1, alignmentEndsOnTig_0BasedExclusive[2*pair],
                    strandedness[2*pair] ? cigars[2*pair] : CigarUtils.invertCigar(cigars[2*pair]),
                    strandedness[2*pair], mapQual[2*pair], mismatches[2*pair], 100, ContigAlignmentsModifier.AlnModType.NONE);
            alignmentIntervalsForSimpleInversion.add(alignmentIntervalLeft);
            final SimpleInterval referenceIntervalRight = new SimpleInterval(refNames.get(0), alignmentStartsOnRef_0Based[2*pair+1]+1, alignmentStartsOnRef_0Based[2*pair+1]+cigars[2*pair+1].getReferenceLength());
            final AlignmentInterval alignmentIntervalRight = new AlignmentInterval(referenceIntervalRight, alignmentStartsOnTig_0BasedInclusive[2*pair+1]+1, alignmentEndsOnTig_0BasedExclusive[2*pair+1],
                    strandedness[2*pair+1] ? cigars[2*pair+1] : CigarUtils.invertCigar(cigars[2*pair+1]),
                    strandedness[2*pair+1], mapQual[2*pair+1], mismatches[2*pair+1], 100, ContigAlignmentsModifier.AlnModType.NONE);
            alignmentIntervalsForSimpleInversion.add(alignmentIntervalRight);

            if (pair == 0) {
                allContigs.add( new AlignedContig(AlignedAssemblyOrExcuse.formatContigName(0, 0), dummySequenceForContigOne, alignmentIntervalsForSimpleInversion) );
            } else if (pair <3) {
                allContigs.add( new AlignedContig(AlignedAssemblyOrExcuse.formatContigName(1, pair-1), pair==1 ? dummySequenceForContigTwo : dummySequenceForContigThree, alignmentIntervalsForSimpleInversion) );
            } else {
                allContigs.add( new AlignedContig(AlignedAssemblyOrExcuse.formatContigName(2, 0), dummySequenceForContigFour, alignmentIntervalsForSimpleInversion) );
            }
        }

        final Object[][] data = new Object[3][];
        data[0] = new Object[]{0, new AlignedAssembly(0, allContigs.subList(0, 1))};
        data[1] = new Object[]{1, new AlignedAssembly(1, allContigs.subList(1, 3))};
        data[2] = new Object[]{2, new AlignedAssembly(2, allContigs.subList(3, 4))};

        return data;
    }

    @Test(dataProvider="AlignedAssemblySerializationTest", groups = "sv")
    public void testAlignedAssemblySerialization(final Integer assemblyID, final AlignedAssembly expectedAssembly) {

        final ByteArrayOutputStream bos = new ByteArrayOutputStream();
        final Output out = new Output(bos);
        final Kryo kryo = new Kryo();
        kryo.writeClassAndObject(out, expectedAssembly);
        out.flush();

        final ByteArrayInputStream bis = new ByteArrayInputStream(bos.toByteArray());
        final Input in = new Input(bis);
        @SuppressWarnings("unchecked")
        final AlignedAssembly roundTrip = (AlignedAssembly) kryo.readClassAndObject(in);
        Assert.assertEquals(roundTrip, expectedAssembly);
    }
}
