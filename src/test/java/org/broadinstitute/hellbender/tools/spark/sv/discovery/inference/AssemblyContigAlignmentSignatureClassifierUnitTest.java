package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.StrandSwitch;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import scala.Tuple3;

import java.util.ArrayList;
import java.util.List;

public class AssemblyContigAlignmentSignatureClassifierUnitTest {

    @DataProvider(name = "forSignatureTest")
    private Object[][] forSignatureTest() {
        final List<Tuple3<AlignmentInterval, AlignmentInterval, SAMSequenceDictionary>> tuple3s = ChimericAlignmentUnitTest.alignmentPairsForSimpleChimeraAndRefSeqDict();
        final List<Object[]> data = new ArrayList<>(tuple3s.size());

        int i=0;
        // simple inversion
        data.add(new Object[]{tuple3s.get(i), StrandSwitch.REVERSE_TO_FORWARD, false, false, false}); ++i;

        data.add(new Object[]{tuple3s.get(i), StrandSwitch.FORWARD_TO_REVERSE, false, false, false}); ++i;

        data.add(new Object[]{tuple3s.get(i), StrandSwitch.FORWARD_TO_REVERSE, false, false, false}); ++i;

        data.add(new Object[]{tuple3s.get(i), StrandSwitch.REVERSE_TO_FORWARD, false, false, false}); ++i;

        data.add(new Object[]{tuple3s.get(i), StrandSwitch.REVERSE_TO_FORWARD, false, false, false}); ++i;

        // simple deletion
        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, false, false, false}); ++i;

        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, false, false, false}); ++i;

        // simple insertion
        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, false, false, false}); ++i;

        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, false, false, false}); ++i;

        // long range substitution
        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, false, false, false}); ++i;

        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, false, false, false}); ++i;

        // simple deletion with homology
        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, false, false, false}); ++i;

        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, false, false, false}); ++i;

        // tandem duplication simple contraction
        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, false, false, false}); ++i;

        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, false, false, false}); ++i;

        // tandem duplication simple expansion
        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, false, false, false}); ++i;

        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, false, false, false}); ++i;

        // tandem duplication simple expansion with novel insertion
        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, false, false, false}); ++i;

        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, false, false, false}); ++i;



        // first test (the original observed event, but assigned to a different chromosome): expansion from 1 unit to 2 units with pseudo-homology
        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, false, false, false}); ++i;

        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, false, false, false}); ++i;


        // second test: contraction from 2 units to 1 unit with pseudo-homology
        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, false, false, false}); ++i;

        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, false, false, false}); ++i;


        // third test: contraction from 3 units to 2 units without pseudo-homology
        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, false, false, false}); ++i;

        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, false, false, false}); ++i;


        // fourth test: expansion from 2 units to 3 units without pseudo-homology
        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, false, false, false}); ++i;

        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, false, false, false}); ++i;


        data.add(new Object[]{tuple3s.get(i), StrandSwitch.REVERSE_TO_FORWARD, false, true, true}); ++i;

        data.add(new Object[]{tuple3s.get(i), StrandSwitch.FORWARD_TO_REVERSE, false, true, true}); ++i;


        // same-chr translocation suspect, forward and reverse representation
        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, true, false, false}); ++i;

        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, true, false, false}); ++i;

        // diff-chr translocation suspect without SS
        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, true, false, false}); ++i;

        // diff-chr translocation suspect with SS
        data.add(new Object[]{tuple3s.get(i), StrandSwitch.FORWARD_TO_REVERSE, true, false, false}); ++i;

        // same-chr reference order switch, but overlaps (hence incomplete picture)
        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, false, false, true}); ++i;

        return data.toArray(new Object[data.size()][]);
    }

    @Test(dataProvider = "forSignatureTest", groups = "sv")
    public void testSignature(Tuple3<AlignmentInterval, AlignmentInterval, SAMSequenceDictionary> chimericPairsAndRefSeqDict,
                              final StrandSwitch expectedStrandSwitch,
                              final boolean expectedIsLikelySimpleTranslocation,
                              final boolean expectedIsLikelyInvDup,
                              final boolean expectedIsIncompletePicture) {

        final AlignmentInterval region1 = chimericPairsAndRefSeqDict._1();
        final AlignmentInterval region2 = chimericPairsAndRefSeqDict._2();

        Assert.assertEquals(AssemblyContigAlignmentSignatureClassifier.isCandidateSimpleTranslocation(region1, region2, expectedStrandSwitch), expectedIsLikelySimpleTranslocation);
        Assert.assertEquals(AssemblyContigAlignmentSignatureClassifier.isCandidateInvertedDuplication(region1, region2), expectedIsLikelyInvDup);
        Assert.assertEquals(AssemblyContigAlignmentSignatureClassifier.hasIncompletePictureFromTwoAlignments(region1, region2), expectedIsIncompletePicture);
    }
}
