package org.broadinstitute.hellbender.tools.walkers.softcliprealignment;

import htsjdk.samtools.*;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.read.CigarBuilder;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

public class RealignSoftClippedReadsIntegrationTest extends CommandLineProgramTest {

    // If true, update the expected outputs in tests that assert an exact or approximate match vs. prior output,
    // instead of actually running the tests.
    public static final boolean UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS = false;

    private static final String BWA_IMAGE_PATH = SubsettingRealignmentEngineTest.BWA_IMAGE_PATH;
    private final String INPUT_BAM_FILE = SubsettingRealignmentEngineTest.BAM_FILE;
    private final String EXPECTED_BAM_FILE = SubsettingRealignmentEngineTest.TEST_DATA_DIR + "/test.expected.bam";

    @Test
    public void testSoftClipRealignment() throws IOException {
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addInput(INPUT_BAM_FILE);
        args.add(SubsettingRealignmentArgumentCollection.BWA_IMAGE_LONG_NAME, BWA_IMAGE_PATH);
        if (UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS) {
            args.addOutput(EXPECTED_BAM_FILE);
            runCommandLine(args);
        } else {
            args.addOutput("%s");
            IntegrationTestSpec testSpec = new IntegrationTestSpec(args.toString(), Arrays.asList(EXPECTED_BAM_FILE));
            testSpec.executeTest("testSoftClipRealignment", this);
        }
    }

    /**
     * Make sure that someone didn't leave the {@value UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS} toggle turned on
     */
    @Test
    public void assertThatExpectedOutputUpdateToggleIsDisabled() {
        Assert.assertFalse(UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS, "The toggle to update expected outputs should not be left enabled");
    }

    @DataProvider(name = "testCheckIfClippedData")
    public Object[][] testCheckIfClippedData() {
        return new Object[][]{
                // 150M
                {new CigarBuilder().add(new CigarElement(150, CigarOperator.MATCH_OR_MISMATCH)), 1, false},
                {new CigarBuilder().add(new CigarElement(150, CigarOperator.MATCH_OR_MISMATCH)), 0, true},
                // 149M/1S
                {new CigarBuilder()
                        .add(new CigarElement(149, CigarOperator.MATCH_OR_MISMATCH))
                        .add(new CigarElement(1, CigarOperator.SOFT_CLIP)), 0, true},
                {new CigarBuilder()
                        .add(new CigarElement(149, CigarOperator.MATCH_OR_MISMATCH))
                        .add(new CigarElement(1, CigarOperator.SOFT_CLIP)), 1, true},
                {new CigarBuilder()
                        .add(new CigarElement(149, CigarOperator.MATCH_OR_MISMATCH))
                        .add(new CigarElement(1, CigarOperator.SOFT_CLIP)), 2, false},
                // 1S/149M
                {new CigarBuilder()
                        .add(new CigarElement(1, CigarOperator.SOFT_CLIP))
                        .add(new CigarElement(149, CigarOperator.MATCH_OR_MISMATCH)), 1, true},
                {new CigarBuilder()
                        .add(new CigarElement(1, CigarOperator.SOFT_CLIP))
                        .add(new CigarElement(149, CigarOperator.MATCH_OR_MISMATCH)), 2, false},
                // 50M/10I/40M/5D/40M/5S
                {new CigarBuilder()
                        .add(new CigarElement(50, CigarOperator.MATCH_OR_MISMATCH))
                        .add(new CigarElement(10, CigarOperator.INSERTION))
                        .add(new CigarElement(40, CigarOperator.MATCH_OR_MISMATCH))
                        .add(new CigarElement(5, CigarOperator.DELETION))
                        .add(new CigarElement(40, CigarOperator.MATCH_OR_MISMATCH))
                        .add(new CigarElement(5, CigarOperator.SOFT_CLIP)), 5, true},
                {new CigarBuilder()
                        .add(new CigarElement(50, CigarOperator.MATCH_OR_MISMATCH))
                        .add(new CigarElement(10, CigarOperator.INSERTION))
                        .add(new CigarElement(40, CigarOperator.MATCH_OR_MISMATCH))
                        .add(new CigarElement(5, CigarOperator.DELETION))
                        .add(new CigarElement(40, CigarOperator.MATCH_OR_MISMATCH))
                        .add(new CigarElement(5, CigarOperator.SOFT_CLIP)), 10, false},
        };
    }

    @Test(dataProvider= "testCheckIfClippedData")
    public void testCheckIfClipped(final CigarBuilder cigar, final int minSoftClipLength, final boolean expected) {
        final GATKRead read = new SAMRecordToGATKReadAdapter(new SAMRecord(new SAMFileHeader()));
        read.setName("test_read");
        read.setCigar(cigar.make());

        // Check soft-clip detection is working
        final boolean result = RealignSoftClippedReads.isValidSoftClip(read, minSoftClipLength);
        Assert.assertEquals(result, expected);

        // Check that the read name gets added to the set if it's a valid clip
        final Set<String> names = new HashSet<>();
        RealignSoftClippedReads.checkIfClipped(read, names, minSoftClipLength);
        if (result) {
            Assert.assertEquals(names.size(), 1);
            Assert.assertTrue(names.contains(read.getName()));
        } else {
            Assert.assertTrue(names.isEmpty());
        }

        // Supplementary read always false
        final GATKRead supplementary = read.copy();
        supplementary.setIsSupplementaryAlignment(true);
        names.clear();
        Assert.assertFalse(RealignSoftClippedReads.isValidSoftClip(supplementary, minSoftClipLength));
        RealignSoftClippedReads.checkIfClipped(supplementary, names, minSoftClipLength);
        Assert.assertTrue(names.isEmpty());

        // Secondary read always false
        final GATKRead secondary = read.copy();
        secondary.setIsSecondaryAlignment(true);
        Assert.assertFalse(RealignSoftClippedReads.isValidSoftClip(secondary, minSoftClipLength));
        RealignSoftClippedReads.checkIfClipped(secondary, names, minSoftClipLength);
        Assert.assertTrue(names.isEmpty());
    }
}