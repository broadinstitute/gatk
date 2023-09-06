package org.broadinstitute.hellbender.tools.walkers.fasta;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.BaseTest;
import org.broadinstitute.hellbender.testutils.FastaTestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

public class ShiftFastaIntegrationTest extends CommandLineProgramTest {

    private static final File MITO_REF = new File(toolsTestDir, "mutect/mito/Homo_sapiens_assembly38.mt_only.fasta");
    private static final File SHIFTED_MITO_REF = new File(largeFileTestDir + "mitochondria_references/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta");
    private static final File EXAMPLE_REF = new File(publicTestDir + "exampleFASTA.fasta");
    private static final File EXPECTED_SHIFTED_INTERVALS = new File(toolsTestDir, "shifted8000.shifted.intervals");
    private static final File EXPECTED_INTERVALS = new File(toolsTestDir, "shifted8000.intervals");

    @Test
    public void testShift8000() throws IOException {
        final File out = BaseTest.createTempFile("shifted8000", ".fasta");
        final File intervals_dir = BaseTest.createTempDir("shifted8000");
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(MITO_REF)
                .addOutput(out)
                .add("shift-back-output", "shiftback.chain")
                .add("shift-offset-list", 8000)
                .add("interval-file-name", intervals_dir.getAbsolutePath() + "/mito");

        runCommandLine(args);
        FastaTestUtils.assertFastaFilesContainTheSameSequence(out.toPath(), SHIFTED_MITO_REF.toPath());
        final List<String> intervals = Files.readAllLines(Paths.get(intervals_dir.getAbsolutePath() + "/mito.intervals"));
        final List<String> expected_intervals = Files.readAllLines(Paths.get(EXPECTED_INTERVALS.getAbsolutePath()));
        final List<String> shifted_intervals = Files.readAllLines(Paths.get(intervals_dir.getAbsolutePath() + "/mito.shifted.intervals"));
        final List<String> shifted_expected_intervals = Files.readAllLines(Paths.get(EXPECTED_SHIFTED_INTERVALS.getAbsolutePath()));
        Assert.assertEquals(intervals.size(), expected_intervals.size());
        Assert.assertEquals(shifted_intervals.size(), shifted_expected_intervals.size());
        for (int i=0; i<intervals.size(); i++) {
            Assert.assertEquals(intervals.get(i), expected_intervals.get(i));
            Assert.assertEquals(shifted_intervals.get(i), shifted_expected_intervals.get(i));
        }
    }

    @Test
    public void testShift() {
        final File out1 = BaseTest.createTempFile("shifted_example", ".fasta");
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(EXAMPLE_REF)
                .addOutput(out1)
                .add("shift-back-output", "shiftback.chain");

        runCommandLine(args);

        final File out2 = BaseTest.createTempFile("reshifted_example", ".fasta");
        ArgumentsBuilder args2 = new ArgumentsBuilder();
        args2.addReference(out1)
                .addOutput(out2)
                .add("shift-back-output", "shiftback.chain");

        runCommandLine(args2);
        FastaTestUtils.assertFastaFilesContainTheSameSequenceCaseInsensitive(out2.toPath(), EXAMPLE_REF.toPath());
    }
}
