package org.broadinstitute.hellbender.tools.walkers.fasta;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.BaseTest;
import org.broadinstitute.hellbender.testutils.FastaTestUtils;
import org.testng.annotations.Test;

import java.io.File;

public class ShiftFastaIntegrationTest extends CommandLineProgramTest {

    private static final File MITO_REF = new File(toolsTestDir, "mutect/mito/Homo_sapiens_assembly38.mt_only.fasta");
    private static final File SHIFTED_MITO_REF = new File(largeFileTestDir + "mitochondria_references/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta");
    private static final File EXAMPLE_REF = new File(publicTestDir + "exampleFASTA.fasta");

    @Test
    public void testShift8000() {
        final File out = BaseTest.createTempFile("shifted8000", ".fasta");
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(MITO_REF)
                .addOutput(out)
                .add("shift-back-output", "shiftback.chain")
                .add("shift-offset-list", 8000);

        runCommandLine(args);
        FastaTestUtils.assertFastaFilesContainTheSameSequence(out.toPath(), SHIFTED_MITO_REF.toPath());
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
