package org.broadinstitute.hellbender.tools.picard.sam;

import htsjdk.samtools.BamFileIoUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.sam.SamAssertionUtils;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;

public class GatherBamFilesIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File(getTestDataDir(), "picard/sam/GatherBamFiles");
    private static final File ORIG_BAM = new File(TEST_DATA_DIR, "orig.bam");
    private static final List<File> SPLIT_BAMS = Arrays.asList(
            new File(TEST_DATA_DIR, "indUnknownChrom.bam"),
            new File(TEST_DATA_DIR, "indchr1.bam"),
            new File(TEST_DATA_DIR, "indchr2.bam"),
            new File(TEST_DATA_DIR, "indchr3.bam"),
            new File(TEST_DATA_DIR, "indchr4.bam"),
            new File(TEST_DATA_DIR, "indchr5.bam"),
            new File(TEST_DATA_DIR, "indchr6.bam"),
            new File(TEST_DATA_DIR, "indchr7.bam"),
            new File(TEST_DATA_DIR, "indchr8.bam")
    );

    public String getCommandLineProgramName() {
        return GatherBamFiles.class.getSimpleName();
    }

    @Test
    public void testTheGathering() throws Exception {
        final File outputFile = File.createTempFile("gatherBamFilesTest.samFile.", BamFileIoUtils.BAM_FILE_EXTENSION);
        final ArgumentsBuilder args = new ArgumentsBuilder();
        for (final File splitBam : SPLIT_BAMS) {
            args.add("INPUT=" + splitBam.getAbsolutePath());
        }
        args.add("OUTPUT=" + outputFile);
        runCommandLine(args.getArgsList());
        SamAssertionUtils.assertSamsEqual(ORIG_BAM, outputFile);
    }

    @Test
    public void sanityCheckTheGathering() throws Exception {
        final File outputFile = File.createTempFile("gatherBamFilesTest.samFile.", BamFileIoUtils.BAM_FILE_EXTENSION);
        final ArgumentsBuilder args = new ArgumentsBuilder();
        for (final File splitBam : SPLIT_BAMS) {
            args.add("INPUT=" + splitBam.getAbsolutePath());
        }
        args.add("OUTPUT=" + outputFile);
        runCommandLine(args.getArgsList());
        SamAssertionUtils.assertSamsEqual(ORIG_BAM, outputFile);
    }
}
