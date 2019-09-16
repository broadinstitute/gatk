package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramIntegrationTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.testutils.SamAssertionUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public final class PrintBAMHeaderIntegrationTest extends CommandLineProgramIntegrationTest {

    protected static final File TEST_DATA_DIR = getTestDataDir();

    @Test
    public void testPrintBAMHeader() throws IOException {
        final File inFile = new File(TEST_DATA_DIR, "print_reads.sorted.bam");
        final File outFile = GATKBaseTest.createTempFile("NA12878.chr17_69k_70k.dictFix.bam.header", ".txt");
        final String[] args = new String[] {
                "--input" , inFile.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_SAM_PROGRAM_RECORD,
                "--output", outFile.getAbsolutePath()
        };
        runCommandLine(args);

        //Make sure contents are the same
        IOUtil.assertFilesEqual(new File(TEST_DATA_DIR, "print_reads.sorted.bam.header.txt"), outFile);
    }
}