package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.PrintReads;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class WalkerOutputIntegrationTest extends CommandLineProgramTest {

    private static final File TEST_DATA_DIR = new File("src/test/resources/org/broadinstitute/hellbender/engine/");

    @Override
    public String getTestedClassName() {
        return PrintReads.class.getSimpleName();
    }

    @Test
    public void testPGOnByDefault() throws IOException {
        final File inFile = new File(TEST_DATA_DIR, "testPGOnByDefault.bam");
        final File outFile = GATKBaseTest.createTempFile("testNoConflictRG", ".bam");
        final String[] args = new String[] {
                //NOTE: no addOutputSAMProgramRecord argument
                "--input" , inFile.getAbsolutePath(),
                "--output", outFile.getAbsolutePath()
        };
        runCommandLine(args);

        //input has GATK PrintReads not NOT GATK PrintReads.1 in headers

        try(final SamReader reader1 = SamReaderFactory.makeDefault().open(inFile)){
            Assert.assertNotNull(reader1.getFileHeader().getProgramRecord("GATK PrintReads"));
        }

        try(final SamReader reader2 = SamReaderFactory.makeDefault().open(inFile)) {
            Assert.assertNull(reader2.getFileHeader().getProgramRecord("GATK PrintReads.1"));
        }
        //output has both GATK PrintReads and GATK PrintReads.1 in headers
        try(final SamReader reader3 = SamReaderFactory.makeDefault().open(outFile)) {
            Assert.assertNotNull(reader3.getFileHeader().getProgramRecord("GATK PrintReads"));
        }

        try(final SamReader reader4 = SamReaderFactory.makeDefault().open(outFile)) {
            Assert.assertNotNull(reader4.getFileHeader().getProgramRecord("GATK PrintReads.1"));
        }
    }
}
