package org.broadinstitute.hellbender.tools.picard.sam;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.read.SamAssertionUtils;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class SamFormatConverterIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File(getTestDataDir(), "picard/sam/SamFormatConverterTest");
    private static final File unmappedSam = new File(TEST_DATA_DIR, "unmapped.sam");
    private static final File unmappedBam = new File(TEST_DATA_DIR, "unmapped.bam");
    private static final File unmappedCram = new File(TEST_DATA_DIR, "unmapped.cram");

    public String getTestedClassName() {
        return SamFormatConverter.class.getSimpleName();
    }

    @Test
    public void testSAMToBAM() throws IOException {
        convertFile(unmappedSam, unmappedBam, ".bam");
    }

    @Test
    public void testSAMToCRAM() throws IOException {
        convertFile(unmappedSam, unmappedCram, ".cram");
    }

    @Test
    public void testBAMToCRAM() throws IOException {
        convertFile(unmappedBam, unmappedCram, ".cram");
    }

    @Test
    public void testBAMToSAM() throws IOException {
        convertFile(unmappedBam, unmappedSam, ".sam");
    }

    @Test
    public void testCRAMToBAM() throws IOException {
        convertFile(unmappedCram, unmappedBam, ".bam");
    }

    @Test
    public void testCRAMToSAM() throws IOException {
        convertFile(unmappedCram, unmappedSam, ".sam");
    }

    private void convertFile(final File inputFile, final File fileToCompare, final String extension) throws IOException {
        final List<String> samFileConverterArgs = new ArrayList<>();
        samFileConverterArgs.add("--INPUT");
        samFileConverterArgs.add(inputFile.getAbsolutePath());
        final File converterOutput = File.createTempFile("SamFileConverterTest." + inputFile.getName(), extension);
        samFileConverterArgs.add("--OUTPUT");
        samFileConverterArgs.add(converterOutput.getAbsolutePath());
        runCommandLine(samFileConverterArgs);
        SamAssertionUtils.assertSamValid(converterOutput);
        SamAssertionUtils.assertSamsEqual(converterOutput, fileToCompare);
    }
}
