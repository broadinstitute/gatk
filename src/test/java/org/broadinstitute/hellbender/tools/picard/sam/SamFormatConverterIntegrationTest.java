package org.broadinstitute.hellbender.tools.picard.sam;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.SamAssertionUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public final class SamFormatConverterIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File(getTestDataDir(), "picard/sam/SamFormatConverterTest");
    private static final File unmappedSam = new File(TEST_DATA_DIR, "unmapped.sam");
    private static final File unmappedBam = new File(TEST_DATA_DIR, "unmapped.bam");
    private static final File unmappedCram = new File(TEST_DATA_DIR, "unmapped.cram");
    private static final File referenceFile = new File(TEST_DATA_DIR, "basic.fasta");

    public String getTestedClassName() {
        return SamFormatConverter.class.getSimpleName();
    }

    @Test
    public void testSAMToBAM() throws IOException {
        convertFile(unmappedSam, unmappedBam, ".bam", null);
    }

    @Test
    public void testSAMToCRAM() throws IOException {
        convertFile(unmappedSam, unmappedCram, ".cram", referenceFile);
    }

    @Test
    public void testBAMToCRAM() throws IOException {
        convertFile(unmappedBam, unmappedCram, ".cram", referenceFile);
    }

    @Test
    public void testBAMToSAM() throws IOException {
        convertFile(unmappedBam, unmappedSam, ".sam", null);
    }

    @Test
    public void testCRAMToBAM() throws IOException {
        convertFile(unmappedCram, unmappedBam, ".bam", null);
    }

    @Test
    public void testCRAMToSAM() throws IOException {
        convertFile(unmappedCram, unmappedSam, ".sam", null);
    }

    private void convertFile(final File inputFile, final File fileToCompare, final String extension, final File referenceFile) throws IOException {
        final List<String> samFileConverterArgs = new ArrayList<>();
        samFileConverterArgs.add("--input");
        samFileConverterArgs.add(inputFile.getAbsolutePath());
        final File converterOutput = BaseTest.createTempFile("SamFileConverterTest." + inputFile.getName(), extension);
        samFileConverterArgs.add("--output");
        samFileConverterArgs.add(converterOutput.getAbsolutePath());
        if (null != referenceFile) {
            samFileConverterArgs.add("--R");
            samFileConverterArgs.add(referenceFile.getAbsolutePath());
        }
        runCommandLine(samFileConverterArgs);
        SamAssertionUtils.assertSamValid(converterOutput);
        SamAssertionUtils.assertSamsEqual(converterOutput, fileToCompare);
    }
}
