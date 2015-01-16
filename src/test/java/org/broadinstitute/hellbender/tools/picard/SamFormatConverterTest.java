package org.broadinstitute.hellbender.tools.picard;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.CompareSAMs;
import org.broadinstitute.hellbender.tools.picard.SamFormatConverter;
import org.broadinstitute.hellbender.tools.picard.ValidateSamFile;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertTrue;

public class SamFormatConverterTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File("src/test/resources/org/broadinstitute/hellbender/tools/SamFormatConverterTest");
    private static final File unmappedSam = new File(TEST_DATA_DIR, "unmapped.sam");
    private static final File unmappedBam = new File(TEST_DATA_DIR, "unmapped.bam");
    private static final File unmappedCram = new File(TEST_DATA_DIR, "unmapped.cram");

    public String getCommandLineProgramName() {
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
        final List<String> samFileConverterArgs = new ArrayList<String>();
        samFileConverterArgs.add("INPUT=" + inputFile);
        final File converterOutput = File.createTempFile("SamFileConverterTest." + inputFile.getName(), extension);
        samFileConverterArgs.add("OUTPUT=" + converterOutput);
        Assert.assertEquals(runCommandLine(samFileConverterArgs), null);

        ValidateSamFile validator = new ValidateSamFile();
        String[] validatorArgs = new String[]{"INPUT=" + converterOutput};
        Assert.assertEquals(validator.instanceMain(validatorArgs), true);

        // TODO this is a bit silly - since doWork is package-protected, we have to call instanceMain;
        // since instanceMain returns null, we rely on the object being mutated and call a public getter.
        CompareSAMs compareSams = new CompareSAMs();
        String[] compareSamsArgs = new String[]{converterOutput.toString(), fileToCompare.toString()};
        compareSams.instanceMain(compareSamsArgs);
        Assert.assertTrue(compareSams.areEqual());
    }
}
