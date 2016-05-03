package org.broadinstitute.hellbender.tools.spark.sv;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.SkipException;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

/**
 * Tests RunBWAMEMViaCommandLine.
 * To make this test runnable, "-DPATHTOBWA" should be set properly.
 */
public final class RunBWAMEMViaCommandLineTest extends CommandLineProgramTest {

    private static final File TEST_DATA_DIR = new File(getTestDataDir(), "spark/sv/RunBWAMEMViaCommandLine");
    private static final File expectedSam = new File(TEST_DATA_DIR, "expected.sam");
    private static final Path ref = Paths.get(b37_reference_20_21);
    private static final Path bwaPath = getBwaPath();
    private static Path getBwaPath(){
        String s = System.getProperty("PATHTOBWA");
        if(null==s){
            s = System.getenv("PATHTOBWA");
        }
        return (null==s) ? null : Paths.get(s);
    }

    @BeforeMethod
    public void checkBWAPATHAvailability(){
        if(null==bwaPath){
            throw new SkipException("Skipping test because \"PATHTOBWA\" is set neither in system property nor as an environment variable.");
        }
    }

    @Test(groups="sv")
    public void testSeparate() throws IOException {

        final ArgumentsBuilder args = new ArgumentsBuilder();
        final File samOutput = boilerPlate(args);

        // input arguments
        final File input = new File(TEST_DATA_DIR, "input_1.fastq");
        final File secondInput = new File(TEST_DATA_DIR, "input_2.fastq");
        args.add(input.getAbsolutePath());
        args.add(secondInput.getAbsolutePath());

        this.runCommandLine(args.getArgsArray());

        Assert.assertEquals(getSamRecords(expectedSam), getSamRecords(samOutput));
    }

    @Test(groups="sv")
    public void testInterLeaved() throws IOException {

        final ArgumentsBuilder args = new ArgumentsBuilder();
        final File samOutput = boilerPlate(args);

        // input arguments
        args.add("-p");
        final File input = new File(TEST_DATA_DIR, "interleaved.fastq");
        args.add(input.getAbsolutePath());

        this.runCommandLine(args.getArgsArray());

        Assert.assertEquals(getSamRecords(expectedSam), getSamRecords(samOutput));
    }

    private static File boilerPlate(final ArgumentsBuilder args) throws IOException{

        args.add("-" + "bwaPath");
        args.add(bwaPath.toString());

        final File output = Files.createTempFile("test", "sam").toFile();
        output.deleteOnExit();
        args.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        args.add(output.getAbsolutePath());

        final File REF = ref.toFile();
        args.add("-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME);
        args.add(REF.getAbsolutePath());

        return output;
    }

    // utility function to extract just alignment records
    // header information, in particular the @PG line, depends on environment in which the test is run and
    // in which expected sam is generated, so is un-reproducible across different testing hosts
    private static List<String> getSamRecords (final File samFile) throws IOException{
        final List<String> samRecords = new ArrayList<>();
        try(final BufferedReader reader = new BufferedReader(new FileReader(samFile))){
            String line;
            while ( (line = reader.readLine()) != null) {
                if(!line.startsWith("@")){
                    samRecords.add(line);
                }
            }
            return samRecords;
        }
    }
}
