package org.broadinstitute.hellbender.tools;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.nio.file.Files;

public final class CountBasesIntegrationTest extends CommandLineProgramTest {

    @DataProvider(name="filenames")
    public Object[][] filenames() {
        return new String[][]{
                {"count_bases.sam", null},
                {"count_bases.bam", null},
                {"count_bases.cram", "count_bases.fasta"}
        };
    }

    @Test(dataProvider = "filenames")
    public void testCountBases(final String fileIn, final String referenceName) throws Exception {
        final File input = new File(getTestDataDir(), fileIn);
        final ArgumentsBuilder args = new ArgumentsBuilder();

        args.addInput(input);
        if (null != referenceName) {
            final File ref = new File(getTestDataDir(), referenceName);
            args.addReference(ref);
        }

        final Object res = runCommandLine(args);
        Assert.assertEquals(res, 808l);
    }

    @Test(dataProvider = "filenames")
    public void testCountBasesWithOutputFile(final String fileIn, final String referenceName) throws Exception {
        final File input = new File(getTestDataDir(), fileIn);
        final File output = createTempFile("testCountBasesWithOutputFile", ".txt");

        final ArgumentsBuilder args = new ArgumentsBuilder();

        args.addInput(input);
        args.addOutput(output);
        if (null != referenceName) {
            final File ref = new File(getTestDataDir(), referenceName);
            args.addReference(ref);
        }

        final Object res = runCommandLine(args);
        Assert.assertEquals(res, 808l);

        Assert.assertEquals(Files.readAllBytes(output.toPath()), "808".getBytes());
    }
}
