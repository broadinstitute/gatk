package org.broadinstitute.hellbender.tools.walkers;

import org.broadinstitute.hellbender.cmdline.CommandLineProgramIntegrationTest;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.SamAssertionUtils;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class ReadSanitizerIntegrationTest extends CommandLineProgramIntegrationTest {

    @DataProvider
    Object[][] provideForTestReadSanitizer() {
        return new Object[][] {
                {
                    publicTestDir + "read_sanitizer_test.sam",
                    publicTestDir + "read_sanitizer_expected.sam"
                }
        };
    }

    @Test(dataProvider = "provideForTestReadSanitizer")
    public void testReadSanitizer(final String inputFilePath, final String expectedOutputFilePath) {
        final ArgumentsBuilder arguments = new ArgumentsBuilder();

        final File outputFile = createTempFile("testReadSanitizer.tmpOut", ".sam");

        arguments.addInput(inputFilePath);
        arguments.addOutput(outputFile);
        arguments.addReference(hg38Reference);

        runCommandLine(arguments);

        // ====
        // Now we validate that the results are correct:
        try {
            SamAssertionUtils.assertSamsEqual(outputFile, new File(expectedOutputFilePath), new File(hg38Reference));
        }
        catch (final IOException ex ) {
            throw new GATKException("Cannot complete test!", ex);
        }

    }

}
