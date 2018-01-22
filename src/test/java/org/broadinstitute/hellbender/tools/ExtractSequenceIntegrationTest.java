package org.broadinstitute.hellbender.tools;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * Integration test class for {@link ExtractSequence}.
 * Created by jonn on 1/19/18.
 */
public class ExtractSequenceIntegrationTest extends CommandLineProgramTest {

    //==================================================================================================================
    // Private Static Members:

    private static final String TEST_FASTA_FILE = GATKBaseTest.largeFileTestDir + File.separator + "human_g1k_v37.20.21.fasta";

    //==================================================================================================================
    // Private Members:

    //==================================================================================================================
    // Helper Methods:

    private void runExtractSequence(final String refFilePath, final String contig, final Integer start, final Integer end) {
        final StringBuilder outFileNameStringBuilder = new StringBuilder();
        outFileNameStringBuilder.append("ExtractSequenceTest_");
        outFileNameStringBuilder.append(contig);
        if ( start != null ) {
            outFileNameStringBuilder.append("s");
            outFileNameStringBuilder.append(start);
        }
        if ( end != null ) {
            outFileNameStringBuilder.append("e");
            outFileNameStringBuilder.append(end);
        }

        final File outputFile = createTempFile(outFileNameStringBuilder.toString(), ".fasta");
        final List<String> arguments = new ArrayList<>();

        arguments.add("-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME);
        arguments.add(refFilePath);
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());

        arguments.add("--" + ExtractSequence.CONTIG_ARG_LONG_NAME);
        arguments.add(contig);
        if ( start != null ) {
            arguments.add("--" + ExtractSequence.START_ARG_LONG_NAME);
            arguments.add(start.toString());
        }
        if ( end != null ) {
            arguments.add("--" + ExtractSequence.END_ARG_LONG_NAME);
            arguments.add(end.toString());
        }

        runCommandLine(arguments);

        // TODO: FILL THIS OUT!
//        IntegrationTestSpec.assertEqualTextFiles();
    }

    //==================================================================================================================
    // Data Providers:

    @DataProvider
    Object[][] provideForTestExtractSequence() {
        return new Object[][] {
            { TEST_FASTA_FILE, "20" , null , null },
            { TEST_FASTA_FILE, "20" , null , 7855 },
            { TEST_FASTA_FILE, "20" , 1 , 7855 },
            { TEST_FASTA_FILE, "20" , 1 , 10 },
            { TEST_FASTA_FILE, "20" , 10 , 10 },

            { TEST_FASTA_FILE, "21" , null , null },
            { TEST_FASTA_FILE, "21" , null , 7855 },
            { TEST_FASTA_FILE, "21" , 1 , 7855 },
            { TEST_FASTA_FILE, "21" , 1 , 10 },
            { TEST_FASTA_FILE, "21" , 10 , 10 },
        };
    }

    @DataProvider
    Object[][] provideForTestExtractSequenceFailures() {
        return new Object[][] {
            { TEST_FASTA_FILE, "22" , null , null },
            { TEST_FASTA_FILE, "20" , null , Integer.MAX_VALUE },
            { TEST_FASTA_FILE, "20" , -1 , 7855 },
            { TEST_FASTA_FILE, "20" , 10 , 1 },
            { TEST_FASTA_FILE, "20" , 10 , -1 },
            { TEST_FASTA_FILE, "20" , 0 , 10 },
        };
    }

    //==================================================================================================================
    // Tests:

    @Test(dataProvider = "provideForTestExtractSequence")
    public void testExtractSequence(final String refFilePath, final String contig, final Integer start, final Integer end) {
        runExtractSequence(refFilePath, contig, start, end);
    }

    @Test(
            dataProvider = "provideForTestExtractSequenceFailures",
            expectedExceptions = { UserException.class, IllegalArgumentException.class }
    )
    public void testExtractSequenceFailures(final String refFilePath, final String contig, final Integer start, final Integer end) {
        runExtractSequence(refFilePath, contig, start, end);
    }

}
