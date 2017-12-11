package org.broadinstitute.hellbender.tools.funcotator;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * An integration test for the {@link Funcotator} tool.
 * Created by jonn on 8/29/17.
 */
public class FuncotatorIntegrationTest extends CommandLineProgramTest {

    //TODO: Add checks on the output file.

    //==================================================================================================================

    @DataProvider
    Object[][] provideDataForIntegrationTest() {
        return new Object[][] {
                {FuncotatorTestConstants.GTF_CHR3_FILE_NAME, FuncotatorTestConstants.HG19_CHR3_REFERENCE_FILE_NAME, FuncotatorTestConstants.GENCODE_TRANSCRIPT_FASTA_FILE_NAME, FuncotatorTestConstants.VARIANT_FILE_HG19_CHR3},
        };
    }

    //==================================================================================================================

    // TODO: REFACTOR ALL THIS:

    @Test(dataProvider = "provideDataForIntegrationTest")
    public void testRun(final String gtfFileName, final String referenceFileName, final String fastaFileName, final String variantFileName) throws IOException {
        final File outputFile = createTempFile("funcotator_tmp_out", ".vcf");
        final List<String> arguments = new ArrayList<>();

        arguments.add("-" + Funcotator.GTF_FILE_ARG_SHORT_NAME);
        arguments.add(gtfFileName);
        arguments.add("-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME);
        arguments.add(referenceFileName);
        arguments.add("-" + Funcotator.GENCODE_FASTA_ARG_NAME);
        arguments.add(fastaFileName);
        arguments.add("-" + StandardArgumentDefinitions.VARIANT_SHORT_NAME);
        arguments.add(variantFileName);
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());

        runCommandLine(arguments);
    }
}
