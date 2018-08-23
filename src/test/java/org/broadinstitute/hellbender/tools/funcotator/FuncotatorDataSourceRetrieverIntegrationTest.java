package org.broadinstitute.hellbender.tools.funcotator;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.annotations.Test;

/**
 * Class to test the {@link FuncotatorDataSourceRetriever}.
 * Created by jonn on 8/24/18.
 */
public class FuncotatorDataSourceRetrieverIntegrationTest extends CommandLineProgramTest {

    //==================================================================================================================
    // Private Static Members:

    private static final boolean doDebugTests = false;

    //==================================================================================================================
    // Private Members:

    //==================================================================================================================
    // Helper Methods:

    //==================================================================================================================
    // Data Providers:

    //==================================================================================================================
    // Tests:

    @Test(enabled = doDebugTests,
            groups = {"funcotatorValidation", "bucket"}
    )
    void testGermlineDownloadWithValidation() {
        final ArgumentsBuilder arguments = new ArgumentsBuilder();

        arguments.addBooleanArgument(FuncotatorDataSourceRetriever.GERMLINE_ARG_LONG_NAME, true);
        arguments.addBooleanArgument(FuncotatorDataSourceRetriever.OVERWRITE_ARG_LONG_NAME, true);
        arguments.addBooleanArgument(FuncotatorDataSourceRetriever.VALIDATE_INTEGRITY_ARG_LONG_NAME, true);
        arguments.addArgument("verbosity", "INFO");

        runCommandLine(arguments);
    }

}
