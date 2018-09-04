package org.broadinstitute.hellbender.tools.funcotator;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

/**
 * Class to test the {@link FuncotatorDataSourceDownloader}.
 * Created by jonn on 8/24/18.
 */
public class FuncotatorDataSourceDownloaderIntegrationTest extends CommandLineProgramTest {

    //==================================================================================================================
    // Private Static Members:

    private static final boolean doDebugTests = false;

    //==================================================================================================================
    // Data Providers:

    @DataProvider
    private Object[][] provideForTestDownload() {
        return new Object[][] {
                {
                        FuncotatorDataSourceDownloader.SOMATIC_ARG_LONG_NAME,
                        true,
                        true,
                        false
                },
                {
                        FuncotatorDataSourceDownloader.SOMATIC_ARG_LONG_NAME,
                        true,
                        false,
                        false
                },
                {
                        FuncotatorDataSourceDownloader.GERMLINE_ARG_LONG_NAME,
                        true,
                        true,
                        false
                },
                {
                        FuncotatorDataSourceDownloader.GERMLINE_ARG_LONG_NAME,
                        true,
                        false,
                        false
                },
                {
                        FuncotatorDataSourceDownloader.SOMATIC_ARG_LONG_NAME,
                        true,
                        false,
                        true
                }
        };
    }

    //==================================================================================================================
    // Tests:

    @Test(enabled = doDebugTests,
            dataProvider = "provideForTestDownload",
            groups = {"funcotatorValidation", "bucket"}
    )
    void testDownload(final String dsTypeArg, final boolean doOverwrite, final boolean doValidate, final boolean doExtract) {
        final ArgumentsBuilder arguments = new ArgumentsBuilder();

        arguments.addBooleanArgument(dsTypeArg, true);
        arguments.addBooleanArgument(FuncotatorDataSourceDownloader.OVERWRITE_ARG_LONG_NAME, doOverwrite);
        arguments.addBooleanArgument(FuncotatorDataSourceDownloader.VALIDATE_INTEGRITY_ARG_LONG_NAME, doValidate);
        arguments.addBooleanArgument(FuncotatorDataSourceDownloader.EXTRACT_AFTER_DOWNLOAD, doExtract);
        arguments.addArgument("verbosity", "INFO");

        runCommandLine(arguments);
    }


}
