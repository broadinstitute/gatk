package org.broadinstitute.hellbender.tools.funcotator;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.nio.file.Path;

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
    private Object[][] provideForTestDownloadSmallDummyDataSources() {
        return new Object[][] {
                {
                        true,
                        true,
                        false
                },
                {
                        true,
                        false,
                        false
                },
                {
                        true,
                        true,
                        true
                },
                {
                        true,
                        false,
                        true
                },
        };
    }

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
    void testDownloadRealDataSources(final String dsTypeArg, final boolean doOverwrite, final boolean doValidate, final boolean doExtract) {
        final ArgumentsBuilder arguments = new ArgumentsBuilder();

        arguments.addBooleanArgument(dsTypeArg, true);
        arguments.addBooleanArgument(FuncotatorDataSourceDownloader.OVERWRITE_ARG_LONG_NAME, doOverwrite);
        arguments.addBooleanArgument(FuncotatorDataSourceDownloader.VALIDATE_INTEGRITY_ARG_LONG_NAME, doValidate);
        arguments.addBooleanArgument(FuncotatorDataSourceDownloader.EXTRACT_AFTER_DOWNLOAD, doExtract);
        arguments.addArgument("verbosity", "INFO");

        runCommandLine(arguments);
    }

    @Test(dataProvider = "provideForTestDownloadSmallDummyDataSources",
            groups = {"bucket"})
    void testDownloadDummySmallDataSources(final boolean doOverwrite, final boolean doValidate, final boolean doExtract) {

        // Create an output location for the data sources to go:
        final File tmpDir = createTempDir("FuncotatorDataSourceDownloaderIntegrationTest_testDownloadDummySmallDataSources");
        tmpDir.deleteOnExit();
        final Path tmpDirPath = tmpDir.toPath();

        final Path outputDataSourcesPath = tmpDirPath.resolve(IOUtils.getPath(FuncotatorTestConstants.DUMMY_DATA_SOURCES_TAR_GZ).getFileName());

        final ArgumentsBuilder arguments = new ArgumentsBuilder();

        arguments.addArgument(FuncotatorDataSourceDownloader.TESTING_OVERRIDE_PATH_FOR_DATA_SOURCES_ARG, FuncotatorTestConstants.DUMMY_DATA_SOURCES_TAR_GZ);
        arguments.addArgument(FuncotatorDataSourceDownloader.TESTING_OVERRIDE_PATH_FOR_DATA_SOURCES_SHA256_ARG, FuncotatorTestConstants.DUMMY_DATA_SOURCES_TAR_GZ_SHA256_FILE);

        arguments.addArgument(StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputDataSourcesPath.toFile().getAbsolutePath());
        arguments.addBooleanArgument(FuncotatorDataSourceDownloader.OVERWRITE_ARG_LONG_NAME, doOverwrite);
        arguments.addBooleanArgument(FuncotatorDataSourceDownloader.VALIDATE_INTEGRITY_ARG_LONG_NAME, doValidate);
        arguments.addBooleanArgument(FuncotatorDataSourceDownloader.EXTRACT_AFTER_DOWNLOAD, doExtract);
        arguments.addArgument("verbosity", "INFO");

        runCommandLine(arguments);
    }


}
