package org.broadinstitute.hellbender.tools.funcotator;

import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

/**
 * Class to test the {@link FuncotatorDataSourceDownloader}.
 * Created by jonn on 8/24/18.
 */
public class FuncotatorDataSourceDownloaderIntegrationTest extends CommandLineProgramTest {

    //==================================================================================================================
    // Private Static Members:

    // Off by default because each test case takes ~1 hour to run:
    private static final boolean doFullScaleTests = false;

    //==================================================================================================================
    // Helper Methods:

    private Path getDataSourceRemotePath(final String dsTypeArg) {
        switch (dsTypeArg) {
            case FuncotatorDataSourceDownloader.SOMATIC_ARG_LONG_NAME:
                return FuncotatorDataSourceDownloader.SOMATIC_GCLOUD_DATASOURCES_PATH;
            case FuncotatorDataSourceDownloader.GERMLINE_ARG_LONG_NAME:
                return FuncotatorDataSourceDownloader.GERMLINE_GCLOUD_DATASOURCES_PATH;
            default: throw new GATKException("Data source type does not exist: " + dsTypeArg);
        }
    }

    private void verifyDataSourcesExistThenDeleteThem(final String dsTypeArg, final boolean doExtract) {
        // Get the path to our files:
        final Path currentPath          = IOUtils.getPath(".");
        final Path remoteDataSourcePath = getDataSourceRemotePath(dsTypeArg);
        final Path expectedDownloadedDataSourcePath = currentPath.resolve(remoteDataSourcePath.getFileName().toString());

        // Verify it exists and delete it:
        verifyDataSourcesExistThenDeleteThem(expectedDownloadedDataSourcePath, doExtract);
    }

    private void verifyDataSourcesExistThenDeleteThem(final Path expectedDownloadedDataSourcePath, final boolean doExtract) {

        // Make sure our file exists:
        Assert.assertTrue( Files.exists(expectedDownloadedDataSourcePath) );

        // Clean up the downloaded files:
        try {
            Files.delete(expectedDownloadedDataSourcePath);
            if ( doExtract ) {
                // Get the base name for our folder.
                // (this way we get rid of all extensions (i.e. both `tar` and `gz`):
                final String baseName = expectedDownloadedDataSourcePath.toFile().getName().replace(".tar.gz", "");
                final Path   extractedDataSourceFolder = expectedDownloadedDataSourcePath.resolveSibling(baseName);
                FileUtils.deleteDirectory(extractedDataSourceFolder.toFile());
            }
        }
        catch ( final IOException ex ) {
            throw new GATKException("Could not clean up downloaded data sources for testDownloadRealDataSources: " + expectedDownloadedDataSourcePath);
        }
    }

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

    @Test(enabled = doFullScaleTests,
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

        // Now verify we got the data sources and clean up the files
        // so we don't have up to 30 gigs of stuff lying around:
        verifyDataSourcesExistThenDeleteThem(dsTypeArg, doExtract);
    }

    @Test(dataProvider = "provideForTestDownloadSmallDummyDataSources",
            groups = {"bucket"})
    void testDownloadDummySmallDataSources(final boolean doOverwrite, final boolean doValidate, final boolean doExtract) {

        // Create an output location for the data sources to go:
        final File tmpDir = createTempDir("FuncotatorDataSourceDownloaderIntegrationTest_testDownloadDummySmallDataSources");
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

        // Now verify we got the data sources and clean up the files
        // so we don't have up to 30 gigs of stuff lying around:
        verifyDataSourcesExistThenDeleteThem(outputDataSourcesPath, doExtract);
    }


}
