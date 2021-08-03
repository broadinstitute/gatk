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
 * Class to test the {@link FuncotatorDataSourceBundler}.
 * Created by Hailey on 8/2/21.
 */

public class FuncotatorDataSourceBundlerIntegrationTest extends CommandLineProgramTest {

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
    private Object[][] provideForTestDownload() {
        return new Object[][]{
                {
                        FuncotatorDataSourceBundler.BACTERIA_ARG_LONG_NAME,
                        "absiella_dolichum_dsm_3991_gca_000154285/",
                        true,
                        true
                }
//                {
//                        FuncotatorDataSourceBundler.BACTERIA_ARG_LONG_NAME,
//                        "species-name acinebactor_baumannii_aye",
//                        true
//                },
//                {
//                        FuncotatorDataSourceBundler.BACTERIA_ARG_LONG_NAME,
//                        "acinebactor_baumannii_aye",
//                        false
//                },
//                {
//                        FuncotatorDataSourceBundler.FUNGI_ARG_LONG_NAME,
//                        "acinebactor_baumannii_aye",
//                        true
//                },
//                {
//                        FuncotatorDataSourceBundler.FUNGI_ARG_LONG_NAME,
//                        "acinebactor_baumannii_aye",
//                        true
//                },
//                {
//                        FuncotatorDataSourceBundler.FUNGI_ARG_LONG_NAME,
//                        "acinebactor_baumannii_aye",
//                        false
//                },
//                {
//                        FuncotatorDataSourceBundler.METAZOA_ARG_LONG_NAME,
//                        "acinebactor_baumannii_aye",
//                        true
//                },
//                {
//                        FuncotatorDataSourceBundler.METAZOA_ARG_LONG_NAME,
//                        "acinebactor_baumannii_aye",
//                        true
//                },
//                {
//                        FuncotatorDataSourceBundler.METAZOA_ARG_LONG_NAME,
//                        "acinebactor_baumannii_aye",
//                        false
//                },
//                {
//                        FuncotatorDataSourceBundler.PLANTS_ARG_LONG_NAME,
//                        "acinebactor_baumannii_aye",
//                        true
//                },
//                {
//                        FuncotatorDataSourceBundler.PLANTS_ARG_LONG_NAME,
//                        "acinebactor_baumannii_aye",
//                        true
//                },
//                {
//                        FuncotatorDataSourceBundler.PLANTS_ARG_LONG_NAME,
//                        "acinebactor_baumannii_aye",
//                        false
//                },
//                {
//                        FuncotatorDataSourceBundler.PROTISTS_ARG_LONG_NAME,
//                        "acinebactor_baumannii_aye",
//                        true
//                },
//                {
//                        FuncotatorDataSourceBundler.PROTISTS_ARG_LONG_NAME,
//                        "acinebactor_baumannii_aye",
//                        true
//                },
//                {
//                        FuncotatorDataSourceBundler.PROTISTS_ARG_LONG_NAME,
//                        "acinebactor_baumannii_aye",
//                        false
//                }
        };
    }

    //==================================================================================================================
    // Tests:

    @Test( dataProvider = "provideForTestDownload")
    void testDownloadRealDataSources (final String dsTypeArg, final String speciesArg, final boolean doOverwrite, final boolean doExtract) {
        final ArgumentsBuilder arguments = new ArgumentsBuilder();

        File outputFile = new File("Absiella_dolichum_dsm_3991_gca_000154285.ASM15428v1.51.gtf");
        arguments.add(dsTypeArg, true);
        arguments.add("species-name", speciesArg);
        arguments.add("output", outputFile);
        arguments.add(FuncotatorDataSourceBundler.OVERWRITE_ARG_LONG_NAME, doOverwrite);
        arguments.add(FuncotatorDataSourceBundler.EXTRACT_AFTER_DOWNLOAD, doExtract);

        runCommandLine(arguments);

        // Now verify we got the data sources and clean up the files
        // so we don't have up to 30 gigs of stuff lying around:
        verifyDataSourcesExistThenDeleteThem(dsTypeArg, doExtract);
    }

}
