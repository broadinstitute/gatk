package org.broadinstitute.hellbender.tools.funcotator;

import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.DataSourceUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;


import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;


/**
 * Class to test the {@link FuncotatorDataSourceBundler} class.
 * Created by Hailey on 8/2/21.
 */

public class FuncotatorDataSourceBundlerIntegrationTest extends CommandLineProgramTest {

    //==================================================================================================================
    // Private Static Members:

    //==================================================================================================================
    // Helper Methods:

//    private String getDataSourceRemoteURL(final String dsTypeArg) {
//        return switch ( dsTypeArg ) {
//            case FuncotatorDataSourceBundler.BACTERIA_ARG_LONG_NAME -> FuncotatorDataSourceBundler.BACTERIA_BASE_URL;
//            case FuncotatorDataSourceBundler.FUNGI_ARG_LONG_NAME -> FuncotatorDataSourceBundler.FUNGI_BASE_URL;
//            case FuncotatorDataSourceBundler.METAZOA_ARG_LONG_NAME -> FuncotatorDataSourceBundler.METAZOA_BASE_URL;
//            case FuncotatorDataSourceBundler.PLANTS_ARG_LONG_NAME -> FuncotatorDataSourceBundler.PLANTS_BASE_URL;
//            case FuncotatorDataSourceBundler.PROTISTS_ARG_LONG_NAME -> FuncotatorDataSourceBundler.PROTISTS_BASE_URL;
//            default -> throw new GATKException("Data source type does not exist: " + dsTypeArg);
//        };
//    }

//    private void verifyDataSourcesExistThenDeleteThem(final String dsOrgArg, final String dsSpeciesArg, final boolean doExtract) {
//        // Get the path to our files:
//        final Path currentPath          = IOUtils.getPath(".");
//        final Path remoteDataSourcePath = IOUtils.getPath(getDataSourceRemoteURL(dsOrgArg) + "/" + FuncotatorDataSourceBundlerUtils.getDSFileName(dsOrgArg, dsSpeciesArg) + DataSourceUtils.GTF_GZ_EXTENSION);
//        final Path expectedDownloadedDataSourcePath = currentPath.resolve(remoteDataSourcePath.getFileName().toString());
//
//        // Verify it exists and delete it:
//        verifyDataSourcesExistThenDeleteThem(expectedDownloadedDataSourcePath, doExtract);
//    }
//    private void verifyDataSourcesExistThenDeleteThem(final Path expectedDownloadedDataSourcePath, final boolean doExtract) {
//
//        // Make sure our file exists:
//        Assert.assertTrue( Files.exists(expectedDownloadedDataSourcePath) );
//
//        // Clean up the downloaded files:
//        try {
//            Files.delete(expectedDownloadedDataSourcePath);
//            if ( doExtract ) {
//                // Get the base name for our folder.
//                // (this way we get rid of all extensions (i.e. both `tar` and `gz`):
//                final String baseName = expectedDownloadedDataSourcePath.toFile().getName().replace(".gtf.gz", "");
//                final Path   extractedDataSourceFolder = expectedDownloadedDataSourcePath.resolveSibling(baseName);
//                FileUtils.deleteDirectory(extractedDataSourceFolder.toFile());
//            }
//        }
//        catch ( final IOException ex ) {
//            throw new GATKException("Could not clean up downloaded data sources for testDownloadRealDataSources: " + expectedDownloadedDataSourcePath);
//        }
//    }

    //==================================================================================================================
    // Data Providers:

    @DataProvider
    private Object[][] provideForTestDownload() {
        return new Object[][]{
                {
                    FuncotatorDataSourceBundler.METAZOA_ARG_LONG_NAME,
                        "octopus_bimaculoides",
                        true,
                }
//                {
//                        FuncotatorDataSourceBundler.BACTERIA_ARG_LONG_NAME,
//                        "absiella_dolichum_dsm_3991_gca_000154285/",
//                        false,
//                        true
//                },
//                {
//                        FuncotatorDataSourceBundler.PLANTS_ARG_LONG_NAME,
//                        "actinidia_chinensis",
//                        true,
//                        true
//                }
        };
    }

    @DataProvider
    private Object[][] provideForTestDownloadWrong() {
        return new Object[][] {
                {
                    FuncotatorDataSourceBundler.METAZOA_ARG_LONG_NAME,
                        "absiella_dolichum_dsm_3991_gca_000154285",
                        false,
                }
        };
    }

    //==================================================================================================================
    // Tests:

    @Test( dataProvider = "provideForTestDownload")
    void testDownloadRealDataSources (final String dsOrgArg, final String dsSpeciesArg, final boolean doOverwrite) {
        final ArgumentsBuilder arguments = new ArgumentsBuilder();
        arguments.add(dsOrgArg, true);
        arguments.add("species-name", dsSpeciesArg);
        arguments.add(FuncotatorDataSourceBundler.OVERWRITE_ARG_LONG_NAME, doOverwrite);

        final Path basePath = getSafeNonExistentPath(dsSpeciesArg);
        arguments.add(FuncotatorDataSourceBundler.OUTPUT_DATASOURCES_FOLDER_ARG_NAME, basePath.toUri().toString());

        runCommandLine(arguments);

        // Now verify we got the data sources and clean up the files
        // so we don't have up to 30 gigs of stuff lying around:
//        verifyDataSourcesExistThenDeleteThem(dsOrgArg, dsSpeciesArg, doExtract);
    }

    // To do: need to make some integration tests for running with incorrect input arguments

    @Test( dataProvider = "provideForTestDownloadWrong", expectedExceptions = UserException.BadInput.class)
    void testDownloadWrongDataSources (final String dsOrgArg, final String dsSpeciesArg, final boolean doOverwrite) {
        final ArgumentsBuilder arguments = new ArgumentsBuilder();
        arguments.add(dsOrgArg, true);
        arguments.add("species-name", dsSpeciesArg);
        arguments.add(FuncotatorDataSourceBundler.OVERWRITE_ARG_LONG_NAME, doOverwrite);

        final Path basePath = getSafeNonExistentPath(dsSpeciesArg);
        arguments.add(FuncotatorDataSourceBundler.OUTPUT_DATASOURCES_FOLDER_ARG_NAME, basePath.toUri().toString());

        runCommandLine(arguments);
    }

}
