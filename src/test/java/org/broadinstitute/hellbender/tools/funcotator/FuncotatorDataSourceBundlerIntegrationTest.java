package org.broadinstitute.hellbender.tools.funcotator;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;


import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;


/**
 * Class to test the {@link FuncotatorDataSourceBundler} class.
 * Created by Hailey on 8/2/21.
 */

public class FuncotatorDataSourceBundlerIntegrationTest extends CommandLineProgramTest {

    //==================================================================================================================
    // Private Static Members:

    //==================================================================================================================
    // Helper Methods:

    //==================================================================================================================
    // Data Providers:

    @DataProvider
    private Object[][] provideForTestDownload() {
        return new Object[][]{
                {
                        FuncotatorDataSourceBundler.OrganismKingdom.METAZOA,
                        "octopus_bimaculoides",
                },
                {
                        FuncotatorDataSourceBundler.OrganismKingdom.BACTERIA,
                        "absiella_dolichum_dsm_3991_gca_000154285",
                },
                {
                        FuncotatorDataSourceBundler.OrganismKingdom.PLANTS,
                        "actinidia_chinensis",
                }
        };
    }

    @DataProvider
    private Object[][] provideForTestDownloadWrong() {
        return new Object[][]{
                {
                        FuncotatorDataSourceBundler.OrganismKingdom.METAZOA,
                        "absiella_dolichum_dsm_3991_gca_000154285",
                }
        };
    }

    //==================================================================================================================
    // Tests:

    @Test(dataProvider = "provideForTestDownload", groups = { "cloud" })
    public void testDownloadRealDataSources(final FuncotatorDataSourceBundler.OrganismKingdom kingdom, final String dsSpeciesArg) {
        final ArgumentsBuilder arguments = new ArgumentsBuilder();
        arguments.add(FuncotatorDataSourceBundler.ORGANISM_KINGDOM_ARG_LONG_NAME, kingdom);
        arguments.add(FuncotatorDataSourceBundler.SPECIES_ARG_LONG_NAME, dsSpeciesArg);

        final Path basePath = getSafeNonExistentPath(dsSpeciesArg);
        arguments.add(StandardArgumentDefinitions.OUTPUT_LONG_NAME, basePath.toUri().toString());

        runCommandLine(arguments);

        // Make sure we clean up:
        IOUtils.deleteOnExit(basePath);
    }

    // To do: need to make some integration tests for running with incorrect input arguments

    @Test(dataProvider = "provideForTestDownloadWrong", groups = { "cloud" }, expectedExceptions = UserException.BadInput.class)
    public void testDownloadWrongDataSources(final FuncotatorDataSourceBundler.OrganismKingdom kingdom, final String dsSpeciesArg) {
        final ArgumentsBuilder arguments = new ArgumentsBuilder();
        arguments.add(FuncotatorDataSourceBundler.ORGANISM_KINGDOM_ARG_LONG_NAME, kingdom);
        arguments.add(FuncotatorDataSourceBundler.SPECIES_ARG_LONG_NAME, dsSpeciesArg);

        final Path basePath = getSafeNonExistentPath(dsSpeciesArg);
        arguments.add(StandardArgumentDefinitions.OUTPUT_LONG_NAME, basePath.toUri().toString());

        runCommandLine(arguments);
    }

    @Test(groups = { "cloud" })
    public void testDownloadDatasourcesWithValidation() {

        final FuncotatorDataSourceBundler.OrganismKingdom kingdom      = FuncotatorDataSourceBundler.OrganismKingdom.BACTERIA;
        final String                                      dsSpeciesArg = "absiella_dolichum_dsm_3991_gca_000154285";

        final ArgumentsBuilder arguments = new ArgumentsBuilder();
        arguments.add(FuncotatorDataSourceBundler.ORGANISM_KINGDOM_ARG_LONG_NAME, kingdom);
        arguments.add(FuncotatorDataSourceBundler.SPECIES_ARG_LONG_NAME, dsSpeciesArg);

        final Path basePath = getSafeNonExistentPath(dsSpeciesArg);
        arguments.add(StandardArgumentDefinitions.OUTPUT_LONG_NAME, basePath.toUri().toString());

        runCommandLine(arguments);

        // Make sure we clean up:
        IOUtils.deleteOnExit(basePath);

        // Make sure our data sources exists:
        Assert.assertTrue(Files.exists(basePath));

        // Now make sure the files are all the same between the actual and expected datasources:
        final Path expectedPath = IOUtils.getPath(largeFileTestDir + "funcotator/absiella_dolichum_dsm_3991_gca_000154285_test_datasources");

        try {
            final Map<String, Path> expectedFiles = Files.walk(expectedPath).filter(Files::isRegularFile).collect(Collectors.toMap(
                    p -> p.toString().replace(expectedPath.toAbsolutePath().toString(), ""), p -> p
            ));
            final Map<String, Path> actualFiles = Files.walk(basePath).filter(Files::isRegularFile).collect(Collectors.toMap(
                    p -> p.toString().replace(basePath.toAbsolutePath().toString(), ""), p -> p
            ));

            // Assert that we have the same number of files in each list:
            Assert.assertEquals(expectedFiles.size(), actualFiles.size(), "Expected and actual datasources have a different number of files!");

            // Now we just check to see that our files are the same:
            int fileComparisonCount = 0;
            for (final String k : expectedFiles.keySet()) {

                final String fileName = expectedFiles.get(k).getFileName().toString();

                // We assume our mechanisms for indexing GTF files and creating sequence dictionaries are correct:
                if (fileName.endsWith(".gtf.idx") || fileName.endsWith(".dict")) {
                    fileComparisonCount++;
                    continue;
                }

                final String commentPrefix;
                if (fileName.equals(FuncotatorDataSourceBundlerHttpClient.MANIFEST_FILE_NAME) ||
                        fileName.equals(FuncotatorDataSourceBundlerHttpClient.TEMPLATE_CONFIG_FILE_NAME) ||
                        fileName.equals(FuncotatorDataSourceBundlerHttpClient.README_FILE_NAME)) {
                    commentPrefix = "Version:";
                }
                else {
                    commentPrefix = null;
                }

                // We just need to make sure that the generated dates / timestamps are ignored.
                // The easiest way to do this is by ignoring "Version:" from the text files:
                IntegrationTestSpec.assertEqualTextFiles(
                    actualFiles.get(k).toFile(),
                    expectedFiles.get(k).toFile(),
                        commentPrefix
                );

                fileComparisonCount++;
            }

            Assert.assertEquals(fileComparisonCount, expectedFiles.size(), "Did not compare all files!  Some files have unexpected names!");
        }
        catch ( final IOException ex ) {
            throw new GATKException("Error: IOException encountered during test execution!", ex);
        }
    }
}
