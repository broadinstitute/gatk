package org.broadinstitute.hellbender.tools.funcotator;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
import java.io.IOException;


/**
 * A unit test suite for the {@link FuncotatorDataSourceBundlerUtils} class.
 * Created by Hailey on 8/5/21.
 */
public class FuncotatorDataSourceBundlerUtilsUnitTest extends CommandLineProgramTest {

    //==================================================================================================================
    // Data Providers:

    @DataProvider
    Object[][] provideForTestGetDatasourceBaseName() {
        return new Object[][] {
                { "plants", "actinidia_chinensis", "Actinidia_chinensis.Red5_PS1_1.69.0.51" }
        };
    }

    @DataProvider
    Object[][] provideForTestGetFastaFileName() {
        return new Object[][] {
                { "plants", "actinidia_chinensis", "Actinidia_chinensis.Red5_PS1_1.69.0.cdna.all" }
        };
    }

    @DataProvider
    Object[][] provideForTestBuildMapWrong() {
        return new Object[][]{
                {
                        "plants",
                        "Aphanomyces astaci GCA_002197585.2"
                }
        };
    }

    @DataProvider
    Object[][] provideForTestExtractGtfGz() {
        return new Object[][] {
            {
                new File(largeFileTestDir + "funcotator/Acyrthosiphon_pisum.Acyr_2.0.51.PARTIAL.gtf.gz").getAbsolutePath(),
                new File(largeFileTestDir + "funcotator/Acyrthosiphon_pisum.Acyr_2.0.51.PARTIAL.gtf").getAbsolutePath(),
            },
        };
    }

    //==================================================================================================================
    // Tests:

    @Test(dataProvider = "provideForTestGetDatasourceBaseName")
    void testGetDatasourceBaseName(final String orgName, final String speciesName, final String expectedBaseName) {

        Assert.assertEquals(
                expectedBaseName,
                FuncotatorDataSourceBundlerUtils.getDatasourceBaseName(orgName, speciesName)
        );
    }

    @Test(dataProvider = "provideForTestGetFastaFileName")
    void testGetFastaFileName(final String orgName, final String speciesName, final String expectedFastaName) {
        Assert.assertEquals(
                FuncotatorDataSourceBundlerUtils.getFastaFileName(orgName, speciesName),
                expectedFastaName
        );
    }

    @Test(dataProvider = "provideForTestBuildMapWrong", expectedExceptions = UserException.BadInput.class)
    void testBuildMapWrong(final String orgName, final String speciesName) {
        String testIncorrectFileName = FuncotatorDataSourceBundlerUtils.getDatasourceBaseName(orgName, speciesName);
        String incorrectFastaFileName = FuncotatorDataSourceBundlerUtils.getDatasourceBaseName(orgName, speciesName);
    }

    @Test(dataProvider = "provideForTestExtractGtfGz")
    void testExtractGtfGz(final String gtfGzFilePath, final String expectedFilePath) {

        final File destFile = getSafeNonExistentFile("testExtractGtfGz");
        FuncotatorDataSourceBundlerUtils.extractGzFile(gtfGzFilePath, destFile.getAbsolutePath(), true);

        try {
            IntegrationTestSpec.assertEqualTextFiles(destFile, new File(expectedFilePath), "#", false);
        }
        catch (final IOException ex) {
            throw new GATKException("Error: IOException encountered during test execution!", ex);
        }
    }

    //To do: write unit tests for testing with an invalid organism name and for testing with an invalid species name
}
