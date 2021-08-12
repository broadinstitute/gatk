package org.broadinstitute.hellbender.tools.funcotator;

import org.apache.hadoop.yarn.webapp.hamlet2.Hamlet;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import org.broadinstitute.hellbender.utils.io.IOUtils;


/**
 * A unit test suite for the {@link FuncotatorDataSourceBundlerUtils} class.
 * Created by Hailey on 8/5/21.
 */
public class FuncotatorDataSourceBundlerUtilsUnitTest extends CommandLineProgramTest{

    //==================================================================================================================
    // Static Variables:
    private static final String TEST_WRONG_ORGNAME      = "virus";
    private static final String TEST_WRONG_SPECIESNAME  = "ashbya_gossypii";
    private static final String PLANTS_NAME             = "plants";
    private static final String PLANTS_SPECIES          = "actinidia_chinensis";
    private static final String PLANTS_FILE_NAME        = "Actinidia_chinensis.Red5_PS1_1.69.0.51.gtf.gz";
    private static final String PLANTS_FASTA_FILE_NAME  = "Actinidia_chinensis.Red5_PS1_1.69.0.cdna.all.fa.gz";
    private static final String PROTISTS_SPECIES        = "Aphanomyces astaci GCA_002197585.2";

    //==================================================================================================================
    // Data Providers:

    @DataProvider
    Object[][] provideForTestBuildMap() {
        return new Object[][] {
                {
                        PLANTS_NAME,
                        PLANTS_SPECIES
                }
        };
    }

    @DataProvider
    Object[][] provideForTestBuildMapWrong() {
        return new Object[][]{
                {
                        PLANTS_NAME,
                        PROTISTS_SPECIES
                }
        };
    }

    @DataProvider
    Object[][] provideForTestExtractGtfGz() {
        return new Object[][] {
                {
                        IOUtils.getPath("Acyrthosiphon_pisum.Acyr_2.0.51.gtf.gz").toAbsolutePath(),
                        IOUtils.getPath("Acyrthosiphon_pisum.Acyr_2.0.51.gtf.gtf").toAbsolutePath()
                },
                {
                        "/Users/hfox/Workspace/gatk/Acyrthosiphon_pisum.Acyr_2.0.51.gtf.gz",
                        "/Users/hfox/Workspace/gatk/Acyrthosiphon_pisum.Acyr_2.0.51.gtf"
                }
        };
    }

    //==================================================================================================================
    // Tests:

    @Test(dataProvider = "provideForTestBuildMap")
    void testBuildMap(String orgName, String speciesName) {

        String testCorrectFileName = FuncotatorDataSourceBundlerUtils.buildMapGetFileName(orgName, speciesName, false);
        if (!testCorrectFileName.equals(PLANTS_FILE_NAME)) {
            throw new UserException("Incorrect gtf file name. File name should be: " + PLANTS_FILE_NAME + ".");
        }

        String testCorrectFastaFileName = FuncotatorDataSourceBundlerUtils.buildMapGetFileName(orgName, speciesName, true);
        if (!testCorrectFastaFileName.equals(PLANTS_FASTA_FILE_NAME)) {
            throw new UserException("Incorrect fasta file name. File name should be: " + PLANTS_FASTA_FILE_NAME + ".");
        }
    }

    @Test(dataProvider = "provideForTestBuildMapWrong", expectedExceptions = UserException.BadInput.class)
    void testBuildMapWrong(String orgName, String speciesName) {
        String testIncorrectFileName = FuncotatorDataSourceBundlerUtils.buildMapGetFileName(orgName, speciesName, false);
        String incorrectFastaFileName = FuncotatorDataSourceBundlerUtils.buildMapGetFileName(orgName, speciesName, false);
    }


    @Test(dataProvider = "provideForTestExtractGtfGz")
    void testExtractGtfGz(String gtfGzFilePath, String decompressedFilePath) {
        FuncotatorDataSourceBundlerUtils.extractGzFile(gtfGzFilePath, decompressedFilePath, true);
    }

    //To do: write unit tests for testing with an invalid organism name and for testing with an invalid species name
}
