package org.broadinstitute.hellbender.tools.funcotator;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;


/**
 * A unit test suite for the {@link FuncotatorDataSourceBundlerUtils} class.
 * Created by Hailey on 8/5/21.
 */
public class FuncotatorDataSourceBundlerUtilsUnitTest extends CommandLineProgramTest{

    //==================================================================================================================
    // Static Variables:
    private static final String TEST_WRONG_ORGNAME = "virus";
    private static final String TEST_WRONG_SPECIESNAME = "ashbya_gossypii";
    private static final String PLANTS_NAME = "plants";
    private static final String PLANTS_SPECIES = "actinidia_chinensis";
    private static final String PROTISTS_SPECIES = "Aphanomyces astaci GCA_002197585.2";

    //==================================================================================================================
    // Data Providers:

    @DataProvider
    Object[][] provideForTestBuildMap() {
        return new Object[][] {
                {
                        PLANTS_NAME,
                        PLANTS_SPECIES,
                        PROTISTS_SPECIES
                }
        };
    }
    //==================================================================================================================
    // Tests:

    @Test(dataProvider = "provideForTestBuildMap")
    void testBuildMap(String orgName, String correctSpeciesName, String incorrectSpeciesName) {
        //String correctFileName = FuncotatorDataSourceBundlerUtils.buildMapGetFileName(orgName, correctSpeciesName, false);
        String correctFastaFileName = FuncotatorDataSourceBundlerUtils.buildMapGetFileName(orgName, correctSpeciesName, true);
        String test = "pause";
        String incorrectFileName = FuncotatorDataSourceBundlerUtils.buildMapGetFileName(orgName, incorrectSpeciesName, false);
        String incorrectFastaFileName = FuncotatorDataSourceBundlerUtils.buildMapGetFileName(orgName, incorrectSpeciesName, false);
        String test2 = "pause2";
        //String results = "Correct species name: " + correctFileName + ", Incorrect species name: " + incorrectSpeciesName;
    }

    @DataProvider
    Object[][] provideForTestExtractGtfGz() {
        return new Object[][] {
//                {
//                    IOUtils.getPath("Acyrthosiphon_pisum.Acyr_2.0.51.gtf.gz").toAbsolutePath(),
//                        IOUtils.getPath("Acyrthosiphon_pisum.Acyr_2.0.51.gtf.txt").toAbsolutePath()
//                }
                {
                        "/Users/hfox/Workspace/gatk/Acyrthosiphon_pisum.Acyr_2.0.51.gtf.gz",
                        "/Users/hfox/Workspace/gatk/Acyrthosiphon_pisum.Acyr_2.0.51.gtf.txt"
                }
        };
    }

    @Test(dataProvider = "provideForTestExtractGtfGz")
    void testExtractGtfGz(String gtfGzFilePath, String decompressedFilePath) {
        FuncotatorDataSourceBundlerUtils.extractGtfGz(gtfGzFilePath, decompressedFilePath, true);
    }
}
