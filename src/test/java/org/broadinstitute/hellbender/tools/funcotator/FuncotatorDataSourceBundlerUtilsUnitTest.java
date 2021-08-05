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
import org.broadinstitute.hellbender.tools.funcotator.FuncotatorDataSourceBundlerUtils;
import org.broadinstitute.hellbender.tools.funcotator.FuncotatorDataSourceBundler;

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
    private static final String PLANTS_SPECIES = "Actinidia chinensis";
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
        String correctFileName = FuncotatorDataSourceBundlerUtils.buildMapGetFileName(orgName, correctSpeciesName);
        String test = "pause";
        String incorrectFileName = FuncotatorDataSourceBundlerUtils.buildMapGetFileName(orgName, incorrectSpeciesName);
        String test2 = "pause2";
        String results = "Correct species name: " + correctFileName + ", Incorrect species name: " + incorrectSpeciesName;
    }
}
