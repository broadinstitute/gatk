package org.broadinstitute.hellbender.tools.funcotator;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;


/**
 * Class to test the {@link FuncotatorDataSourceBundler}.
 * Created by Hailey on 8/3/21.
 */

public class FuncotatorDataSourceBundlerUnitTest extends CommandLineProgramTest{

    //==================================================================================================================
    // Data Providers:

    @DataProvider
    Object[][] provideForTestMakeDataSourcesFolderStructure() {
        return new Object[][] {
                { "speciesName" },
                { "test1" },
                { "test2" },
                { "an_octopus" },
        };
    }

    @Test(dataProvider = "provideForTestMakeDataSourcesFolderStructure")
    public void testMakeDataSourcesFolderStructure(final String speciesName) {

        try {
            final Path   tempDir    = Files.createTempDirectory(speciesName);
            final String folderName = tempDir.toString() + "data_sources_" + speciesName;

            FuncotatorDataSourceBundler.makeDataSourcesFolderStructure(folderName, speciesName);
        }
        catch ( final IOException ex ) {
            throw new UserException("Could not complete test!", ex);
        }
    }
}
