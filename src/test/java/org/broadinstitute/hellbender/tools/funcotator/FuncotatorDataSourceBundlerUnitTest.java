package org.broadinstitute.hellbender.tools.funcotator;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.BaseTest;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;


/**
 * Class to test the {@link FuncotatorDataSourceBundler}.
 * Created by Hailey on 8/3/21.
 */

public class FuncotatorDataSourceBundlerUnitTest extends GATKBaseTest {

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
        final Path   tempDir    = BaseTest.getSafeNonExistentPath(speciesName);
        IOUtils.deleteOnExit(tempDir);

        final String folderName = tempDir.toString() + "/" + "data_sources_" + speciesName;
        FuncotatorDataSourceBundler.makeDataSourcesFolderStructure(folderName, speciesName, true);

        Assert.assertTrue(Files.exists(IOUtils.getPath(folderName)));
    }
}
