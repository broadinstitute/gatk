package org.broadinstitute.hellbender.tools.funcotator;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.nio.file.Path;

/**
 * Class to test the {@link FuncotatorDataSourceBundler}.
 * Created by Hailey on 8/3/21.
 */

public class FuncotatorDataSourceBundlerUnitTest extends CommandLineProgramTest{

    //==================================================================================================================
    // Static Variables:
    private static final String speciesName     = "speciesName";
    private static final String speciesName2    = "speciesName2";
    private static final String speciesName3    = "speciesName3";

    //==================================================================================================================
    // Data Providers:

    @DataProvider
    Object[][] provideForTestMakeFolder1() {
        return new Object[][] {
                {
                    IOUtils.getPath(speciesName)
                },
                {
                    IOUtils.getPath("")
                }
        };
    }

    @DataProvider
    Object[][] provideForTestMakeFolder2() {
        return new Object[][] {
                {
                    IOUtils.getPath(speciesName2).toAbsolutePath().toString(),
                        speciesName2
                }
        };
    }

    @DataProvider
    Object[][] provideForTestMakeFolder3() {
        return new Object[][] {
                {
                    IOUtils.getPath(speciesName3).toAbsolutePath().toString(),
                        speciesName3
                }
        };
    }
    @Test(dataProvider = "provideForTestMakeFolder1")
    public void testMakeFolder1(String folderName) {
        FuncotatorDataSourceBundler.makeFolders(folderName);
    }

    @Test(dataProvider = "provideForTestMakeFolder2", expectedExceptions = UserException.BadInput.class)
    public void testMakeFolder2(String pathName, String speciesName) {
        // Assuming we haven't run testMakeFolder1 yet, this test should cause an exception to be thrown:
        FuncotatorDataSourceBundler.makeFolder2(pathName, speciesName);
    }

    @Test(dataProvider = "provideForTestMakeFolder3", expectedExceptions = UserException.BadInput.class)
    public void testMakeFolder3(String pathName, String speciesName) {
        // Assuming we haven't run testMakeFolder1 or testMakeFolder2 yet, this test should cause an exception to be thrown:
        FuncotatorDataSourceBundler.makeFolder3(pathName, speciesName);
    }

}
