package org.broadinstitute.hellbender.utils.R;

import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

/**
 * Basic unit test for RScriptExecutor in reduced reads
 */
public final class RScriptExecutorUnitTest extends GATKBaseTest {

    private static final String HELLO_WORLD_SCRIPT = "print('hello, world')";
    private static final String GSALIB_LOADED_SCRIPT = "if (!'package:gsalib' %in% search()) stop('gsalib not loaded')";

    @Test(groups = {"R"})
    public void testRscriptExists() {
        Assert.assertTrue(new RScriptExecutor().externalExecutableExists(), "Rscript not found in environment ${PATH}");
    }

    @Test(groups = {"R"})
    public void testRscriptEnsureExists() {
        Assert.assertNotNull(new RScriptExecutor(true), "Rscript not found in environment ${PATH}");
    }

    @Test(groups = {"R"}, dependsOnMethods = "testRscriptExists")
    public void testExistingScript() {
        File script = writeScript(HELLO_WORLD_SCRIPT);
        try {
            RScriptExecutor executor = new RScriptExecutor();
            executor.addScript(script);
            Assert.assertTrue(executor.exec(), "Exec failed");
        } finally {
            FileUtils.deleteQuietly(script);
        }
    }

    @Test(groups = {"R"}, dependsOnMethods = "testRscriptExists", expectedExceptions = RScriptExecutorException.class)
    public void testNonExistentScriptException() {
        RScriptExecutor executor = new RScriptExecutor();
        executor.addScript(BaseTest.getSafeNonExistentFile("does_not_exists.R"));
        executor.exec();
    }

    @Test(groups = {"R"}, dependsOnMethods = "testRscriptExists")
    public void testNonExistentScriptNoException() {
        logger.warn("Testing that warning is printed an no exception thrown for missing script.");
        RScriptExecutor executor = new RScriptExecutor();
        executor.addScript(BaseTest.getSafeNonExistentFile("does_not_exists.R"));
        executor.setIgnoreExceptions(true);
        Assert.assertFalse(executor.exec(), "Exec should have returned false when the job failed");
    }

    @Test(groups = {"R"}, dependsOnMethods = "testRscriptExists")
    public void testLibrary() {
        File script = writeScript(GSALIB_LOADED_SCRIPT);
        try {
            RScriptExecutor executor = new RScriptExecutor();
            executor.addScript(script);
            executor.addLibrary(RScriptLibrary.GSALIB);
            Assert.assertTrue(executor.exec(), "Exec failed");
        } finally {
            FileUtils.deleteQuietly(script);
        }
    }

    @Test(groups = {"R"}, dependsOnMethods = "testRscriptExists", expectedExceptions = RScriptExecutorException.class)
    public void testLibraryMissing() {
        File script = writeScript(GSALIB_LOADED_SCRIPT);
        try {
            RScriptExecutor executor = new RScriptExecutor();
            executor.addScript(script);
            // GSALIB is not added nor imported in the script
            executor.exec();
        } finally {
            FileUtils.deleteQuietly(script);
        }
    }

    private File writeScript(String content) {
        return IOUtils.writeTempFile(content, "myTestScript", ".R");
    }
}
