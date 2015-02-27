/*
* Copyright (c) 2012 The Broad Institute
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.hellbender.utils.R;

import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

/**
 * Basic unit test for RScriptExecutor in reduced reads
 */
public class RScriptExecutorUnitTest extends BaseTest {

    private static final String HELLO_WORLD_SCRIPT = "print('hello, world')";
    private static final String GSALIB_LOADED_SCRIPT = "if (!'package:gsalib' %in% search()) stop('gsalib not loaded')";

    @Test
    public void testRscriptExists() {
        Assert.assertTrue(RScriptExecutor.RSCRIPT_EXISTS, "Rscript not found in environment ${PATH}");
    }

    @Test(dependsOnMethods = "testRscriptExists")
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

    @Test(dependsOnMethods = "testRscriptExists", expectedExceptions = RScriptExecutorException.class)
    public void testNonExistantScriptException() {
        RScriptExecutor executor = new RScriptExecutor();
        executor.addScript(new File("does_not_exists.R"));
        executor.exec();
    }

    @Test(dependsOnMethods = "testRscriptExists")
    public void testNonExistantScriptNoException() {
        logger.warn("Testing that warning is printed an no exception thrown for missing script.");
        RScriptExecutor executor = new RScriptExecutor();
        executor.addScript(new File("does_not_exists.R"));
        executor.setIgnoreExceptions(true);
        Assert.assertFalse(executor.exec(), "Exec should have returned false when the job failed");
    }

    @Test(dependsOnMethods = "testRscriptExists")
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

    @Test(dependsOnMethods = "testRscriptExists", expectedExceptions = RScriptExecutorException.class)
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
