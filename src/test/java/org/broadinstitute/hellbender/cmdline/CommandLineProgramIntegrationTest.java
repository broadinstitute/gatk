/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2015 Daniel Gómez-Sánchez
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

package org.broadinstitute.hellbender.cmdline;

import com.google.common.jimfs.Configuration;
import com.google.common.jimfs.Jimfs;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.nio.file.FileSystem;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.UUID;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class CommandLineProgramIntegrationTest extends GATKBaseTest {

    private static FileSystem jimfs;

    @BeforeClass
    public void setUp() {
        jimfs = Jimfs.newFileSystem(Configuration.unix());
    }

    @AfterClass
    public void tearDown() {
        try {
            jimfs.close();
        } catch (final IOException e) {
            log("Unable to close JimFS");
        }
    }

    @CommandLineProgramProperties(
            summary = "TestCreateTempPathClp",
            oneLineSummary = "TestCreateTempPathClp",
            programGroup = TestProgramGroup.class
    )
    private static final class TestCreateTempPathClp extends CommandLineProgram {

        @Override
        public Path doWork() {
            return createTempPath(UUID.randomUUID().toString(), ".txt");
        }
    }

    @DataProvider
    public Object[][] tmpDirs() throws Exception {
        return new Object[][]{
                {createTempDir("local").toPath()},
                {Files.createDirectory(jimfs.getPath("tmp"))}
        };
    }

    @Test(dataProvider = "tmpDirs", singleThreaded = true)
    public void testCreateTempPath(final Path tmpDir) {
        // store the previous tmp.dir to check that it is working
        final String previousTmpDir = System.getProperty("java.io.tmpdir");

        final Path p = (Path) new TestCreateTempPathClp().instanceMain(new String[]{
                "--" + StandardArgumentDefinitions.TMP_DIR_NAME, tmpDir.toUri().toString()});
        Assert.assertTrue(Files.exists(p));

        // get back to the previous tmp dir
        System.setProperty("java.io.tmpdir", previousTmpDir);
    }

}
