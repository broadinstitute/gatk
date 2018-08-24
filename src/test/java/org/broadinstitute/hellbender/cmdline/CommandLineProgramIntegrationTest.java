package org.broadinstitute.hellbender.cmdline;

import com.google.common.jimfs.Configuration;
import com.google.common.jimfs.Jimfs;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.io.IOUtils;
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
public class CommandLineProgramIntegrationTest extends CommandLineProgramTest {

    private FileSystem jimfs;

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
                // local without the file:// prefix
                {createTempDir("local1").getAbsolutePath()},
                // local with the file:// prefix
                {createTempDir("local2").toPath().toUri().toString()},
                // jimfs:// path
                {Files.createDirectory(jimfs.getPath("tmp")).toUri().toString()}
        };
    }

    @Test(dataProvider = "tmpDirs", singleThreaded = true)
    public void testTmpDirArgument(final String tmpDirArg) throws Exception {
        // store the previous tmp.dir to check that it is working
        final String previousTmpDir = System.getProperty("java.io.tmpdir");
        // expected property string and parent
        final String expectedTmDirArgString = IOUtils.getAbsolutePathWithoutFileProtocol(IOUtils.getPath(tmpDirArg));

        try {
            final Path p = (Path) new TestCreateTempPathClp().instanceMain(new String[] {
                    "--" + StandardArgumentDefinitions.TMP_DIR_NAME, tmpDirArg});
            Assert.assertTrue(Files.exists(p));
            Assert.assertEquals(IOUtils.getAbsolutePathWithoutFileProtocol(p.getParent()), expectedTmDirArgString);
            Assert.assertNotEquals(System.getProperty("java.io.tmpdir"), previousTmpDir);
            Assert.assertEquals(System.getProperty("java.io.tmpdir"), expectedTmDirArgString);
        } finally {
            // get back to the previous tmp dir
            System.setProperty("java.io.tmpdir", previousTmpDir);
        }
    }

}
