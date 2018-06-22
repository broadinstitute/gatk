package org.broadinstitute.hellbender.utils.R;

import org.broadinstitute.hellbender.testutils.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

public final class RScriptLibraryUnitTest extends BaseTest {
    @Test(groups = {"R"})
    public void testProperties() {
        Assert.assertEquals(RScriptLibrary.GSALIB.getLibraryName(), "gsalib");
        Assert.assertEquals(RScriptLibrary.GSALIB.getResourcePath(), "gsalib.tar.gz");
    }

    @Test(groups = {"R"})
    public void testWriteTemp() {
        final File file = RScriptLibrary.GSALIB.writeTemp();
        Assert.assertTrue(file.exists(), "R library was not written to temp file: " + file);
    }

    @Test(groups = {"R"})
    public void testWriteLibraryToTempFileInDir() {
        final RScriptLibrary library = RScriptLibrary.GSALIB;
        final File tempLibSourceDir = createTempDir("testWriteLibraryToTempFileInDir");
        final File tempLibraryFile = library.writeLibraryToTempFile(tempLibSourceDir);
        Assert.assertTrue(tempLibraryFile.exists(), "R library was not written to temp file: " + tempLibraryFile);
        Assert.assertEquals(tempLibraryFile.getParentFile().getAbsolutePath(), (tempLibSourceDir.getAbsolutePath()),
                "R library was not written to temp file: " + tempLibraryFile + " in dir: " + tempLibSourceDir);
    }

}
