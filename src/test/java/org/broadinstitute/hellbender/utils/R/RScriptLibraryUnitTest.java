package org.broadinstitute.hellbender.utils.R;

import org.apache.commons.io.FileUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

public final class RScriptLibraryUnitTest {
    @Test(groups = {"R"})
    public void testProperties() {
        Assert.assertEquals(RScriptLibrary.GSALIB.getLibraryName(), "gsalib");
        Assert.assertEquals(RScriptLibrary.GSALIB.getResourcePath(), "gsalib.tar.gz");
    }

    @Test(groups = {"R"})
    public void testWriteTemp() {
        File file = RScriptLibrary.GSALIB.writeTemp();
        Assert.assertTrue(file.exists(), "R library was not written to temp file: " + file);
        FileUtils.deleteQuietly(file);
    }
}
