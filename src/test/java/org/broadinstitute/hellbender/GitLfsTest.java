package org.broadinstitute.hellbender;

import org.broadinstitute.hellbender.utils.text.XReadLines;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

/**
 * Test that git-lfs is properly pulling down files.
 */
public class GitLfsTest extends GATKBaseTest {

    @Test(groups = "large")
    public void testLargeFilesAreDownloaded() throws IOException {
        File exampleLargeFile = new File(publicTestDir, "large/exampleLargeFile.txt");
        Assert.assertTrue(exampleLargeFile.exists(), "Expected file:" + exampleLargeFile.getAbsoluteFile() + " to exist but it didn't");
        try( XReadLines lines = new XReadLines(exampleLargeFile)){
            Assert.assertEquals(lines.next(), "This is a file to test if git-lfs is working.  Please don't change this text.");
        }
    }
}
