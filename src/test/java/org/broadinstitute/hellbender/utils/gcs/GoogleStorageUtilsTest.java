package org.broadinstitute.hellbender.utils.gcs;

import com.google.cloud.storage.contrib.nio.CloudStorageConfiguration;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.config.ConfigFactory;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.net.URI;
import java.nio.file.Path;
import java.nio.file.Paths;

public final class GoogleStorageUtilsTest extends GATKBaseTest {

    static {
        GoogleStorageUtils.setGlobalNIODefaultOptions(
            ConfigFactory.getInstance().getGATKConfig().gcsMaxRetries(),
            ConfigFactory.getInstance().getGATKConfig().gcsProjectForRequesterPays());
    }

    @Test(groups={"bucket"})
    public void testIsCloudStorageURL(){
        Assert.assertTrue(GoogleStorageUtils.isCloudStorageUrl("gs://abucket/bucket"));
        Assert.assertFalse(GoogleStorageUtils.isCloudStorageUrl("hdfs://namenode/path/to/file"));
        Assert.assertFalse(GoogleStorageUtils.isCloudStorageUrl("localFile"));

        Assert.assertTrue(GoogleStorageUtils.isCloudStorageUrl(Paths.get(URI.create("gs://abucket/bucket"))));
        // We cannot run this one because the HDFS provider looks for the "namenode" machine
        // and throws an exception.
        //Assert.assertFalse(GoogleStorageUtils.isCloudStorageUrl(Paths.get(URI.create("hdfs://namenode/path/to/file"))));
        Assert.assertFalse(GoogleStorageUtils.isCloudStorageUrl(Paths.get("localFile")));

        // this does not throw NullPointerException.
        String x = "" + null + "://";
    }

    @Test
    public void testGetCloudStorageConfiguration() {
        String mockProject = "yes";
        int mockReopens = 100;
        CloudStorageConfiguration config = GoogleStorageUtils.getCloudStorageConfiguration(mockReopens, mockProject);
        Assert.assertEquals(config.maxChannelReopens(), mockReopens);
        Assert.assertEquals(config.userProject(), mockProject);
    }

    @Test(groups={"bucket"})
    public void testGetPathOnGcsDirectory() throws Exception {
        final String dirPath = "gs://bucket/my/dir/";
        final Path pathOnGcs = GoogleStorageUtils.getPathOnGcs(dirPath);
        Assert.assertEquals(pathOnGcs.toUri().toString(), dirPath);
    }
}
