package org.broadinstitute.hellbender.utils.gcs;

import com.google.cloud.storage.contrib.nio.CloudStorageConfiguration;
import htsjdk.samtools.util.IOUtil;
import java.net.URI;
import java.nio.file.Path;
import java.nio.file.Paths;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.MiniClusterUtils;
import org.broadinstitute.hellbender.utils.config.ConfigFactory;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.security.GeneralSecurityException;

public final class BucketUtilsTest extends GATKBaseTest {

    static {
        BucketUtils.setGlobalNIODefaultOptions(
            ConfigFactory.getInstance().getGATKConfig().gcsMaxRetries(),
            ConfigFactory.getInstance().getGATKConfig().gcsProjectForRequesterPays());
    }

    @Test(groups={"bucket"})
    public void testIsCloudStorageURL(){
        Assert.assertTrue(BucketUtils.isCloudStorageUrl("gs://abucket/bucket"));
        Assert.assertFalse(BucketUtils.isCloudStorageUrl("hdfs://namenode/path/to/file"));
        Assert.assertFalse(BucketUtils.isCloudStorageUrl("localFile"));

        Assert.assertTrue(BucketUtils.isCloudStorageUrl(Paths.get(URI.create("gs://abucket/bucket"))));
        // We cannot run this one because the HDFS provider looks for the "namenode" machine
        // and throws an exception.
        //Assert.assertFalse(BucketUtils.isCloudStorageUrl(Paths.get(URI.create("hdfs://namenode/path/to/file"))));
        Assert.assertFalse(BucketUtils.isCloudStorageUrl(Paths.get("localFile")));

        // this does not throw NullPointerException.
        String x = "" + null + "://";
    }

    @Test
    public void testGetCloudStorageConfiguration() {
        String mockProject = "yes";
        int mockReopens = 100;
        CloudStorageConfiguration config = BucketUtils.getCloudStorageConfiguration(mockReopens, mockProject);
        Assert.assertEquals(config.maxChannelReopens(), mockReopens);
        Assert.assertEquals(config.userProject(), mockProject);
    }

    @Test
    public void testIsHadoopURL(){
        Assert.assertFalse(BucketUtils.isHadoopUrl("gs://abucket/bucket"));
        Assert.assertTrue(BucketUtils.isHadoopUrl("hdfs://namenode/path/to/file"));
        Assert.assertFalse(BucketUtils.isHadoopUrl("localFile"));
    }

    @Test
    public void testIsRemoteStorageURL(){
        Assert.assertTrue(BucketUtils.isRemoteStorageUrl("gs://abucket/bucket"));
        Assert.assertTrue(BucketUtils.isRemoteStorageUrl("hdfs://namenode/path/to/file"));
        Assert.assertFalse(BucketUtils.isRemoteStorageUrl("localFile"));
    }

    @Test
    public void testIsFileURL(){
        Assert.assertTrue(BucketUtils.isFileUrl("file:///somefile/something"));
        Assert.assertTrue(BucketUtils.isFileUrl("file:/something"));
        Assert.assertFalse(BucketUtils.isFileUrl("gs://abucket"));
    }

    @Test
    public void testCopyLocal() throws IOException {
        final String src = publicTestDir+"empty.vcf";
        File dest = createTempFile("copy-empty",".vcf");

        BucketUtils.copyFile(src, dest.getPath());
        IOUtil.assertFilesEqual(new File(src), dest);
    }

    @Test
    public void testDeleteLocal() throws IOException, GeneralSecurityException {
        File dest = createTempFile("temp-fortest",".txt");
        try (FileWriter fw = new FileWriter(dest)){
            fw.write("Goodbye, cruel world!");
        }
        BucketUtils.deleteFile(dest.getPath());
        Assert.assertFalse(dest.exists(), "File '"+dest.getPath()+"' was not properly deleted as it should.");
    }

    @Test(groups={"bucket"})
    public void testCopyAndDeleteGCS() throws IOException, GeneralSecurityException {
        final String src = publicTestDir + "empty.vcf";
        File dest = createTempFile("copy-empty", ".vcf");
        final String intermediate = BucketUtils.randomRemotePath(getGCPTestStaging(), "test-copy-empty", ".vcf");
        Assert.assertTrue(BucketUtils.isCloudStorageUrl(intermediate), "!BucketUtils.isCloudStorageUrl(intermediate)");
        BucketUtils.copyFile(src, intermediate);
        BucketUtils.copyFile(intermediate, dest.getPath());
        IOUtil.assertFilesEqual(new File(src), dest);
        Assert.assertTrue(BucketUtils.fileExists(intermediate));
        BucketUtils.deleteFile(intermediate);
        Assert.assertFalse(BucketUtils.fileExists(intermediate));
    }

    @Test(groups={"bucket"})
    public void testGetPathOnGcsDirectory() throws Exception {
        final String dirPath = "gs://bucket/my/dir/";
        final Path pathOnGcs = BucketUtils.getPathOnGcs(dirPath);
        Assert.assertEquals(pathOnGcs.toUri().toString(), dirPath);
    }

    @Test
    public void testCopyAndDeleteHDFS() throws Exception {
        final String src = publicTestDir + "empty.vcf";
        File dest = createTempFile("copy-empty", ".vcf");

        MiniClusterUtils.runOnIsolatedMiniCluster( cluster -> {
            final String intermediate = BucketUtils.randomRemotePath(MiniClusterUtils.getWorkingDir(cluster).toString(), "test-copy-empty", ".vcf");
            Assert.assertTrue(BucketUtils.isHadoopUrl(intermediate), "!BucketUtils.isHadoopUrl(intermediate)");

            BucketUtils.copyFile(src, intermediate);
            BucketUtils.copyFile(intermediate, dest.getPath());
            IOUtil.assertFilesEqual(new File(src), dest);
            Assert.assertTrue(BucketUtils.fileExists(intermediate));
            BucketUtils.deleteFile(intermediate);
            Assert.assertFalse(BucketUtils.fileExists(intermediate));
        });
    }

    @Test
    public void testDirSize() throws IOException {
        File dir = createTempDir("dir");
        File file1 = new File(dir, "file1.txt");
        File file2 = new File(dir, "file2.txt");
        File subdir = new File(dir, "sub");
        subdir.mkdir();
        File file3 = new File(subdir, "file3.txt");

        for (File file : new File[] { file1, file2, file3 }) {
            try (FileWriter fw = new FileWriter(file)){
                fw.write("Hello!");
            }
        }

        long fileSize = BucketUtils.fileSize(file1.getAbsolutePath());
        Assert.assertTrue(fileSize > 0);
        long dirSize = BucketUtils.dirSize(dir.getAbsolutePath());
        Assert.assertEquals(dirSize, fileSize * 2);
    }

    @Test(groups={"bucket"})
    public void testDirSizeGCS() throws IOException, GeneralSecurityException {
        final String src = publicTestDir + "empty.vcf";
        final String gcsSubDir = BucketUtils.randomRemotePath(getGCPTestStaging(), "dir-", "/");
        final String intermediate = BucketUtils.randomRemotePath(gcsSubDir, "test-copy-empty", ".vcf");
        BucketUtils.copyFile(src, intermediate);
        Assert.assertTrue(BucketUtils.fileExists(intermediate));

        long srcFileSize = BucketUtils.fileSize(src);
        Assert.assertTrue(srcFileSize > 0);
        long intermediateFileSize = BucketUtils.fileSize(intermediate);
        Assert.assertEquals(intermediateFileSize, srcFileSize);
        long intermediateDirSize = BucketUtils.dirSize(intermediate);
        Assert.assertEquals(intermediateDirSize, srcFileSize);
        long intermediateParentDirSize = BucketUtils.dirSize(gcsSubDir);
        Assert.assertEquals(intermediateParentDirSize, srcFileSize);

        BucketUtils.deleteFile(intermediate);
        Assert.assertFalse(BucketUtils.fileExists(intermediate));
    }

}
